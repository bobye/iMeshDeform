#include "petscvec.h"
#include "petscmat.h"
#include "petscksp.h"

#include "subspace/subspace.hh"

#include <assert.h>
//#include "pthread.h"

#include <time.h>

#include <mkl.h>
#define _SS_SINGLE // [SINGLE] precision is faster than [DOUBLE]

#ifdef  _SS_DOUBLE
#define _SS_SCALAR          double
#define _SS_FUNC(x)         d ## x
#define _SS_CBLAS_FUNC(x)   cblas_d ## x
#define _SS_LAPACKE_FUNC(x) LAPACKE_d ## x
#elif defined  _SS_SINGLE
#define _SS_SCALAR          float
#define _SS_FUNC(x)         s ## x
#define _SS_CBLAS_FUNC(x)   cblas_s ## x
#define _SS_LAPACKE_FUNC(x) LAPACKE_s ## x
#endif

// MKL is for floating-point op by 16-byte boundary ?
#define _SS_MALLOC_SCALAR(x)       (_SS_SCALAR *) mkl_malloc( (x) *sizeof(_SS_SCALAR), 16) 
#define _SS_MALLOC_INT(x)       (int *) mkl_malloc( (x) *sizeof(int), 16)
#define _SS_FREE(x)         mkl_free(x)

// Timing, count in seconds.
struct timespec start, end;
#define BILLION  1000000000L
  void clock_start(std::string description) {
    std::cout << description <<" ... " << std::flush; 
    clock_gettime(CLOCK_MONOTONIC, &start);
  }
  void clock_end() {
    clock_gettime(CLOCK_MONOTONIC, &end);
    printf("\t[done] %.3f seconds\n", (( end.tv_sec - start.tv_sec )+ (double)( end.tv_nsec - start.tv_nsec ) / (double)BILLION ));
  }

// Timing, count in nano seconds.
#ifdef _SS_SHOW_DEBUG
  struct timespec nstart, nend;
  inline void nclock_start() {clock_gettime(CLOCK_MONOTONIC, &nstart);}
  inline void nclock_end() {clock_gettime(CLOCK_MONOTONIC, &nend);     
    printf("\t%ld nanosec", ((long) ( nend.tv_sec - nstart.tv_sec ) * BILLION + ( nend.tv_nsec - nstart.tv_nsec ) ));}
#endif 

#ifdef _SS_SHOW_DEBUG
#define _SS_PROFILE(x) nclock_start(); x nclock_end();
#else
#define _SS_PROFILE(x) x
#endif

/* Computing the Singular Value Decomposition of 3x3 matrices 
 * with minimal branching and elementary floating point operations 
 * A. McAdams, A. Selle, R. Tamstorf, J. Teran and E. Sifakis
 * http://pages.cs.wisc.edu/~sifakis/project_pages/svd.html
 */
#define NUM_OF_SVD_THREAD 4 // parallel 3x3 svd
#include "fastsvd.hh"

#if 0
inline void apply_rot(float * const y, const _SS_SCALAR *x, const _SS_SCALAR *M, const char Order) {// M in row major
  if (Order == 'R') {
    y[0] = M[0]*x[0] + M[1]*x[1] + M[2]*x[2];
    y[1] = M[3]*x[0] + M[4]*x[1] + M[5]*x[2];
    y[2] = M[6]*x[0] + M[7]*x[1] + M[8]*x[2];
  } else if (Order == 'C') {
    y[0] = M[0]*x[0] + M[3]*x[1] + M[6]*x[2];
    y[1] = M[1]*x[0] + M[4]*x[1] + M[7]*x[2];
    y[2] = M[2]*x[0] + M[5]*x[1] + M[8]*x[2];
  } else {
    y[0] = x[0]; y[1] = x[1]; y[2] = x[2];
  }
}
#endif

namespace subspace {
  /* This section concerns the general variational reduced model for deformable energy.
   * Three parts: 
   * 1. Control layer, control parameters for geometry, default as model geometry
   * 2. Physical layer, elastic model in terms of geometry (triangles tetrahedrons)
   * 3. Display layer, surface mesh displayed
   */

  /*****************************************************************************/
  // Control layer
  static Mat Ctrl2Geom, Ctrl2GeomT; // default identity
  static int cn; // number of artificial controls

  /*****************************************************************************/
  // Physical layer
  
  static int vn, vn3; // vertices number, and geometric dimension
  static int en, en4; // element number
  static _SS_SCALAR totarea, avgarea; // total area and average vertex supporting area

  // number and dimension of linear proxies, rotational proxies 
  static int ln, rn, ln3, rn9; 

  //linear proxies
  // static PetscScalar*  linear_proxies; // column major vn3 x ln3

  //rotational proxies, clusters of vertices
  // static std::vector<int> rotational_proxies;

  /* variational subspace solver
   * assembly: VS, VS_L, VS_R
   * sparse solving: SL, SR
   * vectorization: SL_V, SR_V
   * reduced model: LVS, MVS, LVSR, MVSR, (NR, CR)
   */
  static Mat  VS;//sparse matrix to LU: dimension 7*vn + 3*ln
  static Mat  RE;//linearization artifacts regulazier

#define COEFF_REG_L2      (0.001) // dumping for nonrigid distortion
#define COEFF_REG_RADIO   (0.5)   // radio: normal / (overall + normal)
#define EPSILON           (1E-5)  // N/A
#define NUM_OF_ITERATION  (8)     // number of iterations of reduced model per frame

  static Vec* VS_L, *SL;//solved variational subspace: VS * SL = VS_L
  static Vec* VS_R, *SR;//row major index for each rotation matrix: VS * SR = VS_R

  static _SS_SCALAR *SL_V, *SR_V; // for mesh reconstruction, vectorization of SL and SR.
  static _SS_SCALAR *LVS, *MVS; // reduced model for linear variational subspace
  /*
  static _SS_SCALAR *LVS_DP, *MVS_DP; //reduced model for linear variational subspace (with dumping)
  static _SS_SCALAR *LVS_ND, *MVS_ND; //reduced model for linear variational subspace (without dumping)
  */

  const bool switch_dump = false; // N/A

  static _SS_SCALAR *LVSR, *MVSR;
#ifdef _SS_USE_CONFORMAL
  static _SS_SCALAR *NR, *CR;
#endif 

  static const _SS_SCALAR CFM[9]={1,1,1,1,1,1,1,1,1}; //reduced model for conformal fitting

  static _SS_SCALAR *Lin, *Rot, *Rot_b; //reduced variable, 3x3 rotation matrices are of row major
  static _SS_SCALAR GRot[9], GRot_b[9], GRot_bb[9]; // global rotation estimation

  static _SS_SCALAR *LSYS, *RHS, *RHS_hp; // dense matrix, rhs and rotation 
  static int *LSYS_piv;
  static int hn, hn3, nsys;

  /*****************************************************************************/
  // Display layer

  static _SS_SCALAR *points;
  //static _SS_SCALAR *RotNorm;
  //static float *vertices_f;
  //#define MAX_CONSTRAINT_NUMBER   100
  //  static pthread_t iterate_lin, iterate_rot;


  /*****************************************************************************/
  Subspace::Subspace(int argc, char **argv, vMesh* pm) {
    on_the_fly = true;
    PetscInitialize(&argc,&argv,(char *)0,PETSC_NULL);
    mesh = pm;
    mesh->initialize_subspace_solver();
  };

  void Subspace::set_off_fly() {
    on_the_fly = false;
  }


  void TriangleMesh::initialize_subspace_solver() {

    need_normals();
    need_neighbors();
    need_adjacentfaces();
    need_pointareas();
    //mesh->need_curvatures();

    //    std::cout << mesh->normals[10] << mesh->pdir1[10] << mesh->pdir2[10] << std::endl; exit(0);

    vn = vertices.size(); vn3 = 3*vn;
    en = vertices.size(); en4 = 4*en;

    rotational_proxies.resize(vn);
    points = _SS_MALLOC_SCALAR(vn3);
    //    vertices_f = new float[vn3];

    int count=0;
    for (int i=0; i<vn; ++i) 
      if (!isinf(pointareas[i]))
	{ totarea += pointareas[i]; ++count ;}
    avgarea = totarea/count;
#ifdef _SS_SHOW_DEBUG
    printf("Total area estimation: %e\n", totarea);
#endif
  }

  
  void TetrahedronMesh::initialize_subspace_solver() {

    need_neighbors();
    need_tetravolumes();
    need_facetareas();

    vn = nodes.size(); vn3 = 3*vn;
    en = nodes.size(); en4 = 4*en;
    rotational_proxies.resize(vn);
    points = _SS_MALLOC_SCALAR(vn3);
    
  }

  void vMesh::load_linear_proxies_vg(std::vector<int> &group_ids) {
    assert(vn == group_ids.size());
    ln = *std::max_element(group_ids.begin(), group_ids.end()) + 1; ln3 = 3*ln;
    std::vector<double> count_vertices; count_vertices.resize(ln);
    PetscScalar *linear = new PetscScalar[vn3*ln3]; 
    std::fill(linear, linear+vn3*ln3, 0);

    for (int i=0, j=0; i<vn; ++i, j+=3) {
      linear[3*group_ids[i] + j*ln3] = 1.; // row major
      linear[3*group_ids[i] +1 + (j+1)*ln3] = 1.; 
      linear[3*group_ids[i] +2 + (j+2)*ln3] = 1.; 
      ++count_vertices[group_ids[i]];
    }
    for (int i=0, j=0; i<vn; ++i, j+=3) {
      linear[3*group_ids[i] + j*ln3] /count_vertices[group_ids[i]]; // normalize
      linear[3*group_ids[i] +1 + (j+1)*ln3] /= count_vertices[group_ids[i]]; 
      linear[3*group_ids[i] +2 + (j+2)*ln3] /= count_vertices[group_ids[i]]; 
    }
    for (int i=0; i<vn3*ln3; ++i) linear_proxies.push_back(linear[i]);
  }

  void vMesh::load_rotational_proxies(std::vector<int> &group_ids) {
    assert(vn == group_ids.size());

    int rn_hard =0, rn_soft =0;

    for (int i=0; i<vn; ++i) {
      if (rn_hard <= -group_ids[i]) rn_hard = -group_ids[i];
      if (rn_soft <= group_ids[i]) rn_soft = group_ids[i] + 1;
    }    

    for (int i=0; i<vn; ++i)
      if (group_ids[i] >= 0) rotational_proxies[i] = group_ids[i];
      else rotational_proxies[i] = - group_ids[i] + rn_soft -1; 

    rn = rn_hard + rn_soft;
    rn9 = 9*rn;
    
  }

  void vMesh::load_controls(std::vector<int> &group_ids) {
    assert(vn == group_ids.size());
    clock_start("Preparing control layer");

    int M = 0;
    cn = 0;
    for (int i=0; i<vn; ++i)
      if (cn <= group_ids[i] && group_ids[i] != 0) cn = group_ids[i];
      else if(group_ids[i] == 0) ++M;
    M = 3*M; cn = M + 12*cn;

    int N = vn3 + en4 + ln3;
    int *nnz = new int [N];
    for (int i=0; i<vn; ++i) {
      if (group_ids[i] == 0) nnz[3*i] = nnz[3*i+1] = nnz[3*i+2] = 1;
      else nnz[3*i] = nnz[3*i+1] = nnz[3*i+2] = 4;
    }
    for (int i=0; i<en4; ++i) nnz[vn3 + i] = 1;    
    for (int i=0; i<ln3; ++i)  nnz[vn3 + en4 + i] = 1;


    MatCreateSeqAIJ(PETSC_COMM_SELF, N, cn + en4 + ln3, 0, nnz, &Ctrl2Geom); //delete [] nnz;

    const PetscScalar one = 1;
    for (int i=0, j=0; i<vn; ++i) 
      if (group_ids[i] == 0) {
	int ix=3*i, iy=3*i+1, iz=3*i+2;
	int jx=3*j, jy=3*j+1, jz=3*j+2;
	MatSetValues(Ctrl2Geom, 1, &ix, 1, &jx, &one, INSERT_VALUES);
	MatSetValues(Ctrl2Geom, 1, &iy, 1, &jy, &one, INSERT_VALUES);
	MatSetValues(Ctrl2Geom, 1, &iz, 1, &jz, &one, INSERT_VALUES);	
	++j;
      }
      else {
	int m = M + 12* (group_ids[i] -1);
	int ix=3*i, iy=3*i+1, iz=3*i+2;
	int mx[4] = {m, m+3, m+4, m+5};
	int my[4] = {m+1, m+6, m+7, m+8};
	int mz[4] = {m+2, m+9, m+10, m+11};	
	PetscScalar mv[4] = {1, vertices_tpd[3*i], vertices_tpd[3*i+1], vertices_tpd[3*i+2]};

	MatSetValues(Ctrl2Geom, 1, &ix, 4, mx, mv, INSERT_VALUES);
	MatSetValues(Ctrl2Geom, 1, &iy, 4, my, mv, INSERT_VALUES);
	MatSetValues(Ctrl2Geom, 1, &iz, 4, mz, mv, INSERT_VALUES);

      }

    for (int i=3*vn, j=cn; i< N; ++i, ++j)
      MatSetValues(Ctrl2Geom, 1, &i, 1, &j, &one, INSERT_VALUES);

    MatAssemblyBegin(Ctrl2Geom,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(Ctrl2Geom,MAT_FINAL_ASSEMBLY);

    //MatView(Ctrl2Geom, PETSC_VIEWER_STDOUT_SELF);
    MatTranspose(Ctrl2Geom, MAT_INITIAL_MATRIX, &Ctrl2GeomT);
    clock_end();
  }


  void Subspace::allocate() {
#ifdef _SS_USE_CONFORMAL
    if (!NR) NR = _SS_MALLOC_SCALAR(rn); 
    if (!CR) CR = _SS_MALLOC_SCALAR(rn);
#endif

    if (!SL_V) SL_V = _SS_MALLOC_SCALAR(vn3*ln3); 
    if (!SR_V) SR_V = _SS_MALLOC_SCALAR(vn3*rn9);
    if (!LVS) LVS = _SS_MALLOC_SCALAR (ln3*ln3); 
    if (!MVS) MVS = _SS_MALLOC_SCALAR(ln3*rn9);
    if (!LVSR) LVSR = _SS_MALLOC_SCALAR (rn9*ln3); 
    if (!MVSR) MVSR = _SS_MALLOC_SCALAR(rn9*rn9);
    if (!Rot) Rot = _SS_MALLOC_SCALAR(rn9); 
    if (!Rot_b) Rot_b = _SS_MALLOC_SCALAR(rn9); //RotNorm = _SS_MALLOC_SCALAR(2*rn);

    std::fill(Rot, Rot+rn9, 0);
    for (int i=0; i<rn; ++i) Rot[i] = Rot[i+4*rn] = Rot[i+8*rn] = 1.;
    init_svd(Rot, rn, NUM_OF_SVD_THREAD);

    GRot[0]=GRot[4]=GRot[8]=1;
#ifdef _SS_USE_CONFORMAL
    std::fill(CR, CR+rn, 1);    std::fill(NR, NR+rn, 0);    
#endif
  }

#define MULTIPLY(v,n,w) cblas_dscal(n, w, v, 1);

  inline PetscErrorCode mat_edge_assembly_VS(const PetscInt &v0, const PetscInt &v1, const PetscInt &k, const PetscInt &rk, const PetscScalar &weight, const Vector &v) {
    PetscErrorCode ierr;
    const PetscInt idv[3][2] = {{3*v0, 3*v1}, {3*v0+1, 3*v1+1}, {3*v0+2, 3*v1+2}};    
    const PetscInt idq[4] = {vn3+ 4*k, vn3 + 4*k+1, vn3 + 4*k+2, vn3 + 4*k+3}; 

    PetscScalar vvs[4] = {1, -1, -1, 1}; 
    MULTIPLY(vvs, 4, weight)

    ierr = MatSetValues(VS, 2, idv[0], 2, idv[0], vvs, ADD_VALUES); CHKERRQ(ierr);
    ierr = MatSetValues(VS, 2, idv[1], 2, idv[1], vvs, ADD_VALUES); CHKERRQ(ierr);
    ierr = MatSetValues(VS, 2, idv[2], 2, idv[2], vvs, ADD_VALUES); CHKERRQ(ierr);
    
    PetscScalar vqs[3][8] = {{v[0], 0, v[2], -v[1], -v[0], 0, -v[2], v[1]},
			     {v[1], -v[2], 0, v[0], -v[1], v[2], 0, -v[0]},
			     {v[2], v[1], -v[0], 0, -v[2], -v[1], v[0], 0}};
    MULTIPLY(vqs[0], 8, -weight)
    MULTIPLY(vqs[1], 8, -weight)
    MULTIPLY(vqs[2], 8, -weight)

    ierr = MatSetValues(VS, 2, idv[0], 4, idq, vqs[0], ADD_VALUES); CHKERRQ(ierr);
    ierr = MatSetValues(VS, 2, idv[1], 4, idq, vqs[1], ADD_VALUES); CHKERRQ(ierr);
    ierr = MatSetValues(VS, 2, idv[2], 4, idq, vqs[2], ADD_VALUES); CHKERRQ(ierr);

    PetscScalar qqs[16] = {v[0]*v[0] + v[1]*v[1] + v[2]*v[2], 0, 0, 0,
			   0, v[2]*v[2]+v[1]*v[1], -v[0]*v[1], -v[0]*v[2],
			   0, -v[1]*v[0], v[0]*v[0]+v[2]*v[2], -v[1]*v[2],
			   0, -v[2]*v[0], -v[2]*v[0], v[1]*v[1]+v[0]*v[0]};
    MULTIPLY(qqs, 16, weight)
    ierr = MatSetValues(VS, 4, idq, 4, idq, qqs, ADD_VALUES); CHKERRQ(ierr);


    for (int i=0; i<3; ++i) {
      PetscScalar vs[2] ={v[i], -v[i]};     
      MULTIPLY(vs, 2, weight)
      for (int j=0; j<3; ++j)
	VecSetValues(VS_R[rk+rn*(i+3*j)], 2, idv[j], vs, ADD_VALUES);
      }

    PetscScalar vs[4];
    for (int i=0; i<3; ++i) {
      vs[0] = -v[i]*v[0]; vs[1] = 0; vs[2] = -v[i]*v[2]; vs[3]=v[i]*v[1];
      MULTIPLY(vs, 4, weight)
      VecSetValues(VS_R[rk+rn*i], 4, idq, vs, ADD_VALUES);
    }
    for (int i=0; i<3; ++i) {
      vs[0] = -v[i]*v[1]; vs[1] = v[i]*v[2]; vs[2] = 0; vs[3]=-v[i]*v[0];
      MULTIPLY(vs, 4, weight)
      VecSetValues(VS_R[rk+rn*(i+3)], 4, idq, vs, ADD_VALUES);
    }
    for (int i=0; i<3; ++i) {
      vs[0] = -v[i]*v[2]; vs[1] = -v[i]*v[1]; vs[2] = v[i]*v[0]; vs[3]=0;
      MULTIPLY(vs, 4, weight)
      VecSetValues(VS_R[rk+rn*(i+6)], 4, idq, vs, ADD_VALUES);
    }
#ifdef _SS_USE_CONFORMAL
    NR[rk] += weight *(v DOT v);
#endif
    return ierr;   
  }

  void TriangleMesh::compute_ARAP_approx() {
    for (int i = 0; i<vn; ++i) {
      //trimesh::vec &normal = mesh->normals[i], &pdir1 = mesh->pdir1[i], &pdir2 = mesh->pdir2[2];

      std::vector<int> &nfaces = adjacentfaces[i];
      int fn = nfaces.size();

      PetscScalar area = pointareas[i];
      if (isinf(area) || area < 1E-3 * avgarea) area = 1E-3 * avgarea;

      for (int j= 0; j<fn; ++j) {
	int v0 = faces[nfaces[j]][0], v1 = faces[nfaces[j]][1], v2 = faces[nfaces[j]][2];
	Vector v01 = vertices[v0] - vertices[v1];
	Vector v12 = vertices[v1] - vertices[v2];
	Vector v20 = vertices[v2] - vertices[v0];
	
	
	double tan2 = std::tan(_SS_PI/2 - cornerangle(nfaces[j], 2)); 	
	double tan0 = std::tan(_SS_PI/2 - cornerangle(nfaces[j], 0));	
	double tan1 = std::tan(_SS_PI/2 - cornerangle(nfaces[j], 1));

	mat_edge_assembly_VS(v0, v1, i, rotational_proxies[i], std::fabs(tan2)/avgarea, v01);
	mat_edge_assembly_VS(v1, v2, i, rotational_proxies[i], std::fabs(tan0)/avgarea, v12);
	mat_edge_assembly_VS(v2, v0, i, rotational_proxies[i], std::fabs(tan1)/avgarea, v20);

      }

      /* Dumping radio: overall distortion vs normal directional distortion
       * which only makes sense to surface meshes 
       */

      const PetscInt idq[4] = {vn3+ 4*i, vn3 + 4*i+1, vn3 + 4*i+2, vn3 + 4*i+3}; 
      Vector n = normals[i];
      PetscScalar vqs[16] = {1, 0, 0, 0,
			     0, n[2]*n[2]+n[1]*n[1], -n[0]*n[1], -n[0]*n[2],
			     0, -n[1]*n[0], n[0]*n[0]+n[2]*n[2], -n[1]*n[2],
			     0, -n[2]*n[0], -n[2]*n[0], n[1]*n[1]+n[0]*n[0]};
      MULTIPLY(vqs, 16, COEFF_REG_RADIO)
      vqs[0] += (1-COEFF_REG_RADIO);
      for (int j=1; j<4; ++j) vqs[5*j] += 2*(1-COEFF_REG_RADIO)/3.;

      MULTIPLY(vqs, 16, COEFF_REG_L2 * area/avgarea)
      MatSetValues(VS, 4, idq, 4, idq, vqs, ADD_VALUES);

    }

  }

  void TetrahedronMesh::compute_ARAP_approx() {
    for (int i = 0; i<vn; ++i) {
      //std::vector<int> &nelements = adjacentelements[i];
    }
  }

  void vMesh::compute_subspace_assembly() {
    int N = vn3 + en4 + ln3; 
    int *nnz = new int[N];

    for (int i=0; i<vn3; ++i)  nnz[i] = 5 * (get_numneighbors(i/3)+1) + 3*ln;
    for (int i=0; i<en4; ++i)  nnz[vn3 + i] = 4;// + mesh->neighbors[i/4].size();
    for (int i=0; i<ln3; ++i)  nnz[vn3+en4 + i] = 1;

    MatCreateSeqSBAIJ(PETSC_COMM_SELF, 1, N, N, 0, nnz, &VS); delete [] nnz;
    MatSetOption(VS, MAT_IGNORE_LOWER_TRIANGULAR, PETSC_TRUE);

    MatCreateSeqSBAIJ(PETSC_COMM_SELF, 1, N, N, 1, PETSC_NULL, &RE);

    // create LHS of subspace problem
    VS_L = new Vec[ln3];
    VS_R = new Vec[rn9];
    for (int i=0; i<ln3; ++i) MatGetVecs(VS, &VS_L[i], PETSC_NULL);
    for (int i=0; i<rn9; ++i) MatGetVecs(VS, &VS_R[i], PETSC_NULL);

    for (int i=0; i< ln3; ++i) { 
      PetscInt index = i+vn3+en4; PetscScalar one = 1.; 
      VecSet(VS_L[i], 0.);
      VecSetValues(VS_L[i], 1, &index, &one, INSERT_VALUES);
      VecAssemblyBegin(VS_L[i]);
      VecAssemblyEnd(VS_L[i]);
    }
    for (int i=0; i< rn9; ++i) {
      VecSet(VS_R[i], 0.);
    }
    // assembly VS
    compute_ARAP_approx();

    // assembly VH
    int *irow = new int[vn3], *icol = new int[ln3];
    for (int i=0; i<vn3; ++i) irow[i] = i;
    for (int i=0; i<ln3; ++i) icol[i] = vn3 + en4 + i;

    MatSetValues(VS, vn3, irow, ln3, icol, &linear_proxies[0], ADD_VALUES);
    delete [] irow; delete [] icol;

    PetscScalar zero = 0;
    for (int i=vn3 + en4; i<N; ++i) MatSetValues(VS, 1, &i, 1, &i, &zero, ADD_VALUES);



    MatAssemblyBegin(VS,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(VS,MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(RE,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(RE,MAT_FINAL_ASSEMBLY);


    for (int i=0; i< rn9; ++i) {
      VecAssemblyBegin(VS_R[i]);
      VecAssemblyEnd(VS_R[i]);
    }
  }



  void Subspace::solve() {
    //clock_start("Solving reduced model");
    /**************************************************/
    // assembly matrix
    clock_start("Assembly Physical layer");
    allocate();
    //assembly();
    mesh->compute_subspace_assembly();
    clock_end();
    //    MatAXPY(VS, COEFF_REG, RE, DIFFERENT_NONZERO_PATTERN);
    /**************************************************/
    // solve sparse system

    // convert to control layer
    clock_start("Preparing offline sparse system");
    Mat L, Ltmp, Ltmp2; 
    MatConvert(VS, MATSEQAIJ, MAT_INITIAL_MATRIX, &L);    

    MatMatMult(Ctrl2GeomT, L, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &Ltmp);
    MatTranspose(Ltmp, MAT_INITIAL_MATRIX, &Ltmp2);

    MatMatMult(Ctrl2GeomT, Ltmp2, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &L);

    MatDestroy(&Ltmp); MatDestroy(&Ltmp2);

    KSP ksp;
    PC direct_solver;
    Vec * cVS_L = new Vec[ln3], * cVS_R = new Vec[rn9];
    

    for (int i=0; i<ln3; ++i) MatGetVecs(L, &cVS_L[i], PETSC_NULL);
    for (int i=0; i<rn9; ++i) MatGetVecs(L, &cVS_R[i], PETSC_NULL);
    for (int i=0; i<ln3; ++i) MatMult(Ctrl2GeomT, VS_L[i], cVS_L[i]);
    for (int i=0; i<rn9; ++i) MatMult(Ctrl2GeomT, VS_R[i], cVS_R[i]);

    SL = new Vec[ln3]; SR = new Vec[rn9];
    

    VecDuplicateVecs(cVS_L[0], ln3, &SL); VecDuplicateVecs(cVS_R[0], rn9, &SR);
    clock_end();

    clock_start("Solving offline sparse system");
    KSPCreate(PETSC_COMM_SELF, &ksp);
    KSPSetOperators(ksp, L, L, SAME_PRECONDITIONER);
    KSPSetType(ksp, KSPPREONLY);
    KSPGetPC(ksp, &direct_solver);
    PCSetType(direct_solver, PCLU); // use LU facterization to solve the subspace
    //KSPSetType(ksp, KSPPREONLY);
    KSPSetFromOptions(ksp);

    KSPSetUp(ksp);

    for (int i=0; i<ln3; ++i) KSPSolve(ksp, cVS_L[i], SL[i]); 
    for (int i=0; i<rn9; ++i) KSPSolve(ksp, cVS_R[i], SR[i]); 


    //MatDestroy(&L);
    for (int i=0; i<ln3; ++i) VecDestroy(&VS_L[i]); delete [] VS_L;
    for (int i=0; i<rn9; ++i) VecDestroy(&VS_R[i]); delete [] VS_R;
    KSPDestroy(&ksp);

    clock_end();

    /**************************************************/
    // copy subspace solution data to global array
    clock_start("Preparing display layer");

    Vec tmp; MatGetVecs(VS, &tmp, PETSC_NULL);

     for (int i=0; i<ln3; ++i) {
      PetscScalar *buffer;
      MatMult(Ctrl2Geom, SL[i], tmp);
      VecGetArray(tmp, &buffer);
      for (int j=0; j<vn3; ++j) SL_V[vn3*i+j] = buffer[j];
    }
    for (int i=0; i<rn9; ++i) {
      PetscScalar *buffer;
      MatMult(Ctrl2Geom, SR[i], tmp);
      VecGetArray(tmp, &buffer);
      for (int j=0; j<vn3; ++j) SR_V[vn3*i+j] = buffer[j];
    }
    VecDestroy(&tmp);

    clock_end();
    /**************************************************/
    clock_start("Precomputing online dense system");
    // precompute online dense linear system (with/without dump)
    PetscInt *indices = new PetscInt[ln3];
    PetscScalar *zeros = new PetscScalar[ln3]; std::fill(zeros, zeros + ln3, 0.);
    for (int i=0; i<ln3; ++i) indices[i] = cn + en4 + i;
    for (int i=0; i<ln3; ++i) {
      VecSetValues(SL[i], ln3, indices, zeros, INSERT_VALUES);
      VecAssemblyBegin(SL[i]);
      VecAssemblyEnd(SL[i]);
    }
    for (int i=0; i<rn9; ++i) {
      VecSetValues(SR[i], ln3, indices, zeros, INSERT_VALUES);
      VecAssemblyBegin(SR[i]);
      VecAssemblyEnd(SR[i]);
    }

    delete [] indices;
    delete [] zeros;    

    for (int i=0; i<ln3; ++i) {
      Vec tmp; VecDuplicate(SL[i], &tmp);
      MatMult(L, SL[i], tmp);
      
      for (int j=i; j<ln3; ++j) {
	PetscScalar buffer;
	VecDot(SL[j], tmp, &buffer);
	LVS[ln3*j+i] = buffer;
	LVS[ln3*i+j] = buffer;
      }
      for (int j=0; j<rn9; ++j) {
	PetscScalar buffer;
	VecDot(SL[i], cVS_R[j], &buffer);
	PetscScalar tominus;
	VecDot(SR[j], tmp, &tominus);
	MVS[ln3*j+i] = buffer - tominus;
      }
      VecDestroy(&tmp);
    }
    /**************************************************/
    // precompute rotational fitting system

    indices = new PetscInt[en4];
    zeros = new PetscScalar[en4]; std::fill(zeros, zeros + en4, 0.);
    for (int i=0; i<en4; ++i) indices[i] = cn + i;
    for (int i=0; i<ln3; ++i) {
      VecSetValues(SL[i], en4, indices, zeros, INSERT_VALUES);
      VecAssemblyBegin(SL[i]);
      VecAssemblyEnd(SL[i]);
    }
    for (int i=0; i<rn9; ++i) {
      VecSetValues(SR[i], en4, indices, zeros, INSERT_VALUES);
      VecAssemblyBegin(SR[i]);
      VecAssemblyEnd(SR[i]);
    }

    delete [] indices;
    delete [] zeros;


    for (int i=0; i<rn9; ++i) {
      for (int j=0; j<ln3; ++j) {
	PetscScalar buffer;
	VecDot(cVS_R[i], SL[j], &buffer);
	LVSR[i+j*rn9] = buffer;
      }
      for (int j=0; j<rn9; ++j) {
	PetscScalar buffer;
	VecDot(cVS_R[i], SR[j], &buffer);
	MVSR[i+j*rn9] = buffer;
      }
    }


    MatDestroy(&VS);
    MatDestroy(&L);
    for (int i=0; i<ln3; ++i) VecDestroy(&cVS_L[i]); delete [] cVS_L;
    for (int i=0; i<rn9; ++i) VecDestroy(&cVS_R[i]); delete [] cVS_R;

    clock_end();
    PetscFinalize();

    //clock_end();
    std::cout << "  linear Proxies: " << ln << "; rotational Proxies: " << rn << std::endl;

  }

  Subspace::~Subspace() {        
    //_SS_FREE(linear_proxies);
    
    _SS_FREE(SL_V); 
    _SS_FREE(SR_V);
    _SS_FREE(LVS); _SS_FREE(MVS);
    _SS_FREE(LVSR); _SS_FREE(MVSR);
#ifdef _SS_USE_CONFORMAL
    _SS_FREE(NR); _SS_FREE(CR);
#endif

    _SS_FREE(Rot); _SS_FREE(Rot_b); //_SS_FREE(RotNorm);

    _SS_FREE(points);
    PetscFinalize();
    //delete [] vertices_f;
  }

  void Subspace::prepare(std::vector< std::vector<float> > & constraints,
			 std::vector< Point > & constraint_points) { 
    // precompute LU for dense direct solver, initialize Rot and Lin
    hn = constraints.size(); if (hn<3) {std::cout << "Need more constraints!\n" << std::endl; return;} hn3 =3*hn;

    clock_start("Prepare reduced model");
    nsys = ln3 + hn3;

    _SS_SCALAR *constraints_matrix = _SS_MALLOC_SCALAR(hn3*vn3);
    _SS_SCALAR *one = _SS_MALLOC_SCALAR(hn3*hn3);

    std::fill(constraints_matrix, constraints_matrix + hn3*vn3, 0);
    for (int j=0; j < vn3; j+=3)
      for (int i=0; i < hn3; i+=3)
	constraints_matrix[(i+2) + (j+2)*hn3] = constraints_matrix[i+1 + (j+1)*hn3] = constraints_matrix[i + j*hn3] = constraints[i/3][j/3]; 
      


    LSYS =  _SS_MALLOC_SCALAR (nsys * nsys); RHS =   _SS_MALLOC_SCALAR (nsys * rn9);    
    LSYS_piv = _SS_MALLOC_INT(nsys);

    std::fill(LSYS, LSYS+nsys*nsys, 0);
    std::fill(RHS, RHS+nsys*rn9, 0);
    std::fill(one, one+hn3*hn3, 0); for (int i=0;i<hn3*hn3;i+=hn3+1) one[i] =1;

    _SS_FUNC(lacpy)("F", &ln3, &ln3, LVS, &ln3, LSYS, &nsys); // copy left-top ln3 x ln3 block from LVS
    _SS_CBLAS_FUNC(gemm)(CblasColMajor, CblasNoTrans, CblasNoTrans, hn3, ln3, vn3, 1, constraints_matrix, hn3, SL_V, vn3, 0, LSYS+ln3, nsys); // compute left-bottom hn3 x ln3 block from online constraints

    _SS_CBLAS_FUNC(gemm)(CblasColMajor, CblasTrans, CblasNoTrans, ln3, hn3, hn3, 1, LSYS+ln3, nsys, one, hn3, 0, LSYS+ln3*nsys, nsys); // copy left-bottom block to right-top block by transpose

    _SS_FUNC(lacpy)("F", &ln3, &rn9, MVS, &ln3, RHS, &nsys); // copy to first ln3 rows of RHS from MVS

    _SS_CBLAS_FUNC(gemm)(CblasColMajor, CblasNoTrans, CblasNoTrans, hn3, rn9, vn3, -1, constraints_matrix, hn3, SR_V, vn3, 0, RHS+ln3, nsys); // compute last hn3 rows of RHS

    _SS_LAPACKE_FUNC(getrf)(LAPACK_COL_MAJOR, nsys, nsys, LSYS, nsys, LSYS_piv);// compute LU of LSYS

    ready = true;

    RHS_hp = _SS_MALLOC_SCALAR(hn3);
    Lin  = _SS_MALLOC_SCALAR(nsys);

    _SS_FREE(constraints_matrix); _SS_FREE(one);
    update(constraint_points, true);
    clock_end();
  }



  void reduced_linsolve() {//solve reduced linear variables via dense direct solve
    _SS_CBLAS_FUNC(gemv)(CblasColMajor, CblasNoTrans, nsys, rn9, 1, RHS, nsys, Rot, 1, 0, Lin, 1);
    _SS_CBLAS_FUNC(gemm)(CblasColMajor, CblasTrans, CblasNoTrans, 3, hn, 3, 1, GRot, 3, RHS_hp, 3, 1, Lin+ln3, 3);
    //_SS_CBLAS_FUNC(axpy)(hn3, 1, RHS_hp, 1, Lin+ln3, 1);
    _SS_LAPACKE_FUNC(getrs)(LAPACK_COL_MAJOR, 'N', nsys, 1, LSYS, nsys, LSYS_piv, Lin, nsys);
  }

  //#define ROT_STEP_SIZE  10/(totarea * COEFF_REG)
  void reduced_rotsolve() {//solve reduced rotational variables via SVD of gradient
    _SS_CBLAS_FUNC(gemv)(CblasColMajor, CblasNoTrans, rn9, rn9, 1, MVSR, rn9, Rot, 1, 0, Rot_b, 1);
    _SS_CBLAS_FUNC(gemv)(CblasColMajor, CblasNoTrans, rn9, ln3, 1, LVSR, rn9, Lin, 1, 1, Rot_b, 1);
#ifdef _SS_USE_CONFORMAL
    for (int i=0; i<rn; ++i) {
      CR[i] = _SS_CBLAS_FUNC(dot)(9, Rot_b+i, rn, Rot+i, rn) / (CR[i] * NR[i]);
    }
    for (int i=0; i<rn; ++i) 
      if (CR[i]>MAX_CFM_STRECH) CR[i] = MAX_CFM_STRECH;
      else if (CR[i]<1./MAX_CFM_STRECH) CR[i] = 1./MAX_CFM_STRECH;
    //for (int i=0; i<rn; ++i) {printf("%f ", CR[i]);}
#endif
    //apply global rotations
    for (int i=0; i<9; ++i) 
      {_SS_SCALAR sum=0; for (int j=i*rn; j< i*rn+rn; ++j) sum += Rot_b[j]; GRot_b[i] = sum;}
    proj_rot(GRot_b, 1); 

    _SS_CBLAS_FUNC(copy)(9, GRot, 1, GRot_bb, 1);
    _SS_CBLAS_FUNC(gemm)(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3, 3, 3, 1, GRot_b, 3, GRot_bb, 3, 0, GRot, 3);

    proj_rot(GRot, 1); //proj_rot(GRot, 1); 
    std::swap<_SS_SCALAR>(GRot[1], GRot[3]);
    std::swap<_SS_SCALAR>(GRot[2], GRot[6]);
    std::swap<_SS_SCALAR>(GRot[5], GRot[7]);


    _SS_CBLAS_FUNC(gemm)(CblasColMajor, CblasNoTrans, CblasNoTrans, 3*rn, 3, 3, 1, Rot_b, 3*rn, GRot_b, 3, 0, Rot, 3*rn);
    
    proj_rot(Rot, rn);

#ifdef _SS_USE_CONFORMAL
    for (int j=0; j<9; ++j) {
      _SS_SCALAR *ptr=Rot + j*rn;
      for (int i=0; i<rn; ++i) ptr[i]*=CR[i];
    }
#endif
    //    for (int i=0; i<9; ++i) printf("%.3f ", Rot[10+i*rn]); printf("\n");

  }

  void update_mesh(vMesh *mesh) {
    //update mesh vertices
    // reduced variable to mesh vertices, often computational expensive
    //    _SS_PROFILE(
    _SS_CBLAS_FUNC(gemv)(CblasColMajor, CblasNoTrans, 
			 3*mesh->numberofvertices, ln3, 
			 1, 
			 SL_V, vn3, 
			 Lin, 1, 
			 0, 
			 points, 1);
    _SS_CBLAS_FUNC(gemv)(CblasColMajor, CblasNoTrans, 
			 3*mesh->numberofvertices, rn9, 
			 1, 
			 SR_V, vn3, 
			 Rot, 1, 
			 1, 
			 points, 1);
    //		)
    /*
    for (int i=0, j=0; i<vn; ++i, j+=3)
      apply_rot(mesh->vertices[i], &vertices[j], GRot, 'C');
    */
    _SS_CBLAS_FUNC(gemm)(CblasColMajor, CblasNoTrans, CblasNoTrans, 
			 3, mesh->numberofvertices, 3,
			 1, 
			 GRot, 3, 
			 points, 3, 
			 0, 
			 mesh->vertices_tpd, 3);
   
  }


  void Subspace::update(std::vector<Point> & constraint_points, bool inf){
    if (!on_the_fly && !inf) return;
    if (ready) {
      for (int i=0, j=0; j<hn; i+=3, ++j) {
	RHS_hp[i] = constraint_points[j][0];
	RHS_hp[i+1] = constraint_points[j][1];
	RHS_hp[i+2] = constraint_points[j][2];
      }
      
      int N = inf? 42: NUM_OF_ITERATION;
      //int N=1;

      for (int i=0; i<N; ++i) {
	reduced_linsolve();
	reduced_rotsolve();
      }
      reduced_linsolve();

      update_mesh(mesh); 

      //recompute normals

      if (inf) { mesh->recompute_normals(); }
    }
  }
  void Subspace::terminate() {
    _SS_FREE(LSYS); _SS_FREE(LSYS_piv);
    _SS_FREE(RHS); 

    _SS_FREE(RHS_hp); _SS_FREE(Lin);
    ready = false;
  }

  void Subspace::read(std::string ifilename) {
    std::cout << "Import subspace data from " << ifilename << std::endl;
    allocate();
    std::ifstream binfid(ifilename.c_str(), std::ios::in|std::ios::binary);
    // read sizes
    binfid.read((char*) &ln,   sizeof(int));
    binfid.read((char*) &ln3,  sizeof(int));
    binfid.read((char*) &rn,   sizeof(int));
    binfid.read((char*) &rn9,  sizeof(int));
    binfid.read((char*) &hn,   sizeof(int));
    binfid.read((char*) &hn3,  sizeof(int));

    // read reduced variable to vertices mapping
    binfid.read((char*) SL_V,  sizeof(_SS_SCALAR) * vn3 * ln3);
    binfid.read((char*) SR_V,  sizeof(_SS_SCALAR) * vn3 * rn9);
    
    // read reduced system : linear variables
    binfid.read((char*) LVS, sizeof(_SS_SCALAR) * ln3 * ln3);
    binfid.read((char*) MVS, sizeof(_SS_SCALAR) * ln3 * rn9);
    // read reduced system : rotational variables
    binfid.read((char*) LVSR, sizeof(_SS_SCALAR) * rn9 * ln3);
    binfid.read((char*) MVSR, sizeof(_SS_SCALAR) * rn9 * rn9);
#ifdef _SS_USE_CONFORMAL
    // read reduced system : scaling factors
    binfid.read((char*) NR, sizeof(_SS_SCALAR) * rn);
#endif

    binfid.close();

  }
  
  void Subspace::write(std::string ofilename) {
    std::cout << "Export subspace data to " << ofilename << std::endl;
    std::ofstream binfid(ofilename.c_str(), std::ios::out|std::ios::binary);
    // write sizes
    binfid.write((char*) &ln,   sizeof(int));
    binfid.write((char*) &ln3,  sizeof(int));
    binfid.write((char*) &rn,   sizeof(int));
    binfid.write((char*) &rn9,  sizeof(int));
    binfid.write((char*) &hn,   sizeof(int));
    binfid.write((char*) &hn3,  sizeof(int));

    // write reduced variable and full mapping
    binfid.write((char*) SL_V,  sizeof(_SS_SCALAR) * vn3 * ln3);
    binfid.write((char*) SR_V,  sizeof(_SS_SCALAR) * vn3 * rn9);
    
    // write reduced system : linear variables
    binfid.write((char*) LVS, sizeof(_SS_SCALAR) * ln3 * ln3);
    binfid.write((char*) MVS, sizeof(_SS_SCALAR) * ln3 * rn9);
    // write reduced system : rotational variables
    binfid.write((char*) LVSR, sizeof(_SS_SCALAR) * rn9 * ln3);
    binfid.write((char*) MVSR, sizeof(_SS_SCALAR) * rn9 * rn9);
    // write reduced system : global coordinate rotation
#ifdef _SS_USE_CONFORMAL
    // write reduced system : scaling factors
    binfid.write((char*) NR, sizeof(_SS_SCALAR) * rn);
#endif

    binfid.close();
  }

#ifdef _SS_SHOW_DEBUG
  void Subspace::show_debug() {
    //for (int i=0; i<9; ++i) printf("%.3f ", GRot[i]); printf("\n");
    /*
    for (int j=0; j<rn; ++j) 
      {printf("%E ", RotNorm[j]);} printf("\n");
    for (int j=0; j<rn; ++j) 
      {printf("%E ", RotNorm[j+rn]);} printf("\n");
    */
  }
#endif
}



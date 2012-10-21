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

#define _SS_MALLOC_SCALAR(x)       (_SS_SCALAR *) mkl_malloc( (x) *sizeof(_SS_SCALAR), 64)
#define _SS_MALLOC_INT(x)       (int *) mkl_malloc( (x) *sizeof(int), 64)
#define _SS_FREE(x)         mkl_free(x)

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

#define NUM_OF_SVD_THREAD 3 // parallel 3x3 svd
#include "fastsvd.hh"

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

namespace subspace {

  static int vn, vn3;
  static _SS_SCALAR totarea, avgarea;

  static int ln, rn, ln3, rn9;
  //linear proxies
  static PetscScalar*  linear_proxies; // column major
  //rotational proxies
  static std::vector<int> rotational_proxies;

  //variational subspace solver

  static Mat  VS;//sparse matrix to LU
  static Mat  RE;//linearization artifacts regulazier
#define COEFF_REG_L2      (0.001)
#define COEFF_REG_RADIO   (0.5)
#define EPSILON           1E-5
#define NUM_OF_ITERATION 8


  static Vec* VS_L, *SL;//solved variational subspace
  static Vec* VS_R, *SR;//row major index for each rotation matrix

  static _SS_SCALAR *SL_V, *SR_V; // for mesh reconstruction
  static _SS_SCALAR *LVS, *MVS; //reduced model for linear variational subspace
  /*
  static _SS_SCALAR *LVS_DP, *MVS_DP; //reduced model for linear variational subspace (with dumping)
  static _SS_SCALAR *LVS_ND, *MVS_ND; //reduced model for linear variational subspace (without dumping)
  */
  const bool switch_dump = false;

  static _SS_SCALAR *LVSR, *MVSR; //reduced model for rotational fitting

  static _SS_SCALAR *Lin, *Rot, *Rot_b; //reduced variable, 3x3 rotation matrices are of row major
  static _SS_SCALAR GRot[9], GRot_b[9], GRot_bb[9]; // global rotation estimation

  static _SS_SCALAR *LSYS, *RHS, *RHS_hp; // dense matrix, rhs and rotation 
  static int *LSYS_piv;
  static int hn, hn3, nsys;

  static _SS_SCALAR *vertices;
  //static _SS_SCALAR *RotNorm;
  //static float *vertices_f;
  //#define MAX_CONSTRAINT_NUMBER   100
  //  static pthread_t iterate_lin, iterate_rot;


  Subspace::Subspace(int argc, char **argv) {
    PetscInitialize(&argc,&argv,(char *)0,PETSC_NULL);
  };


  void Subspace::init(Mesh * pm) {
    mesh = pm;

    mesh->need_normals();
    mesh->need_neighbors();
    mesh->need_adjacentfaces();
    mesh->need_pointareas();
    //mesh->need_curvatures();

    //    std::cout << mesh->normals[10] << mesh->pdir1[10] << mesh->pdir2[10] << std::endl; exit(0);

    vn = mesh->vertices.size(); vn3 = 3*vn;
    rotational_proxies.resize(vn);
    vertices = _SS_MALLOC_SCALAR(vn3);
    //    vertices_f = new float[vn3];

    int count=0;
    for (int i=0; i<vn; ++i) 
      if (!isinf(mesh->pointareas[i]))
	{ totarea += mesh->pointareas[i]; ++count ;}
    avgarea = totarea/count;
#ifdef _SS_SHOW_DEBUG
    printf("Total area estimation: %e\n", totarea);
#endif
  }



  void Subspace::load_linear_proxies_vg(std::vector<int> &group_ids) {
    assert(vn == group_ids.size());
    ln = *std::max_element(group_ids.begin(), group_ids.end()) + 1; ln3 = 3*ln;
    std::vector<double> count_vertices; count_vertices.resize(ln);
    linear_proxies = new PetscScalar[vn3*ln3]; 
    std::fill(linear_proxies, linear_proxies+vn3*ln3, 0);

    for (int i=0, j=0; i<vn; ++i, j+=3) {
      linear_proxies[3*group_ids[i] + j*ln3] = 1.; // row major
      linear_proxies[3*group_ids[i] +1 + (j+1)*ln3] = 1.; 
      linear_proxies[3*group_ids[i] +2 + (j+2)*ln3] = 1.; 
      ++count_vertices[group_ids[i]];
    }
    for (int i=0, j=0; i<vn; ++i, j+=3) {
      linear_proxies[3*group_ids[i] + j*ln3] /= count_vertices[group_ids[i]]; // normalize
      linear_proxies[3*group_ids[i] +1 + (j+1)*ln3] /= count_vertices[group_ids[i]]; 
      linear_proxies[3*group_ids[i] +2 + (j+2)*ln3] /= count_vertices[group_ids[i]]; 
    }
  }

  void Subspace::load_rotational_proxies(std::vector<int> &group_ids) {
    assert(vn == group_ids.size());
    
    for (int i=0; i<vn; ++i)
      if (rn <= (rotational_proxies[i] = group_ids[i])) rn = group_ids[i]+1;
    rn9 = 9*rn;
    
  }

#define MULTIPLY(v,n,w) cblas_dscal(n, w, v, 1);

  inline PetscErrorCode mat_edge_assembly_VS(const PetscInt &v0, const PetscInt &v1, const PetscInt &k, const PetscScalar &weight, const Vector &v) {
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
	VecSetValues(VS_R[rotational_proxies[k]+rn*(i+3*j)], 2, idv[j], vs, ADD_VALUES);
      }

    PetscScalar vs[4];
    for (int i=0; i<3; ++i) {
      vs[0] = -v[i]*v[0]; vs[1] = 0; vs[2] = -v[i]*v[2]; vs[3]=v[i]*v[1];
      MULTIPLY(vs, 4, weight)
      VecSetValues(VS_R[rotational_proxies[k]+rn*i], 4, idq, vs, ADD_VALUES);
    }
    for (int i=0; i<3; ++i) {
      vs[0] = -v[i]*v[1]; vs[1] = v[i]*v[2]; vs[2] = 0; vs[3]=-v[i]*v[0];
      MULTIPLY(vs, 4, weight)
      VecSetValues(VS_R[rotational_proxies[k]+rn*(i+3)], 4, idq, vs, ADD_VALUES);
    }
    for (int i=0; i<3; ++i) {
      vs[0] = -v[i]*v[2]; vs[1] = -v[i]*v[1]; vs[2] = v[i]*v[0]; vs[3]=0;
      MULTIPLY(vs, 4, weight)
      VecSetValues(VS_R[rotational_proxies[k]+rn*(i+6)], 4, idq, vs, ADD_VALUES);
    }

    return ierr;   
  }

  void Subspace::assembly() {
    int N = 7*vn + 3*ln; 
    int *nnz = new int[N];
    for (int i=0; i<vn3; ++i)  nnz[i] = 5 * (mesh->neighbors[i/3].size()+1) + 3*ln;
    for (int i=0; i<4*vn; ++i) nnz[vn3 + i] = 4;// + mesh->neighbors[i/4].size();
    for (int i=0; i<ln3; ++i)  nnz[7*vn + i] = 1;

    MatCreateSeqSBAIJ(PETSC_COMM_SELF, 1, N, N, 0, nnz, &VS); delete [] nnz;
    MatSetOption(VS, MAT_IGNORE_LOWER_TRIANGULAR, PETSC_TRUE);

    MatCreateSeqSBAIJ(PETSC_COMM_SELF, 1, N, N, 1, PETSC_NULL, &RE);

    // create LHS of subspace problem
    VS_L = new Vec[ln3];
    VS_R = new Vec[rn9];
    for (int i=0; i<ln3; ++i) MatGetVecs(VS, &VS_L[i], PETSC_NULL);
    for (int i=0; i<rn9; ++i) MatGetVecs(VS, &VS_R[i], PETSC_NULL);

    for (int i=0; i< ln3; ++i) { 
      PetscInt index = i+7*vn; PetscScalar one = 1.; 
      VecSet(VS_L[i], 0.);
      VecSetValues(VS_L[i], 1, &index, &one, INSERT_VALUES);
      VecAssemblyBegin(VS_L[i]);
      VecAssemblyEnd(VS_L[i]);
    }
    for (int i=0; i< rn9; ++i) {
      VecSet(VS_R[i], 0.);
    }

    // assembly VS

    for (int i = 0; i<vn; ++i) {
      //trimesh::vec &normal = mesh->normals[i], &pdir1 = mesh->pdir1[i], &pdir2 = mesh->pdir2[2];

      std::vector<int> &faces = mesh->adjacentfaces[i];
      int fn = faces.size();
      PetscScalar area = mesh->pointareas[i];
      if (isinf(area) || area < 1E-3 * avgarea) area = 1E-3 * avgarea;

      for (int j= 0; j<fn; ++j) {
	int v0 = mesh->faces[faces[j]][0], v1 = mesh->faces[faces[j]][1], v2 = mesh->faces[faces[j]][2];
	Vector v01 = mesh->vertices[v0] - mesh->vertices[v1];
	Vector v12 = mesh->vertices[v1] - mesh->vertices[v2];
	Vector v20 = mesh->vertices[v2] - mesh->vertices[v0];
	
	
	double tan2 = std::tan(_SS_PI/2 - mesh->cornerangle(2,j)); 
	mat_edge_assembly_VS(v0, v1, i, std::fabs(tan2)/avgarea, v01);
	
	double tan0 = std::tan(_SS_PI/2 - mesh->cornerangle(0,j));
	mat_edge_assembly_VS(v1, v2, i, std::fabs(tan0)/avgarea, v12);
	
	double tan1 = std::tan(_SS_PI/2 - mesh->cornerangle(1,j));
	mat_edge_assembly_VS(v2, v0, i, std::fabs(tan1)/avgarea, v20);

      }

      const PetscInt idq[4] = {vn3+ 4*i, vn3 + 4*i+1, vn3 + 4*i+2, vn3 + 4*i+3}; 
      Vector n = mesh->normals[i];
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

    // assembly VH
    int *irow = new int[vn3], *icol = new int[ln3];
    for (int i=0; i<vn3; ++i) irow[i] = i;
    for (int i=0; i<ln3; ++i) icol[i] = 7*vn + i;

    MatSetValues(VS, vn3, irow, ln3, icol, linear_proxies, ADD_VALUES);
    delete [] irow; delete [] icol;

    PetscScalar zero = 0;
    for (int i=7*vn; i<N; ++i) MatSetValues(VS, 1, &i, 1, &i, &zero, ADD_VALUES);



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
    clock_start("Solving reduced model");
    /**************************************************/
    // assembly matrix
    assembly();
    //    MatAXPY(VS, COEFF_REG, RE, DIFFERENT_NONZERO_PATTERN);
    /**************************************************/
    // solve sparse system
    Mat L; MatConvert(VS, MATSEQAIJ, MAT_INITIAL_MATRIX, &L);
    KSP ksp;
    PC direct_solver;
    SL = new Vec[ln3]; SR = new Vec[rn9];
    VecDuplicateVecs(VS_L[0], ln3, &SL); VecDuplicateVecs(VS_R[0], rn9, &SR);
    
    KSPCreate(PETSC_COMM_SELF, &ksp);
    KSPSetOperators(ksp, L, L, SAME_PRECONDITIONER);
    KSPSetType(ksp, KSPPREONLY);
    KSPGetPC(ksp, &direct_solver);
    PCSetType(direct_solver, PCLU); // use LU facterization to solve the subspace
    //KSPSetType(ksp, KSPPREONLY);
    KSPSetFromOptions(ksp);
    KSPSetUp(ksp);
    
    for (int i=0; i<ln3; ++i) KSPSolve(ksp, VS_L[i], SL[i]); 
    for (int i=0; i<rn9; ++i) KSPSolve(ksp, VS_R[i], SR[i]); 


    MatDestroy(&L);
    KSPDestroy(&ksp);
    /**************************************************/
    // copy subspace solution data to global array

    SL_V = _SS_MALLOC_SCALAR(vn3*ln3); SR_V = _SS_MALLOC_SCALAR(vn3*rn9);
    for (int i=0; i<ln3; ++i) {
      PetscScalar *buffer;
      VecGetArray(SL[i], &buffer);
      for (int j=0; j<vn3; ++j) SL_V[vn3*i+j] = buffer[j];
    }
    for (int i=0; i<rn9; ++i) {
      PetscScalar *buffer;
      VecGetArray(SR[i], &buffer);
      for (int j=0; j<vn3; ++j) SR_V[vn3*i+j] = buffer[j];
    }

    /**************************************************/
    // precompute online dense linear system (with/without dump)
    LVS = _SS_MALLOC_SCALAR (ln3*ln3); MVS = _SS_MALLOC_SCALAR(ln3*rn9);

    PetscInt *indices = new PetscInt[ln3];
    PetscScalar *zeros = new PetscScalar[ln3]; std::fill(zeros, zeros + ln3, 0.);
    for (int i=0; i<ln3; ++i) indices[i] = 7*vn + i;
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
      MatMult(VS, SL[i], tmp);
      
      for (int j=i; j<ln3; ++j) {
	PetscScalar buffer;
	VecDot(SL[j], tmp, &buffer);
	LVS[ln3*j+i] = buffer;
	LVS[ln3*i+j] = buffer;
      }
      for (int j=0; j<rn9; ++j) {
	PetscScalar buffer;
	VecDot(SL[i], VS_R[j], &buffer);
	PetscScalar tominus;
	VecDot(SR[j], tmp, &tominus);
	MVS[ln3*j+i] = buffer - tominus;
      }
      VecDestroy(&tmp);
    }


    /**************************************************/
    // precompute rotational fitting system
    LVSR = _SS_MALLOC_SCALAR (rn9*ln3); MVSR = _SS_MALLOC_SCALAR(rn9*rn9);

    indices = new PetscInt[4*vn];
    zeros = new PetscScalar[4*vn]; std::fill(zeros, zeros + 4*vn, 0.);
    for (int i=0; i<4*vn; ++i) indices[i] = 3*vn + i;
    for (int i=0; i<ln3; ++i) {
      VecSetValues(SL[i], 4*vn, indices, zeros, INSERT_VALUES);
      VecAssemblyBegin(SL[i]);
      VecAssemblyEnd(SL[i]);
    }
    for (int i=0; i<rn9; ++i) {
      VecSetValues(SR[i], 4*vn, indices, zeros, INSERT_VALUES);
      VecAssemblyBegin(SR[i]);
      VecAssemblyEnd(SR[i]);
    }

    delete [] indices;
    delete [] zeros;


    for (int i=0; i<rn9; ++i) {
      for (int j=0; j<ln3; ++j) {
	PetscScalar buffer;
	VecDot(VS_R[i], SL[j], &buffer);
	LVSR[i+j*rn9] = buffer;
      }
      for (int j=0; j<rn9; ++j) {
	PetscScalar buffer;
	VecDot(VS_R[i], SR[j], &buffer);
	MVSR[i+j*rn9] = buffer;
      }
    }

    Rot = _SS_MALLOC_SCALAR(rn9); Rot_b = _SS_MALLOC_SCALAR(rn9); //RotNorm = _SS_MALLOC_SCALAR(2*rn);

    std::fill(Rot, Rot+rn9, 0);
    for (int i=0; i<rn; ++i) Rot[i] = Rot[i+4*rn] = Rot[i+8*rn] = 1.;
    init_svd(Rot, rn, NUM_OF_SVD_THREAD);

    GRot[0]=GRot[4]=GRot[8]=1;


    clock_end();
    std::cout << "  linear Proxies: " << ln << "; rotational Proxies: " << rn << std::endl;

  }

  Subspace::~Subspace() {        
    _SS_FREE(linear_proxies);
    
    MatDestroy(&VS);
    for (int i=0; i<ln3; ++i) VecDestroy(&VS_L[i]); delete [] VS_L;
    for (int i=0; i<rn9; ++i) VecDestroy(&VS_R[i]); delete [] VS_R;

    _SS_FREE(SL_V); 
    _SS_FREE(SR_V);
    _SS_FREE(LVS); _SS_FREE(MVS);
    _SS_FREE(LVSR); _SS_FREE(MVSR);

    _SS_FREE(Rot); _SS_FREE(Rot_b); //_SS_FREE(RotNorm);

    _SS_FREE(vertices);
    //delete [] vertices_f;
    PetscFinalize();
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
    //    for (int i=0; i<9; ++i) printf("%.3f ", Rot[10+i*rn]); printf("\n");

  }

  void update_mesh(Mesh *mesh) {
    //update mesh vertices
    // reduced variable to mesh vertices, often computational expensive
    _SS_CBLAS_FUNC(gemv)(CblasColMajor, CblasNoTrans, vn3, ln3, 1, SL_V, vn3, Lin, 1, 0, vertices, 1);
    _SS_CBLAS_FUNC(gemv)(CblasColMajor, CblasNoTrans, vn3, rn9, 1, SR_V, vn3, Rot, 1, 1, vertices, 1);


    for (int i=0, j=0; i<vn; ++i, j+=3)
      apply_rot(mesh->vertices[i], &vertices[j], GRot, 'C');
  }


  void Subspace::update(std::vector<Point> & constraint_points, bool inf){
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

      if (inf) { mesh->normals.clear(); mesh->need_normals();}
    }
  }
  void Subspace::terminate() {
    _SS_FREE(LSYS); _SS_FREE(LSYS_piv);
    _SS_FREE(RHS); 

    _SS_FREE(RHS_hp); _SS_FREE(Lin);
    ready = false;
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



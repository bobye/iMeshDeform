#include "petscvec.h"
#include "petscmat.h"
#include "petscksp.h"

#include "subspace/subspace.hh"

#include <assert.h>
//#include "pthread.h"

#include <time.h>

#include <mkl.h>
#define _SS_SCALAR          double
#define _SS_MALLOC_SCALAR(x)       (_SS_SCALAR *) mkl_malloc( (x) *sizeof(_SS_SCALAR), 64)
#define _SS_MALLOC_INT(x)       (int *) mkl_malloc( (x) *sizeof(int), 64)
#define _SS_FREE(x)         mkl_free(x)

#include "fsvd.hh"

namespace subspace {

  static int vn, vn3;

  static int ln, rn, ln3, rn9;
  //linear proxies
  static PetscScalar*  linear_proxies; // column major
  //rotational proxies
  static std::vector<int> rotational_proxies;

  //variational subspace solver

  static Mat  VS;//sparse matrix to LU
  static Mat  RE;//linearization artifacts regulazier
#define COEFF_REG   (0)

  static Vec* VS_L, *SL;//solved variational subspace
  static Vec* VS_R, *SR;//row major index for each rotation matrix

  static _SS_SCALAR *SL_V, *SR_V; // for mesh reconstruction
  static _SS_SCALAR *LVS, *MVS; //reduced model for linear variational subspace
  static _SS_SCALAR *LVSR, *MVSR; //reduced model for rotational fitting

  static _SS_SCALAR *Lin, *Rot, *Rot_b; //reduced variable, 3x3 rotation matrices are of row major

  static _SS_SCALAR *LSYS, *RHS, *RHS_hp; // dense matrix, rhs and rotation 
  static int *LSYS_piv;
  static int hn, hn3, nsys;

  static _SS_SCALAR *vertices;
  //static float *vertices_f;
  //#define MAX_CONSTRAINT_NUMBER   100
  //  static pthread_t iterate_lin, iterate_rot;


  struct timespec start, end;
#define BILLION  1000000000L
  void clock_start(std::string description) {
    clock_gettime(CLOCK_MONOTONIC, &start);
    std::cout << description <<" ... " << std::flush; 
  }
  void clock_end() {
    clock_gettime(CLOCK_MONOTONIC, &end);
    printf("\t[done] %.3f seconds\n",
	   (( end.tv_sec - start.tv_sec )+ (double)( end.tv_nsec - start.tv_nsec ) / (double)BILLION ));
  }


  Subspace::Subspace(int argc, char **argv) {
    PetscInitialize(&argc,&argv,(char *)0,PETSC_NULL);
  };


  void Subspace::init(trimesh::TriMesh * pm) {
    mesh = pm;

    mesh->need_neighbors();
    mesh->need_adjacentfaces();
    mesh->need_pointareas();
    //mesh->need_curvatures();

    //    std::cout << mesh->normals[10] << mesh->pdir1[10] << mesh->pdir2[10] << std::endl; exit(0);

    vn = mesh->vertices.size(); vn3 = 3*vn;
    rotational_proxies.resize(vn);
    vertices = _SS_MALLOC_SCALAR(vn3);
    //    vertices_f = new float[vn3];
  }



  void Subspace::load_linear_proxies_vg(std::vector<int> &group_ids) {
    assert(vn == group_ids.size());
    ln = *std::max_element(group_ids.begin(), group_ids.end()) + 1; ln3 = 3*ln;
    std::vector<double> count_vertices; count_vertices.resize(ln);
    linear_proxies = new PetscScalar[vn3*ln3]; std::fill(linear_proxies, linear_proxies+vn3*ln3, 0);

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

  inline PetscErrorCode mat_edge_assembly_VS(const PetscInt &v0, const PetscInt &v1, const PetscInt &k, const PetscScalar &weight, const trimesh::vec &v) {
    PetscErrorCode ierr;
    const PetscInt idx[2] = {3*v0, 3*v1}, idy[2] = {3*v0+1, 3*v1+1}, idz[2] = {3*v0+2, 3*v1+2};    
    const PetscInt idq[4] = {vn3+ 4*k, vn3 + 4*k+1, vn3 + 4*k+2, vn3 + 4*k+3}; 

    PetscScalar vvs[4] = {1, -1, -1, 1}; 
    MULTIPLY(vvs, 4, weight)

    ierr = MatSetValues(VS, 2, idx, 2, idx, vvs, ADD_VALUES); CHKERRQ(ierr);
    ierr = MatSetValues(VS, 2, idy, 2, idy, vvs, ADD_VALUES); CHKERRQ(ierr);
    ierr = MatSetValues(VS, 2, idz, 2, idz, vvs, ADD_VALUES); CHKERRQ(ierr);

    PetscScalar vqsx[8] = {v[0], 0, v[2], -v[1], -v[0], 0, -v[2], v[1]};
    PetscScalar vqsy[8] = {v[1], -v[2], 0, v[0], -v[1], v[2], 0, -v[0]};
    PetscScalar vqsz[8] = {v[2], v[1], -v[0], 0, -v[2], -v[1], v[0], 0};
    MULTIPLY(vqsx, 8, -weight)
    MULTIPLY(vqsy, 8, -weight)
    MULTIPLY(vqsz, 8, -weight)

    ierr = MatSetValues(VS, 2, idx, 4, idq, vqsx, ADD_VALUES); CHKERRQ(ierr);
    ierr = MatSetValues(VS, 2, idy, 4, idq, vqsy, ADD_VALUES); CHKERRQ(ierr);
    ierr = MatSetValues(VS, 2, idz, 4, idq, vqsz, ADD_VALUES); CHKERRQ(ierr);

    PetscScalar vqs[16] = {v[0]*v[0] + v[1]*v[1] + v[2]*v[2], 0, 0, 0,
			   0, v[2]*v[2]+v[1]*v[1], -v[0]*v[1], -v[0]*v[2],
			   0, -v[1]*v[0], v[0]*v[0]+v[2]*v[2], -v[1]*v[2],
			   0, -v[2]*v[0], -v[2]*v[0], v[1]*v[1]+v[0]*v[0]};
    MULTIPLY(vqs, 16, weight)
    ierr = MatSetValues(VS, 4, idq, 4, idq, vqs, ADD_VALUES); CHKERRQ(ierr);

    PetscScalar vs[4];
    for (int i=0; i<3; ++i) {
      vs[0] = v[i]; vs[1] = -v[i];     
      MULTIPLY(vs, 2, weight)
      VecSetValues(VS_R[9*rotational_proxies[k]+i], 2, idx, vs, ADD_VALUES);
    }
    for (int i=0; i<3; ++i) {
      vs[0] = v[i]; vs[1] = -v[i];
      MULTIPLY(vs, 2, weight)
      VecSetValues(VS_R[9*rotational_proxies[k]+i+3], 2, idy, vs, ADD_VALUES);
    }
    for (int i=0; i<3; ++i) {
      vs[0] = v[i]; vs[1] = -v[i];
      MULTIPLY(vs, 2, weight)
      VecSetValues(VS_R[9*rotational_proxies[k]+i+6], 2, idz, vs, ADD_VALUES);
    }

    for (int i=0; i<3; ++i) {
      vs[0] = -v[i]*v[0]; vs[1] = 0; vs[2] = -v[i]*v[2]; vs[3]=v[i]*v[1];
      MULTIPLY(vs, 4, weight)
      VecSetValues(VS_R[9*rotational_proxies[k]+i], 4, idq, vs, ADD_VALUES);
    }
    for (int i=0; i<3; ++i) {
      vs[0] = -v[i]*v[1]; vs[1] = v[i]*v[2]; vs[2] = 0; vs[3]=-v[i]*v[0];
      MULTIPLY(vs, 4, weight)
      VecSetValues(VS_R[9*rotational_proxies[k]+i+3], 4, idq, vs, ADD_VALUES);
    }
    for (int i=0; i<3; ++i) {
      vs[0] = -v[i]*v[2]; vs[1] = -v[i]*v[1]; vs[2] = v[i]*v[0]; vs[3]=0;
      MULTIPLY(vs, 4, weight)
      VecSetValues(VS_R[9*rotational_proxies[k]+i+6], 4, idq, vs, ADD_VALUES);
    }
    return ierr;   
  }

  void Subspace::assembly() {
    int N = 7*vn + 3*ln; 
    int *nnz = new int[N];
    for (int i=0; i<vn3; ++i) nnz[i] = 7*mesh->neighbors[i/3].size() + 3*ln + 5;
    for (int i=0; i<4*vn; ++i) nnz[vn3 + i] = 4;
    for (int i=0; i<ln3; ++i) nnz[7*vn + i] = 1;

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
      for (int j= 0; j<fn; ++j) {
	int v0 = mesh->faces[faces[j]][0], v1 = mesh->faces[faces[j]][1], v2 = mesh->faces[faces[j]][2];
	trimesh::vec v01 = mesh->vertices[v0] - mesh->vertices[v1];
	trimesh::vec v12 = mesh->vertices[v1] - mesh->vertices[v2];
	trimesh::vec v20 = mesh->vertices[v2] - mesh->vertices[v0];
	/*
	trimesh::vec u01(v01 DOT normal, v01 DOT pdir1, v01 DOT pdir2);
	trimesh::vec u12(v12 DOT normal, v12 DOT pdir1, v12 DOT pdir2);
	trimesh::vec u20(v20 DOT normal, v20 DOT pdir1, v20 DOT pdir2);
	*/
	mat_edge_assembly_VS(v0, v1, i, std::fabs(1./std::tan(mesh->cornerangle(2,j))), v01);
	mat_edge_assembly_VS(v1, v2, i, std::fabs(1./std::tan(mesh->cornerangle(0,j))), v12);
	mat_edge_assembly_VS(v2, v0, i, std::fabs(1./std::tan(mesh->cornerangle(1,j))), v20);

      }

      PetscScalar parea = mesh->pointareas[i];
      for (PetscInt idq = vn3+4*i; idq < vn3+4*i+4; ++idq)
	MatSetValues(RE, 1, &idq, 1, &idq, &parea, INSERT_VALUES);

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
    /**************************************************/
    // solve sparse system
    Mat L; MatConvert(VS, MATSEQAIJ, MAT_INITIAL_MATRIX, &L);
    KSP ksp;
    SL = new Vec[ln3]; SR = new Vec[rn9];
    VecDuplicateVecs(VS_L[0], ln3, &SL); VecDuplicateVecs(VS_R[0], rn9, &SR);
    
    KSPCreate(PETSC_COMM_SELF, &ksp);
    KSPSetOperators(ksp, L, L, SAME_PRECONDITIONER);
    //KSPSetType(ksp, KSPPREONLY);
    KSPSetFromOptions(ksp);


    for (int i=0; i<ln3; ++i) KSPSolve(ksp, VS_L[i], SL[i]); 
    for (int i=0; i<rn9; ++i) KSPSolve(ksp, VS_R[i], SR[i]); 

    MatDestroy(&L);
    KSPDestroy(&ksp);
    /**************************************************/
    // copy subspace solution data to global array

    SL_V = _SS_MALLOC_SCALAR(vn3*ln3); SR_V = _SS_MALLOC_SCALAR(vn3*rn9);
    PetscInt *indices = new PetscInt[vn3];
    for (int i=0; i<vn3; ++i) indices[i] = i;
    for (int i=0; i<ln3; ++i)
      VecGetValues(SL[i], vn3, indices, &SL_V[vn3*i]);
    for (int i=0; i<rn9; ++i)
      VecGetValues(SR[i], vn3, indices, &SR_V[vn3*i]);
    delete [] indices;

    /**************************************************/
    // precompute online dense linear system
    MatAXPY(VS, COEFF_REG, RE, DIFFERENT_NONZERO_PATTERN);

    LVS = _SS_MALLOC_SCALAR (ln3*ln3); MVS = _SS_MALLOC_SCALAR(ln3*rn9);

    indices = new PetscInt[ln3];
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
	VecDot(SL[j], tmp, &LVS[ln3*j+i]);
	LVS[ln3*i+j] = LVS[ln3*j+i];
      }
      for (int j=0; j<rn9; ++j) {
	VecDot(SL[i], VS_R[j], &MVS[ln3*j+i]);
	PetscScalar tominus;
	VecDot(SR[j], tmp, &tominus);
	MVS[ln3*j+i] -= tominus;
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
      for (int j=0; j<ln3; ++j)
	VecDot(VS_R[i], SL[j], &LVSR[i+j*rn9]);
      for (int j=0; j<rn9; ++j)
	VecDot(VS_R[i], SR[j], &MVSR[i+j*rn9]);
    }

    Rot = _SS_MALLOC_SCALAR(rn9); Rot_b = _SS_MALLOC_SCALAR(rn9);
    std::fill(Rot, Rot+rn9, 0);
    for (int i=0; i<rn9; i+=9) Rot[i] = Rot[i+4] = Rot[i+8] = 1.;


    clock_end();
    std::cout << "  linear Proxies: " << ln << "; rotational Proxies: " << rn << std::endl;

  }

  Subspace::~Subspace() {        
    delete [] linear_proxies;
    
    MatDestroy(&VS);
    for (int i=0; i<ln3; ++i) VecDestroy(&VS_L[i]); delete [] VS_L;
    for (int i=0; i<rn9; ++i) VecDestroy(&VS_R[i]); delete [] VS_R;

    _SS_FREE(SL_V); 
    _SS_FREE(SR_V);
    _SS_FREE(LVS); _SS_FREE(MVS);
    _SS_FREE(LVSR); _SS_FREE(MVSR);

    _SS_FREE(Rot); _SS_FREE(Rot_b);

    _SS_FREE(vertices);
    //delete [] vertices_f;
    PetscFinalize();
  }


  void Subspace::prepare(std::vector< std::vector<float> > & constraints, std::vector<trimesh::point> & constraint_points) { // precompute LU for dense direct solver, initialize Rot and Lin
    hn = constraint_points.size(); if (hn<3) {std::cout << "Need more constraints!\n" << std::endl; return;} hn3 =3*hn;

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
    /*
    for (int i=0; i<ln3; ++i)
      for (int j=0; j<ln3; ++j)
	LSYS[j + i*nsys] = LVS[j + i*ln3];
    */    
    dlacpy("F", &ln3, &ln3, LVS, &ln3, LSYS, &nsys); // copy left-top ln3 x ln3 block from LVS
    /*
    for (int k=0; k<ln3; ++k)
      for (int i=0; i<hn3; ++i) 
	for (int j=0; j<vn3; ++j) 
	  LSYS[k + (ln3 + i)*nsys] += constraints[i/3][j/3] * SL_V[j + vn3*(k)];
    */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, hn3, ln3, vn3, 1, constraints_matrix, hn3, SL_V, vn3, 0, LSYS+ln3, nsys); // compute left-bottom hn3 x ln3 block from online constraints

    /*
    for (int i=ln3; i<nsys; ++i)
      for (int j=0; j<ln3; ++j)
	LSYS[j+i*nsys] = LSYS[i+j*nsys];
    */
    cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, ln3, hn3, hn3, 1, LSYS+ln3, nsys, one, hn3, 0, LSYS+ln3*nsys, nsys); // copy left-bottom block to right-top block by transpose

    /*
    for (int i=0; i<rn9; ++i)
      for (int j=0; j<ln3; ++j)
	RHS[j + i*nsys] = MVS[j + i*ln3];
    */
    dlacpy("F", &ln3, &rn9, MVS, &ln3, RHS, &nsys); // copy to first ln3 rows of RHS from MVS

    /*
    for (int k=0; k<rn9; ++k)
      for (int i=0; i<hn3; ++i) 
	for (int j=0; j<vn3; ++j)
	  RHS[ln3+i + k*nsys] -= constraints[i/3][j/3] * SR_V[j + vn3*k];
    */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, hn3, rn9, vn3, -1, constraints_matrix, hn3, SR_V, vn3, 0, RHS+ln3, nsys); // compute last hn3 rows of RHS

    LAPACKE_dgetrf(LAPACK_COL_MAJOR, nsys, nsys, LSYS, nsys, LSYS_piv);// compute LU of LSYS

    ready = true;

    RHS_hp = _SS_MALLOC_SCALAR(hn3);
    Lin  = _SS_MALLOC_SCALAR(nsys);

    _SS_FREE(constraints_matrix); _SS_FREE(one);
    
    update(constraint_points, true);
    clock_end();
  }



  void reduced_linsolve() {//solve reduced linear variables via dense direct solve
    cblas_dgemv(CblasColMajor, CblasNoTrans, nsys, rn9, 1, RHS, nsys, Rot, 1, 0, Lin, 1);
    cblas_daxpy(hn3, 1, RHS_hp, 1, Lin+ln3, 1);
    LAPACKE_dgetrs(LAPACK_COL_MAJOR, 'N', nsys, 1, LSYS, nsys, LSYS_piv, Lin, nsys);
  }

  void reduced_rotsolve() {//solve reduced rotational variables via SVD of gradient
    cblas_dcopy(rn9, Rot, 1, Rot_b, 1);
    cblas_dgemv(CblasColMajor, CblasNoTrans, rn9, rn9, 1, MVSR, rn9, Rot_b, 1, 0, Rot, 1);
    cblas_dgemv(CblasColMajor, CblasNoTrans, rn9, ln3, 1, LVSR, rn9, Lin, 1, 1, Rot, 1);   

    //cblas_dgemv(CblasColMajor, CblasTrans, ln3, rn9, 1, MVS, ln3, Lin, 1, 0, Rot, 1);       
    dfastsvd(Rot, rn);
    /*
    for (int j=0; j<rn; ++j) 
      {for (int i=0; i<9; ++i) printf("%.3f ", Rot[9*j+i]); printf("\n");}
    */
  }

  void update_mesh(trimesh::TriMesh *mesh) {
    //update mesh vertices
    cblas_dgemv(CblasColMajor, CblasNoTrans, vn3, ln3, 1, SL_V, vn3, Lin, 1, 0, vertices, 1);
    cblas_dgemv(CblasColMajor, CblasNoTrans, vn3, rn9, 1, SR_V, vn3, Rot, 1, 1, vertices, 1);    
    /*
    for (int i=0; i<vn3; ++i) vertices_f[i] = (float) vertices[i];
    cblas_scopy(vn, vertices_f, 3, &mesh->vertices[0][0], sizeof(mesh->vertices[0]));
    cblas_scopy(vn, vertices_f+1, 3, &mesh->vertices[0][1], sizeof(mesh->vertices[0]));
    cblas_scopy(vn, vertices_f+2, 3, &mesh->vertices[0][2], sizeof(mesh->vertices[0]));
    */

    for (int i=0, j=0; i<vn; ++i, j+=3) {
      mesh->vertices[i][0] = vertices[j];
      mesh->vertices[i][1] = vertices[j+1];
      mesh->vertices[i][2] = vertices[j+2];
    }
  }


#define NUM_OF_ITERATION 10
  void Subspace::update(std::vector<trimesh::point> & constraint_points, bool inf){
    if (ready) {
      //      clock_start("Run reduced model");
      for (int i=0, j=0; j<hn; i+=3, ++j) {
	RHS_hp[i] = constraint_points[j][0];
	RHS_hp[i+1] = constraint_points[j][1];
	RHS_hp[i+2] = constraint_points[j][2];
      }
      
      int N = inf? 100: NUM_OF_ITERATION;
      for (int i=0; i<N; ++i) {
	reduced_linsolve();
	reduced_rotsolve();
      }
      update_mesh(mesh);
      //recompute normals

      if (inf) {mesh->normals.clear(); mesh->need_normals();}

      //      clock_end();
    }
  }
  void Subspace::terminate() {
    _SS_FREE(LSYS); _SS_FREE(LSYS_piv);
    _SS_FREE(RHS); 

    _SS_FREE(RHS_hp); _SS_FREE(Lin);
    ready = false;
  }
}



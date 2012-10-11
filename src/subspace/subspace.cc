#include "petscvec.h"
#include "petscmat.h"
#include "petscksp.h"

#include "subspace/subspace.hh"

#include "assert.h"
namespace subspace {

  static int vn;

  static int ln, rn;
  //linear proxies
  static PetscScalar*  linear_proxies; // row major
  //rotational proxies
  static std::vector<int> rotational_proxies;

  //variational subspace solver

  static Mat  VS;//sparse matrix to LU

  static Vec* VS_L, *SL;//solved variational subspace
  static Vec* VS_R, *SR;//row major index for each rotation matrix

  static double *SL_V, *SR_V; // for mesh reconstruction
  static double *LVS, *MVS; //reduced model

  clock_t start, end;
  void clock_start(std::string description) {
    start = clock();
    std::cout << description <<" ... " << std::flush; 
  }

  void clock_end() {
    end = clock();
    std::cout << "\t[done] " << (static_cast<double> (end) - static_cast<double> (start)) / CLOCKS_PER_SEC <<" seconds" << std::endl;
  }


  Subspace::Subspace(int argc, char **argv) {
    PetscInitialize(&argc,&argv,(char *)0,PETSC_NULL);
  };


  void Subspace::init(trimesh::TriMesh * pm) {
    mesh = pm;

    mesh->need_neighbors();
    mesh->need_adjacentfaces();

    vn = mesh->vertices.size();
    rotational_proxies.resize(vn);
  }



  void Subspace::load_linear_proxies_vg(std::vector<int> &group_ids) {
    assert(vn == group_ids.size());
    ln = *std::max_element(group_ids.begin(), group_ids.end()) + 1;
    std::vector<double> count_vertices; count_vertices.resize(ln);
    linear_proxies = new PetscScalar[9*vn*ln]; std::fill(linear_proxies, linear_proxies+9*vn*ln, 0);

    for (int i=0; i<vn; ++i) {
      linear_proxies[3*group_ids[i] + 3*i*3*ln] = 1.; // row major
      linear_proxies[3*group_ids[i] +1 + (3*i+1)*3*ln] = 1.; 
      linear_proxies[3*group_ids[i] +2 + (3*i+2)*3*ln] = 1.; 
      ++count_vertices[group_ids[i]];
    }
    for (int i=0; i<vn; ++i) {
      linear_proxies[3*group_ids[i] + 3*i*3*ln] /= count_vertices[group_ids[i]]; // normalize
      linear_proxies[3*group_ids[i] +1 + (3*i+1)*3*ln] /= count_vertices[group_ids[i]]; 
      linear_proxies[3*group_ids[i] +2 + (3*i+2)*3*ln] /= count_vertices[group_ids[i]]; 
    }
  }

  void Subspace::load_rotational_proxies(std::vector<int> &group_ids) {
    assert(vn == group_ids.size());
    
    for (int i=0; i<vn; ++i)
      if (rn <= (rotational_proxies[i] = group_ids[i])) rn = group_ids[i]+1;
    
  }

#define MULTIPLY(v,n,w) for (int i=0; i<n; ++i) v[i] *= w;

  inline PetscErrorCode mat_edge_assembly_VS(const PetscInt &v0, const PetscInt &v1, const PetscInt &k, const PetscScalar &weight, const trimesh::vec &v) {
    PetscErrorCode ierr;
    const PetscInt idx[2] = {3*v0, 3*v1}, idy[2] = {3*v0+1, 3*v1+1}, idz[2] = {3*v0+2, 3*v1+2};    
    const PetscInt idq[4] = {3*vn+ 4*k, 3*vn + 4*k+1, 3*vn + 4*k+2, 3*vn + 4*k+3}; 

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
      VecSetValues(VS_R[9*rotational_proxies[k]+i], 2, idx, vs, ADD_VALUES);
    }
    for (int i=0; i<3; ++i) {
      vs[0] = v[i]; vs[1] = -v[i];
      VecSetValues(VS_R[9*rotational_proxies[k]+i+3], 2, idy, vs, ADD_VALUES);
    }
    for (int i=0; i<3; ++i) {
      vs[0] = v[i]; vs[1] = -v[i];
      VecSetValues(VS_R[9*rotational_proxies[k]+i+6], 2, idz, vs, ADD_VALUES);
    }

    for (int i=0; i<3; ++i) {
      vs[0] = -v[i]*v[0]; vs[1] = 0; vs[2] = -v[i]*v[2]; vs[3]=v[i]*v[1];
      VecSetValues(VS_R[9*rotational_proxies[k]+i], 4, idq, vs, ADD_VALUES);
    }
    for (int i=0; i<3; ++i) {
      vs[0] = -v[i]*v[1]; vs[1] = v[i]*v[2]; vs[2] = 0; vs[3]=-v[i]*v[0];
      VecSetValues(VS_R[9*rotational_proxies[k]+i+3], 4, idq, vs, ADD_VALUES);
    }
    for (int i=0; i<3; ++i) {
      vs[0] = -v[i]*v[2]; vs[1] = -v[i]*v[1]; vs[2] = v[i]*v[0]; vs[3]=0;
      VecSetValues(VS_R[9*rotational_proxies[k]+i+6], 4, idq, vs, ADD_VALUES);
    }
    return ierr;   
  }

  void Subspace::assembly() {
    int N = 7*vn + 3*ln; 
    int *nnz = new int[N];
    for (int i=0; i<3*vn; ++i) nnz[i] = 7*mesh->neighbors[i/3].size() + 3*ln + 5;
    for (int i=0; i<4*vn; ++i) nnz[3*vn + i] = 4;
    for (int i=0; i<3*ln; ++i) nnz[7*vn + i] = 1;

    MatCreateSeqSBAIJ(PETSC_COMM_SELF, 1, N, N, 0, nnz, &VS); delete [] nnz;
    MatSetOption(VS, MAT_IGNORE_LOWER_TRIANGULAR, PETSC_TRUE);

    // create LHS of subspace problem
    VS_L = new Vec[3*ln];
    VS_R = new Vec[9*rn];
    for (int i=0; i<3*ln; ++i) MatGetVecs(VS, &VS_L[i], PETSC_NULL);
    for (int i=0; i<9*rn; ++i) MatGetVecs(VS, &VS_R[i], PETSC_NULL);

    for (int i=0; i< 3*ln; ++i) { 
      PetscInt index = i+7*vn; PetscScalar one = 1.; 
      VecSet(VS_L[i], 0.);
      VecSetValues(VS_L[i], 1, &index, &one, INSERT_VALUES);
      VecAssemblyBegin(VS_L[i]);
      VecAssemblyEnd(VS_L[i]);
    }
    for (int i=0; i< 9*rn; ++i) {
      VecSet(VS_R[i], 0.);
    }

    // assembly VS

    for (int i = 0; i<vn; ++i) {
      std::vector<int> &faces = mesh->adjacentfaces[i];
      int fn = faces.size();
      for (int j= 0; j<fn; ++j) {
	int v0 = mesh->faces[faces[j]][0], v1 = mesh->faces[faces[j]][1], v2 = mesh->faces[faces[j]][2];
	trimesh::vec v01 = mesh->vertices[v0] - mesh->vertices[v1];
	trimesh::vec v12 = mesh->vertices[v1] - mesh->vertices[v2];
	trimesh::vec v20 = mesh->vertices[v2] - mesh->vertices[v0];
	
	mat_edge_assembly_VS(v0, v1, i, std::fabs(1./std::tan(mesh->cornerangle(2,j))), v01);
	mat_edge_assembly_VS(v1, v2, i, std::fabs(1./std::tan(mesh->cornerangle(0,j))), v12);
	mat_edge_assembly_VS(v2, v0, i, std::fabs(1./std::tan(mesh->cornerangle(1,j))), v20);
      }
    }

    // assembly VH
    int *irow = new int[3*vn], *icol = new int[3*ln];
    for (int i=0; i<3*vn; ++i) irow[i] = i;
    for (int i=0; i<3*ln; ++i) icol[i] = 7*vn + i;

    MatSetValues(VS, 3*vn, irow, 3*ln, icol, linear_proxies, ADD_VALUES);
    delete [] irow; delete [] icol;

    PetscScalar zero = 0;
    for (int i=7*vn; i<N; ++i) MatSetValues(VS, 1, &i, 1, &i, &zero, ADD_VALUES);



    MatAssemblyBegin(VS,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(VS,MAT_FINAL_ASSEMBLY);


    for (int i=0; i< 9*rn; ++i) {
      VecAssemblyBegin(VS_R[i]);
      VecAssemblyEnd(VS_R[i]);
    }
  }



  void Subspace::solve() {
    clock_start("Solving reduced model");
    assembly();

    Mat L; MatConvert(VS, MATSEQAIJ, MAT_INITIAL_MATRIX, &L);
    KSP ksp;
    SL = new Vec[3*ln]; SR = new Vec[9*rn];
    VecDuplicateVecs(VS_L[0], 3*ln, &SL); VecDuplicateVecs(VS_R[0], 9*rn, &SR);
    
    KSPCreate(PETSC_COMM_SELF, &ksp);
    KSPSetOperators(ksp, L, L, SAME_PRECONDITIONER);
    //KSPSetType(ksp, KSPPREONLY);
    KSPSetFromOptions(ksp);


    for (int i=0; i<3*ln; ++i) KSPSolve(ksp, VS_L[i], SL[i]); 
    for (int i=0; i<9*rn; ++i) KSPSolve(ksp, VS_R[i], SR[i]); 

    MatDestroy(&L);
    KSPDestroy(&ksp);

    SL_V = new double [3*vn*3*ln]; SR_V = new double [3*vn*9*rn];
    PetscInt *indices = new PetscInt[3*vn];
    for (int i=0; i<3*vn; ++i) indices[i] = i;
    for (int i=0; i<3*ln; ++i)
      VecGetValues(SL[i], 3*vn, indices, &SL_V[3*vn*i]);
    for (int i=0; i<9*rn; ++i)
      VecGetValues(SR[i], 3*vn, indices, &SR_V[3*vn*i]);
    delete [] indices;

    LVS = new double [3*ln*3*ln]; MVS = new double [3*ln*9*rn];

    indices = new PetscInt[3*ln];
    PetscScalar *zeros = new PetscScalar[3*ln]; std::fill(zeros, zeros + 3*ln, 0.);
    for (int i=0; i<3*ln; ++i) indices[i] = 7*vn + i;
    for (int i=0; i<3*ln; ++i) {
      VecSetValues(SL[i], 3*ln, indices, zeros, INSERT_VALUES);
      VecAssemblyBegin(SL[i]);
      VecAssemblyEnd(SL[i]);
    }
    for (int i=0; i<9*rn; ++i) {
      VecSetValues(SR[i], 3*ln, indices, zeros, INSERT_VALUES);
      VecAssemblyBegin(SR[i]);
      VecAssemblyEnd(SR[i]);
    }

    for (int i=0; i<3*ln; ++i) {
      Vec tmp; VecDuplicate(SL[i], &tmp);
      MatMult(VS, SL[i], tmp);
      for (int j=i; j<3*ln; ++j) {
	VecDot(SL[j], tmp, &LVS[3*ln*j+i]);
	LVS[3*ln*i+j] = LVS[3*ln*j+i];
      }
      for (int j=0; j<9*rn; ++j) {
	VecDot(SR[j], tmp, &MVS[3*ln*j+i]);
	PetscScalar toadd;
	VecDot(SL[i], VS_R[j], &toadd);
	MVS[3*ln*j+i] += toadd;
      }
      VecDestroy(&tmp);
    }

    clock_end();
  }

  Subspace::~Subspace() {        
    delete [] linear_proxies;
    
    MatDestroy(&VS);
    for (int i=0; i<3*ln; ++i) VecDestroy(&VS_L[i]); delete [] VS_L;
    for (int i=0; i<9*rn; ++i) VecDestroy(&VS_R[i]); delete [] VS_R;

    delete [] SL_V; 
    delete [] SR_V;
    delete [] LVS;
    delete [] MVS;

    PetscFinalize();
  }

  
}

#include "subspace/LB.hh"
#include <petscmat.h>
#include "slepceps.h"


namespace subspace {

  Mat mass_mat, stiff_mat;
  PetscInt mat_size;

  EPS         	 eps;		  /* eigenproblem solver context */
  ST          	 st;		  /* spectral transformation context */
  Vec            x;

  std::vector<double> eigvalues;
  std::vector< PetscScalar *> eigvectors;
  std::vector<Vec> eigVecs;

  LB::LB(int argc, char** argv) {
    SlepcInitialize(&argc,&argv,(char *)0,PETSC_NULL);
  }

  LB::~LB(){
    MatDestroy(&mass_mat); MatDestroy(&stiff_mat);    

    EPSDestroy(&eps); VecDestroy(&x);
    for (int i=0; i<eigVecs.size(); ++i) VecDestroy(&eigVecs[i]);

    SlepcFinalize();
  }

  void LB::init(Mesh* pm) {
    mesh = pm;
    mesh->need_neighbors();
    mesh->need_faceareas();
    mesh->need_edgelengths();
  }

const PetscInt J1[9]
= {2, -1, -1, 
   -1, 0, 1, 
   -1, 1, 0};
const PetscInt J2[9]
= {0, -1, 1, 
   -1, 2, -1, 
   1, -1, 0};
const PetscInt J3[9]
= {0, 1, -1,
   1, 0, -1, 
   -1, -1, 2};
const PetscInt J4[9]
= {2, 1, 1, 1, 2, 1, 1, 1, 2};

  void LB::compute_operator(){
    mat_size = mesh->vertices.size();
    PetscInt *nnz = new PetscInt[mat_size];
    for (PetscInt i=0;i<mat_size;++i)
      nnz[i]=mesh->neighbors[i].size() +1;

    MatCreateSeqAIJ(PETSC_COMM_SELF, mat_size, mat_size, 0, nnz, &mass_mat);
    MatCreateSeqAIJ(PETSC_COMM_SELF, mat_size, mat_size, 0, nnz, &stiff_mat);
    delete [] nnz;
    //MatSetOption(mass_mat, MAT_IGNORE_LOWER_TRIANGULAR, PETSC_TRUE);
    //MatSetOption(stiff_mat, MAT_IGNORE_LOWER_TRIANGULAR, PETSC_TRUE);

    for (PetscInt i=0; i<mesh->faces.size(); ++i) {
      float &l1 = mesh->edgelengths[i][0],
	&l2 = mesh->edgelengths[i][1],
	&l3 = mesh->edgelengths[i][2];

      int *idx = mesh->faces[i];
      PetscScalar mass[9], stiff[9];

      for (PetscInt j = 0; j < 9; ++j) {
	mass[j]= mesh->faceareas[i]*J4[j]/12.;
	stiff[j]= (l1*l1*J1[j]+l2*l2*J2[j]+l3*l3*J3[j])/(8.* mesh->faceareas[i]);
      }
      MatSetValues(stiff_mat, 3, idx, 3, idx, stiff, ADD_VALUES);
      MatSetValues(mass_mat, 3, idx, 3, idx, mass, ADD_VALUES);	       

    }

    MatAssemblyBegin(stiff_mat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(stiff_mat, MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(mass_mat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(mass_mat, MAT_FINAL_ASSEMBLY);
    
  }

  void LB::solve_eigen(){
    /* 
       Create eigensolver context
    */    
    EPSCreate(PETSC_COMM_WORLD,&eps);
    /* 
       Set operators. In this case, it is a generalized eigenvalue problem
    */
    EPSSetOperators(eps,stiff_mat,mass_mat);
    EPSSetProblemType(eps, EPS_GHEP); 

    EPSGetST(eps,&st);
    EPSSetTarget(eps,-0.005);
    EPSSetWhichEigenpairs(eps,EPS_TARGET_MAGNITUDE);
    STSetType(st,STSINVERT);

    MatGetVecs(stiff_mat,&x,PETSC_NULL);
    VecSet(x,1.0);
    EPSSetDeflationSpace(eps,1,&x);

    EPSSetFromOptions(eps);

    EPSSolve(eps);

    int nconv;
    PetscScalar 	 kr, ki, re, im, error;

    EPSGetConverged(eps,&nconv);
    eigvalues.resize(nconv);
    eigvectors.resize(nconv);
    eigVecs.resize(nconv);
    
    for(int i=0; i<nconv; i++ ) {
      VecDuplicate(x, &eigVecs[i]);
      EPSGetEigenpair(eps,i,&kr,&ki,eigVecs[i],PETSC_NULL);
      EPSComputeRelativeError(eps,i,&error);
#if defined(PETSC_USE_COMPLEX)
      re = PetscRealPart(kr);
      im = PetscImaginaryPart(kr);
#else
      re = kr;
      im = ki;
#endif

      if ( im != 0.0 ) {
	PetscPrintf(PETSC_COMM_WORLD," % 6g %+6g i",re,im);
      } else {
	PetscPrintf(PETSC_COMM_WORLD,"  % 12.12e      ",re);
      }
      PetscPrintf(PETSC_COMM_WORLD," % 12g\n",error);

      eigvalues[i] = re;
      VecGetArray(eigVecs[i], &eigvectors[i]);

    }

    std::ofstream fid("default_examples.embd");
    fid << mat_size << " " << nconv << std::endl;
    for (int i = 0; i < mat_size; ++i) {
      for (int j = 0; j < nconv; ++j) 
	fid << eigvectors[j][i]/std::sqrt(eigvalues[j]) << " ";
      fid << std::endl;
    }
    fid.close();


  }
}






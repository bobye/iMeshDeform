//#include "petscvec.h"
//#include "petscmat.h"

#include "subspace/subspace.hh"

namespace subspace {

  std::vector<bool*> linear_constraint_handler;  
  bool* is_vertex_rigid;
  std::vector<bool*> rigid_transformer;

  //proxies
  //Mat  linear_proxy_handler;

  //variational subspace solver
  /*  
  Mat  VS_mat;//sparse matrix for LU
  Vec* VS_N;//solved variational subspace
  Vec* VS_U;
  */

  Subspace::Subspace(int argc, char **argv) {
    //PetscInitialize(&argc,&argv,(char *)0,PETSC_NULL);
  };


  void Subspace::init(TriMesh * pm) {
    mesh = pm;
    is_vertex_rigid = new bool[mesh->vertices.size()];
  }

  void Subspace::add_rigid_transformer(bool * selected) {
    int vn = mesh->vertices.size(), new_rigid_count=0, rt_count = rigid_transformer.size();
    rigid_transformer.push_back(new bool[vn]); 
    for (int i=0; i<vn; ++i) 
      if (selected[i] && !is_vertex_rigid[i]) 
	new_rigid_count += (is_vertex_rigid[i] = rigid_transformer[rt_count][i] = true);
    std::cout << "Add rigid transformer (#vert " << new_rigid_count << ")" << std::endl;
  }

  void Subspace::add_linear_constraint_handler(bool * selected) {
    int vn = mesh->vertices.size(), vertex_count=0, rt_count = linear_constraint_handler.size();
    linear_constraint_handler.push_back(new bool[vn]); 
    for (int i=0; i<vn; ++i)       
      vertex_count += (linear_constraint_handler[rt_count][i] = selected[i]);
    std::cout << "Add linear constraint handler (#vert " << vertex_count << ")" << std::endl;
  }

  Subspace::~Subspace() {    
    delete [] is_vertex_rigid;
    for (unsigned int i=0; i<rigid_transformer.size(); ++i)
      delete [] rigid_transformer[i];
    for (unsigned int i=0; i<linear_constraint_handler.size(); ++i)
      delete [] linear_constraint_handler[i];
  }
}

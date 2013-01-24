#include "TriMesh.h"
#include "subspace/gui.hh"

using namespace trimesh;
using namespace subspace;


#define UI_DEBUG

int main(int argc, char *argv[])
{
  // read mesh  
  TriMesh *mesh = TriMesh::read(argv[1]);
  if (!mesh) exit(1);
  int vn = mesh->vertices.size();
  /*
  float *vertices = new float[3*vn];
  for (int i=0; i<vn; ++i) {
    vertices[3*i] = mesh->vertices[i][0];
    vertices[3*i+1] = mesh->vertices[i][1];
    vertices[3*i+2] = mesh->vertices[i][2];
  }
  */
  // load linear proxies
  /*
  std::fstream fid(argv[2]); std::vector<int> group_ids1; group_ids1.resize(vn);
  if (!fid) exit(1); 
  for (int i=0; i<vn; ++i) fid >> group_ids1[i];
  fid.close();  

  // load rotational proxies
  fid.open(argv[3]); std::vector<int> group_ids2; group_ids2.resize(vn);
  if (!fid) exit(1); 
  for (int i=0; i<vn; ++i) fid >> group_ids2[i];
  fid.close();  

  */

  /*
  Subspace ss_solver = Subspace(argc,argv);
  ss_solver.init(mesh);
  ss_solver.load_linear_proxies_vg(group_ids1);
  ss_solver.load_rotational_proxies(group_ids2);
  */

#ifndef UI_DEBUG
  //  ss_solver.solve();
#else
  //  ss_solver.read("scene.ss");
#endif
  
  Scene scene(argc, argv);


  Object object = Object(mesh);
  object.name = std::string(argv[1]);
  //  object.register_mesh(vertices);
  object.register_mesh();
  scene.bind(&object);
  //  scene.bind(&ss_solver); // bind mesh to subspace solver


  scene.view();
  return 0;
}

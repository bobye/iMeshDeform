#include "TriMesh.h"
#include "subspace/gui.hh"

using namespace trimesh;
using namespace subspace;
int main(int argc, char *argv[])
{
  TriMesh *mesh = TriMesh::read(argv[1]);
  if (!mesh) exit(1);
  int vn = mesh->vertices.size();
  std::fstream fid(argv[2]); std::vector<int> group_ids; group_ids.resize(vn);
  if (!fid) exit(1); 
  for (int i=0; i<vn; ++i) fid >> group_ids[i];
  fid.close();  

  Subspace ss_solver = Subspace(argc,argv);
  ss_solver.init(mesh);
  ss_solver.load_linear_proxies_vg(group_ids);
  ss_solver.load_rotational_proxies(group_ids);
  ss_solver.solve();

  
  Scene scene(argc, argv);

  Object object = Object(mesh);
  object.register_mesh();
  scene.bind(&object);

  scene.bind(&ss_solver); // bind mesh to subspace solver

  scene.view();
  return 0;
}

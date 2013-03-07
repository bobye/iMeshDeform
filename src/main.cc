#include "subspace/geometry.hh"
#include "subspace/gui.hh"

//using namespace trimesh;
using namespace subspace;


//#define UI_DEBUG
//#define OFF_THE_FLY

int main(int argc, char *argv[])
{
  // read mesh
  Mesh *mesh = (Mesh *) Mesh::read(argv[1]);
  if (!mesh) exit(1);
  int vn = mesh->vertices.size();

  // load linear proxies
  std::fstream fid(argv[2]); std::vector<int> group_ids1(vn);
  if (!fid) exit(1); 
  for (int i=0; i<vn; ++i) fid >> group_ids1[i];
  fid.close();  

  // load rotational proxies
  fid.open(argv[3]); std::vector<int> group_ids2(vn);
  if (!fid) exit(1); 
  for (int i=0; i<vn; ++i) fid >> group_ids2[i];
  fid.close();  

  //std::vector<int> group_ids3(vn, 0);

  fid.open(argv[4]); std::vector<int> group_ids3(vn);
  if (!fid) exit(1); 
  for (int i=0; i<vn; ++i) fid >> group_ids3[i];
  fid.close();  


  Subspace ss_solver(argc,argv);
  ss_solver.init(mesh);
  ss_solver.load_linear_proxies_vg(group_ids1);
  ss_solver.load_rotational_proxies(group_ids2);
  ss_solver.load_controls(group_ids3);

#ifndef UI_DEBUG
  ss_solver.solve();
#else
  ss_solver.read("scene.ss");
#endif

#ifdef OFF_THE_FLY  
  ss_solver.set_off_fly();//deform mesh off the fly
#endif

  Scene scene(argc, argv);

  Object object = Object(mesh);
  object.register_mesh();
  scene.bind(&object);
  scene.bind(&ss_solver); // bind mesh to subspace solver


  scene.view();
  return 0;
}







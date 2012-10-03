#include "TriMesh.h"
#include "subspace/gui.hh"

using namespace trimesh;
using namespace subspace;
int main(int argc, char *argv[])
{
  TriMesh *mesh = TriMesh::read(argv[1]);
  if (!mesh) exit(1);
  
  Scene scene(argc, argv);

  Object object = Object(mesh);
  object.register_mesh();
  scene.init(&object);

  Subspace ss_solver = Subspace(argc,argv);
  scene.init(&ss_solver);

  scene.view();
  return 0;
}

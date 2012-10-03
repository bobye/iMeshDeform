#include "TriMesh.h"
#include "subspace/gui.hh"

int main(int argc, char *argv[])
{
  TriMesh *mesh = TriMesh::read(argv[1]);
  if (!mesh) exit(1);
  
  subspace::Scene scene(argc, argv);

  subspace::Object object = subspace::Object(mesh);
  object.register_mesh();
  scene.init(&object);

  subspace::Subspace ss_solver = subspace::Subspace(argc,argv);
  scene.init(&ss_solver);

  scene.view();
  return 0;
}

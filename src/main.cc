#include "subspace/mesh.hh"
#include "subspace/gui.hh"

int main(int argc, char *argv[])
{
  subspace::TriMesh mesh;
  mesh.read("test","off");
  mesh.update_internal();

  subspace::Scene scene(argc, argv);
  subspace::Object object = subspace::Object(&mesh);
  object.register_mesh();
  scene.init(&object);
  scene.view();
  return 0;
}

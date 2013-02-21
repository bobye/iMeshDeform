#include "subspace/geometry.hh"
#include "subspace/LB.hh"

using namespace subspace;

int main(int argc, char* argv[]) {
  Mesh *mesh = (Mesh *) Mesh::read(argv[1]);
  if (!mesh) exit(1);

  LB op(argc, argv);
  op.init(mesh);
  op.compute_operator();
  op.solve_eigen();

  return 0;
}

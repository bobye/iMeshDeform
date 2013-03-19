#include "subspace/geometry.hh"
#include "subspace/LB.hh"

using namespace subspace;

typedef TriangleMesh Mesh;
    

int main(int argc, char* argv[]) {
  Mesh *mesh = (Mesh *) Mesh::read(argv[1]);
  if (!mesh) exit(1);

  LB op(argc, argv, mesh);
  op.solve_eigen();
  std::system("matlab -r \"cluster_data('default', 25);exit;\" -nojvm -nodesktop -nodisplay");

  return 0;
}

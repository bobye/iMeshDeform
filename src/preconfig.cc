#include "subspace/geometry.hh"
#include "subspace/LB.hh"

using namespace subspace;

typedef TriangleMesh Mesh;
//typedef TetrahedronMesh Mesh;
    

int main(int argc, char* argv[]) {
  Mesh *mesh = (Mesh *) Mesh::read(argv[1]);
  if (!mesh) exit(1);

  LB op(argc, argv, mesh);
  op.solve_eigen();
  char cmd[1024]; int n = 25;
  n = atoi(argv[2]);
  sprintf(cmd, "matlab -r \"cluster_data('default', %d);exit;\" -nojvm -nodesktop -nodisplay", n);
  std::system(cmd);


  return 0;
}

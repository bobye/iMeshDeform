#include "subspace/geometry.hh"
#include "subspace/LB.hh"

using namespace subspace;


//#define USE_TRIANGLE
#define USE_TETRAHEDRON

#ifdef USE_TRIANGLE
typedef TriangleMesh Mesh;
#elif defined USE_TETRAHEDRON
typedef TetrahedronMesh Mesh;    
#endif

    

int main(int argc, char* argv[]) {
  Mesh *mesh = (Mesh *) Mesh::read(argv[1]);
  if (!mesh) exit(1);

#ifdef USE_TETRAHEDRON
  mesh->write("default.node");// write re-rank nodes;
  mesh->surface.write("default.off");
#endif

  LB op(argc, argv, mesh);
  op.solve_eigen();
  char cmd[1024]; int n = 25;
  n = atoi(argv[2]);
  sprintf(cmd, "matlab -r \"cluster_data('default', %d);exit;\" -nojvm -nodesktop -nodisplay", n);
  std::system(cmd);


  return 0;
}

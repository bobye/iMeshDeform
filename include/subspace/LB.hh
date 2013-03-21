#ifndef _LB_HH_
#define _LB_HH_

#include "geometry.hh"
namespace subspace {
  class LB {
  protected:
    int dimension;
  public:
    LB (int, char**, vMesh*);
    ~LB();
    void init(TriangleMesh*);
    void init(TetrahedronMesh*);

    // solve eigenpairs of LB
    void solve_eigen();
  };
}

#endif /* _LB_HH_ */

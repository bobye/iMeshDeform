#ifndef _LB_HH_
#define _LB_HH_

#include "geometry.hh"
namespace subspace {
  class LB {
  protected:

  public:
    LB (int, char**, vMesh*);
    ~LB();
    void init(Mesh*);

    // solve eigenpairs of LB
    void solve_eigen();
  };
}

#endif /* _LB_HH_ */

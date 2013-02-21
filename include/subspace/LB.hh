#ifndef _LB_HH_
#define _LB_HH_

#include "geometry.hh"
namespace subspace {
  class LB {
  protected:
    // geometry
    Mesh* mesh;
  public:
    LB (int, char**);
    ~LB();
    void init(Mesh*);

    // compute LB operator over mesh: stiff and mass matrices
    void compute_operator();
    // solve eigenpairs of LB
    void solve_eigen();
    // landmark based clustering
    void clustering();
  };
}

#endif /* _LB_HH_ */

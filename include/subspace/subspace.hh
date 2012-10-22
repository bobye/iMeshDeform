#ifndef _SUBSPACE_H_
#define _SUBSPACE_H_

#define _SS_USE_CONFORMAL
#define _SS_SHOW_DEBUG

#include <vector>
#include "geometry.hh"

#define _SS_PI                       (3.141592653589793238)

namespace subspace {
  
  class Subspace {
  protected:
    Mesh* mesh;       
    void assembly();
    
    bool ready;
  public:
    //void add_rigid_constraint(int );

    Subspace(int, char**);
    ~Subspace();
    void init(Mesh *);
    

    // vg: vertices group
    void load_linear_proxies_vg(std::vector<int> &);

    void load_rotational_proxies(std::vector<int> &);

    void solve();


    // precompute LU and start updating thread
    void prepare(std::vector< std::vector<float> > &, std::vector<Point> &);

    // online solve and update vertices
    void update(std::vector<Point> &, bool );

    // terminate online updating thread
    void terminate();
    void toggle_dump();
#ifdef _SS_SHOW_DEBUG
    void show_debug();
#endif
  };


}
#endif /* _SUBSPACE_H_ */

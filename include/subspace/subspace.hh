#ifndef _SUBSPACE_H_
#define _SUBSPACE_H_


#include <vector>
#include "TriMesh.h"

namespace subspace {
  
  class Subspace {
  protected:
    trimesh::TriMesh* mesh;       
    void assembly();
    
    bool ready;
  public:
    //void add_rigid_constraint(int );

    Subspace(int, char**);
    ~Subspace();
    void init(trimesh::TriMesh *);
    

    // vg: vertices group
    void load_linear_proxies_vg(std::vector<int> &);

    void load_rotational_proxies(std::vector<int> &);

    void solve();


    // precompute LU and start updating thread
    void prepare(std::vector< std::vector<float> > &, std::vector<trimesh::point> &);

    // online solve and update vertices
    void update(std::vector<trimesh::point> &);

    // terminate online updating thread
    void terminate();
  };


}
#endif /* _SUBSPACE_H_ */

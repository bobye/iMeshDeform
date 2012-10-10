#ifndef _SUBSPACE_H_
#define _SUBSPACE_H_


#include <vector>
#include "TriMesh.h"

namespace subspace {
  
  class Subspace {
  protected:
    trimesh::TriMesh* mesh;    

    void assembly();

  public:
    Subspace(int, char**);
    ~Subspace();
    void init(trimesh::TriMesh *);
    
    void add_rigid_transformer(bool*);
    void add_linear_constraint_handler(bool*);

    // vg: vertices group
    void load_linear_proxies_vg(std::vector<int> &);

    void load_rotational_proxies(std::vector<int> &);

    void solve();


    // online solve and update vertices
    void update(std::vector<double> );
  };

}
#endif /* _SUBSPACE_H_ */

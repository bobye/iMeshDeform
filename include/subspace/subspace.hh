#ifndef _SUBSPACE_H_
#define _SUBSPACE_H_


#include <vector>
#include "TriMesh.h"

namespace subspace {
  
  class Subspace {
  protected:
    TriMesh* mesh;    

  public:
    Subspace(int, char**);
    ~Subspace();
    void init(TriMesh *);
    
    void add_rigid_transformer(bool*);
    void add_linear_constraint_handler(bool*);

    
  };

}
#endif /* _SUBSPACE_H_ */

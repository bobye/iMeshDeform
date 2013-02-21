#ifndef _GEOMETRY_H_
#define _GEOMETRY_H_

#include "XForm.h"
#include "TriMesh.h"
#include <vector>

namespace subspace {

  typedef trimesh::point Point;
  typedef trimesh::vec   Vector;
  typedef trimesh::TriMesh Mesh2d;

  typedef trimesh::fxform XForm;
  

  typedef Mesh2d        Mesh;
}
#endif /* _GEOMETRY_H_ */

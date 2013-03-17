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
  

  //  typedef Mesh2d        Mesh;

  class vMesh {
    // generic class for triangle and tetrahedral mesh 
    // which serves as a communicator with functional package, e.g. GUI, Subspace, LB
  public:
    // geometry of surface mesh
    float *vertices_tpd;//positions of vertices to be displayed
    float *normals_tpd;
    int numberofvertices;    
    
    // geometry of elements
    float *vertices_ele;

    // virtual functions:

    // draw mesh visiable
    virtual void draw() = 0;
    // compute surface mesh normals
    virtual void recompute_normals() = 0;
    // write mesh to files
    virtual void write(const char*) = 0;

    // compute Laplacian-Beltrami operator
    virtual void compute_LB_operator() = 0;

    virtual void initialize_subspace_solver() = 0;
    // assembly large sparse matrix for solving subspaces
    virtual void compute_subspace_assembly() = 0;
    virtual void extract_display_layer() = 0;
  };
  
  class TriangleMesh : public vMesh, public Mesh2d {
    // TriangleMesh is inherited from both vMesh and Mesh2d
  public:
    static TriangleMesh *read(const char* filename) {
      TriangleMesh *mesh = new TriangleMesh();
      if (read_helper(filename, mesh)) {
	mesh->numberofvertices = mesh->vertices.size();
	mesh->allocate_data_tightpacked();// to be called only once!
	mesh->vertices_tpd = mesh->vertices_tightpacked;//bind pointer to geometric data
	mesh->normals_tpd  = mesh->normals_tightpacked;
	return mesh;
      }
      delete mesh;
      return NULL;
    }

    void draw();

    void recompute_normals() { recompute_normals_tightpacked(); }

    void write(const char* filename) {
      ((Mesh2d*) this)->write(filename);
    }

    void compute_LB_operator();

    void initialize_subspace_solver();
    void compute_subspace_assembly();
    void extract_display_layer();
  };

  /*
  class TetrahedralMesh : public vMesh {
  };
  */
  typedef TriangleMesh Mesh;
    
}
#endif /* _GEOMETRY_H_ */

#ifndef _GEOMETRY_H_
#define _GEOMETRY_H_

#include "XForm.h"
#include "TriMesh.h"
#include "TetMesh.h"
#include <vector>

namespace subspace {

  typedef trimesh::point Point;
  typedef trimesh::vec   Vector;
  typedef trimesh::TriMesh Mesh2d;
  typedef trimesh::TetMesh Mesh3d;
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
    
    // geometry of elements (optional)
    float *vertices_ele;


    std::vector<double> linear_proxies;
    // load linear_proxies via vertices group
    void load_linear_proxies_vg(std::vector<int> &);

    //rotational proxies, clusters of vertices
    std::vector<int> rotational_proxies;
    void load_rotational_proxies(std::vector<int> &);
    // virtual functions:

    virtual void load_controls(std::vector<int> &);

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
  protected:
    virtual int get_numneighbors(int ) = 0;
    virtual void compute_ARAP_approx() = 0;
  public:
    void compute_subspace_assembly();


    ~vMesh() {
      free(vertices_tpd);
      free(normals_tpd);
    }
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

    //    void load_controls(std::vector<int> &);
    void initialize_subspace_solver();
  protected:
    int get_numneighbors(int i) {return neighbors[i].size();}
    void compute_ARAP_approx();
  };


  class TetrahedronMesh : public vMesh, public Mesh3d {
  public:
    static TetrahedronMesh *read(const char* filename) {
      TetrahedronMesh *mesh = new TetrahedronMesh();
      if (read_helper(filename,mesh)) {
	mesh->numberofvertices = mesh->surface.vertices.size();
	mesh->surface.allocate_data_tightpacked();// to be called only once!
	mesh->vertices_tpd = mesh->surface.vertices_tightpacked;//bind pointer to geometric data
	mesh->normals_tpd  = mesh->surface.normals_tightpacked;
	return mesh;
      }
      delete mesh;
      return NULL;
    }

    void draw();
    void recompute_normals() { surface.recompute_normals_tightpacked();}
    void write(const char* filename) {};

    void compute_LB_operator();

    //    void load_controls(std::vector<int> &){};
    void initialize_subspace_solver();
  protected:
    int get_numneighbors(int i) {return neighbors[i].size();}
    void compute_ARAP_approx();

  };

}
#endif /* _GEOMETRY_H_ */

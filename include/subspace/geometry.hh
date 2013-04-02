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
  typedef std::vector< Point > ConstraintPointList;

  //  typedef Mesh2d        Mesh;

  class vMesh {
    // generic class for triangle and tetrahedral mesh 
    // which serves as a communicator with functional package, e.g. GUI, Subspace, LB
  public:
    int dimension;
    // geometry of surface mesh
    float *vertices_tpd;//positions of vertices to be displayed
    float *normals_tpd;
    int numberofvertices, numberofpoints;    
    
    // geometry of elements (optional)
    float *vertices_ele;
    


    // proxy I/O

    // load linear_proxies via vertices group
    std::vector<double> linear_proxies;
    void add_linear_proxies_vg(std::vector<int> &); // vg: vertices group
    void add_linear_proxies_ss(std::vector<int> &); // ss: sparse sampling
    void add_linear_proxies_custom(std::vector<double> &);// customed linear proxies

    //rotational proxies, clusters of vertices
    std::vector<int> rotational_proxies;
    void load_rotational_proxies(std::vector<int> &);

    // virtual functions:
    // load affine controller if applicable, finalize input
    std::vector<int> is_rigid;
    virtual void load_controls(std::vector<int> &);

    // draw mesh visiable
    virtual void draw() = 0;
    // compute surface mesh normals
    virtual void recompute_normals() = 0;
    // write mesh to files
    virtual void write(const char*) = 0;

    //    virtual void write_deformed_surface_mesh(const char*) = 0;
    virtual void write_deformed_surface_mesh(const char*, const unsigned char * colors =NULL) = 0;

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
	mesh->numberofpoints = mesh->numberofvertices;
	mesh->allocate_data_tightpacked();// to be called only once!
	mesh->vertices_tpd = mesh->vertices_tightpacked;//bind pointer to geometric data
	mesh->normals_tpd  = mesh->normals_tightpacked;
	mesh->dimension = 2;
	return mesh;
      }
      delete mesh;
      return NULL;
    }

    void draw();

    void recompute_normals() { recompute_normals_tightpacked(); }

    void write(const char* filename) {((Mesh2d*) this)->write(filename);}
    void write_deformed_surface_mesh(const char* filename, const unsigned char * colors) {
      FILE * f = fopen(filename, "wb");
      if (colors) fprintf(f,"C"); fprintf(f,"OFF\n");
      fprintf(f,"%d\t%d\t0\n", numberofvertices, (int) faces.size());
      for (int i=0; i<numberofvertices; ++i)  {	
	fprintf(f,"%f %f %f", vertices_tpd[3*i], vertices_tpd[3*i+1], vertices_tpd[3*i+2]);
	if (colors) 
	  fprintf(f, " %f %f %f %f\n", colors[4*i]/255., colors[4*i+1]/255., colors[4*i+2]/255., colors[4*i+3]/255.);
      }
      for (int i=0; i<faces.size(); ++i)
	fprintf(f,"3 %d %d %d\n", faces[i][0], faces[i][1], faces[i][2]);
      fclose(f);
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
	mesh->numberofpoints = mesh->nodes.size();
	mesh->surface.allocate_data_tightpacked();// to be called only once!
	mesh->vertices_tpd = mesh->surface.vertices_tightpacked;//bind pointer to geometric data
	mesh->normals_tpd  = mesh->surface.normals_tightpacked;
	mesh->dimension = 3;
	return mesh;
      }
      delete mesh;
      return NULL;
    }

    void draw();
    void recompute_normals() { surface.recompute_normals_tightpacked();}
    void write(const char* filename) { ((Mesh3d*) this)->write(filename);};
    void write_deformed_surface_mesh(const char* filename, const unsigned char * colors) {
      FILE * f = fopen(filename, "wb");
      if(colors) fprintf(f,"C"); fprintf(f,"OFF\n");
      fprintf(f,"%d\t%d\t0\n", numberofvertices, (int) surface.faces.size());
      for (int i=0; i<numberofvertices; ++i) {
	fprintf(f,"%f %f %f", vertices_tpd[3*i], vertices_tpd[3*i+1], vertices_tpd[3*i+2]);
	if (colors) 
	  fprintf(f, " %f %f %f %f\n", colors[4*i]/255., colors[4*i+1]/255., colors[4*i+2]/255., colors[4*i+3]/255.);
      }
      for (int i=0; i<surface.faces.size(); ++i)
	fprintf(f,"3 %d %d %d\n", surface.faces[i][0], surface.faces[i][1], surface.faces[i][2]);
      fclose(f);
    }

    void compute_LB_operator();

    //    void load_controls(std::vector<int> &){};
    void initialize_subspace_solver();
  protected:
    int get_numneighbors(int i) {return neighbors[i].size();}
    void compute_ARAP_approx();

  };

}
#endif /* _GEOMETRY_H_ */

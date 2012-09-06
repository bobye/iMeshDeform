#ifndef _MESH_H_
#define _MESH_H_

#include <vector>
#include "mesh_internal.hh"

namespace subspace {
  class TriMesh {
  protected:
    // CGAL Polyhedron data
    Polyhedron P;

    // index system provide random access to object handles
    ISHalfedgeList IH; // 
    ISVertexList IV;
    ISFacetList IF;

    VectorFunction vertex_norm;
    VectorFunction facet_norm;
    VectorFunction halfedge_vec;    

    ScalarFunction facet_area;
    ScalarFunction vertex_area;
    ScalarFunction vertex_avg_len;

    double update_halfedge();// to update: halfedge_vec, avg_edge_len    
    double update_facet();// to update: facet_norm, facet_area
    void update_vertex();// to update: vertex_norm, vertex_area, vertex_avg_len

    // build connection with lower CGAL layer, should be called
    // immediatelly after loading the mesh, update: IH, IV, IF
    // and [vertex, facet, halfedge]_handle->index 
    virtual void init_index();

    template <class T>
    void facet2vertex_point_average(std::vector<T> &f, std::vector<T> &v, T zero){
      double sigma;

      for (int i=0;i<vertex_num;i++){
	Vector tmp;
	double scale, total_scale=0;
	HV_circulator hv=IV[i]->vertex_begin();
	T vec = zero;

	sigma = vertex_avg_len[i] *2. /3.;
	
	do{
	  if (hv->facet()==NULL) continue;

	  tmp = (-halfedge_vec[hv->index]+halfedge_vec[hv->next()->index]);
	  scale = CGAL::sqrt(tmp * tmp) /3.;		    
	  scale = facet_area[hv->facet()->index] * std::exp(- scale * scale / (2 * sigma * sigma));	    

	  total_scale += scale;
	  vec = vec + scale * f[hv->facet()->index];

	} while (++ hv!=IV[i]->vertex_begin());


	v[i] = vec/total_scale;
      }
    };

  public:
    int halfedge_num;
    int vertex_num;
    int facet_num;

    // data accessed on run-time stage
    double*          vertex_array;
    double*          normal_array;
    unsigned int*    facet_array;


    // read mesh with format specification, for example 
    //   M.read("input", "off"); 
    // would load input.off into the TriMesh instance M
    void read(std::string , std::string);
    // like read( , ), just output mesh with specific format
    void write(std::string, std::string);
    
    // update data internal from vertices[] and facets[]
    void update_internal();
    void destroy();

  };
}

#endif /* _MESH_H_ */

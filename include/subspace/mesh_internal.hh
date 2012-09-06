#ifndef _MESH_INTERNAL_H_
#define _MESH_INTERNAL_H_

#include "mesh_precompile.hh"
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>

namespace subspace {
  typedef CGAL::Simple_cartesian<SS_SCALAR_TYPE>                       Kernel;
  //typedef CGAL::Cartesian<double>                              Kernel;
  typedef Kernel::Vector_3                                     Vector;
  //typedef Kerenl::Vector_2                                     Vector2D;
  typedef Kernel::Point_3                                      Point;
  //typedef Kernel::Point_2                                      Point2D;

  template <class Refs, class T, class P>
  class Perel_vertex : public CGAL::HalfedgeDS_vertex_base<Refs, T, P> {
  public:
    Perel_vertex() {} // repeat mandatory constructors
    Perel_vertex( const P& pt) : CGAL::HalfedgeDS_vertex_base<Refs, T, P>(pt) {}
    //CGAL::Color color;//
    int index;
  };


  template <class Refs>
  class Perel_face : public CGAL::HalfedgeDS_face_base<Refs> {
  public:
    int index;
  };

  template <class Refs>
  class Perel_halfedge : public CGAL::HalfedgeDS_halfedge_base<Refs> {
  public:
    int index;
  };

  struct Perel_items : public CGAL::Polyhedron_items_3 {
    template <class Refs, class Traits>
    struct Vertex_wrapper {
      typedef typename Traits::Point_3  Point;
      typedef Perel_vertex<Refs, CGAL::Tag_true, Point> Vertex;
    };

    template <class Refs, class Traits>
    struct Face_wrapper {
      typedef Perel_face<Refs> Face;
    };

    template <class Refs, class Traits>
    struct Halfedge_wrapper {
      typedef Perel_halfedge<Refs> Halfedge;
    };

  };






  typedef CGAL::Polyhedron_3<Kernel, Perel_items>              Polyhedron;

  typedef Polyhedron::Vertex                                   Vertex;
  typedef Polyhedron::Vertex_handle                            Vertex_handle;
  typedef Polyhedron::Vertex_iterator                          Vertex_iterator;

  typedef Polyhedron::Halfedge                                 Halfedge;
  typedef Polyhedron::Halfedge_handle                          Halfedge_handle;
  typedef Polyhedron::Halfedge_iterator                        Halfedge_iterator;

  typedef Polyhedron::Facet                                    Facet;
  typedef Polyhedron::Facet_handle                             Facet_handle;
  typedef Polyhedron::Facet_iterator                           Facet_iterator;

  typedef Polyhedron::Edge_iterator                            Edge_iterator;


  typedef Polyhedron::Halfedge_around_vertex_circulator        HV_circulator;
  typedef Polyhedron::Halfedge_around_facet_circulator         HF_circulator;


  typedef std::vector<Halfedge_handle>                         ISHalfedgeList;
  typedef std::vector<Vertex_handle>                           ISVertexList;
  typedef std::vector<Facet_handle>                            ISFacetList;


  // Feature type defined from mesh, VectorFunction represents 3D vector function
  // define over mesh domain, with respect to halfedges, vertices or facets.
  // ScalarFunction and BooleanFunction correspond to scalar and boolean
  // function defined over mesh domain.
  typedef std::vector<Vector>              VectorFunction; // 
  typedef std::vector<SS_SCALAR_TYPE>      ScalarFunction; // displayed by color ramper
  typedef std::vector<bool>                BooleanFunction;// displayed by point marker.

}
#endif

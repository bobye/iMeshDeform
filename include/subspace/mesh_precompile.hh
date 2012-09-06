#ifndef _MESH_PRECOMPILE_HH_
#define _MESH_PRECOMPILE_HH_

#define SS_VERSION                 ("0.5.0")



#define SS_USER_ATTRIBUTE_START (128)
///////////////////////////////////////////////////////////////
//Following are used in function TriMesh::attribute_allocate()
// 1st argument
#define SS_VERTEX                   0
#define SS_FACET                    1
#define SS_HALFEDGE                 2
// 2nd argument
#define SS_SCALAR                   3
#define SS_VECTOR                   4
#define SS_BOOLEAN                  5
#define SS_POINT                    6
/////////////////////////////////////////////////////////////////
// Following are registration number of inherent mesh attributes
// Scalar section: 0-31
#define SS_VERTEX_PC0               0
#define SS_VERTEX_PC1               1
#define SS_VERTEX_HCURV             2
#define SS_VERTEX_KCURV             3

#define SS_FACET_PC0                4
#define SS_FACET_PC1                5
#define SS_FACET_HCURV              6
#define SS_FACET_KCURV              7
// Vector section: 32-63
#define SS_VERTEX_NORM              32

#define SS_FACET_NORM               33
// Boolean section: 64 - 95             
#define SS_VERTEX_SALIENT           64
#define SS_VERTEX_SALIENT_SUP       65
#define SS_VERTEX_SALIENT_INF       66

// Point section: 96 - 127
#define SS_VERTEX_COORD             96






#define SS_PI                       (3.1415926535898)
#define SS_SCALAR_TYPE              double
#define SS_COLOR_HIST               (0b00000001) // histogram equalization
#define SS_COLOR_CONT               (0b00000010) // contour plot
#endif /* _MESH_PRECOMPILE_HH_ */

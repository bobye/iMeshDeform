#ifndef _GUI_H_
#define _GUI_H_

#include <vector>
#include "TriMesh.h"
#include "subspace.hh"
#include <GL/glut.h>

namespace subspace {

  class Geometry {
  public:
    // center point
    trimesh::point center;
    // transformer
    GLfloat transMat[16], transMat_buffer[16];    
    virtual void destroy();
  };


  class Object : public Geometry {
  public:
    // mesh data
    trimesh::TriMesh *mesh;    
    // bounding box
    GLfloat bbox[6], size;
    
    Object(trimesh::TriMesh*);


    // render!
    void register_mesh();
    void back_draw();
    void draw();
    virtual void destroy();
  };


  class VertSelect : public Geometry {
  protected:
    bool *buffer_selected;
    bool buffered;
    void update_color();
  public:
    Object  *object;
    int vn;//    trimesh::TriMesh *mesh;
    bool *selected;

    GLubyte *black;
    GLubyte *index;
    GLubyte *color_solid, *color_wire;

    VertSelect(Object*);
    void register_selected(int,int,int,int,bool, bool onlyone=false);
    void toggle_selected();
    void destroy();
  };

  class HandlerSelect : public Geometry {
  protected:
    std::vector< std::vector<float> > constraints; 
    std::vector< trimesh::point > constraint_points;
    /*
    std::vector< std::vector<bool> > rigids;
    std::vector< float[12] > rigid_transforms; 
    std::vector< bool > is_vertex_rigid;
    */

  public:    
    int vn;
    HandlerSelect(Object*);
    Object  *object;

    void add_rigid(bool*);
    void add_constraint(bool*);

    GLubyte *index;
    std::vector< bool > selected;
    bool register_selected(int, int, bool);
    void delete_selected();

    void set_buffer();
    void restore_buffer();

    void update();

    void draw(double);
    void destroy();
  };

  class Scene {
  protected:
    Object *object;
    VertSelect *vertsel;
    HandlerSelect *handsel;
    Geometry *context;//default set to object

    Subspace *ss_solver;//subspace solver
    void get_window_world_radio();
    void add_lights();

    void set_buffer();
    void restore_buffer();
  public:
    int width,height;
    trimesh::point cursor;
    // transformer
    static Scene * currentScene;
    Scene(int , char**);
    

    static void display();
    static void reshape(int, int);
    static void keyboard(unsigned char, int, int);
    static void skeyboard(int, int, int);
    static void motion(int, int);
    static void pmotion(int, int);
    static void mouse(int, int, int, int);
    void bind(Object*);
    void bind(Subspace*);
    void view();
  };


}
#endif /* _GUI_H_ */

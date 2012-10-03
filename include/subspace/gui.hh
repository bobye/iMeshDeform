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
    point center;
    // transformer
    GLfloat transMat[16];
    virtual void destroy();
  };


  class Object : public Geometry {
  public:
    // mesh data
    TriMesh *mesh;    
    // bounding box
    GLfloat bbox[6], size;
    
    Object(TriMesh*);


    // render!
    void register_mesh();
    void back_draw();
    void draw();
    virtual void destroy();
  };


  class VertSelect : public Geometry {
  public:
    Object  *object;
    TriMesh *mesh;
    bool *selected;

    GLubyte *black;
    GLubyte *index;
    GLubyte *color_solid, *color_wire;

    VertSelect(Object*);
    void register_selected(int,int,int,int,bool, bool onlyone=false);
    void destroy();
  };


  class Scene {
  protected:
    Object *object;
    VertSelect *vertsel;
    Geometry *context;//default set to object

    Subspace *ss_solver;//subspace solver
    void get_window_world_radio();
    void add_lights();
  public:
    int width,height;
    point cursor;
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
    void init(Object*);
    void init(Subspace*);
    void view();
  };


}
#endif /* _GUI_H_ */

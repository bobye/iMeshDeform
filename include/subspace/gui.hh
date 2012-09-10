#ifndef _GUI_H_
#define _GUI_H_

#include <vector>
#include "mesh.hh"
#include <GL/glut.h>

namespace subspace {

  class Geometry {
  public:
    // center point
    Point center;
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

    GLubyte *index;
    GLubyte *color;

    VertSelect(Object*);
    void register_selected(int,int,int,int,bool);
    void destroy();
  };


  class Scene {
  protected:
    Object *object;
    Geometry *context;//default set to object
    void get_window_world_radio();
    void add_lights();
  public:
    int width,height;
    Point cursor;
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
    void view();
  };


}
#endif /* _GUI_H_ */

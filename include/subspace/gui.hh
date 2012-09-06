#ifndef _GUI_H_
#define _GUI_H_

#include <vector>
#include "mesh.hh"
#include <GL/glut.h>

namespace subspace {

  class Object {
  public:
    GLuint LIST_NAME; //name of display list
    // mesh data
    TriMesh *mesh;
    // center point
    Point center;
    // bounding box
    GLfloat bbox[6];
    // transformer
    GLfloat transMat[16];
    // selected vertices
    bool selected_status;

    std::vector<int> selected_vertices;

    Object(TriMesh*);
    // render!
    void register_mesh();
    void draw();
  };

  class Scene {
  protected:
    Object *object;
    void get_window_world_radio();
  public:
    int width,height;
    Point cursor;
    // transformer
    static Scene * currentScene;
    Scene(int , char**);
    GLfloat transMat[16];

    static void display();
    static void reshape(int, int);
    static void keyboard(unsigned char, int, int);
    static void motion(int, int);
    static void mouse(int, int, int, int);
    void init(Object*);
    void add_lights();
    void view();
  };


}
#endif /* _GUI_H_ */

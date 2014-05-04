#ifndef _GUI_H_
#define _GUI_H_

#include <vector>
//#include "subspace.hh"
#include "geometry.hh"
#define _SS_PI                       (3.141592653589793238)
#define GL_GLEXT_PROTOTYPES




#include <GL/glut.h>

namespace subspace {

  class Geometry {
  public:
    // center point
    Point center;
    // transformer
    //GLfloat transMat[16], transMat_buffer[16];    
    XForm xf, xf_buf;
  };


  class Object : public Geometry {
  public:
    // mesh data
    Mesh *mesh;    
    // bounding box
    GLfloat bbox[6], size;
    //float *vbo;
    GLuint vbo_reg, ebo_reg;
    std::string name;

    Object();
    Object(Mesh*);


    // render!
    void register_mesh();
    void register_mesh(float *);//vbo register
    void back_draw();
    void draw();
    //    friend class Subspace;
  };


  class VertSelect : public Geometry {
  protected:
    bool *buffer_selected;
    bool buffered;
    void update_color();
  public:
    Object  *object;
    int vn;//    
    bool *selected;

    GLubyte *black;
    GLubyte *index;
    GLubyte *color_solid, *color_wire;

    VertSelect(Object*);
    ~VertSelect();
    int pick_vertex(int, int);
    void register_selected(int,int,int,int,bool, bool onlyone=false);
    void toggle_selected();
  };


  class HandlerSelect : public Geometry {
  protected:
    std::vector< std::vector<float> > constraints; 
    std::vector< Point > constraint_points;
  
  public:    
    int vn;
    HandlerSelect(Object*);
    Object  *object;
    //Subspace *ss_solver;

    void add_rigid(bool*);
    void add_constraint(bool*);

    GLubyte *index;
    std::vector< bool > selected;
    bool register_selected(int, int, bool);
    void delete_selected();
    void export_selected();
    void import_selected();
    int set_editing_selected();
    void set_selected_pos(int, int, bool reset = false);

    void set_buffer();
    void restore_buffer();

    void update(bool inf=false);

    void draw(double);


    //void set_solver(Subspace *);
    //void unset_solver();

    Point set_focus();
  };

  class Scene {
  public:
    Object *object;
    VertSelect *vertsel;
    HandlerSelect *handsel;
    Geometry *context;//default set to object

    const GLvoid *rightImage;

    int (*animator)();
    //Subspace *ss_solver;//subspace solver

    void set_buffer();
    void restore_buffer();

    Point cursor;
    // transformer
    Scene();
    Scene(int , char**);
    
    /*
    static void display();
    static void reshape(int, int);
    static void keyboard(unsigned char, int, int);
    static void skeyboard(int, int, int);
    static void motion(int, int);
    static void pmotion(int, int);
    static void mouse(int, int, int, int);
    */
    void bind(Object*);
    //void bind(Subspace*);
    void view();

    
    void set_animator( int (*func)());
    void animate();

    void read(std::string ); 
    void write(std::string );

  };


}
#endif /* _GUI_H_ */

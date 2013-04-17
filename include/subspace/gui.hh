#ifndef _GUI_H_
#define _GUI_H_

#include <vector>
#include "subspace.hh"

#define GL_GLEXT_PROTOTYPES
#include <GL/gl.h>
#include <GL/glut.h>

namespace subspace {
  class Scene;
  class Animator {
  public:    
    std::vector< std::vector<float> > constraints; 
    std::vector< ConstraintPointList > constraint_points_list;
    std::vector< XForm > projectionmatrixes;
    std::vector< XForm > modelmatrixes;
    std::vector< XForm > objectmatrixes;
    int numberofframes;
    int currentframeid;
    std::vector<float> original_constraints;
    ConstraintPointList original_constraint_points;
    XForm original_projection;
    XForm original_model;
    XForm original_object;
  public:
    Animator() {numberofframes = 0;currentframeid=0;};
    Animator(const Animator& other);
    void reset(Scene* cs=NULL);
    void clear();
    Animator merge(const Animator& other);

    bool set_constraints(Subspace* ss_solver,const std::vector< std::vector<float> >& con);
    bool add_frame(ConstraintPointList* cp, XForm* proj, XForm* model, XForm* obj);
    void remove_constrains() { constraints.clear(); constraint_points_list.clear(); }
    void remove_constrain_points() { constraint_points_list.clear(); }
    void remove_scene_matrixes() { projectionmatrixes.clear(); modelmatrixes.clear(); }
    void remove_object_matrixes() { objectmatrixes.clear(); }
    int run(Scene* currentScene);

    bool read(const char* filename);

    bool write(const char* filename);
  }; 

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
    vMesh *mesh;    
    GLuint vboId_vertices;
    bool is_vbo_updated;
    // bounding box
    GLfloat bbox[6], size;
    GLubyte *color_base, *color_render;
    Object(TriangleMesh*);
    Object(TetrahedronMesh*);
    ~Object();
    // render!
    //void register_mesh();
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
    void register_selected(int,int,int,int,bool, bool onlyone=false);
    void toggle_selected();
  };

  class HandlerSelect : public Geometry {
  protected:
    std::vector< std::vector<float> > constraints; 
    std::vector< Point > constraint_points;
    /*
    std::vector< std::vector<bool> > rigids;
    std::vector< float[12] > rigid_transforms; 
    std::vector< bool > is_vertex_rigid;
    */

  public:    
    int vn, pn;
    HandlerSelect(Object*);
    Object  *object;
    Subspace *ss_solver;

    Animator *animator;
    void add_rigid(bool*);
    void add_constraint(bool*);

    GLubyte *index;
    std::vector< bool > selected;
    bool register_selected(int, int, bool);
    void delete_selected();

    void set_buffer();
    void restore_buffer();
    void record();
    void write_record(const char* str);
    void update(bool inf=false);

    void draw(double);
    

    void set_solver(Subspace *);
    void unset_solver();

    Point set_focus();
  };

  class Scene {
  public:
    Object *object;
    VertSelect *vertsel;
    HandlerSelect *handsel;
    Geometry *context;//default set to object
    Animator *animator;
    Subspace *ss_solver;//subspace solver
    int render_mode; //0: default black scene; 1: white scene without grid, text and constraint point

    void set_buffer();
    void restore_buffer();
    
    int (*animatorfunc)();
    Point cursor;
    // transformer

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
    void bind(Subspace*);
    void view();
    void set_animator( int (*func)());
    void read(std::string ); 
    void write(std::string );

  };


}
#endif /* _GUI_H_ */

#ifndef _GUI_H_
#define _GUI_H_

#include <vector>
#include "subspace.hh"
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
  public:
    Animator() {numberofframes = 0;};
    Animator(const Animator& other) {
      this->numberofframes = other.numberofframes;
      this->constraints = other.constraints;
      this->constraint_points_list = other.constraint_points_list;
      this->projectionmatrixes = other.projectionmatrixes;
      this->modelmatrixes = other.modelmatrixes;
      this->objectmatrixes = other.objectmatrixes;
    }
    void reset() { numberofframes = 0; constraints.clear(); constraint_points_list.clear(); projectionmatrixes.clear(); modelmatrixes.clear(); objectmatrixes.clear(); }
    Animator merge(const Animator& other) {
      if(this->numberofframes>0) {
	Animator newAnimator(*this);
	if(newAnimator.constraints.empty())
	  newAnimator.constraints = other.constraints;
	if(newAnimator.constraint_points_list.empty())
	  newAnimator.constraint_points_list = other.constraint_points_list;
	if(newAnimator.projectionmatrixes.empty())
	  newAnimator.projectionmatrixes = other.projectionmatrixes;
	if(newAnimator.modelmatrixes.empty())
	  newAnimator.modelmatrixes = other.modelmatrixes;
	if(newAnimator.objectmatrixes.empty())
	  newAnimator.objectmatrixes = other.objectmatrixes;
      }
      else {
	Animator newAnimator(other);
      }
      return newAnimator;
    }
    bool set_constraints(const std::vector< std::vector<float> >& con) {
      if(con.empty()) {
	printf("No Constrains while calling set_constraints.\n");
	return false;
      }
      this->constraints = con;
      return true;
    }
    bool add_frame(ConstraintPointList* cp, XForm* proj, XForm* model, XForm* obj) {
      if(cp==NULL && proj==NULL && model==NULL && obj==NULL ){
	printf("Nothing to be added to current animator frame.\n");
	return false;
      }
      if(cp!=NULL)
	constraint_points_list.push_back(cp);
      if(proj!=NULL)
	projectionmatrixes.push_back(proj);
      if(model!=NULL)
	modelmatrixes.push_back(model);
      return true;
    }
    
    int run(Scene* currentScene) {
      return 0;
    }

    bool read(const char* filename) {
      return true;
    }

    bool write(const char* filename) {
      return true;
    }
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
    // bounding box
    GLfloat bbox[6], size;
    GLubyte *color_base;
    Object(TriangleMesh*);
    Object(TetrahedronMesh*);

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
    int animate_mode;  //0 normal mode; 1 recording mode; 2 play mode;
    void add_rigid(bool*);
    void add_constraint(bool*);

    GLubyte *index;
    std::vector< bool > selected;
    bool register_selected(int, int, bool);
    void delete_selected();

    void set_buffer();
    void restore_buffer();

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

    Subspace *ss_solver;//subspace solver
    int render_mode; //0: default black scene; 1: white scene without grid, text and constraint point

    void set_buffer();
    void restore_buffer();
    
    int (*animator)();
    
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
    int animate();
    void read(std::string ); 
    void write(std::string );

  };


}
#endif /* _GUI_H_ */

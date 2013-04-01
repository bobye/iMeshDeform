#include "subspace/gui.hh"
#include <fstream>
namespace subspace{
  Animator::Animator(const Animator& other) {
    this->currentframeid = other.currentframeid;
    this->numberofframes = other.numberofframes;
    this->constraints = other.constraints;
    this->constraint_points_list = other.constraint_points_list;
    this->projectionmatrixes = other.projectionmatrixes;
    this->modelmatrixes = other.modelmatrixes;
    this->objectmatrixes = other.objectmatrixes;
  }
  void Animator::clear() { 
    numberofframes = 0; 
    constraints.clear(); 
    constraint_points_list.clear(); 
    projectionmatrixes.clear(); 
    modelmatrixes.clear(); 
    objectmatrixes.clear(); 
    //reset();
  }
  void Animator::reset(Scene* cs) {
    currentframeid = 0;
    // reset Scene 
    
  }
  Animator Animator::merge(const Animator& other) {
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
      newAnimator.reset();
      return newAnimator;
    }
    else {
      Animator newAnimator(other);
      newAnimator.reset();
    return newAnimator;
    }
    
  }
  bool Animator::set_constraints(Subspace* ss_solver, const std::vector< std::vector<float> >& con) {
    if(con.empty()) {
      printf("No Constrains while calling set_constraints.\n");
      return false;
    }
    this->constraints = con;
    if(constraint_points_list.empty()) {
      printf("Please add one frame first, then set constraints!");
      return false;
    }
    ss_solver->prepare(constraints, constraint_points_list[0]);
    return true;
  }
  bool Animator::add_frame(ConstraintPointList* cp, XForm* proj, XForm* model, XForm* obj) {
    if(cp==NULL && proj==NULL && model==NULL && obj==NULL ){
      printf("Nothing to be added to current animator frame.\n");
      return false;
    }
    if(cp!=NULL)
      constraint_points_list.push_back(*cp);
    if(proj!=NULL)
      projectionmatrixes.push_back(*proj);
    if(model!=NULL)
      modelmatrixes.push_back(*model);
    if(obj!=NULL)
      objectmatrixes.push_back(*obj);
    numberofframes++;
    return true;
  }
    
  int Animator::run(Scene* currentScene) {
    if(!constraint_points_list.empty()){
      //set constraints' positions
      bool inf = false;
      if(currentframeid==0||currentframeid==numberofframes-1)
	inf = true;
      if (currentScene->ss_solver) currentScene->ss_solver->update(constraint_points_list[currentframeid], inf);
    }
    if(!projectionmatrixes.empty()) {
      glMatrixMode(GL_PROJECTION);
      glLoadIdentity();
      glMultMatrixf(projectionmatrixes[currentframeid]);
    }
    if(!modelmatrixes.empty()) {
      glMatrixMode(GL_MODELVIEW);    
      glLoadIdentity();
      glMultMatrixf(modelmatrixes[currentframeid]);
    }
    if(!objectmatrixes.empty()) {
      currentScene->object->xf = objectmatrixes[currentframeid];
    }
    
    currentframeid++;
    if(currentframeid>=numberofframes)
      return 0;
    return 1;
  }

  bool Animator::read(const char* filename) {
    std::ifstream ifs(filename);
    if(!ifs.is_open()){
      printf("No Such File...\n");
      return false;
    }
    //Line 1
    ifs >> numberofframes;
    int pn; //constraints degree;
    int csize,cplsize,projsize,modelsize,objsize; //constraint number, comstraint points list, projection matrixes, model matrixes, object matrixes;
    ifs >> csize >> pn >> cplsize >> projsize >> modelsize >> objsize;
    for(int i=0; i<csize; i++) {
      std::vector<float> cp;
      for(int j=0; j<pn; j++) {
	float tmpcpdata;
	ifs >> tmpcpdata ;
	cp.push_back(tmpcpdata);
      }
      constraints.push_back(cp);
    }
   
    for(int i=0; i<cplsize; i++) {
      std::vector< Point > temppointlist;
      for(int j=0; j<csize; j++) {
	Point temppoint;
	ifs >> temppoint[0] ;
	ifs >> temppoint[1] ;
	ifs >> temppoint[2] ;
	temppointlist.push_back(temppoint);
      }
      constraint_points_list.push_back(temppointlist);
    }
   
    for(int i=0; i<projsize; i++) {
      XForm tempproj;
      ifs >> tempproj;
      projectionmatrixes.push_back(tempproj);
    }
    
    for(int i=0; i<modelsize; i++) {
      XForm tempmodel;
      ifs >> tempmodel;
      modelmatrixes.push_back(tempmodel);
    }
    
    for(int i=0; i<objsize; i++) {
      XForm tempobj;
      ifs >> tempobj;
      objectmatrixes.push_back(tempobj);
    }
   
    ifs.close();
    return true;
  }

  bool Animator::write(const char* filename) {
    std::ofstream ofs(filename);
    /* 
     * L1:  number_of_frames constraint_number constraint_degree constraint_point_list_length projection_matrix_length model_matrix_length object_matrix_length
     * e.g. 10 4 10 10 0 0 10
     * L2-L(constraint_number+1):  all constraints
     * the rest lines: constraint points list and the three xforms output one by one, and one frame a line for all the for above. 
     */
    //Line 1
    if(!ofs.is_open()){
      printf("No Such File...\n");
    }
    ofs << numberofframes << " ";
    ofs << (int)constraints.size() << " ";
    if(constraints.empty())
      ofs << "0 ";
    else
      ofs << (int)constraints[0].size() << " ";
    ofs << (int)constraint_points_list.size() << " ";
    ofs << (int)projectionmatrixes.size() << " ";
    ofs << (int)modelmatrixes.size() << " ";
    ofs << (int)objectmatrixes.size() << "\n";
    //constraints
    for(int i=0; i<constraints.size(); i++) {
      for(int j=0; j<constraints[i].size(); j++) {
	ofs << constraints[i][j] << " ";
      }
      ofs << "\n";
    }
    for(int i=0; i<constraint_points_list.size(); i++) {
      for(int j=0; j<constraint_points_list[i].size(); j++) {
	ofs << constraint_points_list[i][j][0] << " ";
	ofs << constraint_points_list[i][j][1] << " ";
	ofs << constraint_points_list[i][j][2] << " ";
      }
      ofs << "\n";
    }
    for(int i=0; i<projectionmatrixes.size(); i++) {
      ofs << projectionmatrixes[i];
    }
    for(int i=0; i<modelmatrixes.size(); i++) {
      ofs << modelmatrixes[i];
    }
    for(int i=0; i<objectmatrixes.size(); i++) {
      ofs << objectmatrixes[i];
    }
    ofs.close();
    return true;
  }
}

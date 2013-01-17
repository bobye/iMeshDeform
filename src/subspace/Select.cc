#include "subspace/gui.hh"

namespace subspace {  
  VertSelect::VertSelect(Object * obj) {
    object = obj; vn = obj->mesh->vertices.size();

    selected = new bool [vn];
    buffer_selected = new bool [vn];
    color_solid = new GLubyte [4*vn];
    color_wire  = new GLubyte [4*vn];
    for (int i=0; i < vn; ++i) {
      selected[i] = false;
      color_solid[4*i] = 77; color_solid[4*i+1] = 128; color_solid[4*i+2] = 154; color_solid[4*i+3] = 192;
      color_wire[4*i] = color_wire[4*i+1] = color_wire[4*i+2] = color_wire[4*i+3] = 0;
    }
    xf = XForm::identity();//for (int i=0; i< 16; ++i) xf[i] = (i%5 ==0);


    index = new GLubyte [3*vn],
    black = new GLubyte [3*vn];    

    for (unsigned int i =0; i< vn; ++i) {
      index[3*i] = (i+1) % 256; 
      index[3*i+1] = ((i+1)>>8) % 256; 
      index[3*i+2] = ((i+1)>>16) % 256; 
      //std::cout << (int) index[3*i] << " " << (int) (index[3*i+1]<<8) << " " << (int) (index[3*i+2]<<16) << std::endl;
      black[3*i] = black[3*i+1] = black[3*i+2] = 0;
    }
    //    index = (GLubyte *) &iIndex[0];
    //    black = (GLubyte *) &iBlack[0];

  }
  void VertSelect::register_selected(int winX, int winY, int nWidth, int nHeight, bool toselect, bool onlyone){

    glDrawBuffer(GL_BACK);
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);//(NEW) setup our buffers
    /** CullFace enable */
    glEnable(GL_CULL_FACE);
    /*
    glPushMatrix();
    glPolygonMode(GL_FRONT, GL_FILL);
    glColorPointer(3, GL_UNSIGNED_BYTE, 0, black);
    object->back_draw();
    glPopMatrix();
    */
    glPushMatrix();
    glPolygonMode(GL_FRONT, GL_POINT);
    glColorPointer(3, GL_UNSIGNED_BYTE, 0, index);
    object->back_draw();
    glDisable(GL_CULL_FACE);
    glPopMatrix();

    unsigned int nTotal = nWidth * nHeight;
    GLubyte *pRGBA = new GLubyte [4*nTotal];
    glReadBuffer(GL_BACK);
    glReadPixels(winX, winY, nWidth, nHeight, GL_RGBA, GL_UNSIGNED_BYTE, &pRGBA[0]);
    int *ptr = new int[nTotal];


    for (int i=0; i<nTotal; ++i)
      ptr[i] = pRGBA[4*i] + (pRGBA[4*i+1]<<8) + (pRGBA[4*i+2]<<16);
    
    if (onlyone) {
      int mindis=(nWidth*nHeight)/4, minindex = 0;
      int center_x = nWidth/2, center_y=nHeight/2;
      for (int i=0; i<nWidth; ++i)
	for (int j=0; j<nHeight; ++j) 
	  if (ptr[nHeight*i+j] > 0){
	    int dis = (i-center_x) * (i-center_x) + (j-center_y) * (j-center_y);
	    if (mindis > dis) {
	      mindis = dis;
	      minindex = ptr[nHeight*i+j];
	    }
	  }
      if (toselect) 
	{
	  for (int i=0; i< vn; ++i) selected[i] = false;
	  if (minindex>0) selected[minindex-1] = true;
	}
      else if (minindex>0) selected[minindex-1] = !selected[minindex-1];//toggle
    } else {
      for (unsigned int i=0; i< nTotal; ++i) 
	if (ptr[i]>0) selected[ptr[i]-1] = toselect;
    }

    delete [] pRGBA;
    delete [] ptr;

    update_color();
    buffered = false;
  }


  void VertSelect::toggle_selected() {
    if (buffered) {
      for (int i=0; i<vn; ++i) selected[i] = buffer_selected[i];
      buffered = false;
    } else {
      for (int i=0; i<vn; ++i) 
	{ buffer_selected[i] = selected[i]; selected[i] = false;}    
      buffered = true;
    }
    update_color();
  }

  void VertSelect::update_color() {
    for (int i=0; i < vn; ++i) 
      if (selected[i]){
	color_solid[4*i] = 126; color_solid[4*i+1] = 123; color_solid[4*i+2] = 57; color_solid[4*i+3] = 192;	
	color_wire[4*i] = 255; color_wire[4*i+1] = 118; color_wire[4*i+2] = 0; color_wire[4*i+3] = 192;
      } else {
	color_solid[4*i] = 77; color_solid[4*i+1] = 128; color_solid[4*i+2] = 154; color_solid[4*i+3] = 192;
	color_wire[4*i] = 0; color_wire[4*i+1] = 0; color_wire[4*i+2] = 0; color_wire[4*i+3] = 0;
      }    
  }
  VertSelect::~VertSelect() {
    delete [] selected;
    delete [] black;
    delete [] index;
    delete [] color_solid;
    delete [] color_wire;
  };

  HandlerSelect::HandlerSelect(Object * obj) {
    object = obj; 
    vn = obj->mesh->vertices.size();
    //is_vertex_rigid.resize(vn);

    xf = XForm::identity();//for (int i=0; i< 16; ++i) xf[i] = (i%5 ==0);
    ss_solver = NULL;
  };

  bool HandlerSelect::register_selected(int winX, int winY, bool toselect) {
    int hn = constraint_points.size();

    index = new GLubyte[hn*3];
    for (unsigned int i =0; i< hn; ++i) { 
      index[i*3] = (i+1) % 256;
      index[i*3+1] = ((i+1)>>8) % 256;
      index[i*3+2] = ((i+1)>>16) % 256;
    }

    int nWidth = 21, nHeight = 21; winX -=10; winY -=10;

    glDrawBuffer(GL_BACK);
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);//(NEW) setup our buffers
    glPushMatrix();
    glDisable(GL_LIGHTING);
    glMultMatrixf(object->xf);
    glBegin(GL_POINTS);
    for (int i=0; i<hn; ++i) {
      glColor3ub(index[3*i], index[3*i+1], index[3*i+2]);
      glVertex3f(constraint_points[i][0], constraint_points[i][1], constraint_points[i][2]);
    }
    glEnd();
    glEnable(GL_LIGHTING);
    glPopMatrix();

    unsigned int nTotal=nWidth*nHeight;
    unsigned char *pRGBA = new unsigned char [4*nTotal];
    glReadBuffer(GL_BACK);
    glReadPixels(winX, winY, nWidth, nHeight, GL_RGBA, GL_UNSIGNED_BYTE, &pRGBA[0]);
    unsigned int *ptr = new unsigned int[nTotal]; 
    for (int i=0; i<nTotal; ++i) 
      ptr[i] = pRGBA[4*i] + (pRGBA[4*i+1]<<8) + (pRGBA[4*i+2]<<16);


    int mindis=(nWidth*nHeight)/4, minindex = 0;
    int center_x = nWidth/2, center_y=nHeight/2;

    for (int i=0; i<nWidth; ++i)
      for (int j=0; j<nHeight; ++j) 
	if (ptr[nHeight*i+j] > 0){
	  int dis = (i-center_x) * (i-center_x) + (j-center_y) * (j-center_y);
	  if (mindis > dis) {
	    mindis = dis;
	    minindex = ptr[nHeight*i+j];
	  }
	}

    if (toselect) 
      {
	for (int i=0; i< hn; ++i) selected[i] = false;
	if (minindex>0) selected[minindex-1] = true;
      }
    else if (minindex>0) selected[minindex-1] = !selected[minindex-1];//toggle
    delete [] pRGBA;     delete [] index; delete [] ptr;

    for (int i=0; i<hn; ++i) if (selected[i]) return true;
    return false;
  }

  /*
  void HandlerSelect::add_rigid(bool * selected) {
    int new_rigid_count=0, rt_count = rigids.size();
    rigids.push_back(std::vector<bool>(vn)); 
    for (int i=0; i<vn; ++i) 
      if (selected[i] && !is_vertex_rigid[i])
	new_rigid_count += (is_vertex_rigid[i] = rigids[rt_count][i] = true);

    std::cout << "Add rigid transformer (#vert " << new_rigid_count << ")" << std::endl;
  }
  */

  static std::vector<Point> constraint_points_buffer;

  void HandlerSelect::add_constraint(bool * selected) {
    int vertex_count=0, rt_count = constraints.size();
    constraints.push_back(std::vector<float>(vn)); 
    for (int i=0; i<vn; ++i)       
      vertex_count += (constraints[rt_count][i] = selected[i]);

    Point constraint_center;
    for (int i=0; i<vn; ++i) {
      constraints[rt_count][i]/=vertex_count;
      constraint_center += constraints[rt_count][i]*object->mesh->vertices[i];
    }
    constraint_points.push_back(constraint_center);
    this->selected.push_back(false);
    std::cout << "Add linear constraint (#vert " << vertex_count << ")" << std::endl;

    constraint_points_buffer = constraint_points;
    if (ss_solver) unset_solver();
  }

  void HandlerSelect::delete_selected(){
    std::vector<bool>::iterator iter = selected.begin();
    std::vector<Point>::iterator piter = constraint_points.begin();
    std::vector< std::vector<float> >::iterator viter=constraints.begin();

    while (iter < selected.end()) {
      if (*iter) {
	selected.erase(iter);
	constraint_points.erase(piter);
	constraints.erase(viter);
      } else {
	++iter; ++piter; ++viter;
      }
    }

    constraint_points_buffer = constraint_points;
    if (ss_solver) unset_solver();
  }

  /*
  void HandlerSelect::manipulate_constraint_points(GLfloat *xf){
    int hn = constraint_points.size();
    //    for (int i=0; i<hn; ++i)
    
  }
  */

  void HandlerSelect::set_buffer() {
    constraint_points_buffer = constraint_points;
  }
  void HandlerSelect::restore_buffer() {
    constraint_points = constraint_points_buffer;
    if (ss_solver) ss_solver->update(constraint_points, true);
  }

  void HandlerSelect::update(bool inf) {
    int hn = constraint_points.size();
    for (int i = 0; i < hn; ++i)
      if(selected[i])
	constraint_points[i] = xf * constraint_points_buffer[i];      
    
    if (ss_solver) ss_solver->update(constraint_points, inf);
  }

  void HandlerSelect::draw(double win_world_radio) {
    int hn = constraint_points.size();
    for (int i=0; i<hn; ++i) {

      glPushMatrix();       
      if (selected[i]) 
	glColor3f(.5, .2, .2); 
      else glColor3f(0, .5, .5);
      glTranslatef(constraint_points[i][0], constraint_points[i][1], constraint_points[i][2]);
      glutSolidSphere(7*win_world_radio, 20, 20);
      glPopMatrix();
    }
  }

  void HandlerSelect::set_solver(Subspace * ss){    
    ss_solver = ss;
    ss->prepare(constraints, constraint_points);
  }
  void HandlerSelect::unset_solver() {
    ss_solver->terminate();
    ss_solver = NULL;
  }
  /*
  void HandlerSelect::toggle_dump(Subspace *ss) {
    ss_solver = ss;
    ss_solver->toggle_dump();
    ss_solver->prepare(constraints, constraint_points);
  }
  */

  Point HandlerSelect::set_focus(){
    int hn = constraint_points.size(), count =0;
    Point pbuf(0,0,0);
    
    for (int i= 0; i<hn; ++i)
      if (selected[i]) {
	pbuf += constraint_points[i]; ++count;
      }
    if (count) 
      return  (float) (1./(float) count) * pbuf;
    else return object->center;
  }
}

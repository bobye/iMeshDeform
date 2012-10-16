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
    for (int i=0; i< 16; ++i) transMat[i] = (i%5 ==0);

    unsigned int *iIndex = new unsigned int [vn],
      *iBlack = new unsigned int [vn];    
    for (unsigned int i =0; i< vn; ++i) {iIndex[i] = i+1; iBlack[i] = 0;}
    index = (GLubyte *) &iIndex[0];
    black = (GLubyte *) &iBlack[0];

  }
  void VertSelect::register_selected(int winX, int winY, int nWidth, int nHeight, bool toselect, bool onlyone){
    glDrawBuffer(GL_BACK);
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);//(NEW) setup our buffers
    glPushMatrix();
    glPolygonMode(GL_FRONT, GL_FILL);
    glColorPointer(4, GL_UNSIGNED_BYTE, 0, black);
    object->back_draw();
    glPopMatrix();

    glPushMatrix();
    glPolygonMode(GL_FRONT, GL_POINT);
    glColorPointer(4, GL_UNSIGNED_BYTE, 0, index);
    object->back_draw();
    glPopMatrix();

    unsigned char *pRGBA = new unsigned char [4*nWidth*nHeight];
    glReadBuffer(GL_BACK);
    glReadPixels(winX, winY, nWidth, nHeight, GL_RGBA, GL_UNSIGNED_BYTE, &pRGBA[0]);
    unsigned int *ptr = (unsigned int*) & pRGBA[0], nTotal=nWidth*nHeight;

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
  void VertSelect::destroy() {
    delete [] selected;
    delete [] black;
    delete [] index;
    delete [] color_solid;
    delete [] color_wire;
  }

  HandlerSelect::HandlerSelect(Object * obj) {
    object = obj; 
    vn = obj->mesh->vertices.size();
    //is_vertex_rigid.resize(vn);

    for (int i=0; i< 16; ++i) transMat[i] = (i%5 ==0);
    ss_solver = NULL;
  };

  bool HandlerSelect::register_selected(int winX, int winY, bool toselect) {
    int hn = constraint_points.size();
    unsigned int *iIndex = new unsigned int [hn];  
    for (unsigned int i =0; i< hn; ++i) {iIndex[i] = i+1;}
    index = (GLubyte *) &iIndex[0];

    int nWidth = 21, nHeight = 21; winX -=10; winY -=10;

    glDrawBuffer(GL_BACK);
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);//(NEW) setup our buffers
    glPushMatrix();
    glDisable(GL_LIGHTING);
    glMultMatrixf(object->transMat);
    glBegin(GL_POINTS);
    for (int i=0; i<hn; ++i) {
      glColor4ub(index[4*i], index[4*i+1], index[4*i+2], index[4*i+3]);
      glVertex3f(constraint_points[i][0], constraint_points[i][1], constraint_points[i][2]);
    }
    glEnd();
    glEnable(GL_LIGHTING);
    glPopMatrix();

    unsigned char *pRGBA = new unsigned char [4*nWidth*nHeight];
    glReadBuffer(GL_BACK);
    glReadPixels(winX, winY, nWidth, nHeight, GL_RGBA, GL_UNSIGNED_BYTE, &pRGBA[0]);
    unsigned int *ptr = (unsigned int*) & pRGBA[0], nTotal=nWidth*nHeight;

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
    delete [] pRGBA;     delete [] index;

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

  static std::vector<trimesh::point> constraint_points_buffer;

  void HandlerSelect::add_constraint(bool * selected) {
    int vertex_count=0, rt_count = constraints.size();
    constraints.push_back(std::vector<float>(vn)); 
    for (int i=0; i<vn; ++i)       
      vertex_count += (constraints[rt_count][i] = selected[i]);

    trimesh::point constraint_center;
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
    std::vector<trimesh::point>::iterator piter = constraint_points.begin();
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
  void HandlerSelect::manipulate_constraint_points(GLfloat *transMat){
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
      if(selected[i]) {
	trimesh::point &x = constraint_points_buffer[i], &y = constraint_points[i];
	y[0] = transMat[12] + transMat[0] * x[0] + transMat[4] * x[1] + transMat[8] * x[2];
	y[1] = transMat[13] + transMat[1] * x[0] + transMat[5] * x[1] + transMat[9] * x[2];
	y[2] = transMat[14] + transMat[2] * x[0] + transMat[6] * x[1] + transMat[10]* x[2];
      }
    
    if (ss_solver) ss_solver->update(constraint_points, inf);
  }

  void HandlerSelect::draw(double win_world_radio) {
    int hn = constraint_points.size();
    for (int i=0; i<hn; ++i) {

      glPushMatrix();       
      if (selected[i]) 
	glColor4f(.5, .2, .2, .75); 
      else glColor4f(0, .5, .5, .75);
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

  void HandlerSelect::destroy() {
  }

}

#include "subspace/gui.hh"
#include "subspace/mesh.hh"

namespace subspace {
  VertSelect::VertSelect(Object * obj) {
    object = obj; mesh= obj->mesh;

    selected = new bool [mesh->vertex_num];
    color_solid = new GLubyte [4*mesh->vertex_num];
    color_wire  = new GLubyte [4*mesh->vertex_num];
    for (int i=0; i < mesh->vertex_num; ++i) {
      selected[i] = false;
      color_solid[4*i] = 77; color_solid[4*i+1] = 128; color_solid[4*i+2] = 154; color_solid[4*i+3] = 192;
      color_wire[4*i] = color_wire[4*i+1] = color_wire[4*i+2] = color_wire[4*i+3] = 0;
    }
    for (int i=0; i< 16; ++i) transMat[i] = (i%5 ==0);

    unsigned int *iIndex = new unsigned int [mesh->vertex_num], vn = mesh->vertex_num,
      *iBlack = new unsigned int [mesh->vertex_num];    
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
	  for (int i=0; i< mesh->vertex_num; ++i) selected[i] = false;
	  if (minindex>0) selected[minindex-1] = true;
	}
      else if (minindex>0) selected[minindex-1] = !selected[minindex-1];//toggle
    } else {
      for (unsigned int i=0; i< nTotal; ++i) 
	if (ptr[i]>0) selected[ptr[i]-1] = toselect;
    }

    delete [] pRGBA;

    for (int i=0; i < mesh->vertex_num; ++i) 
      if (selected[i]){
	color_solid[4*i] = 126; color_solid[4*i+1] = 123; color_solid[4*i+2] = 57; color_solid[4*i+3] = 192;	
	color_wire[4*i] = 255; color_wire[4*i+1] = 118; color_wire[4*i+2] = 0; color_wire[4*i+3] = 192;
      } else {
	color_solid[4*i] = 77; color_solid[4*i+1] = 128; color_solid[4*i+2] = 154; color_solid[4*i+3] = 192;
	color_wire[4*i] = 0; color_wire[4*i+1] = 0; color_wire[4*i+2] = 0; color_wire[4*i+3] = 0;
      }

    glColorPointer(4, GL_UNSIGNED_BYTE, 0, color_solid);

  }

  void VertSelect::destroy() {
    delete [] selected;
    delete [] black;
    delete [] index;
    delete [] color_solid;
    delete [] color_wire;
  }
}

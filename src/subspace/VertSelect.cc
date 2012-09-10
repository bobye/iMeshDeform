#include "subspace/gui.hh"
#include "subspace/mesh.hh"

namespace subspace {
  VertSelect::VertSelect(Object * obj) {
    object = obj; mesh= obj->mesh;

    selected = new bool [mesh->vertex_num];
    color = new GLubyte [4*mesh->vertex_num];
    for (int i=0; i < mesh->vertex_num; ++i) {
      selected[i] = false;
      color[4*i] = 77; color[4*i+1] = 128; color[4*i+2] = 154; color[4*i+3] = 192;
    }
    for (int i=0; i< 16; ++i) transMat[i] = (i%5 ==0);

    unsigned int *iIndex = new unsigned int [mesh->vertex_num], vn = mesh->vertex_num;    
    for (unsigned int i =0; i< vn; ++i) iIndex[i] = i+1;
    index = (GLubyte *) &iIndex[0];

  }
  void VertSelect::register_selected(int winX, int winY, int nWidth, int nHeight, bool toselect){
    glDrawBuffer(GL_BACK);
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);//(NEW) setup our buffers
    glPolygonMode(GL_FRONT_AND_BACK, GL_POINT);
    glColorPointer(4, GL_UNSIGNED_BYTE, 0, index);
    glPushMatrix();
    object->back_draw();
    glPopMatrix();

    unsigned char *pRGBA = new unsigned char [4*nWidth*nHeight];
    glReadBuffer(GL_BACK);
    glReadPixels(winX, winY, nWidth, nHeight, GL_RGBA, GL_UNSIGNED_BYTE, &pRGBA[0]);
    unsigned int *ptr = (unsigned int*) & pRGBA[0], nTotal=nWidth*nHeight;

    for (unsigned int i=0; i< nTotal; ++i) 
      if (ptr[i]>0) selected[ptr[i]-1] = toselect;
    
    delete [] pRGBA;

    for (int i=0; i < mesh->vertex_num; ++i) 
      if (selected[i]){
	color[4*i] = 255; color[4*i+1] = 255; color[4*i+2] = 255; color[4*i+3] = 192;	
      } else {
	color[4*i] = 77; color[4*i+1] = 128; color[4*i+2] = 154; color[4*i+3] = 192;
      }

    glColorPointer(4, GL_UNSIGNED_BYTE, 0, color);

  }

  void VertSelect::destroy() {
    delete [] selected;
    delete [] index;
    delete [] color;
  }
}

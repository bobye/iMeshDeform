#include "subspace/gui.hh"
#include "subspace/mesh.hh"

namespace subspace {
  Object::Object(TriMesh *pmesh) : mesh(pmesh){

    //compute bounding box
    for (int i=0;i < mesh->vertex_num; ++i) {
      if (i==0) {
	bbox[0] = bbox[1] = mesh->vertex_array[0];
	bbox[2] = bbox[3] = mesh->vertex_array[1];
	bbox[4] = bbox[5] = mesh->vertex_array[2];	
      }
      else {
	if (bbox[0] > mesh->vertex_array[3*i])      bbox[0]=mesh->vertex_array[3*i];
	else if (bbox[1] < mesh->vertex_array[3*i]) bbox[1]=mesh->vertex_array[3*i];
	if (bbox[2] > mesh->vertex_array[3*i+1])      bbox[2]=mesh->vertex_array[3*i+1];
	else if (bbox[3] < mesh->vertex_array[3*i+1]) bbox[3]=mesh->vertex_array[3*i+1];
	if (bbox[4] > mesh->vertex_array[3*i+2])      bbox[4]=mesh->vertex_array[3*i+2];
	else if (bbox[5] < mesh->vertex_array[3*i+2]) bbox[5]=mesh->vertex_array[3*i+2];
      }
    }

    size = bbox[1] - bbox[0];
    if (size < bbox[3] -bbox[2]) size = bbox[3] - bbox[2];
    if (size < bbox[5]- bbox[4]) size = bbox[5] - bbox[4];

    center = Point((bbox[0]+bbox[1])/2,(bbox[2]+bbox[3])/2,(bbox[4]+bbox[5])/2);

    for (int i=0; i< 16; ++i) transMat[i] = (i%5 ==0);


    //    for (int i=0; i<16;++i) std::cout << (int) index[i] << " ";
    //    std::cout<< std::endl;
  }


  void Object::register_mesh() {
    glEnableClientState(GL_NORMAL_ARRAY);
    glEnableClientState(GL_VERTEX_ARRAY);
    //load geometric information 
    glNormalPointer(GL_DOUBLE, 0, mesh->normal_array);
    glVertexPointer(3, GL_DOUBLE, 0, mesh->vertex_array);    
  }

  void Object::draw() {

    glEnable(GL_COLOR_MATERIAL);
    glColor4f(0.3, 0.5, 0.6, 0.75);

    glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT);

    //glCullFace(GL_BACK);

    //glDepthMask(GL_FALSE);
    //glEnable(GL_BLEND);

    //glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  
    glMultMatrixf(transMat);	  
    glDrawElements(GL_TRIANGLES, 3*mesh->facet_num, GL_UNSIGNED_INT, mesh->facet_array);

    //glDisable(GL_BLEND);
    //glDepthMask(GL_TRUE);
    glDisable(GL_COLOR_MATERIAL);
  }

  void Object::back_draw(){
    glDisable(GL_LIGHTING);
    glMultMatrixf(transMat);	  
    glDrawElements(GL_TRIANGLES, 3*mesh->facet_num, GL_UNSIGNED_INT, mesh->facet_array);
    glEnable(GL_LIGHTING);
  }
  
  void Object::destroy(){
  }
}


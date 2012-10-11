#include "subspace/gui.hh"
#include <string.h>

namespace subspace {
  Object::Object(trimesh::TriMesh *pmesh) : mesh(pmesh){

    //compute bounding box
    mesh->need_tstrips();
    mesh->need_bsphere();
    size = 2*mesh->bsphere.r;
    center = mesh->bsphere.center;

    for (int i=0; i< 16; ++i) transMat[i] = (i%5 ==0);

  }


  void Object::register_mesh() {
    //load geometric information 
    mesh->need_normals();
    glNormalPointer(GL_FLOAT, sizeof(mesh->normals[0]), &mesh->normals[0][0]);
    glVertexPointer(3, GL_FLOAT, sizeof(mesh->vertices[0]), &mesh->vertices[0][0]);    
  }

  void draw_tstrips(const trimesh::TriMesh *themesh)
  {
    static bool use_glArrayElement = false;
    static bool tested_renderer = false;
    if (!tested_renderer) {
      use_glArrayElement = !!strstr(
				    (const char *) glGetString(GL_RENDERER), "Intel");
      tested_renderer = true;
    }

    const int *t = &themesh->tstrips[0];
    const int *end = t + themesh->tstrips.size();
    if (use_glArrayElement) {
      while (likely(t < end)) {
	glBegin(GL_TRIANGLE_STRIP);
	int striplen = *t++;
	for (int i = 0; i < striplen; i++)
	  glArrayElement(*t++);
	glEnd();
      }
    } else {
      while (likely(t < end)) {
	int striplen = *t++;
	glDrawElements(GL_TRIANGLE_STRIP, striplen, GL_UNSIGNED_INT, t);
	t += striplen;
      }
    }
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
    //glDrawElements(GL_TRIANGLES, 3*mesh->faces.size(), GL_UNSIGNED_INT, &mesh->faces[0][0]);
    draw_tstrips(mesh);
    //glDisable(GL_BLEND);
    //glDepthMask(GL_TRUE);
    glDisable(GL_COLOR_MATERIAL);
  }

  void Object::back_draw(){
    glDisable(GL_LIGHTING);
    glMultMatrixf(transMat);	  
    //glDrawElements(GL_TRIANGLES, 3*mesh->faces.size(), GL_UNSIGNED_INT, &mesh->faces[0][0]);
    draw_tstrips(mesh);
    glEnable(GL_LIGHTING);
  }
  
  void Object::destroy(){
  }
}


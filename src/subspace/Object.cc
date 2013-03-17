#include "subspace/gui.hh"
#include <string.h>

namespace subspace {
  Object::Object(TriangleMesh *pmesh) : mesh(pmesh){

    //compute bounding box
    pmesh->need_tstrips();
    pmesh->need_bsphere();
    size = 2*pmesh->bsphere.r;
    center = pmesh->bsphere.center;

    xf = XForm::identity();//for (int i=0; i< 16; ++i) xf[i] = (i%5 ==0);

    glNormalPointer(GL_FLOAT, 0, pmesh->normals_tpd);
    glVertexPointer(3, GL_FLOAT, 0, pmesh->vertices_tpd);    

  }


  //  void draw_tstrips(const Mesh *themesh)
  void TriangleMesh::draw()
  {
    static bool use_glArrayElement = false;
    static bool tested_renderer = false;
    if (!tested_renderer) {
      use_glArrayElement = !!strstr(
				    (const char *) glGetString(GL_RENDERER), "Intel");
      tested_renderer = true;
    }

    const int *t = &tstrips[0];
    const int *end = t + tstrips.size();
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
  
    glMultMatrixf(xf);	  
    //glDrawElements(GL_TRIANGLES, 3*mesh->faces.size(), GL_UNSIGNED_INT, &mesh->faces[0][0]);
    //draw_tstrips(mesh);
    mesh->draw();
    //glDisable(GL_BLEND);
    //glDepthMask(GL_TRUE);
    glDisable(GL_COLOR_MATERIAL);
  }

  void Object::back_draw(){
    glDisable(GL_LIGHTING);
    glMultMatrixf(xf);	  
    //glDrawElements(GL_TRIANGLES, 3*mesh->faces.size(), GL_UNSIGNED_INT, &mesh->faces[0][0]);
    //draw_tstrips(mesh);
    mesh->draw();
    glEnable(GL_LIGHTING);
  }
  
}


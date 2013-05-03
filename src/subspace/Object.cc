#include "subspace/gui.hh"
#include <string.h>

namespace subspace {
  Object::Object(){
  }
  Object::Object(Mesh *pmesh) : mesh(pmesh){

    //compute bounding box
    mesh->need_tstrips();
    mesh->need_bsphere();
    size = 2*mesh->bsphere.r;
    center = mesh->bsphere.center;

    xf = XForm::identity();//for (int i=0; i< 16; ++i) xf[i] = (i%5 ==0);

  }


  void Object::register_mesh() {
    //load geometric information 
    /*
    mesh->need_normals();
    glNormalPointer(GL_FLOAT, sizeof(mesh->normals[0]), &mesh->normals[0][0]);
    glVertexPointer(3, GL_FLOAT, sizeof(mesh->vertices[0]), &mesh->vertices[0][0]);    
    */
    
    glNormalPointer(GL_FLOAT, 0, mesh->normals_tightpacked);
    glVertexPointer(3, GL_FLOAT, 0, mesh->vertices_tightpacked);
  }

  void Object::register_mesh(float *vbo) {
    //    mesh->need_normals();
    //glNormalPointer(GL_FLOAT, sizeof(mesh->normals[0]), &mesh->normals[0][0]);
    //glVertexPointer(3, GL_FLOAT, 0 , vbo);    

    glNormalPointer(GL_FLOAT, 0, mesh->normals_tightpacked);
    glVertexPointer(3, GL_FLOAT, 0, vbo); 

    /*
    glGenBuffers(1, &vbo_reg);
    glBindBuffer(GL_ARRAY_BUFFER, vbo_reg);
    glBufferData(GL_ARRAY_BUFFER, sizeof(float)*3*mesh->vertices.size(), vbo, GL_DYNAMIC_DRAW);
    */
  }


  void draw_tstrips(const Mesh *themesh)
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


  void draw_triangles(const Mesh*themesh) {

    glDrawElements(GL_TRIANGLES, 3*themesh->faces.size(), GL_UNSIGNED_INT, 
		   themesh->face_indices_tightpacked);


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
    /*
    glBindBuffer(GL_ARRAY_BUFFER, vbo_reg);
    glVertexPointer(3, GL_FLOAT, 0, 0);
    glEnableClientState(GL_VERTEX_ARRAY);
    */
    //glDrawElements(GL_TRIANGLES, 3*mesh->faces.size(), GL_UNSIGNED_INT, &mesh->faces[0][0]);
    draw_triangles(mesh);
    //draw_tstrips(mesh);
    
    
    //glDisable(GL_BLEND);
    //glDepthMask(GL_TRUE);
    glDisable(GL_COLOR_MATERIAL);
  }

  void Object::back_draw(){
    glDisable(GL_LIGHTING);
    glMultMatrixf(xf);	  
    //glDrawElements(GL_TRIANGLES, 3*mesh->faces.size(), GL_UNSIGNED_INT, &mesh->faces[0][0]);
    draw_triangles(mesh);
    //draw_tstrips(mesh);
    glEnable(GL_LIGHTING);
  }
  
}





#include "subspace/gui.hh"
#include <string.h>

namespace subspace {
  Object::Object(TriangleMesh *pmesh) : mesh(pmesh){

    //compute bounding box
    pmesh->need_tstrips();
    pmesh->need_bsphere();
    size = 2*pmesh->bsphere.r;
    center = pmesh->bsphere.center;

    xf = XForm::identity();
    int vn = pmesh->numberofvertices;


    glGenBuffers(1, &vboId_vertices);
    glBindBuffer(GL_ARRAY_BUFFER, vboId_vertices);// bind VBO in order to use for vertices
    glBufferData(GL_ARRAY_BUFFER, 3*vn*sizeof(float), pmesh->vertices_tpd, GL_DYNAMIC_DRAW);

    //glVertexPointer(3, GL_FLOAT, 0, pmesh->vertices_tpd);    
    glVertexPointer(3, GL_FLOAT, 0, 0); // last param is offset, not ptr
    glBindBuffer(GL_ARRAY_BUFFER, 0); // bind with 0, so, switch back to normal pointer operations

    glNormalPointer(GL_FLOAT, 0, pmesh->normals_tpd);


    color_base = new GLubyte[4*vn];
    color_render = new GLubyte[4*vn];
    
    for( int i=0 ; i<vn ; ++i ) {
      if(pmesh->is_rigid[i]!=0) {
	color_base[4*i]=128;color_base[4*i+1]=128;color_base[4*i+2]=128;color_base[4*i+3]=192;
	color_render[4*i]=0;color_render[4*i+1]=0;color_render[4*i+2]=255;color_render[4*i+3]=0;
      }
      else {
	color_base[4*i]=77; color_base[4*i+1]=128; color_base[4*i+2]=154; color_base[4*i+3]=192;
	color_render[4*i]=255;color_render[4*i+1]=215;color_render[4*i+2]=0; color_render[4*i+3]=0;
      }
    }
  }

  Object::Object(TetrahedronMesh *pmesh) : mesh(pmesh) {
    //compute bounding box
    pmesh->surface.need_tstrips();
    pmesh->surface.need_bsphere();
    size = 2*pmesh->surface.bsphere.r;
    center = pmesh->surface.bsphere.center;

    xf = XForm::identity();
    
    glNormalPointer(GL_FLOAT, 0, pmesh->normals_tpd);
    glVertexPointer(3, GL_FLOAT, 0, pmesh->vertices_tpd);    
    int vn = pmesh->numberofvertices;
    color_base = new GLubyte[4*vn];
    color_render = new GLubyte[4*vn];
    
    for( int i=0 ; i<vn ; ++i ) {
      if(pmesh->is_rigid[i]!=0) {
	color_base[4*i]=128;color_base[4*i+1]=128;color_base[4*i+2]=128;color_base[4*i+3]=192;
	color_render[4*i]=0;color_render[4*i+1]=0;color_render[4*i+2]=255;color_render[4*i+3]=0;
      }
      else {
	color_base[4*i]=77; color_base[4*i+1]=128; color_base[4*i+2]=154; color_base[4*i+3]=192;
	color_render[4*i]=255;color_render[4*i+1]=215;color_render[4*i+2]=0; color_render[4*i+3]=0;
      }
    }
  }

  Object::~Object() {
    delete color_base;
    delete color_render;
  }

  void draw_tstrips(const Mesh2d *mesh) {
    static bool use_glArrayElement = false;
    static bool tested_renderer = false;
    if (!tested_renderer) {
      use_glArrayElement = !!strstr((const char *) glGetString(GL_RENDERER), "Intel");
      tested_renderer = true;
    }

    const int *t = &mesh->tstrips[0];
    const int *end = t + mesh->tstrips.size();
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

  void draw_simple(const Mesh2d *mesh) {
    glDrawElements(GL_TRIANGLES, 
		   3*mesh->faces.size(), 
		   GL_UNSIGNED_INT, 
		   mesh->face_indices_tightpacked);
  }


  void TriangleMesh::draw()
  {
    draw_tstrips(this);
  }

  void TetrahedronMesh::draw()
  {
    draw_tstrips(&surface);
  }

  void Object::draw() {

    glEnable(GL_COLOR_MATERIAL);
    glColor4f(0.3, 0.5, 0.6, 0.75);

    //glCullFace(GL_BACK);

    //glDepthMask(GL_FALSE);
    //glEnable(GL_BLEND);

    //glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    
    glMultMatrixf(xf);	  

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


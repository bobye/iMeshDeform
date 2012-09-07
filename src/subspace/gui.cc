#include "subspace/gui.hh"
#include "subspace/mesh.hh"

namespace subspace {

  Scene* Scene::currentScene;

  GLfloat scale = 1., win_world_radio;
  GLdouble origin_x, origin_y, origin_z, depth, d2_x, d2_y,
    axis_x, axis_y, axis_z;
  int transform_x=0., transform_y=0.;
  GLfloat ground_wire = 30;

  bool view_translate = false , view_rotate = false, 
    object_translate = false, object_rotate = false,
    orthOrNot=false,
    wireOrNot=false;

  GLfloat CTM[16], transMat_buffer[16];

  GLfloat light_ambient[] = { .4, .4, .4, 1.0 };
  GLfloat light_diffuse[] = { .8, .8, .8, 1.0 };
  GLfloat light_specular[] = { .5, .5, .5, 1.0 };



  GLfloat light_position0[] = { 1.0, 1.0, 1.0, 0.0 };
  GLfloat light_position1[] = { -1.0, -1.0, -1.0, .0};

  GLfloat mat_ambient[] = { .5, .5, .5, 1.0 };
  GLfloat mat_emission[] = { 0, 0, 0, 0.6 };
  GLfloat mat_diffuse[] = { .5, .5, .5, .6 };
  GLfloat mat_specular[] = { .1, .1, .1, .6 };
  GLfloat mat_shininess[] = {100};

  const GLfloat perfect_factor = 1.414;

  Object::Object(TriMesh *pmesh) : mesh(pmesh){
    LIST_NAME =glGenLists(1);

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


  }

  void Object::register_mesh() {
    glEnableClientState(GL_NORMAL_ARRAY);
    glEnableClientState(GL_VERTEX_ARRAY);
    //load geometric information 
    glNormalPointer(GL_DOUBLE, 0, mesh->normal_array);
    glVertexPointer(3, GL_DOUBLE, 0, mesh->vertex_array); 
  }

  void Object::draw() {
    int fn = mesh->facet_num;

    glEnable(GL_COLOR_MATERIAL);
    glColor4f(0.3, 0.5, 0.6, 0.75);

    glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT);

    //glCullFace(GL_BACK);

    //glDepthMask(GL_FALSE);
    //glEnable(GL_BLEND);

    //glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  

    if (wireOrNot) 
      glPolygonMode(GL_FRONT, GL_LINE);
    else
      glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);



    glDrawElements(GL_TRIANGLES, 3*fn, GL_UNSIGNED_INT, mesh->facet_array);

    //glDisable(GL_BLEND);
    //glDepthMask(GL_TRUE);
    glDisable(GL_COLOR_MATERIAL);
  }



  Scene::Scene(int argc, char** argv) 
    :width(800), height(800) {
    currentScene = this; //static member need definition

    glutInit(&argc, argv);

    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
    glutInitWindowPosition(200.0, 0.0);
    glutInitWindowSize(width, height);
    glutCreateWindow("iMeshDeform - Viewer");

    glShadeModel(GL_SMOOTH);// Enable Smooth Shading
    

    glEnable(GL_LINE_SMOOTH);
    glEnable(GL_POLYGON_SMOOTH);
    //glHint(GL_LINE_SMOOTH_HINT, GL_DONT_CARE);
    //glHint(GL_POLYGON_SMOOTH_HINT, GL_DONT_CARE);
    //
    //glEnable(GL_CULL_FACE);
    //glBlendFunc(GL_SRC_ALPHA_SATURATE, GL_ONE);
    
    //glutIdleFunc(idle);
    

    glEnable(GL_LIGHTING);

    //glLightModelf(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);

    glEnable(GL_DEPTH_TEST);


  }

  void Scene::display(){

    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);//(NEW) setup our buffers



    // draw ground
    glEnable(GL_COLOR_MATERIAL);
    glColor4f(0.3, 0.3, 0.3, 0.75);


    glBegin(GL_LINES);    
    for (int i = -20; i<=20; ++i) {
      glVertex3f( i * ground_wire, 0,  (-20) * ground_wire);
      glVertex3f( i * ground_wire, 0,  (+20) * ground_wire);
      glVertex3f( (-20) * ground_wire, 0,  i * ground_wire);
      glVertex3f( (+20) * ground_wire, 0,  i * ground_wire);
    }
    glEnd();
    

    // draw object
    glPushMatrix();//push i-th matrix
    //glCallList(currentScene->object->LIST_NAME);       
    glMultMatrixf(currentScene->object->transMat);
    currentScene->object->draw();		    



    glDisable(GL_DEPTH_TEST);
    glPushMatrix();
    glEnable(GL_COLOR_MATERIAL);
    glColor4f(.5, .5, 0., 0.75);
    glTranslatef(currentScene->cursor.x(), currentScene->cursor.y(), currentScene->cursor.z());
    glutSolidSphere(5*win_world_radio , 20, 20);
    glPopMatrix();    
    glEnable(GL_DEPTH_TEST);

    glPopMatrix();//pop i-th matrix




	

    glutSwapBuffers();

    
      
  }


  void Scene::get_window_world_radio() {
    GLdouble corner1[3], corner2[3];
    GLdouble modelview[16], projection[16];
    int viewport[4];

    glGetDoublev(GL_PROJECTION_MATRIX, projection);
    glGetDoublev(GL_MODELVIEW_MATRIX, modelview);
    glGetIntegerv( GL_VIEWPORT, viewport );

    gluProject( cursor.x(), cursor.y(), cursor.z(), modelview, projection, viewport, &corner1[0], &corner1[1], &depth);

    gluUnProject( viewport[0], viewport[1], depth, modelview, 
		  projection, viewport, &corner1[0], &corner1[1], &corner1[2]);

    gluUnProject( viewport[2], viewport[3], depth, modelview, 
		  projection, viewport, &corner2[0], &corner2[1], &corner2[2]);
	
    GLdouble vect_world[3];
    vect_world[0] = corner1[0] - corner2[0];
    vect_world[1] = corner1[1] - corner2[1];
    vect_world[2] = corner1[2] - corner2[2];
    win_world_radio =std::sqrt((vect_world[0] * vect_world[0] + vect_world[1] * vect_world[1] + vect_world[2] * vect_world[2])  / (width * width + height * height));
  }

  void Scene::init(Object* obj) {
    object = obj; cursor = object->center;

    glClearColor(0.0, 0.0, 0.0, 0.0);
    glClearDepth(1.0);

    //GLfloat center_x = (object->bbox[0] + object->bbox[1]) /2.;
    //GLfloat center_y = (object->bbox[2] + object->bbox[3]) /2.;
    //GLfloat center_z = (object->bbox[4] + object->bbox[5]) /2.;
    //GLfloat length_z = (object->bbox[5] - object->bbox[4]) /2.;
    

    //std::cout << center_x << " " << center_y << " " << center_z << " " << length_z << std::endl;

    glViewport( 0, 0, width, height );

    /*
    glMatrixMode(GL_PROJECTION);

    glLoadIdentity();
    gluPerspective(30, (float)width/(float)height, 5*object->size, 20*object->size);
    glTranslated(-center_x, -center_y, -center_z - 10 * object->size );

    glMatrixMode(GL_MODELVIEW);
    */
    /*    glLoadIdentity();

    gluLookAt(center_x, center_y, center_z - 20*length_z,
	      center_x, center_y, center_z,
	      0,1,0);
    */
    //glOrtho(-1,1,-1,1,-10,10);


    scale =1.;


    add_lights();

    /*
    glNewList(object->LIST_NAME, GL_COMPILE_AND_EXECUTE);
    object->register_mesh();
    object->draw();
    glEndList();
    */
    get_window_world_radio();
    ground_wire = 5 * ( (int) (100000 * currentScene->height * win_world_radio) / 50) * 0.00001;    


    glutReshapeFunc(reshape);
    glutKeyboardFunc(keyboard);
    glutSpecialFunc(skeyboard);
    glutMotionFunc(motion);
    glutPassiveMotionFunc(pmotion);
    glutMouseFunc(mouse);
    glutDisplayFunc(display);

  }

  void Scene::add_lights(){
    /*
      light_position0[0] = coordinate_max_x * perfect_factor;
      light_position0[1] = coordinate_max_y * perfect_factor;
      light_position0[2] = coordinate_max_z * perfect_factor;

      light_position1[0] =  coordinate_min_x * perfect_factor;
      light_position1[1] =  coordinate_min_y * perfect_factor;
      light_position1[2] =  coordinate_min_z * perfect_factor;
    */

    glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
    glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
    glLightfv(GL_LIGHT0, GL_POSITION, light_position0);

    glLightfv(GL_LIGHT1, GL_AMBIENT, light_ambient);
    glLightfv(GL_LIGHT1, GL_DIFFUSE, light_diffuse);
    glLightfv(GL_LIGHT1, GL_SPECULAR, light_specular);
    glLightfv(GL_LIGHT1, GL_POSITION, light_position1);

    glEnable(GL_LIGHT0);
    glEnable(GL_LIGHT1);

    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, mat_ambient);
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mat_diffuse);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular);
    glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, mat_emission);

    glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess);


  }


  void Scene::view(){

    glutMainLoop();
  }

  
  void Scene::reshape(GLsizei w, GLsizei h)
  {
    currentScene-> width = w; currentScene->height =h;



    
    //GLfloat center_x = (currentScene->object->bbox[0] + currentScene->object->bbox[1]) /2.;
    //GLfloat center_y = (currentScene->object->bbox[2] + currentScene->object->bbox[3]) /2.;
    //GLfloat center_z = (currentScene->object->bbox[4] + currentScene->object->bbox[5]) /2.;
    //GLfloat length_z = (currentScene->object->bbox[5] - currentScene->object->bbox[4]) /2.;





    glViewport( 0, 0, w, h );

    glMatrixMode(GL_PROJECTION);

    glLoadIdentity();
    
    if (orthOrNot) {
      glOrtho(  - win_world_radio *  w/2,   + win_world_radio * w/2,
		- win_world_radio *  h/2,  +  win_world_radio * h/2,
		//100000, -100000);
		currentScene->object->size,   20* currentScene->object->size);//(NEW) set up our viewing area
    }

    else            
      gluPerspective(30*scale, (float)w/(float)h, currentScene->object->size, 20*currentScene->object->size);
    
	/*
      glFrustum(  - radio * w/perfect_factor,   radio * w/perfect_factor,
		  - radio * h/perfect_factor,  radio * h/perfect_factor,
		  //100000, -100000);
		  50* length_z,  200* length_z);//(NEW) set up our viewing area
	*/

    // glTranslated(-center_x, -center_y, -center_z - 100 * length_z);
    glTranslated(0, 0 , - 3 * currentScene->object->size);

    glMatrixMode(GL_MODELVIEW);
    keyboard('.',0,0);
    //gluLookAt( -50 * length_z , 70 * length_z ,  70 * length_z, 0, 0,  - 100 * length_z , 0, 1, 0);
    currentScene->get_window_world_radio();

  }



  void Scene::keyboard(unsigned char key, int x, int y) {
    if (key=='q'||key=='Q') exit(0); //quit
    else if (key=='o'||key=='O') { // switch between orth and perspective
      orthOrNot = !orthOrNot;
      reshape(currentScene->width, currentScene->height);
      glutPostRedisplay();
    }
    else if (key=='w' || key == 'W') { // render object as wire or solid
      wireOrNot = !wireOrNot;
      glutPostRedisplay();
    }
    else if (key=='x' || key == 'X' || key=='y' || key == 'Y' || key=='z' || key == 'Z') { // canonical views
      glGetFloatv( GL_MODELVIEW_MATRIX, CTM);
      GLfloat tx, ty, tz;
      tx = CTM[0]*CTM[12] + CTM[1]*CTM[13] + CTM[2]*CTM[14];
      ty = CTM[4]*CTM[12] + CTM[5]*CTM[13] + CTM[6]*CTM[14];
      tz = CTM[8]*CTM[12] + CTM[9]*CTM[13] + CTM[10]*CTM[14];
      glLoadIdentity();
      
      if (key=='x') glRotatef(90, 0,1,0); 
      else if (key == 'X') glRotatef(-90,0,1,0);      
      else if (key == 'y') glRotatef(90, 1,0,0); 
      else if (key == 'Y') glRotatef(-90,1,0,0);
      else if (key == 'Z') glRotatef(180, 0, 1, 0);
      glTranslated(tx,ty,tz);
      glutPostRedisplay();
    }
    else if (key == '.') { // focus to object view
      glGetFloatv( GL_MODELVIEW_MATRIX, CTM);
      GLfloat cx,cy,cz;
      cx = currentScene->cursor.x(); 
      cy = currentScene->cursor.y();
      cz = currentScene->cursor.z();
      CTM[12] = - cx * CTM[0] - cy * CTM[4] - cz * CTM[8];
      CTM[13] = - cx * CTM[1] - cy * CTM[5] - cz * CTM[9];
      CTM[14] = - cx * CTM[2] - cy * CTM[6] - cz * CTM[10];
      glLoadIdentity();
      glMultMatrixf(CTM);
      
      glutPostRedisplay();
    }
    else if (key == 'g') {
      object_translate =true;
      glPushMatrix();
      glMultMatrixf(currentScene->object->transMat);
      
      GLdouble modelview[16], projection[16];
      int viewport[4];

      glGetDoublev(GL_PROJECTION_MATRIX, projection);
      glGetDoublev(GL_MODELVIEW_MATRIX, modelview);
      glGetIntegerv( GL_VIEWPORT, viewport );
      gluUnProject(x, viewport[3] - y, depth, modelview, projection, viewport, &origin_x, &origin_y, &origin_z);
      for (int i =0;i < 16; ++i) transMat_buffer[i] = currentScene->object->transMat[i];
      glPopMatrix();
    }
    else if (key == 'r') {
      object_rotate =true;      
      glPushMatrix();
      glMultMatrixf(currentScene->object->transMat);
      GLdouble modelview[16], projection[16], t;
      int viewport[4];

      d2_x = x; d2_y = y;

      glGetDoublev(GL_PROJECTION_MATRIX, projection);
      glGetDoublev(GL_MODELVIEW_MATRIX, modelview);
      glGetIntegerv( GL_VIEWPORT, viewport );
      gluProject(currentScene->cursor.x(), currentScene->cursor.y(), currentScene->cursor.z(), modelview, projection, viewport, &origin_x, &origin_y, &origin_z);
      gluUnProject(origin_x, origin_y, 2*origin_z+1, modelview, projection, viewport, &axis_x, &axis_y, &axis_z);
      //rotation axis
      axis_x -= currentScene->cursor.x(); axis_y -= currentScene->cursor.y(); axis_z-=currentScene->cursor.z();
      t = std::sqrt(axis_x * axis_x + axis_y * axis_y + axis_z * axis_z);
      axis_x /= t; axis_y /= t; axis_z /= t; //normalize

      for (int i =0;i < 16; ++i) transMat_buffer[i] = currentScene->object->transMat[i];           
      glPopMatrix();
    }
  }

  void Scene::skeyboard(int key, int x, int y) {
    if (key == GLUT_KEY_UP) {
      glGetFloatv(GL_MODELVIEW_MATRIX, CTM);
      glLoadIdentity();
      glRotatef(30, 1, 0, 0);
      glMultMatrixf(CTM);
      glutPostRedisplay();
    }
    else if (key == GLUT_KEY_DOWN) {
      glGetFloatv(GL_MODELVIEW_MATRIX, CTM);
      glLoadIdentity();
      glRotatef(-30, 1, 0, 0);
      glMultMatrixf(CTM);
      glutPostRedisplay();
    }
    else if (key == GLUT_KEY_LEFT) {
      glGetFloatv(GL_MODELVIEW_MATRIX, CTM);
      glLoadIdentity();
      glRotatef(30, 0, 1, 0);
      glMultMatrixf(CTM);
      glutPostRedisplay();
    }
    else if (key == GLUT_KEY_RIGHT) {
      glGetFloatv(GL_MODELVIEW_MATRIX, CTM);
      glLoadIdentity();
      glRotatef(-30, 0, 1, 0);
      glMultMatrixf(CTM);
      glutPostRedisplay();
    }
  }


  void Scene::motion(int x, int y) {    
    if (view_translate) {
      
      glLoadIdentity();
      glTranslated((x-origin_x)*win_world_radio, - (y-origin_y)*win_world_radio, 0);
      glMultMatrixf(CTM);
      glutPostRedisplay();
    }

    else if (view_rotate) {
      glLoadIdentity();
      
      glRotatef(360.0 * (x-origin_x)/currentScene->width/3.14159265, 0.0, 1.0, 0.0);
      glRotatef(360.0 * (y-origin_y)/currentScene->height/3.14159265, 1.0, 0.0, 0.0);

      glMultMatrixf(CTM);	

      glutPostRedisplay();
    }


  }


  inline void MatxTranslate(GLfloat* Mat, GLfloat* BMat, GLdouble x, GLdouble y, GLdouble z) {
    for (int i=0; i<12; ++i) Mat[i] = BMat[i];
    Mat[12] = BMat[12] + BMat[0] * x + BMat[4] * y + BMat[8] * z;
    Mat[13] = BMat[13] + BMat[1] * x + BMat[5] * y + BMat[9] * z;
    Mat[14] = BMat[14] + BMat[2] * x + BMat[6] * y + BMat[10]* z;
  }

  inline void MatxMat(GLfloat*Mat, GLfloat* MMat) {
    GLfloat BMat[16];
    for (int i=0; i<16; ++i) { BMat[i] = Mat[i]; Mat[i] =0;}
    for (int i=0; i<4; ++i)
      for (int j=0; j<4; ++j)
	for (int k=0; k<4; ++k)
	  Mat[4*i+j] += BMat[4*k+j] * MMat[4*i+k];
    
  }

  void Scene::pmotion(int x, int y) {
    if (object_translate) {
      glPushMatrix();
      glMultMatrixf(transMat_buffer);

      GLdouble modelview[16], projection[16], tx, ty, tz;
      GLfloat *transMat = currentScene->object->transMat;
      int viewport[4];


      glGetDoublev(GL_PROJECTION_MATRIX, projection);
      glGetDoublev(GL_MODELVIEW_MATRIX, modelview);
      glGetIntegerv( GL_VIEWPORT, viewport );
      gluUnProject(x, viewport[3] - y, depth, modelview, projection, viewport, 
		   &tx, &ty, &tz);             
      tx -=origin_x; ty-=origin_y; tz-=origin_z;
      MatxTranslate(transMat, transMat_buffer, tx, ty, tz);
      glPopMatrix();
      glutPostRedisplay();
    }    
    else if (object_rotate) {
      GLdouble sin, cos; GLfloat RTM[16];
      GLdouble v1_x, v1_y, v2_x, v2_y, v1, v2;
      v1_x = d2_x - origin_x; v1_y = d2_y - origin_y;
      v2_x = x - origin_x; v2_y = y - origin_y;
      v1 = std::sqrt(v1_x * v1_x + v1_y * v1_y);
      v2 = std::sqrt(v2_x * v2_x + v2_y * v2_y);
      v1_x /= v1; v1_y /= v1;
      v2_x /= v2; v2_y /= v2;

      v1=v1_x - v2_x; v2=v1_y - v2_y;
      
      sin = v1_y * v2_x - v1_x * v2_y;
      cos = (1+1-(v1*v1+v2*v2))/2;

      GLdouble xy, xz, yz;

      xy = axis_x * axis_y * (1-cos);
      xz = axis_x * axis_z * (1-cos);
      yz = axis_y * axis_z * (1-cos);
      RTM[0] = axis_x * axis_x * (1-cos) + cos;
      RTM[1] = xy + axis_z * sin;
      RTM[2] = xz - axis_y * sin;

      RTM[4] = xy - axis_z * sin;
      RTM[5] = axis_y * axis_y * (1-cos) + cos;
      RTM[6] = yz + axis_x * sin;

      RTM[8] = xz + axis_y * sin;
      RTM[9] = yz - axis_x * sin;
      RTM[10]= axis_z * axis_z * (1-cos) + cos;

      RTM[3]=RTM[7]=RTM[11]=RTM[12]=RTM[13]=RTM[14]=0;
      RTM[15] =1;
      /*
      for (int i=0; i< 4; ++i){
	for (int j=0; j<4; ++j)
	  std::cout << RTM[i+4*j] << " ";
	std::cout << std::endl;
      }
      */

      GLfloat *transMat = currentScene->object->transMat;

      MatxTranslate(transMat, transMat_buffer, currentScene->cursor.x(), currentScene->cursor.y(), currentScene->cursor.z());
      MatxMat(transMat, RTM);
      MatxTranslate(transMat, transMat, -currentScene->cursor.x(), -currentScene->cursor.y(), -currentScene->cursor.z());

      glutPostRedisplay();
      
    }
  }

// compatibility with original GLUT

#if !defined(GLUT_WHEEL_UP)
#  define GLUT_WHEEL_UP   3
#  define GLUT_WHEEL_DOWN 4
#endif

#define scale_coeff       (1.1)

  void Scene::mouse(int button, int state, int x, int y) {
    switch(button) {
    case GLUT_MIDDLE_BUTTON:
      if (state == GLUT_DOWN) {
	/*
	glLoadIdentity();
	scale =1.;
	reshape(currentScene->width, currentScene->height);
	glutPostRedisplay();
	*/

	if (glutGetModifiers() == GLUT_ACTIVE_SHIFT) {
	  origin_x = x; origin_y = y; 
	  view_translate =true;
	  glGetFloatv( GL_MODELVIEW_MATRIX, CTM);
	}
	else {
	  origin_x = x; origin_y = y; view_rotate =true;
	  glGetFloatv(GL_MODELVIEW_MATRIX, CTM);	  
	}
      }
      else if (state == GLUT_UP) {
	if (view_translate) {
	  view_translate = false;
	  glutPostRedisplay();
	}
	else if (view_rotate) {
	  view_rotate = false;
	  glutPostRedisplay();	  
	}
      }
      break;
    case GLUT_LEFT_BUTTON:
      if (state == GLUT_DOWN) {	
	if (object_translate) object_translate = false;	
	if (object_rotate) object_rotate = false;
	glutPostRedisplay();
      }
      else if (state == GLUT_UP) {
      }
      break;
    case GLUT_RIGHT_BUTTON:
      if (state == GLUT_DOWN) {
	if (object_translate) {
	  object_translate = false;	
	  for (int i=0; i<16; ++i) currentScene->object->transMat[i] = transMat_buffer[i];
	}
	if (object_rotate) {
	  object_rotate = false;
	  for (int i=0; i<16; ++i) currentScene->object->transMat[i] = transMat_buffer[i];	  
	}
	glutPostRedisplay();
      }
      else if (state == GLUT_UP) {
      }

      break;
    case GLUT_WHEEL_UP:
      if (state == GLUT_UP){
	if (scale * scale_coeff < 6) scale *= scale_coeff;	
	reshape(currentScene->width, currentScene->height);
	glutPostRedisplay();      
      }
      break;
    case GLUT_WHEEL_DOWN:
      if (state == GLUT_DOWN){	
	scale /= scale_coeff;
	reshape(currentScene->width, currentScene->height);
	glutPostRedisplay();      
      }
      break;
    default:
      break;
    }

  }


}



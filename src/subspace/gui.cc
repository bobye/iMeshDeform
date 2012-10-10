#include "subspace/gui.hh"
#define SS_PI                       (3.1415926535898)

namespace subspace {


  void Geometry::destroy(){}

  Scene* Scene::currentScene;

  static GLfloat scale = 1., win_world_radio;
  static GLdouble origin_x, origin_y, origin_z, depth, d2_x, d2_y,
    axis_x, axis_y, axis_z;
  static int transform_x=0, transform_y=0, viewport[4];
  static GLfloat ground_wire = 30;

  static bool object_rotate_switch = false,
    orthOrNot=false,
    wireOrNot=false;

  char current_state = 0x00;
#define LOCK_VIEW_TRANSLATE   0x01
#define LOCK_VIEW_ROTATE      0x02
#define LOCK_OBJECT_TRANSLATE 0x04
#define LOCK_OBJECT_ROTATE    0x08
#define LOCK_BACK_BUFFER_SELECT 0x10
#define LOCK_MODE_SELECT      0x20
#define LOCK_MODE_SPEC        0x40

  std::string spec_info="";


  static GLfloat CTM[16], transMat_buffer[16];

  const GLfloat light_ambient[] = { .4, .4, .4, 1.0 };
  const GLfloat light_diffuse[] = { .8, .8, .8, 1.0 };
  const GLfloat light_specular[] = { .5, .5, .5, 1.0 };



  const GLfloat light_position0[] = { 1.0, 1.0, 1.0, 0.0 };
  const GLfloat light_position1[] = { -1.0, -1.0, -1.0, .0};

  const GLfloat mat_ambient[] = { .5, .5, .5, 1.0 };
  const GLfloat mat_emission[] = { 0, 0, 0, 0.6 };
  const GLfloat mat_diffuse[] = { .5, .5, .5, .6 };
  const GLfloat mat_specular[] = { .1, .1, .1, .6 };
  const GLfloat mat_shininess[] = {100};

  const GLfloat perfect_factor = 1.414;




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
  inline void MatxVec(GLfloat*Mat, GLdouble &x, GLdouble &y, GLdouble &z) {
    GLdouble mx,my,mz;
    mx = Mat[12] + Mat[0] * x + Mat[4] * y + Mat[8] * z;
    my = Mat[13] + Mat[1] * x + Mat[5] * y + Mat[9] * z;
    mz = Mat[14] + Mat[2] * x + Mat[6] * y + Mat[10]* z;
    x = mx; y = my; z = mz;
  }

  inline void MatxRotate(GLfloat* Mat, GLfloat* BMat, GLdouble axis_x, GLdouble axis_y, GLdouble axis_z, GLdouble sin, GLdouble cos, GLdouble center_x, GLdouble center_y, GLdouble center_z) {
    GLfloat RTM[16];
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


    MatxTranslate(Mat, BMat, center_x, center_y, center_z);
    MatxMat(Mat, RTM);
    MatxTranslate(Mat, Mat, -center_x, -center_y, -center_z);
	
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
    

    //glEnable(GL_LINE_SMOOTH);
    glEnable(GL_POLYGON_SMOOTH);
    //glLineWidth(0.);
    //glHint(GL_LINE_SMOOTH_HINT, GL_DONT_CARE);
    //glHint(GL_POLYGON_SMOOTH_HINT, GL_DONT_CARE);
    //
    //glEnable(GL_CULL_FACE);
    //glBlendFunc(GL_SRC_ALPHA_SATURATE, GL_ONE);
    
    //glutIdleFunc(idle);
    

    glEnable(GL_LIGHTING);

    //glLightModelf(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_POLYGON_OFFSET_FILL);
    glPolygonOffset(1.,1.);

  }


  void display_text() {
    glColor4f(0.0, 1.0, 0.0, 0.75); // Green
    std::string text2render;
    if (current_state & LOCK_MODE_SPEC)
      text2render =  spec_info;
    else {
      text2render = "iMeshDeform\t";
      if (current_state & LOCK_MODE_SELECT)
	text2render += "| Selection Mode\t";
      else 
	text2render += "| Normal Mode\t";
    }
    


    glRasterPos2f( 10,10 );
    for (std::string::iterator i=text2render.begin(); i!= text2render.end(); ++i) 
      glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, (char ) *i );
    
  }

  void display_action() {
    glColor4f(1, 1, 1, 0.75);
    if (current_state & LOCK_OBJECT_ROTATE) {
      glBegin(GL_LINES);

      if (object_rotate_switch) {
      } else {
	glVertex2f(origin_x, viewport[3] - origin_y);
	glVertex2f(transform_x, viewport[3] - transform_y);
      }
      glEnd();
    }
    else if (current_state & LOCK_BACK_BUFFER_SELECT) {
      glBegin(GL_LINE_LOOP);
      glVertex2f(origin_x, viewport[3] - origin_y);
      glVertex2f(origin_x, viewport[3] - transform_y);      
      glVertex2f(transform_x, viewport[3] - transform_y);      
      glVertex2f(transform_x, viewport[3] - origin_y);      
      glEnd();
    }
    

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
    


    if (current_state & LOCK_MODE_SELECT) {
      glPushMatrix();//push i-th matrix
      glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
      currentScene->object->draw();
      glPopMatrix();//pop i-th matrix
      glColorPointer(4, GL_UNSIGNED_BYTE, 0, ((VertSelect*) currentScene->context)->color_wire);
    }


    // draw object
    glPushMatrix();//push i-th matrix
    //glCallList(currentScene->object->LIST_NAME);       
    //    glMultMatrixf(currentScene->object->transMat);
    if (wireOrNot) 
      glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    else
      glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

    //    if (current_state & LOCK_MODE_SELECT) 
    //      currentScene->context->draw();
      //    else 
    currentScene->object->draw();		    
    if (current_state & LOCK_MODE_SELECT)
      glColorPointer(4, GL_UNSIGNED_BYTE, 0, ((VertSelect*) currentScene->context)->color_solid);


    glDisable(GL_DEPTH_TEST);
    glPushMatrix();
    glEnable(GL_COLOR_MATERIAL);
    glColor4f(.5, .5, 0., 0.75);
    glTranslatef(currentScene->cursor[0], currentScene->cursor[1], currentScene->cursor[2]);
    glutSolidSphere(5*win_world_radio , 20, 20);
    glPopMatrix();    
    glEnable(GL_DEPTH_TEST);

    glPopMatrix();//pop i-th matrix



    glDisable(GL_DEPTH_TEST);
    glDisable(GL_LIGHTING);
    glPushMatrix();
    glLoadIdentity();
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();

    gluOrtho2D(0.0, viewport[2], 0.0, viewport[3]);
    display_action();
    display_text();

    glEnable(GL_LIGHTING);
    glEnable(GL_DEPTH_TEST);    

    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();


    glutSwapBuffers();

    
      
  }


  void Scene::get_window_world_radio() {
    GLdouble corner1[3], corner2[3];
    GLdouble modelview[16], projection[16];

    glGetDoublev(GL_PROJECTION_MATRIX, projection);
    glGetDoublev(GL_MODELVIEW_MATRIX, modelview);
    glGetIntegerv( GL_VIEWPORT, viewport );

    gluProject( cursor[0], cursor[1], cursor[2], modelview, projection, viewport, &corner1[0], &corner1[1], &depth);

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

  void Scene::bind(Object* obj) {
    context = obj; object = obj; 
    vertsel = new VertSelect(obj);
    cursor = context->center;
    
    glClearColor(0, 0, 0, 0.0);
    glClearDepth(1.0);


    glViewport( 0, 0, width, height );



    scale =1.;


    add_lights();

    /*
    glNewList(object->LIST_NAME, GL_COMPILE_AND_EXECUTE);
    object->register_mesh();
    object->draw();
    glEndList();
    */
    reshape(width, height); keyboard('.',0,0);
    ground_wire = 5 * ( (int) (100000 * currentScene->height * win_world_radio) / 50) * 0.00001;    

    glutReshapeFunc(reshape);
    glutKeyboardFunc(keyboard);
    glutSpecialFunc(skeyboard);
    glutMotionFunc(motion);
    glutPassiveMotionFunc(pmotion);
    glutMouseFunc(mouse);
    glutDisplayFunc(display);

  }

  void Scene::bind(Subspace* ss) {
    ss_solver = ss;
    //    ss_solver->init(object->mesh);
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
		- win_world_radio *  h/2,   + win_world_radio * h/2,
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
    //keyboard('.',0,0);
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
      GLdouble cx,cy,cz;
      cx = currentScene->cursor[0]; 
      cy = currentScene->cursor[1];
      cz = currentScene->cursor[2];
      MatxVec(currentScene->context->transMat, cx, cy, cz);
      CTM[12] = - cx * CTM[0] - cy * CTM[4] - cz * CTM[8];
      CTM[13] = - cx * CTM[1] - cy * CTM[5] - cz * CTM[9];
      CTM[14] = - cx * CTM[2] - cy * CTM[6] - cz * CTM[10];
      glLoadIdentity();
      glMultMatrixf(CTM);
      
      glutPostRedisplay();
    }


    if (key == 9) {// TAB key, switch between selection mode and normal mode
      if (current_state & LOCK_MODE_SELECT) {
	current_state &= ~LOCK_MODE_SELECT;
	wireOrNot = false;	
	currentScene->context = currentScene->object;
	glDisableClientState(GL_COLOR_ARRAY);
      } else {
	current_state |= LOCK_MODE_SELECT;
	wireOrNot = true;
	currentScene->context = currentScene->vertsel;
	glColorPointer(4, GL_UNSIGNED_BYTE, 0, ((VertSelect*) currentScene->context)->color_solid);
	glEnableClientState(GL_COLOR_ARRAY);
      }
      glutPostRedisplay();
    }

    if (!(current_state & 0xf0)) {//Normal mode      
      if (key == 'g' && !(current_state & ~LOCK_OBJECT_TRANSLATE)) {      
	if (glutGetModifiers() == GLUT_ACTIVE_ALT) {
	  GLfloat *transMat = currentScene->context->transMat;
	  transMat[12] = currentScene->cursor[0];
	  transMat[13] = currentScene->cursor[1];
	  transMat[14] = currentScene->cursor[2];
	  MatxTranslate(transMat, transMat, -currentScene->cursor[0], -currentScene->cursor[1], -currentScene->cursor[2]);
	  glutPostRedisplay();
	}else {
	  current_state |= LOCK_OBJECT_TRANSLATE;	
	
	  glPushMatrix();
	  glMultMatrixf(currentScene->context->transMat);
      
	  GLdouble modelview[16], projection[16];

	  glGetDoublev(GL_PROJECTION_MATRIX, projection);
	  glGetDoublev(GL_MODELVIEW_MATRIX, modelview);
	  gluUnProject(x, viewport[3] - y, depth, modelview, projection, viewport, &origin_x, &origin_y, &origin_z);
	  glPopMatrix();

	  for (int i =0;i < 16; ++i) transMat_buffer[i] = currentScene->context->transMat[i];
	}
      }
      else if (key == 'r' && !(current_state & ~LOCK_OBJECT_ROTATE)) {
	if (glutGetModifiers() == GLUT_ACTIVE_ALT) {
	  GLfloat *transMat = currentScene->context->transMat;
	  MatxTranslate(transMat, transMat, currentScene->cursor[0], currentScene->cursor[1], currentScene->cursor[2]);
	  for (int i=0; i < 12; ++i) transMat[i] = (i%5==0);
	  MatxTranslate(transMat, transMat, -currentScene->cursor[0], -currentScene->cursor[1], -currentScene->cursor[2]);
	  glutPostRedisplay(); return;
	}
      
	if (current_state & LOCK_OBJECT_ROTATE){
	  object_rotate_switch =!object_rotate_switch;
	  for (int i = 0; i < 16; ++i) currentScene->context->transMat[i] = transMat_buffer[i];
	  glutPostRedisplay();
	} else {
	  for (int i =0;i < 16; ++i) transMat_buffer[i] = currentScene->context->transMat[i];
	  current_state |= LOCK_OBJECT_ROTATE;
	}

	if (object_rotate_switch) {
	  glPushMatrix();
	  glMultMatrixf(transMat_buffer);
	  GLdouble modelview[16], projection[16];

	  d2_x = x; d2_y = y;

	  glGetDoublev(GL_PROJECTION_MATRIX, projection);
	  glGetDoublev(GL_MODELVIEW_MATRIX, modelview);
	  glPopMatrix();
	  gluProject(currentScene->cursor[0], currentScene->cursor[1], currentScene->cursor[2], modelview, projection, viewport, &origin_x, &origin_y, &origin_z); 


	} else {
	  glPushMatrix();
	  glMultMatrixf(transMat_buffer);
	  GLdouble modelview[16], projection[16], t;

	  d2_x = x; d2_y = y;

	  glGetDoublev(GL_PROJECTION_MATRIX, projection);
	  glGetDoublev(GL_MODELVIEW_MATRIX, modelview);
	  gluProject(currentScene->cursor[0], currentScene->cursor[1], currentScene->cursor[2], modelview, projection, viewport, &origin_x, &origin_y, &origin_z); 
	  glPopMatrix();	
	  gluUnProject(origin_x, origin_y, 0, modelview, projection, viewport, &axis_x, &axis_y, &axis_z);
	  //rotation axis
	  axis_x -= currentScene->cursor[0]; axis_y -= currentScene->cursor[1]; axis_z-=currentScene->cursor[2];
	  t = std::sqrt(axis_x * axis_x + axis_y * axis_y + axis_z * axis_z);
	  axis_x /= t; axis_y /= t; axis_z /= t; //normalize	
	  origin_y = viewport[3] - origin_y;

	}
      }
    } else if (current_state & LOCK_MODE_SELECT) {   
      if (key == 'b') {
	current_state |= LOCK_BACK_BUFFER_SELECT;     
      }
      else if (key == 'A') {
	current_state |= LOCK_MODE_SPEC;
	spec_info = "[r]: Add rigid transformer [h]: Add handler";
	glutPostRedisplay();
      }
      else if (key == 'r' && (current_state & LOCK_MODE_SPEC)) {
	current_state &= ~LOCK_MODE_SPEC;
	spec_info = "";
	// add rigid transformer
	currentScene->ss_solver->add_rigid_transformer(((VertSelect*) currentScene->context)->selected);
	glutPostRedisplay();
      }
      else if (key == 'h' && (current_state & LOCK_MODE_SPEC)) {
	current_state &= ~LOCK_MODE_SPEC;
	spec_info = "";
	// add linear constraint handler
	currentScene->ss_solver->add_linear_constraint_handler(((VertSelect*) currentScene->context)->selected);
	glutPostRedisplay();
      }
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
    if (current_state & LOCK_VIEW_TRANSLATE) {
      
      glLoadIdentity();
      glTranslated((x-origin_x)*win_world_radio, - (y-origin_y)*win_world_radio, 0);
      glMultMatrixf(CTM);
      glutPostRedisplay();
    }

    else if (current_state & LOCK_VIEW_ROTATE) {
      glLoadIdentity();
      
      glRotatef(360.0 * (x-origin_x)/currentScene->width/SS_PI, 0.0, 1.0, 0.0);
      glRotatef(360.0 * (y-origin_y)/currentScene->height/SS_PI, 1.0, 0.0, 0.0);

      glMultMatrixf(CTM);	

      glutPostRedisplay();
    }
    else if (current_state & LOCK_BACK_BUFFER_SELECT) {      
      transform_x = x; transform_y =y;
      glutPostRedisplay();
    }


  }



  void Scene::pmotion(int x, int y) {
    if (current_state & LOCK_OBJECT_TRANSLATE) {
      glPushMatrix();
      glMultMatrixf(transMat_buffer);

      GLdouble modelview[16], projection[16], tx, ty, tz;
      GLfloat *transMat = currentScene->context->transMat;


      glGetDoublev(GL_PROJECTION_MATRIX, projection);
      glGetDoublev(GL_MODELVIEW_MATRIX, modelview);
      glPopMatrix();
      gluUnProject(x, viewport[3] - y, depth, modelview, projection, viewport, &tx, &ty, &tz);
      tx -=origin_x; ty-=origin_y; tz-=origin_z;
      MatxTranslate(transMat, transMat_buffer, tx, ty, tz);

      glutPostRedisplay();
    }    
    else if (current_state & LOCK_OBJECT_ROTATE) {
      transform_x = x; transform_y = y;
      if (object_rotate_switch) {
	glPushMatrix();
	glMultMatrixf(transMat_buffer);
	GLdouble modelview[16], projection[16],t;


	glGetDoublev(GL_PROJECTION_MATRIX, projection);
	glGetDoublev(GL_MODELVIEW_MATRIX, modelview);
	glPopMatrix();

	gluUnProject(origin_x + (y-d2_y), origin_y + (x-d2_x), origin_z, modelview, projection, viewport,  &axis_x, &axis_y, &axis_z);
	axis_x -= currentScene->cursor[0];
	axis_y -= currentScene->cursor[1];
	axis_z -= currentScene->cursor[2];
	t = std::sqrt(axis_x * axis_x + axis_y * axis_y + axis_z * axis_z);
	axis_x /= t; axis_y /= t; axis_z /= t; //normalize	

	t = 2*SS_PI * std::sqrt(((x-d2_x)*(x-d2_x) + (y-d2_y)*(y-d2_y))/(viewport[2]*viewport[2]+viewport[3]*viewport[3])) ;
	GLdouble sin,cos;
	sin = std::sin(t); cos = std::cos(t);
	MatxRotate(currentScene->context->transMat, transMat_buffer, axis_x, axis_y, axis_z, sin, cos, currentScene->cursor[0], currentScene->cursor[1], currentScene->cursor[2]);

      } else {
	GLdouble sin, cos; 
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
	MatxRotate(currentScene->context->transMat, transMat_buffer, axis_x, axis_y, axis_z, sin, cos, currentScene->cursor[0], currentScene->cursor[1], currentScene->cursor[2]);

      }
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
      if (state == GLUT_DOWN) { //&& !(current_state & ~LOCK_VIEW_TRANSLATE)) {
	if (current_state & LOCK_BACK_BUFFER_SELECT) {
	  origin_x = x; origin_y = y;
	} else {
	  if (glutGetModifiers() == GLUT_ACTIVE_SHIFT) {
	    origin_x = x; origin_y = y; 
	    current_state |= LOCK_VIEW_TRANSLATE;

	    glGetFloatv( GL_MODELVIEW_MATRIX, CTM);
	  }
	  else {
	    origin_x = x; origin_y = y; 
	    current_state |= LOCK_VIEW_ROTATE;
	    glGetFloatv(GL_MODELVIEW_MATRIX, CTM);	  
	  }
	}
      }
      else if (state == GLUT_UP) {
	if (current_state & LOCK_VIEW_TRANSLATE) {
	  current_state &= ~LOCK_VIEW_TRANSLATE;
	  glutPostRedisplay();
	}
	else if (current_state & LOCK_VIEW_ROTATE) {
	  current_state &= ~LOCK_VIEW_ROTATE;
	  glutPostRedisplay();	  
	}
	else if (current_state & LOCK_BACK_BUFFER_SELECT) {
	  current_state &= ~LOCK_BACK_BUFFER_SELECT;
	  int mx, Mx, my, My, nHeight, nWidth;
	  if (origin_x > x) { mx = x; Mx = origin_x;} else {mx = origin_x; Mx = x;}
	  if (origin_y > y) { my = y; My = origin_y;} else {my = origin_y; My = y;}

	  nWidth = Mx -mx +1; nHeight = My-my +1;
	  
	  ((VertSelect*) currentScene->context)->register_selected(mx, viewport[3]-My, nWidth, nHeight, false);
	  glutPostRedisplay();	  
	}
      }
      break;
    case GLUT_LEFT_BUTTON:
      if (state == GLUT_DOWN) {	
	if (current_state & LOCK_OBJECT_TRANSLATE) { current_state &= ~LOCK_OBJECT_TRANSLATE; glutPostRedisplay(); }
	else if (current_state & LOCK_OBJECT_ROTATE) { current_state &= ~LOCK_OBJECT_ROTATE; object_rotate_switch = false; glutPostRedisplay();}
	else if (current_state & LOCK_BACK_BUFFER_SELECT) { origin_x = x; origin_y = y;}

      }
      else if (state == GLUT_UP) {
	if (current_state & LOCK_BACK_BUFFER_SELECT) {
	  current_state &= ~LOCK_BACK_BUFFER_SELECT;
	  int mx, Mx, my, My, nHeight, nWidth;
	  if (origin_x > x) { mx = x; Mx = origin_x;} else {mx = origin_x; Mx = x;}
	  if (origin_y > y) { my = y; My = origin_y;} else {my = origin_y; My = y;}

	  nWidth = Mx -mx +1; nHeight = My-my +1;
	  
	  ((VertSelect*) currentScene->context)->register_selected(mx, viewport[3]-My, nWidth, nHeight, true);
	  glutPostRedisplay();
	}
      }
      break;
    case GLUT_RIGHT_BUTTON:
      if (state == GLUT_DOWN) {
	if (current_state & LOCK_OBJECT_TRANSLATE) {
	  current_state &= ~LOCK_OBJECT_TRANSLATE;
	  for (int i=0; i<16; ++i) currentScene->context->transMat[i] = transMat_buffer[i];
	}
	if (current_state & LOCK_OBJECT_ROTATE) {
	  current_state &= ~LOCK_OBJECT_ROTATE;
	  object_rotate_switch = false;
	  for (int i=0; i<16; ++i) currentScene->context->transMat[i] = transMat_buffer[i];	  
	}

	if (current_state & LOCK_MODE_SELECT) {
	  if (glutGetModifiers() == GLUT_ACTIVE_SHIFT) 
	    ((VertSelect*) currentScene->context)->register_selected(x-10, viewport[3]-y-10, 21, 21, false, true);	  
	  else 
	    ((VertSelect*) currentScene->context)->register_selected(x-10, viewport[3]-y-10, 21, 21, true, true);	  
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



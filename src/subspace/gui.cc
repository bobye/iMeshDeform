#include "subspace/gui.hh"
#include <fstream>
#include <iostream>

/*
#include "sys/time.h"
#define GetTime(x) clock_gettime(CLOCK_MONOTONIC, x)
#define TimeType struct timespec
#define BILLION  1000000000L
#define DiffTime(s,e) (( s.tv_sec - e.tv_sec )+ (double)( s.tv_nsec - e.tv_nsec ) / (double)BILLION)
*/

#define GetTime(x) x = glutGet(GLUT_ELAPSED_TIME);
#define TimeType unsigned int
#define DiffTime(s,e) ((double) (s-e)/ (double) 1000)

namespace subspace {

  TimeType time, timebase;
  static unsigned frame =0;
  char s_fps[10];

  Scene* currentScene;

  static GLfloat  win_world_radio;
  static GLfloat camera_displace;
  static GLdouble origin_x, origin_y, origin_z, depth, d2_x, d2_y,
    axis_x, axis_y, axis_z;
  static GLint transform_x=0, transform_y=0, viewport[4], width, height;
  static GLfloat ground_wire = 30;

  static bool object_rotate_switch = false,
    orthOrNot=false,
    wireOrNot=false;
  static Geometry *context_buffer;

  char current_state = 0x00;
#define LOCK_VIEW_TRANSLATE   0x01
#define LOCK_VIEW_ROTATE      0x02
#define LOCK_OBJECT_TRANSLATE 0x04
#define LOCK_OBJECT_ROTATE    0x08
#define LOCK_BACK_BUFFER_SELECT 0x10
#define LOCK_MODE_SELECT      0x20
#define LOCK_MODE_SPEC        0x40
  //#define LOCK_MODE_DEFORM      0x80

  std::string spec_info="";


  static XForm CTM;//static GLfloat CTM[16];  

  const GLfloat light_ambient[][4] = {{ .4, .4, .4, 1.0 }, { .15, .15, .15, 1.0 }};
  const GLfloat light_diffuse[][4] = {{ .8, .8, .8, 1.0 }, { 1, 1, 1, 1.0 }};
  const GLfloat light_specular[][4] = {{ .5, .5, .5, 1.0 }, { .3, .3, .3, 1.0 }};



  const GLfloat light_position0[][4] = {{ 1.0, 1.0, 1.0, 0.0 },{ 0.0, 0.0, 1.0, 0.0 }} ;
  const GLfloat light_position1[][4] = {{ -1.0, -1.0, -1.0, .0}, { -1.0, -1.0, -1.0, .0}};

  const GLfloat mat_ambient[][4] = {{ .5, .5, .5, 1.0 }, { .0, .0, .0, 1.0 }};
  const GLfloat mat_emission[][4] = {{ 0, 0, 0, 0.6 }, { .0, .0, .0, 0.6 }};
  const GLfloat mat_diffuse[][4] = {{ .5, .5, .5, .6 }, { 1., 1., 1., .6 }};
  const GLfloat mat_specular[][4] = {{ .1, .1, .1, .6 }, { .8, .8, .8, .6 }};
  const GLfloat mat_shininess[] = {100, 30};

  const GLfloat perfect_factor = 1.414;

  int record_switch = 0;

  void get_window_world_radio();
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
  {
    width = height = 800;
    currentScene = this; //static member need definition
    render_mode = 0;
    glutInit(&argc, argv);

    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
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
    glPolygonOffset(1.,1.);
    glEnable(GL_POLYGON_OFFSET_FILL);    
    //glDisable(GL_POLYGON_OFFSET_LINE);
    //glDisable(GL_POLYGON_OFFSET_POINT);

    glEnable(GL_DEPTH_TEST);

    glEnableClientState(GL_NORMAL_ARRAY);
    glEnableClientState(GL_VERTEX_ARRAY);

    GetTime(timebase)

  }


  void Scene::set_animator( int (*func)()) {
    animatorfunc = func;
  }

  int record_animate() {

    if(currentScene->animator==NULL)
      currentScene->animator = new Animator();
    XForm proj, model;
    glGetFloatv( GL_PROJECTION_MATRIX, proj);
    glGetFloatv( GL_MODELVIEW_MATRIX, model);
    currentScene->animator->add_frame(NULL, &proj, &model, &currentScene->object->xf);
    currentScene->handsel->record();
    return 1;
  }
  int play_animate() {
    return currentScene->animator->run(currentScene);
  }

  void animate() {

    if ((*currentScene->animatorfunc)() == 0) {
      glutIdleFunc(NULL);
      delete currentScene->animator;
      currentScene->animator = NULL;
    }
    get_window_world_radio();
    glutPostRedisplay();
  }


  void display_text() {
    glColor3f(0.0, 1.0, 0.0); // Green
    std::string text2render;
    if (current_state & LOCK_MODE_SPEC)
      text2render =  spec_info;
    else {
      text2render = "iMeshDeform\t";
      if (current_state & LOCK_MODE_SELECT)
	text2render += "| Select Mode\t";
      else
	text2render += "| Normal Mode\t";
    }
    
    text2render += "| " + std::string(s_fps);


    glRasterPos2f( 10,10 );
    for (std::string::iterator i=text2render.begin(); i!= text2render.end(); ++i) 
      glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, (char ) *i );
    
  }

  void display_action() {
    glColor3f(1, 1, 1);
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

  void display(){

    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);//(NEW) setup our buffers

    if(currentScene->render_mode==0) {
      // draw ground
      glClearColor(0.0,0.0,0.0,0.0);
      glEnable(GL_COLOR_MATERIAL);
      glColor3f(0.3, 0.3, 0.3);

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
	glColorPointer(3, GL_UNSIGNED_BYTE, 4, ((VertSelect*) currentScene->context)->color_solid);
	currentScene->object->draw();
	glPopMatrix();//pop i-th matrix
	glColorPointer(3, GL_UNSIGNED_BYTE, 4, ((VertSelect*) currentScene->context)->color_wire);
      }  


      // draw object
      glPushMatrix();//push i-th matrix
      if (wireOrNot) 
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
      else
	{
	  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	  glColorPointer(3, GL_UNSIGNED_BYTE, 4, currentScene->object->color_base);
	}
      currentScene->object->draw();		    
      glDisable(GL_DEPTH_TEST);
      glEnable(GL_COLOR_MATERIAL);


      currentScene->handsel->draw(win_world_radio);


      glPushMatrix();
      glColor3f(.5, .5, 0.);
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
   }
    else
    {
      glClearColor(1.0,1.0,1.0,0.0);
      glColor3f(0.3, 0.3, 0.3);
      glEnable(GL_COLOR_MATERIAL);
      glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
      glPushMatrix();
      glColorPointer(3, GL_UNSIGNED_BYTE, 4, currentScene->object->color_base);
      currentScene->object->draw();		    
      
      glDisable(GL_DEPTH_TEST);
      glEnable(GL_COLOR_MATERIAL);


      currentScene->handsel->draw(win_world_radio);

      glPushMatrix();
      
      glColor3f(.5, .5, 0.);
      glTranslatef(currentScene->cursor[0], currentScene->cursor[1], currentScene->cursor[2]);
      glutSolidSphere(5*win_world_radio , 20, 20);
      glPopMatrix();    
      glEnable(GL_DEPTH_TEST);

      glPopMatrix();//pop i-th matrix
      

    }    
    ++frame;
    GetTime(time)
    double accum = DiffTime(time,timebase);
    if (accum > 1) {
      sprintf(s_fps,"FPS:%4.1f", frame/accum);
      timebase = time;
      frame = 0;
    }
    glutSwapBuffers();

    
      
  }


  void get_window_world_radio() {
    GLdouble corner1[3], corner2[3];
    GLdouble modelview[16], projection[16];

    glGetDoublev(GL_PROJECTION_MATRIX, projection);
    glGetDoublev(GL_MODELVIEW_MATRIX, modelview);
    glGetIntegerv( GL_VIEWPORT, viewport );

    gluProject( currentScene->cursor[0], currentScene->cursor[1], currentScene->cursor[2], modelview, projection, viewport, &corner1[0], &corner1[1], &depth);

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


  void add_lights(int mode = 0){

    glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient[mode]);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse[mode]);
    glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular[mode]);
    glLightfv(GL_LIGHT0, GL_POSITION, light_position0[mode]);


    glLightfv(GL_LIGHT1, GL_AMBIENT, light_ambient[mode]);
    glLightfv(GL_LIGHT1, GL_DIFFUSE, light_diffuse[mode]);
    glLightfv(GL_LIGHT1, GL_SPECULAR, light_specular[mode]);
    glLightfv(GL_LIGHT1, GL_POSITION, light_position1[mode]);


    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, mat_ambient[mode]);
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mat_diffuse[mode]);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular[mode]);
    glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, mat_emission[mode]);

    glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, &mat_shininess[mode]);


  }


  void Scene::view(){

    glutMainLoop();
  }

  
  void reshape(GLsizei w, GLsizei h)
  {
    width = w; height =h;



    
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

    else {           
      gluPerspective(30, (float)w/(float)h, currentScene->object->size, 20*currentScene->object->size);
    }
	/*
      glFrustum(  - radio * w/perfect_factor,   radio * w/perfect_factor,
		  - radio * h/perfect_factor,  radio * h/perfect_factor,
		  //100000, -100000);
		  50* length_z,  200* length_z);//(NEW) set up our viewing area
	*/

    // glTranslated(-center_x, -center_y, -center_z - 100 * length_z);
    glTranslated(0, 0 , - camera_displace);

    glMatrixMode(GL_MODELVIEW);
    //keyboard('.',0,0);
    //gluLookAt( -50 * length_z , 70 * length_z ,  70 * length_z, 0, 0,  - 100 * length_z , 0, 1, 0);
    get_window_world_radio();
  }

  void Scene::set_buffer() {
    
    object->xf_buf = object->xf; //for (int i =0;i < 16; ++i) currentScene->object->xf_buf[i] = currentScene->object->xf[i];

    if (context == handsel) {
      context->xf_buf = XForm::identity(); //for (int i =0;i < 16; ++i) currentScene->context->xf_buf[i] = (i%5 ==0 );
      handsel->set_buffer();
    }
  }

  void Scene::restore_buffer() {
    context->xf = context->xf_buf; //for (int i =0;i < 16; ++i) currentScene->context->xf[i] = currentScene->context->xf_buf[i];
    if (context == handsel) 
      currentScene->handsel->restore_buffer();
  }


  static const int ctrl = 1 - 'a';
  void keyboard(unsigned char key, int x, int y) {
    if (key=='q'||key=='Q') exit(0); //quit
    else if (key=='s' + ctrl) {
      currentScene->ss_solver->write("scene.ss");
    }
    else if (key=='o'||key=='O') { // switch between orth and perspective
      orthOrNot = !orthOrNot;
      reshape(width, height);
      glutPostRedisplay();
    }
    else if (key=='w' || key == 'W') { // render object as wire or solid
      wireOrNot = !wireOrNot;
      glutPostRedisplay();
    }
    else if (key=='x' || key == 'X' || key=='y' || key == 'Y' || key=='z' || key == 'Z') { // canonical views
      glGetFloatv( GL_MODELVIEW_MATRIX, CTM);
      /*      GLfloat tx, ty, tz;
      tx = CTM[0]*CTM[12] + CTM[1]*CTM[13] + CTM[2]*CTM[14];
      ty = CTM[4]*CTM[12] + CTM[5]*CTM[13] + CTM[6]*CTM[14];
      tz = CTM[8]*CTM[12] + CTM[9]*CTM[13] + CTM[10]*CTM[14];
      */
      invert(CTM);
      glLoadIdentity();
      
      if (key=='x') glRotatef(90, 0,1,0); 
      else if (key == 'X') glRotatef(-90,0,1,0);      
      else if (key == 'y') glRotatef(90, 1,0,0); 
      else if (key == 'Y') glRotatef(-90,1,0,0);
      else if (key == 'Z') glRotatef(180, 0, 1, 0);
      glTranslated(-CTM[12],-CTM[13],-CTM[14]);
      glutPostRedisplay();
    }
    else if (key == '.') { // focus to object view
      glGetFloatv( GL_MODELVIEW_MATRIX, CTM);
      Point c = currentScene->object->xf * currentScene->cursor;
      CTM = rot_only(CTM) * XForm::trans(-c[0],-c[1],-c[2]);
      glLoadIdentity();
      glMultMatrixf(CTM);
      
      glutPostRedisplay();
    }
    else if(key == 'e') {
      if(currentScene->render_mode==0) {
	currentScene->render_mode = 1;
	glClearColor(1.0,1.0,1.0,0.0);
	glutPostRedisplay();
      }
      else
      {
	currentScene->render_mode = 0;
	glClearColor(0.0,0.0,0.0,0.0);
	glutPostRedisplay();
      }
    }
    else if (key == 'l') {//record animate

      if(record_switch==0) {
	//currentScene->animator = new Animator();
	record_switch = 1;
	if(currentScene->handsel->animator!=NULL){
	  currentScene->handsel->animator->clear();
	  delete currentScene->handsel->animator;
	  currentScene->handsel->animator = NULL;
	}
	  
	currentScene->handsel->record();
	//currentScene->set_animator(&record_animate);
	//glutIdleFunc(animate);
	//glutPostRedisplay();
      }
      else {
	char text[256];
	printf("Please Enter File Name: ");
	char cc = getchar();
	if(cc!='\n') {
	  int cci=0;
	  text[cci++] = cc;
	  while(1) {
	    cc = getchar();
	    if(cc=='\n'){
	      text[cci] = '\0';
	      break;
	    }
	    text[cci++] = cc;
	  }
	  currentScene->handsel->write_record(text);
	}
	else
	  printf("\n Cancel Writing...");
	//char fix[32] = "_scene.anim";
	//strcat(text,fix);
	//currentScene->animator->write(text);
	//animator->clear();
	//delete currentScene->animator;
	//currentScene->animator = NULL;
	printf("\n");
	record_switch = 0;
      }
    }
    else if ( key == 'L') {//play animate
      printf("Please Enter Animator File Name: ");
      char filepath[256];
      scanf("%s",filepath);
      currentScene->animator = new Animator();//change to a vector
      bool res = currentScene->animator->read(filepath);
      if(res==false)
	return;
      if(!currentScene->animator->constraint_points_list.empty()) {
	currentScene->ss_solver->prepare(currentScene->animator->constraints, currentScene->animator->constraint_points_list[0]);
	 
      }

      currentScene->set_animator(&play_animate);
      glutIdleFunc(animate);
      glutPostRedisplay();
    }
    if (key == 9) {// TAB key, switch between selection mode and normal mode
      if (current_state & LOCK_MODE_SELECT ) { // return to normal mode
	current_state = 0;
	wireOrNot = false;
	//currentScene->context = currentScene->object;
	currentScene->context = context_buffer;
      } else { // set select mode
	current_state = LOCK_MODE_SELECT;
	wireOrNot = true;
	context_buffer = currentScene->context;
	currentScene->context = currentScene->vertsel;
      }
      glutPostRedisplay();
    }

    if (!(current_state & 0xf0)) {//Normal mode      
      if (key == 'd') {
	currentScene->handsel->delete_selected();
	currentScene->context = currentScene->object;
	glutPostRedisplay();
      }

      if (key == 'P') {
	currentScene->handsel->set_solver(currentScene->ss_solver);
	glutPostRedisplay();
      } else if (key == 'S') {
	currentScene->cursor = currentScene->handsel->set_focus();
	//	currentScene->handsel->toggle_dump(currentScene->ss_solver);
	glutPostRedisplay();
      }
#ifdef _SS_SHOW_DEBUG
      else if (key == '0') {
	currentScene->ss_solver->show_debug();
      }
#endif


      if (key == 'g' && !(current_state & ~LOCK_OBJECT_TRANSLATE)) {      
	if (glutGetModifiers() == GLUT_ACTIVE_ALT) {
	  XForm &xf = currentScene->context->xf;
	  xf[12] = currentScene->cursor[0];
	  xf[13] = currentScene->cursor[1];
	  xf[14] = currentScene->cursor[2];
	  //MatxTranslate(xf, xf, -currentScene->cursor[0], -currentScene->cursor[1], -currentScene->cursor[2]);
	  xf = xf * XForm::trans( -currentScene->cursor[0], -currentScene->cursor[1], -currentScene->cursor[2]);
	  glutPostRedisplay();
	}else {
	  current_state |= LOCK_OBJECT_TRANSLATE;	
	  currentScene->set_buffer();
	  glPushMatrix();
	  glMultMatrixf(currentScene->object->xf);
      
	  GLdouble modelview[16], projection[16];

	  glGetDoublev(GL_PROJECTION_MATRIX, projection);
	  glGetDoublev(GL_MODELVIEW_MATRIX, modelview);
	  gluUnProject(x, viewport[3] - y, depth, modelview, projection, viewport, &origin_x, &origin_y, &origin_z);
	  glPopMatrix();

	  
	}
      }
      else if (key == 'r' && !(current_state & ~LOCK_OBJECT_ROTATE)) {
	if (glutGetModifiers() == GLUT_ACTIVE_ALT) {
	  XForm &xf = currentScene->context->xf;
	  MatxTranslate(xf, xf, currentScene->cursor[0], currentScene->cursor[1], currentScene->cursor[2]);
	  for (int i=0; i < 12; ++i) xf[i] = (i%5==0);
	  MatxTranslate(xf, xf, -currentScene->cursor[0], -currentScene->cursor[1], -currentScene->cursor[2]);
	  glutPostRedisplay(); return;
	}
      
	if (current_state & LOCK_OBJECT_ROTATE){
	  object_rotate_switch =!object_rotate_switch;
	  currentScene->restore_buffer();
	  glutPostRedisplay();
	} else {
	  currentScene->set_buffer();
	  current_state |= LOCK_OBJECT_ROTATE;
	}

	if (object_rotate_switch) {
	  glPushMatrix();
	  glMultMatrixf(currentScene->object->xf);
	  GLdouble modelview[16], projection[16];

	  d2_x = x; d2_y = y;

	  glGetDoublev(GL_PROJECTION_MATRIX, projection);
	  glGetDoublev(GL_MODELVIEW_MATRIX, modelview);
	  glPopMatrix();
	  gluProject(currentScene->cursor[0], currentScene->cursor[1], currentScene->cursor[2], modelview, projection, viewport, &origin_x, &origin_y, &origin_z); 


	} else {
	  glPushMatrix();
	  glMultMatrixf(currentScene->object->xf);
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
      else if (key == 'a') {
	((VertSelect*) currentScene->context)->toggle_selected();	
	glutPostRedisplay();
      }
      else if (key == 'A') {
	current_state |= LOCK_MODE_SPEC;
	spec_info = "[r]: Add rigid transformer [h]: Add contraint handler";
	glutPostRedisplay();
      }
      else if (key == 'r' && (current_state & LOCK_MODE_SPEC)) {
	current_state &= ~LOCK_MODE_SPEC;
	spec_info = "";
	// add rigid transformer
	currentScene->handsel->add_rigid(((VertSelect*) currentScene->context)->selected);
	glutPostRedisplay();
      }
      else if (key == 'h' && (current_state & LOCK_MODE_SPEC)) {
	current_state &= ~LOCK_MODE_SPEC;
	spec_info = "";
	// add linear constraint handler
	currentScene->handsel->add_constraint(((VertSelect*) currentScene->context)->selected);
	glutPostRedisplay();
      }
    }
  }


  void skeyboard(int key, int x, int y) {
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
    else if (key == GLUT_KEY_F1) {
      if (glutGetModifiers() == GLUT_ACTIVE_SHIFT) currentScene->read("scene.1");
      else currentScene->write("scene.1");
    }
    else if (key == GLUT_KEY_F2) {
      if (glutGetModifiers() == GLUT_ACTIVE_SHIFT) currentScene->read("scene.2");
      else currentScene->write("scene.2");
    }
    else if (key == GLUT_KEY_F3) {
      if (glutGetModifiers() == GLUT_ACTIVE_SHIFT) currentScene->read("scene.3");
      else currentScene->write("scene.3");
    }
    else if (key == GLUT_KEY_F4) {
      if (glutGetModifiers() == GLUT_ACTIVE_SHIFT) currentScene->read("scene.4");
      else currentScene->write("scene.4");
    }
    else if (key == GLUT_KEY_F5) {
      if (glutGetModifiers() == GLUT_ACTIVE_SHIFT) currentScene->read("scene.5");
      else currentScene->write("scene.5");
    }
    else if (key == GLUT_KEY_F6) {
      if (glutGetModifiers() == GLUT_ACTIVE_SHIFT) currentScene->read("scene.6");
      else currentScene->write("scene.6");
    }
    else if (key == GLUT_KEY_F7) {
      if (glutGetModifiers() == GLUT_ACTIVE_SHIFT) currentScene->read("scene.7");
      else currentScene->write("scene.7");
    }
    else if (key == GLUT_KEY_F8) {
      if (glutGetModifiers() == GLUT_ACTIVE_SHIFT) currentScene->read("scene.8");
      else currentScene->write("scene.8");
    }
    else if (key == GLUT_KEY_F9) {
      if (glutGetModifiers() == GLUT_ACTIVE_SHIFT) currentScene->read("scene.9");
      else currentScene->write("scene.9");
    }
    else if (key == GLUT_KEY_F10) {
      if (glutGetModifiers() == GLUT_ACTIVE_SHIFT) currentScene->read("scene.10");
      else currentScene->write("scene.10");
    }
    else if (key == GLUT_KEY_F11) {
      if (glutGetModifiers() == GLUT_ACTIVE_SHIFT) currentScene->read("scene.11");
      else currentScene->write("scene.11");
    }
    else if (key == GLUT_KEY_F12) {
      if (glutGetModifiers() == GLUT_ACTIVE_SHIFT) currentScene->read("scene.12");
      else currentScene->write("scene.12");
    }
  }


  void motion(int x, int y) {    
    if (current_state & LOCK_VIEW_TRANSLATE) {
      
      glLoadIdentity();
      glTranslated((x-origin_x)*win_world_radio, - (y-origin_y)*win_world_radio, 0);
      glMultMatrixf(CTM);
      glutPostRedisplay();
    }

    else if (current_state & LOCK_VIEW_ROTATE) {
      glLoadIdentity();
      
      glRotatef(360.0 * (x-origin_x)/width/_SS_PI, 0.0, 1.0, 0.0);
      glRotatef(360.0 * (y-origin_y)/height/_SS_PI, 1.0, 0.0, 0.0);

      glMultMatrixf(CTM);	

      glutPostRedisplay();
    }
    else if (current_state & LOCK_BACK_BUFFER_SELECT) {      
      transform_x = x; transform_y =y;
      glutPostRedisplay();
    }


  }



  void pmotion(int x, int y) {
    if (current_state & LOCK_OBJECT_TRANSLATE) {
      glPushMatrix();
      glMultMatrixf(currentScene->object->xf_buf);

      GLdouble modelview[16], projection[16], tx, ty, tz;

      glGetDoublev(GL_PROJECTION_MATRIX, projection);
      glGetDoublev(GL_MODELVIEW_MATRIX, modelview);
      glPopMatrix();
      gluUnProject(x, viewport[3] - y, depth, modelview, projection, viewport, &tx, &ty, &tz);
      tx -=origin_x; ty-=origin_y; tz-=origin_z;
      MatxTranslate(currentScene->context->xf, currentScene->context->xf_buf, tx, ty, tz);

      glutPostRedisplay();
    }    
    else if (current_state & LOCK_OBJECT_ROTATE) {
      transform_x = x; transform_y = y;
      if (object_rotate_switch) {
	glPushMatrix();
	glMultMatrixf(currentScene->object->xf_buf);
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

	t = 2*_SS_PI * std::sqrt(((x-d2_x)*(x-d2_x) + (y-d2_y)*(y-d2_y))/(viewport[2]*viewport[2]+viewport[3]*viewport[3])) ;
	GLdouble sin,cos;
	sin = std::sin(t); cos = std::cos(t);
	MatxRotate(currentScene->context->xf, currentScene->context->xf_buf, axis_x, axis_y, axis_z, sin, cos, currentScene->cursor[0], currentScene->cursor[1], currentScene->cursor[2]);

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
	MatxRotate(currentScene->context->xf, currentScene->context->xf_buf, axis_x, axis_y, axis_z, sin, cos, currentScene->cursor[0], currentScene->cursor[1], currentScene->cursor[2]);

      }
      
      glutPostRedisplay();
      
    }

    if(currentScene->context == currentScene->handsel &&
       (current_state & (LOCK_OBJECT_ROTATE|LOCK_OBJECT_TRANSLATE))) {
      currentScene->handsel->update();
      if(record_switch==1)
	currentScene->handsel->record();
    }

  }

// compatibility with original GLUT

#if !defined(GLUT_WHEEL_UP)
#  define GLUT_WHEEL_UP   3
#  define GLUT_WHEEL_DOWN 4
#endif

#define scale_coeff       (1.1)

  void mouse(int button, int state, int x, int y) {
    switch(button) {
    case GLUT_MIDDLE_BUTTON:
      if (state == GLUT_DOWN) { //&& !(current_state & ~LOCK_VIEW_TRANSLATE)) {
	if (current_state & LOCK_BACK_BUFFER_SELECT) {
	  origin_x = x; origin_y = y;
	} else if (!(current_state & (LOCK_OBJECT_TRANSLATE | LOCK_OBJECT_ROTATE))){
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
	if(currentScene->context == currentScene->handsel &&
	   (current_state & (LOCK_OBJECT_ROTATE|LOCK_OBJECT_TRANSLATE))) {
	  currentScene->handsel->update(true);
	}

	if (current_state & LOCK_OBJECT_TRANSLATE) { 
	  current_state &= ~LOCK_OBJECT_TRANSLATE; 
	  glutPostRedisplay(); 
	}
	else if (current_state & LOCK_OBJECT_ROTATE) { 
	  current_state &= ~LOCK_OBJECT_ROTATE; 
	  object_rotate_switch = false; 
	  glutPostRedisplay();
	}
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
	  currentScene->restore_buffer();
	}
	if (current_state & LOCK_OBJECT_ROTATE) {
	  current_state &= ~LOCK_OBJECT_ROTATE;
	  object_rotate_switch = false;
	  currentScene->restore_buffer();
	}
	if (current_state & LOCK_MODE_SELECT) {
	  if (glutGetModifiers() == GLUT_ACTIVE_SHIFT) 
	    ((VertSelect*) currentScene->context)->register_selected(x-10, viewport[3]-y-10, 21, 21, false, true);	  
	  else 
	    ((VertSelect*) currentScene->context)->register_selected(x-10, viewport[3]-y-10, 21, 21, true, true);	  
	} else {
	  bool check_if_selected;
	  if (glutGetModifiers() == GLUT_ACTIVE_SHIFT) 
	    check_if_selected = currentScene->handsel->register_selected(x, viewport[3]-y, false);	  
	  else 
	    check_if_selected = currentScene->handsel->register_selected(x, viewport[3]-y, true);	  
	  if (check_if_selected) currentScene->context = currentScene->handsel;
	  else currentScene->context = currentScene->object;
	}

	glutPostRedisplay();
      }
      else if (state == GLUT_UP) {
      }

      break;
    case GLUT_WHEEL_UP:
      if (state == GLUT_UP){
	//if (scale * scale_coeff < 6) scale *= scale_coeff;	
	glMatrixMode(GL_PROJECTION);
	glScalef(scale_coeff, scale_coeff, scale_coeff);
	glMatrixMode(GL_MODELVIEW);
	get_window_world_radio();
	glutPostRedisplay();      
      }
      break;
    case GLUT_WHEEL_DOWN:
      if (state == GLUT_DOWN){	
	glMatrixMode(GL_PROJECTION);
	glScalef(1/scale_coeff, 1/scale_coeff, 1/scale_coeff);
	glMatrixMode(GL_MODELVIEW);
	get_window_world_radio();
	glutPostRedisplay();      
      }
      break;
    default:
      break;
    }

  }

  void Scene::bind(Object* obj) {
    context = obj; object = obj; 
    vertsel = new VertSelect(obj); 
    handsel = new HandlerSelect(obj);

    cursor = context->center;
    
    glClearColor(0, 0, 0, 0);
    glClearDepth(1.0);


    glViewport( 0, 0, width, height );



    add_lights(0);
    glEnable(GL_LIGHT0);
    glEnable(GL_LIGHT1);

    camera_displace = 3 * object->size;
    glEnableClientState(GL_COLOR_ARRAY);
    /*
    glNewList(object->LIST_NAME, GL_COMPILE_AND_EXECUTE);
    object->register_mesh();
    object->draw();
    glEndList();
    */
    reshape(width, height); keyboard('.',0,0);
    ground_wire = 5 * ( (int) (100000 * height * win_world_radio) / 50) * 0.00001;    

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

  void dump_image(std::string filename)
  {
    printf("Saving image %s... ", filename.c_str());
    FILE *f = fopen(filename.c_str(), "wb");

    // Read pixels
    //GLUI_Master.auto_set_viewport();
    GLint V[4];
    glGetIntegerv(GL_VIEWPORT, V);
    GLint width = V[2], height = V[3];
    char *buf = new char[width*height*3];
    glPixelStorei(GL_PACK_ALIGNMENT, 1);
    glReadPixels(V[0], V[1], width, height, GL_RGB, GL_UNSIGNED_BYTE, buf);

    // Flip top-to-bottom
    for (int i = 0; i < height/2; i++) {
      char *row1 = buf + 3 * width * i;
      char *row2 = buf + 3 * width * (height - 1 - i);
      for (int j = 0; j < 3 * width; j++)
	std::swap(row1[j], row2[j]);
    }

    // Write out file
    fprintf(f, "P6\n%d %d\n255\n", width, height);
    fwrite(buf, width*height*3, 1, f);
    fclose(f);
    delete [] buf;

    printf("Done.\n");
  }


  void Scene::read(std::string filename) {
    
    object->xf.read(filename + ".obj.xf");
    //object->mesh = TriangleMesh::read(std::string(filename + ".obj.off").c_str());
    glMatrixMode(GL_PROJECTION);
    CTM.read(filename + ".prj.xf");
    glLoadIdentity();
    glMultMatrixf(CTM);

    glMatrixMode(GL_MODELVIEW);    
    CTM.read(filename + ".mod.xf");
    glLoadIdentity();
    glMultMatrixf(CTM);
    get_window_world_radio();
    std::cout << "Import from " << filename << std::endl;
    glutPostRedisplay();

    
  }
  void Scene::write(std::string filename) {
    std::cout << "Press Enter to confirm ... "; char check = getchar();
    if (check != '\n') {
      std::cout << "canceled" << std::endl; 
      while (getchar()!='\n');
      currentScene->render_mode = 0;
      glutPostRedisplay();
      glClearColor(0.0,0.0,0.0,0.0);
      return;}

    glGetFloatv( GL_PROJECTION_MATRIX, CTM); 
    CTM.write(filename + ".prj.xf");

    XForm xf_rtsc; // XForm to be loaded by rtsc NPR rendering

    glGetFloatv( GL_MODELVIEW_MATRIX, CTM); xf_rtsc = xf_rtsc * CTM;
    CTM.write(filename + ".mod.xf");   

    object->xf.write(filename + ".obj.xf"); xf_rtsc = xf_rtsc * object->xf;

    std::string mesh_export = filename + ".obj.off";
    object->mesh->write_deformed_surface_mesh(mesh_export.c_str());
    
    xf_rtsc[14] -= camera_displace;
    xf_rtsc.write(filename + ".obj.rtsc.xf"); 

    dump_image(filename + ".cap.ppm");
    std::cout << "Export to " << filename << std::endl;
    currentScene->render_mode = 0;
    glutPostRedisplay();
    glClearColor(0.0,0.0,0.0,0.0);
  }


}



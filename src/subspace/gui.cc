#include "subspace/gui.hh"
#include "subspace/mesh.hh"

namespace subspace {

  Scene* Scene::currentScene;

  GLfloat scale = 1., win_world_radio, depth;
  int origin_x, origin_y;
  int transform_x=0., transform_y=0.;


  bool left_button = false , right_button= false, 
    orthOrNot=false,
    wireOrNot=false;

  GLfloat CTM[16], PTM[16];

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

    center = Point(0,0,0);
      

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
    
    //glEnable(GL_LINE_SMOOTH);
    glEnable(GL_LINE_SMOOTH);
    glEnable(GL_POLYGON_SMOOTH);
    glHint(GL_LINE_SMOOTH_HINT, GL_DONT_CARE);
    glHint(GL_POLYGON_SMOOTH_HINT, GL_DONT_CARE);
    glLineWidth(1.0);
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

	glPushMatrix();//push i-th matrix
	//glCallList(currentScene->object->LIST_NAME);
	currentScene->object->draw();
		

      
	//currentMeshViewer->Painters[i]->prepare();    
	//currentMeshViewer->Painters[i]->draw();
	glPushMatrix();//pop i-th matrix

    glutSwapBuffers();

    
      
  }


  void Scene::get_window_world_radio() {
    GLdouble corner1[3], corner2[3], depth;
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

    GLfloat center_x = (object->bbox[0] + object->bbox[1]) /2.;
    GLfloat center_y = (object->bbox[2] + object->bbox[3]) /2.;
    GLfloat center_z = (object->bbox[4] + object->bbox[5]) /2.;
    GLfloat length_z = (object->bbox[5] - object->bbox[4]) /2.;


    //std::cout << center_x << " " << center_y << " " << center_z << " " << length_z << std::endl;

    glViewport( 0, 0, width, height );

    glMatrixMode(GL_PROJECTION);

    glLoadIdentity();
    /*
    glOrtho( center_x - radio * width/perfect_factor, center_x + radio * width/perfect_factor,
	     center_y - radio * height/perfect_factor, center_y + radio * height/perfect_factor,
	     //100000, -100000);
	     1.0 ,  center_z + 1000* length_z);//(NEW) set up our viewing area

    glTranslated(0,0, -center_z - 100* length_z);
    */
    gluPerspective(30, (float)width/(float)height, 50*length_z, 200*length_z);
    /*
    glFrustum(  - radio * width/perfect_factor,   radio * width/perfect_factor,
	        - radio * height/perfect_factor,  radio * height/perfect_factor,
	       //100000, -100000);
		50* length_z,  200* length_z);//(NEW) set up our viewing area
    */
    glTranslated(-center_x, -center_y, -center_z - 100 * length_z);


    


    glMatrixMode(GL_MODELVIEW);

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
    
    glutReshapeFunc(reshape);
    glutKeyboardFunc(keyboard);
    glutMotionFunc(motion);
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




    GLfloat center_x = (currentScene->object->bbox[0] + currentScene->object->bbox[1]) /2.;
    GLfloat center_y = (currentScene->object->bbox[2] + currentScene->object->bbox[3]) /2.;
    GLfloat center_z = (currentScene->object->bbox[4] + currentScene->object->bbox[5]) /2.;
    GLfloat length_z = (currentScene->object->bbox[5] - currentScene->object->bbox[4]) /2.;





    glViewport( 0, 0, w, h );

    glMatrixMode(GL_PROJECTION);

    glLoadIdentity();
    
    if (orthOrNot) {
      glOrtho(  - win_world_radio *  w/2,   + win_world_radio * w/2,
		- win_world_radio *  h/2,  +  win_world_radio * h/2,
		//100000, -100000);
		50*length_z,   200* length_z);//(NEW) set up our viewing area
    }

    else      
      gluPerspective(30*scale, (float)w/(float)h, 50*length_z, 200*length_z);
	/*
      glFrustum(  - radio * w/perfect_factor,   radio * w/perfect_factor,
		  - radio * h/perfect_factor,  radio * h/perfect_factor,
		  //100000, -100000);
		  50* length_z,  200* length_z);//(NEW) set up our viewing area
	*/

    glTranslated(-center_x, -center_y, -center_z - 100 * length_z);



    glMatrixMode(GL_MODELVIEW);
    currentScene->get_window_world_radio();

  }



  void Scene::keyboard(unsigned char key, int x, int y) {
    if (key=='q'||key=='Q') exit(0);
    else if (key=='o'||key=='O') {
      orthOrNot = !orthOrNot;
      reshape(currentScene->width, currentScene->height);
      glutPostRedisplay();
    }
    else if (key=='w' || key == 'W') {
      wireOrNot = !wireOrNot;
      glutPostRedisplay();
    }
  }


  void Scene::motion(int x, int y) {
    if (left_button) {
      
      glLoadIdentity();
      glTranslated((x-origin_x)*win_world_radio, - (y-origin_y)*win_world_radio, 0);
      glMultMatrixf(CTM);
      glutPostRedisplay();
    }

    else if (right_button) {
      glLoadIdentity();
      
      glRotatef(360.0 * (x-origin_x)/currentScene->width/3.14159265, 0.0, 1.0, 0.0);
      glRotatef(360.0 * (y-origin_y)/currentScene->height/3.14159265, 1.0, 0.0, 0.0);

      glMultMatrixf(CTM);	

      glutPostRedisplay();
    }

  }



// compatibility with original GLUT

#if !defined(GLUT_WHEEL_UP)
#  define GLUT_WHEEL_UP   3
#  define GLUT_WHEEL_DOWN 4
#endif

#define scale_coeff       (1.2)

  void Scene::mouse(int button, int state, int x, int y) {
    switch(button) {
    case GLUT_LEFT_BUTTON:
      if (state == GLUT_DOWN) {	
	origin_x = x; origin_y = y; 
	left_button =true;
	glGetFloatv( GL_MODELVIEW_MATRIX, CTM);
      }
      else if (state == GLUT_UP) {
	left_button = false;
	glutPostRedisplay();
      }
      break;
    case GLUT_MIDDLE_BUTTON:
      if (state == GLUT_DOWN) {
	glLoadIdentity();
	scale =1.;
	reshape(currentScene->width, currentScene->height);
	glutPostRedisplay();
      }
      break;
    case GLUT_RIGHT_BUTTON:
      if (state == GLUT_DOWN) {
	origin_x = x; origin_y = y; right_button =true;
	glGetFloatv(GL_MODELVIEW_MATRIX, CTM);
      }
      else if (state == GLUT_UP) {
	right_button = false;
	glutPostRedisplay();
      }

      break;
    case GLUT_WHEEL_UP:
      if (state == GLUT_UP){	
	scale *= scale_coeff;
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



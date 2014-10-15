#include "Voronoi.h"

void display()
{
	glClear(GL_COLOR_BUFFER_BIT);
	glColor3f(0.0, 0.0, 0.0);
	glPointSize(5.0);

	MESH mesh;
	CreatePionts(&mesh);
	IncrementalDelaunay(&mesh);

	//printDelaunay(&mesh);
	printVoronoi(&mesh);
	glFlush();
}
void init()
{
	glClearColor(1.0,1.0,1.0,0.0);
	glColor3f(0.0,0.0,0.0);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(200.0,600.0,200,400);
}
void myKeyboard(unsigned char key,int x,int y)
{
	if(key==27)exit(0);
}



// FUNCTIONS //////////////////////////////////////////////
int main(int argc, char** argv)
{
	glutInit(&argc,argv);
	glutInitDisplayMode(GLUT_SINGLE|GLUT_RGB);
	glutInitWindowSize(600,500);
	glutInitWindowPosition(600,100);
	glutCreateWindow("άŵͼ");
	glutDisplayFunc(display);
	glutKeyboardFunc(myKeyboard);
	init();
	glutMainLoop();

	return 0;
}
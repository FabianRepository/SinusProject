/*
Copyright (c) 2017, Fabian Prada
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of
conditions and the following disclaimer. Redistributions in binary form must reproduce
the above copyright notice, this list of conditions and the following disclaimer
in the documentation and/or other materials provided with the distribution.

Neither the name of the Johns Hopkins University nor the names of its contributors
may be used to endorse or promote products derived from this software without specific
prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO THE IMPLIED WARRANTIES
OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
TO, PROCUREMENT OF SUBSTITUTE  GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.
*/

#include <Util/Timer.h>
#include <Util/Visualization.h>
#include <Util/Image.h>
#include <Util/Camera.h>
#include <Eigen/Dense>
#ifndef M_PI
#define M_PI		3.14159265358979323846
#endif // M_PI

class PathSelectionVisualization : public Visualization{
public:
	PathSelectionVisualization();

	int diskX, diskY;
	bool showDisk;
	bool isBrushActive;
	bool brushAdd;
	float brushRadius;

	std::unordered_set<int> roiTriangles;

	Point3D<double> pathStart;
	Point3D<double> pathEnd;
	std::vector<Point3D<double>> pathCenters;
	Point3D<double> closestPoint;
	Point3D<double> selectedPoint;

	GLuint offscreen_depth_texture = 0;
	GLuint offscreen_color_texture = 0;
	GLuint offscreen_framebuffer_handle = 0;
	int offscreen_frame_width, offscreen_frame_height;
	void SetupOffScreenBuffer();
	void RenderOffScreenBuffer(Image<Point3D<float>> & image);

	double nearPlane;
	double farPlane;

	Camera worldCamera;
	float zoom;
	int oldX, oldY, newX, newY;

	Eigen::Vector3f centroid;
	float scale;
	std::vector< Eigen::Vector3f> colors;
	std::vector< Eigen::Vector3f> normals;
	std::vector< Eigen::Vector3f> vertices;
	std::vector< Eigen::Vector3f> baricenters;
	std::vector< TriangleIndex > triangles;

	GLfloat lightAmbient[4], lightDiffuse[4], lightSpecular[4], lightPosition[4], shapeSpecular[4], shapeSpecularShininess;
	bool showEdges;
	bool showMesh;
	bool useLight;
	bool rotating, scaling, panning;
	bool showSelection;
	int selectedVertex;
	int selectedTriangle;

	GLuint vertexBuffer = 0;
	GLuint normalBuffer = 0;
	GLuint colorBuffer = 0;
	GLuint faceBuffer = 0;

	void UpdateVertexBuffer();
	void UpdateFaceBuffer();


	//void idle(void);
	void keyboardFunc(unsigned char key, int x, int y);
	//void specialFunc(int key, int x, int y);
	void display(void);
	void mouseFunc(int button, int state, int x, int y);
	void motionFunc(int x, int y);

	static void ShowMeshCallBack(Visualization* v, const char*){
		PathSelectionVisualization* av = (PathSelectionVisualization*)v;
		av->showMesh = !av->showMesh;
	}

	static void ShowEdgesCallBack(Visualization* v, const char*){
		PathSelectionVisualization* av = (PathSelectionVisualization*)v;
		av->showEdges = !av->showEdges;
	}

	static void ScreenshotCallBack(Visualization* v, const char* prompt);

	bool select(int x, int  y, Eigen::Vector3f& out);
	void UpdateSelectedVertex(int x, int y);

	void UpdateSelectedTriangle(int x, int y);
};


void PathSelectionVisualization::SetupOffScreenBuffer(){
	// The depth buffer texture
	glGenTextures(1, &offscreen_depth_texture);
	glBindTexture(GL_TEXTURE_2D, offscreen_depth_texture);
	glTexStorage2D(GL_TEXTURE_2D, 1, GL_DEPTH_COMPONENT24, offscreen_frame_width, offscreen_frame_height);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);

	// The color buffer texture
	glGenTextures(1, &offscreen_color_texture);
	glBindTexture(GL_TEXTURE_2D, offscreen_color_texture);
	glTexStorage2D(GL_TEXTURE_2D, 1, GL_RGBA8, offscreen_frame_width, offscreen_frame_height);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);


	// Create and set up the FBO
	glGenFramebuffers(1, &offscreen_framebuffer_handle);
	glBindFramebuffer(GL_FRAMEBUFFER, offscreen_framebuffer_handle);
	glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, offscreen_depth_texture, 0);
	glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, offscreen_color_texture, 0);
	GLenum drawBuffers[] = { GL_COLOR_ATTACHMENT0 };
	glDrawBuffers(1, drawBuffers);

	GLenum result = glCheckFramebufferStatus(GL_FRAMEBUFFER);
	if (result == GL_FRAMEBUFFER_COMPLETE) {
		printf("Framebuffer is complete.\n");
	}
	else {
		printf("Framebuffer is not complete.\n");
	}
	glBindFramebuffer(GL_FRAMEBUFFER, 0);
}

void PathSelectionVisualization::RenderOffScreenBuffer(Image<Point3D<float>> & image){
	if (!offscreen_framebuffer_handle) SetupOffScreenBuffer();
	glViewport(0, 0, offscreen_frame_width, offscreen_frame_height);
	glBindFramebuffer(GL_FRAMEBUFFER, offscreen_framebuffer_handle);
	int windowScreenWidth = screenWidth;
	int windowScreenHeight = screenHeight;
	screenWidth = offscreen_frame_width;
	screenHeight = offscreen_frame_height;
	display();
	screenWidth = windowScreenWidth;
	screenHeight = windowScreenHeight;
	glFlush();

	//Save color buffer to image
	Pointer(float) GLColorBuffer = AllocPointer< float >(sizeof(float)* 3 * offscreen_frame_width * offscreen_frame_height);
	glReadBuffer(GL_COLOR_ATTACHMENT0);
	glReadPixels(0, 0, offscreen_frame_width, offscreen_frame_height, GL_RGB, GL_FLOAT, GLColorBuffer);
	glFinish();
	//colorBuffer.resize(offscreen_frame_width, offscreen_frame_height);
	//image.init(hh2::Ar(offscreen_frame_height, offscreen_frame_width));
	image.resize(offscreen_frame_width, offscreen_frame_height);
	for (int i = 0; i<offscreen_frame_width; i++) for (int j = 0; j<offscreen_frame_height; j++)  for (int c = 0; c<3; c++){
		image(i, j)[c] = GLColorBuffer[c + i * 3 + (offscreen_frame_height - 1 - j) * offscreen_frame_width * 3];
	}
	FreePointer(GLColorBuffer);
	glBindFramebuffer(GL_FRAMEBUFFER, 0);
	glViewport(0, 0, screenWidth, screenHeight);
}


void PathSelectionVisualization::ScreenshotCallBack(Visualization* v, const char* prompt){
	Image<Point3D<float>> image;
	PathSelectionVisualization* av = (PathSelectionVisualization*)v;
	av->RenderOffScreenBuffer(image);
	image.write(prompt);
}


void PathSelectionVisualization::UpdateVertexBuffer(){
	if (!glIsBuffer(vertexBuffer)){
		glGenBuffers(1, &vertexBuffer);
	}
	if (!glIsBuffer(normalBuffer)){
		glGenBuffers(1, &normalBuffer);
	}
	if (!glIsBuffer(colorBuffer)){
		glGenBuffers(1, &colorBuffer);
	}

	int vCount = vertices.size();

	glBindBuffer(GL_ARRAY_BUFFER, vertexBuffer);
	glBufferData(GL_ARRAY_BUFFER, vCount * sizeof(Eigen::Vector3f), &vertices[0], GL_DYNAMIC_DRAW);
	glBindBuffer(GL_ARRAY_BUFFER, normalBuffer);
	glBufferData(GL_ARRAY_BUFFER, vCount * sizeof(Eigen::Vector3f), &normals[0], GL_DYNAMIC_DRAW);
	glBindBuffer(GL_ARRAY_BUFFER, colorBuffer);
	glBufferData(GL_ARRAY_BUFFER, vCount * sizeof(Eigen::Vector3f), &colors[0], GL_DYNAMIC_DRAW);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
}

void PathSelectionVisualization::UpdateFaceBuffer(){
	if (!glIsBuffer(faceBuffer)){
		glGenBuffers(1, &faceBuffer);
	}
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, faceBuffer);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, triangles.size() * sizeof(int)* 3, &triangles[0][0], GL_STATIC_DRAW);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
}

void PathSelectionVisualization::display(){
	if (!vertexBuffer || !normalBuffer || !colorBuffer) UpdateVertexBuffer();
	if (!faceBuffer) UpdateFaceBuffer();

	glClearColor(1, 1, 1, 1);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	//glDisable(GL_CULL_FACE);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	float ar = (float)screenWidth / (float)screenHeight, ar_r = 1.f / ar;
	//if (screenWidth > screenHeight) glOrtho(-ar*zoom, ar*zoom, -zoom, zoom, -10.5, 10.5);
	//else                           glOrtho(-zoom, zoom, -ar_r*zoom, ar_r*zoom, -10.5, 10.5);
	//if (screenWidth > screenHeight) gluPerspective(30, ar, 0.01f, 100.f);
	//else                           glOrtho(-zoom, zoom, -ar_r*zoom, ar_r*zoom, -10.5, 10.5);
	gluPerspective(30, ar, nearPlane, farPlane);
	//gluPerspective(60.f, (float)screenWidth / (float)screenHeight, nearPlane, farPlane);
	glMatrixMode(GL_MODELVIEW);

	glLoadIdentity();


	Point3D<double> camera_position = worldCamera.position;
	Point3D<double> camera_right = worldCamera.right;
	Point3D<double> camera_forward = worldCamera.forward;
	Point3D<double> camera_up = worldCamera.up;

	//Draw Camera
	//camera.draw();
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	//if (useToolCamera){
	//	Point3D<double> position = toolObject->GetWorldPosition(toolCamera.position);
	//	Point3D<double> forward = toolObject->GetWorldDirection(toolCamera.forward);
	//	Point3D<double> up = toolObject->GetWorldDirection(toolCamera.up);
	//	gluLookAt(position[0], position[1], position[2],
	//		position[0] + forward[0], position[1] + forward[1], position[2] + forward[2],
	//		up[0], up[1], up[2]
	//		);
	//}
	//else{
	//	gluLookAt(
	//		worldCamera.position[0], worldCamera.position[1], worldCamera.position[2],
	//		worldCamera.position[0] + worldCamera.forward[0], worldCamera.position[1] + worldCamera.forward[1], worldCamera.position[2] + worldCamera.forward[2],
	//		worldCamera.up[0], worldCamera.up[1], worldCamera.up[2]
	//		);
	//}

	gluLookAt(
		camera_position[0], camera_position[1], camera_position[2],
		camera_position[0] + camera_forward[0], camera_position[1] + camera_forward[1], camera_position[2] + camera_forward[2],
		camera_up[0], camera_up[1], camera_up[2]
		);

	glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_FALSE);
	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
	glLightfv(GL_LIGHT0, GL_AMBIENT, lightAmbient);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, lightDiffuse);
	glLightfv(GL_LIGHT0, GL_SPECULAR, lightSpecular);

	lightPosition[0] = -camera_forward[0];
	lightPosition[1] = -camera_forward[1];
	lightPosition[2] = -camera_forward[2];
	lightPosition[3] = 0.f;

	glLightfv(GL_LIGHT0, GL_POSITION, lightPosition);
	glEnable(GL_LIGHT0);
	if (useLight) glEnable(GL_LIGHTING);
	else           glDisable(GL_LIGHTING);
	glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
	glEnable(GL_COLOR_MATERIAL);

	glEnable(GL_DEPTH_TEST);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, shapeSpecular);
	glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, shapeSpecularShininess);

	glEnable(GL_DEPTH_TEST);

	glBindBuffer(GL_ARRAY_BUFFER, vertexBuffer);
	glEnableClientState(GL_VERTEX_ARRAY);
	glVertexPointer(3, GL_FLOAT, 0, NULL);

	glBindBuffer(GL_ARRAY_BUFFER, normalBuffer);
	glEnableClientState(GL_NORMAL_ARRAY);
	glNormalPointer(GL_FLOAT, 0, NULL);

	glBindBuffer(GL_ARRAY_BUFFER, colorBuffer);
	glEnableClientState(GL_COLOR_ARRAY);
	glColorPointer(3, GL_FLOAT, 0, NULL);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, faceBuffer);

	if (showMesh) glDrawElements(GL_TRIANGLES, (GLsizei)(triangles.size() * 3), GL_UNSIGNED_INT, NULL);
	glDisableClientState(GL_COLOR_ARRAY);


	if (showEdges)
	{
		GLint src, dst;
		glGetIntegerv(GL_BLEND_SRC, &src);
		glGetIntegerv(GL_BLEND_DST, &dst);
		Point3D<float> f = camera_forward / 256;
		glPushMatrix();
		glTranslatef(-f[0], -f[1], -f[2]);
		glColor3f(0.125, 0.125, 0.125);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		glEnable(GL_BLEND);
		glEnable(GL_LINE_SMOOTH);
		//glLineWidth(0.25f);
		glLineWidth(0.75f);
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		glDrawElements(GL_TRIANGLES, (GLsizei)(triangles.size() * 3), GL_UNSIGNED_INT, NULL);
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		glDisable(GL_LINE_SMOOTH);
		glPopMatrix();
		glDisable(GL_BLEND);
		glBlendFunc(src, dst);
	}

	glDisableClientState(GL_VERTEX_ARRAY);
	glDisableClientState(GL_NORMAL_ARRAY);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

	if (pathStart[0] != -1 || pathStart[1] != -1 || pathStart[2] != -1){
		glColor3f(0.0, 0.f, 1.f);
		GLUquadric* quad = gluNewQuadric();
		glPushMatrix();
		glTranslatef(pathStart[0], pathStart[1], pathStart[2]);
		gluSphere(quad, 0.01f, 20, 20);
		glPopMatrix();
		gluDeleteQuadric(quad);
	}

	if (pathEnd[0] != -1 || pathEnd[1] != -1 || pathEnd[2] != -1){
		glColor3f(1.0, 0.f, 0.f);
		GLUquadric* quad = gluNewQuadric();
		glPushMatrix();
		glTranslatef(pathEnd[0], pathEnd[1], pathEnd[2]);
		gluSphere(quad, 0.01f, 20, 20);
		glPopMatrix();
		gluDeleteQuadric(quad);
	}

	if (pathCenters.size()){
		glLineWidth(3);
		glColor3f(0.f, 0.f, 0.f);
		glEnable(GL_LINE_SMOOTH);
		glBegin(GL_LINE_STRIP);
		for (int i = 0; i < pathCenters.size(); i++) glVertex3f(pathCenters[i][0], pathCenters[i][1], pathCenters[i][2]);
		glEnd();
	}

	if (closestPoint[0] != -1 || closestPoint[1] != -1 || closestPoint[2] != -1){
		glColor3f(0.0, 1.f, 0.f);
		GLUquadric* quad = gluNewQuadric();
		glPushMatrix();
		glTranslatef(closestPoint[0], closestPoint[1], closestPoint[2]);
		gluSphere(quad, 0.01f, 20, 20);
		glPopMatrix();

		glPushMatrix();
		glTranslatef(selectedPoint[0], selectedPoint[1], selectedPoint[2]);
		gluSphere(quad, 0.01f, 20, 20);
		glPopMatrix();

		glEnable(GL_LINE_SMOOTH);
		glBegin(GL_LINE_STRIP);
		glVertex3f(closestPoint[0], closestPoint[1], closestPoint[2]);
		glVertex3f(selectedPoint[0], selectedPoint[1], selectedPoint[2]);
		glEnd();

		gluDeleteQuadric(quad);
	}

	if (roiTriangles.size()){
		Point3D<float> f = camera_forward / 256;
		glPushMatrix();
		glTranslatef(-f[0], -f[1], -f[2]);

		glColor3f(1.0, 0.2, 0.2);
		glBegin(GL_TRIANGLES);
		for (auto iter = roiTriangles.begin(); iter != roiTriangles.end(); iter++){
			Eigen::Vector3f normal = (vertices[triangles[*iter][1]] - vertices[triangles[*iter][0]]).cross(vertices[triangles[*iter][2]] - vertices[triangles[*iter][0]]);
			normal /= normal.norm();
			glNormal3f(normal[0], normal[1], normal[2]);
			for (int k = 0; k < 3; k++) glVertex3f(vertices[triangles[*iter][k]][0], vertices[triangles[*iter][k]][1], vertices[triangles[*iter][k]][2]);
		}
		glEnd();

		glPopMatrix();
	}

	if (showDisk && isBrushActive){
		glDisable(GL_DEPTH_TEST);
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		glOrtho(0, screenWidth, 0, screenHeight, -1, 1);

		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();

		if (brushAdd) glColor3f(0, 1.0, 0);
		else glColor3f(1.0, 0.0, 0);

		GLUquadric* quad = gluNewQuadric();
		glTranslatef(diskX, screenHeight - diskY, 0);
		gluDisk(quad, 18, 22, 40, 3);
		gluDeleteQuadric(quad);
		glEnable(GL_DEPTH_TEST);

		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		float ar = (float)screenWidth / (float)screenHeight, ar_r = 1.f / ar;
		gluPerspective(30, ar, nearPlane, farPlane);
	}
}


void PathSelectionVisualization::keyboardFunc(unsigned char key, int x, int y)
{
	switch (key)
	{
		//case 'q': worldCamera.rotateUp(-M_PI / 128); break;
		//case 'Q': worldCamera.rotateUp(-M_PI / 4); break;
		//case 'w': worldCamera.rotateUp(M_PI / 128); break;
		//case 'W': worldCamera.rotateUp(M_PI / 4); break;
		//case 'a': worldCamera.rotateRight(-M_PI / 128); break;
		//case 'A': worldCamera.rotateRight(-M_PI / 4); break;
		//case 'z': worldCamera.rotateRight(M_PI / 128); break;
		//case 'Z': worldCamera.rotateRight(M_PI / 4); break;
		//case 's': worldCamera.rotateForward(-M_PI / 128); break;
		//case 'S': worldCamera.rotateForward(-M_PI / 4); break;
		//case 'x': worldCamera.rotateForward(M_PI / 128); break;
		//case 'X': worldCamera.rotateForward(M_PI / 4); break;
	}
}

void PathSelectionVisualization::mouseFunc(int button, int state, int x, int y)
{
	newX = x; newY = y;

	rotating = scaling = panning = false;
	if (button == GLUT_LEFT_BUTTON)
	if (glutGetModifiers() & GLUT_ACTIVE_CTRL) panning = true;
	else                                        rotating = true;
	else if (button == GLUT_RIGHT_BUTTON) scaling = true;
}

void PathSelectionVisualization::motionFunc(int x, int y)
{
	oldX = newX, oldY = newY, newX = x, newY = y;
	int screenSize = std::min< int >(screenWidth, screenHeight);
	float rel_x = (newX - oldX) / (float)screenSize * 2;
	float rel_y = (newY - oldY) / (float)screenSize * 2;
	float pRight = -rel_x * zoom, pUp = rel_y * zoom;
	float pForward = rel_y * zoom;
	float rRight = rel_y, rUp = rel_x;
	if (rotating) worldCamera.rotateUp(rUp), worldCamera.rotateRight(rRight);
	else if (scaling) worldCamera.translate(worldCamera.forward *pForward);
	else if (panning) worldCamera.translate(worldCamera.right * pRight + worldCamera.up * pUp);
	glutPostRedisplay();
}

PathSelectionVisualization::PathSelectionVisualization(){

	showDisk = true;
	isBrushActive = false;
	brushAdd = true;
	brushRadius = 0.025;

	closestPoint = Point3D<double>(-1, -1, -1);
	selectedPoint = Point3D<double>(-1, -1, -1);
	pathStart = Point3D<double>(-1, -1, -1);
	pathEnd = Point3D<double>(-1, -1, -1);

	showMesh = true;
	useLight = true;
	selectedVertex = -1;
	selectedTriangle = -1; 
	showSelection = true;
	showEdges = false;
	screenHeight = 600;
	screenWidth = 600;
	offscreen_frame_width = 1200;
	offscreen_frame_height = 1200;

	nearPlane = 0.01f;
	farPlane = 100.f;

	zoom = 1.f;
	rotating = scaling = panning = false;

	worldCamera.position = Point3D<float>(0.f, 0.f, -5.f);
	callBacks.push_back(KeyboardCallBack(this, 'e', "show edges", ShowEdgesCallBack));
	callBacks.push_back(KeyboardCallBack(this, 'm', "show mesh", ShowMeshCallBack));
	callBacks.push_back(KeyboardCallBack(this, 'K', "screenshot", "File Name", ScreenshotCallBack));

	lightAmbient[0] = lightAmbient[1] = lightAmbient[2] = 0.25f, lightAmbient[3] = 1.f;
	lightDiffuse[0] = lightDiffuse[1] = lightDiffuse[2] = 0.70f, lightDiffuse[3] = 1.f;
	lightSpecular[0] = lightSpecular[1] = lightSpecular[2] = 1.00f, lightSpecular[3] = 1.f;
	lightPosition[3] = 0.f;
	shapeSpecular[0] = shapeSpecular[1] = shapeSpecular[2] = 0.0f, shapeSpecular[3] = 1.f;
	shapeSpecularShininess = 128;
}
//
////ORTHOGRAPHIC
////bool PathSelectionVisualization::SelectPoint(int x, int  y, Point3D<float>& out){
////	bool ret = false;
////	Pointer(float) depthBuffer = AllocPointer< float >(sizeof(float)* screenWidth * screenHeight);
////	glReadPixels(0, 0, screenWidth, screenHeight, GL_DEPTH_COMPONENT, GL_FLOAT, depthBuffer);
////	float ar = (float)screenWidth / (float)screenHeight;
////	float _screenWidth, _screenHeight;
////	if (screenWidth>screenHeight) _screenWidth = screenWidth * ar, _screenHeight = screenHeight;
////	else                           _screenWidth = screenWidth, _screenHeight = screenHeight / ar;
////	{
////		double _x = (double)x / screenWidth - 0.5, _y = 1. - (double)y / screenHeight - 0.5, _z;
////		if (screenWidth>screenHeight) _x *= zoom*ar, _y *= zoom;
////		else                           _x *= zoom, _y *= zoom / ar;
////		_x *= 2., _y *= 2;
////		int x1 = (int)floor(x), y1 = (int)floor(y), x2 = x1 + 1, y2 = y1 + 1;
////		float dx = x - x1, dy = y - y1;
////		x1 = std::max< int >(0.f, std::min< int >(x1, screenWidth - 1));
////		y1 = std::max< int >(0.f, std::min< int >(y1, screenHeight - 1));
////		x2 = std::max< int >(0.f, std::min< int >(x2, screenWidth - 1));
////		y2 = std::max< int >(0.f, std::min< int >(y2, screenHeight - 1));
////		_z =
////			depthBuffer[(screenHeight - 1 - y1)*screenWidth + x1] * (1.f - dx) * (1.f - dy) +
////			depthBuffer[(screenHeight - 1 - y1)*screenWidth + x2] * (dx)* (1.f - dy) +
////			depthBuffer[(screenHeight - 1 - y2)*screenWidth + x1] * (1.f - dx) * (dy)+
////			depthBuffer[(screenHeight - 1 - y2)*screenWidth + x2] * (dx)* (dy);
////		if (_z<1) out = Point3D<float>(worldCamera.forward * (-1.5 + 3. * _z) + worldCamera.right * _x + worldCamera.up * _y + worldCamera.position), ret = true;
////	}
////	FreePointer(depthBuffer);
////	return ret;
////}
//
////PERSPECTIVE
//bool PathSelectionVisualization::SelectPoint(int x, int  y, Eigen::Vector3f& out)
//{
//	GLfloat projectionMatrix[16];
//	glGetFloatv(GL_PROJECTION_MATRIX, projectionMatrix);
//	if (0) printf("GL Projection Matrix\n");
//	if (0) for (int i = 0; i < 4; i++)printf("%f %f %f %f \n", projectionMatrix[i], projectionMatrix[4 + i], projectionMatrix[8 + i], projectionMatrix[12 + i]);
//	float frustrumWidth = 2.f*nearPlane / projectionMatrix[0];
//	float frustrumHeight = 2.f*nearPlane / projectionMatrix[5];
//	printf("Frustrum Dimensions  %f x %f \n", frustrumWidth, frustrumHeight);
//
//	bool ret = false;
//	Pointer(float) depthBuffer = AllocPointer< float >(sizeof(float)* screenWidth * screenHeight);
//	glReadPixels(0, 0, screenWidth, screenHeight, GL_DEPTH_COMPONENT, GL_FLOAT, depthBuffer);
//	int x1 = (int)floor(x), y1 = (int)floor(y), x2 = x1 + 1, y2 = y1 + 1;
//	float dx = x - x1, dy = y - y1;
//	x1 = std::max< int >(0.f, std::min< int >(x1, screenWidth - 1));
//	y1 = std::max< int >(0.f, std::min< int >(y1, screenHeight - 1));
//	x2 = std::max< int >(0.f, std::min< int >(x2, screenWidth - 1));
//	y2 = std::max< int >(0.f, std::min< int >(y2, screenHeight - 1));
//	float z_depth =
//		depthBuffer[(screenHeight - 1 - y1)*screenWidth + x1] * (1.f - dx) * (1.f - dy) +
//		depthBuffer[(screenHeight - 1 - y1)*screenWidth + x2] * (dx)* (1.f - dy) +
//		depthBuffer[(screenHeight - 1 - y2)*screenWidth + x1] * (1.f - dx) * (dy)+
//		depthBuffer[(screenHeight - 1 - y2)*screenWidth + x2] * (dx)* (dy);
//	if (z_depth < 1.f){
//		float nz = z_depth*2.f - 1.f;
//		float ez = -(2.f*(farPlane*nearPlane) / (farPlane - nearPlane)) / (-(farPlane + nearPlane) / (farPlane - nearPlane) + nz);
//		float ex = (float(x) / float(screenWidth) - 0.5f)*frustrumWidth*(ez / nearPlane);
//		float ey = -(float(y) / float(screenHeight) - 0.5f)*frustrumHeight*(ez / nearPlane);
//		printf("Eye Coordinates %f %f %f \n", ex, ey, ez);
//		Point3D<float> worldCoordinate = worldCamera.position + (worldCamera.right / Point3D<float>::Length(worldCamera.right)) *ex + (worldCamera.up / Point3D<float>::Length(worldCamera.up)) *ey + (worldCamera.forward / Point3D<float>::Length(worldCamera.forward)) *ez;
//		printf("World Coordinates %f %f %f \n", worldCoordinate[0], worldCoordinate[1], worldCoordinate[2]);
//		out = worldCoordinate;
//		ret = true;
//	}
//
//
//	//float ar = (float)screenWidth / (float)screenHeight;
//	//float _screenWidth, _screenHeight;
//	//if (screenWidth>screenHeight) _screenWidth = screenWidth * ar, _screenHeight = screenHeight;
//	//else                           _screenWidth = screenWidth, _screenHeight = screenHeight / ar;
//	//{
//	//	double _x = (double)x / screenWidth - 0.5, _y = 1. - (double)y / screenHeight - 0.5, _z;
//	//	if (screenWidth>screenHeight) _x *= zoom*ar, _y *= zoom;
//	//	else                           _x *= zoom, _y *= zoom / ar;
//	//	_x *= 2., _y *= 2;
//	//	int x1 = (int)floor(x), y1 = (int)floor(y), x2 = x1 + 1, y2 = y1 + 1;
//	//	float dx = x - x1, dy = y - y1;
//	//	x1 = std::max< int >(0.f, std::min< int >(x1, screenWidth - 1));
//	//	y1 = std::max< int >(0.f, std::min< int >(y1, screenHeight - 1));
//	//	x2 = std::max< int >(0.f, std::min< int >(x2, screenWidth - 1));
//	//	y2 = std::max< int >(0.f, std::min< int >(y2, screenHeight - 1));
//	//	_z =
//	//		depthBuffer[(screenHeight - 1 - y1)*screenWidth + x1] * (1.f - dx) * (1.f - dy) +
//	//		depthBuffer[(screenHeight - 1 - y1)*screenWidth + x2] * (dx)* (1.f - dy) +
//	//		depthBuffer[(screenHeight - 1 - y2)*screenWidth + x1] * (1.f - dx) * (dy)+
//	//		depthBuffer[(screenHeight - 1 - y2)*screenWidth + x2] * (dx)* (dy);
//	//	if (_z<1) out = Eigen::Vector3f(worldCamera.forward * (-1.5 + 3. * _z) + worldCamera.right * _x + worldCamera.up * _y + worldCamera.position), ret = true;
//	//}
//	FreePointer(depthBuffer);
//	return ret;
//}
//
//void PathSelectionVisualization::UpdateSelectedTriangle(int x, int y){
//	Point3D<float> point;
//	if (SelectPoint(x, y, point)){
//		int tCount = triangles.size();
//		int minIndex = -1;
//		float minDistance = FLT_MAX;
//		for (int i = 0; i < tCount; i++){
//			Point3D<float> centroid = (vertices[triangles[i][0]] + vertices[triangles[i][1]] + vertices[triangles[i][2]]) / 3.0;
//			float tDistance = Point3D<float>::SquareNorm(centroid - point);
//			if (tDistance < minDistance){
//				minDistance = tDistance;
//				minIndex = i;
//			}
//		}
//		selectedTriangle = minIndex;
//		glutPostRedisplay();
//	}
//	else{
//		selectedTriangle = -1;
//	}
//}



//PathSelectionVisualization::PathSelectionVisualization(void)
//{
//	vbo = 0;
//	ebo = 0;
//	zoom = 1.05f;
//	lightAmbient[0] = lightAmbient[1] = lightAmbient[2] = 0.25f, lightAmbient[3] = 1.f;
//	lightDiffuse[0] = lightDiffuse[1] = lightDiffuse[2] = 0.70f, lightDiffuse[3] = 1.f;
//	lightSpecular[0] = lightSpecular[1] = lightSpecular[2] = 1.00f, lightSpecular[3] = 1.f;
//	shapeSpecular[0] = shapeSpecular[1] = shapeSpecular[2] = 1.00f, shapeSpecular[3] = 1.f;
//	shapeSpecularShininess = 128;
//	oldX, oldY, newX, newY;
//	rotating = scaling = panning = false;
//	useLight = false;
//	showEdges = true;
//	showMesh = true;
//	selectedTriangle = -1;
//}
//
//
void PathSelectionVisualization::UpdateSelectedVertex(int x, int y){
	Eigen::Vector3f point;
	if (select(x, y, point)){
		int vCount = vertices.size();
		int minIndex = -1;
		float minDistance = FLT_MAX;
		for (int i = 0; i < vCount; i++){
			float tDistance = (vertices[i] - point).norm();
			if (tDistance < minDistance){
				minDistance = tDistance;
				minIndex = i;
			}
		}
		selectedVertex = minIndex;
		glutPostRedisplay();
	}
	else{
		selectedVertex = -1;
	}
}

void PathSelectionVisualization::UpdateSelectedTriangle(int x, int y){
	Eigen::Vector3f point;
	if (select(x, y, point)){
		int tCount = triangles.size();
		int minIndex = -1;
		float minDistance = FLT_MAX;
		for (int i = 0; i < tCount; i++){
			float tDistance = (baricenters[i]- point).norm();
			if (tDistance < minDistance){
					minDistance = tDistance;
					minIndex = i;
			}
		}
		selectedTriangle = minIndex;
		glutPostRedisplay();
	}
	else{
		selectedTriangle = -1;
	}
}

//bool PathSelectionVisualization::select(int x, int  y, Eigen::Vector3f& out) //ORTHOGRAPHIC
//{
//	bool ret = false;
//	Pointer(float) depthBuffer = AllocPointer< float >(sizeof(float)* screenWidth * screenHeight);
//	glReadPixels(0, 0, screenWidth, screenHeight, GL_DEPTH_COMPONENT, GL_FLOAT, depthBuffer);
//	float ar = (float)screenWidth / (float)screenHeight;
//	float _screenWidth, _screenHeight;
//	if (screenWidth>screenHeight) _screenWidth = screenWidth * ar, _screenHeight = screenHeight;
//	else                           _screenWidth = screenWidth, _screenHeight = screenHeight / ar;
//	{
//		double _x = (double)x / screenWidth - 0.5, _y = 1. - (double)y / screenHeight - 0.5, _z;
//		if (screenWidth>screenHeight) _x *= zoom*ar, _y *= zoom;
//		else                           _x *= zoom, _y *= zoom / ar;
//		_x *= 2., _y *= 2;
//		int x1 = (int)floor(x), y1 = (int)floor(y), x2 = x1 + 1, y2 = y1 + 1;
//		float dx = x - x1, dy = y - y1;
//		x1 = std::max< int >(0.f, std::min< int >(x1, screenWidth - 1));
//		y1 = std::max< int >(0.f, std::min< int >(y1, screenHeight - 1));
//		x2 = std::max< int >(0.f, std::min< int >(x2, screenWidth - 1));
//		y2 = std::max< int >(0.f, std::min< int >(y2, screenHeight - 1));
//		_z =
//			depthBuffer[(screenHeight - 1 - y1)*screenWidth + x1] * (1.f - dx) * (1.f - dy) +
//			depthBuffer[(screenHeight - 1 - y1)*screenWidth + x2] * (dx)* (1.f - dy) +
//			depthBuffer[(screenHeight - 1 - y2)*screenWidth + x1] * (1.f - dx) * (dy)+
//			depthBuffer[(screenHeight - 1 - y2)*screenWidth + x2] * (dx)* (dy);
//		if (_z<1) out = Eigen::Vector3f(worldCamera.forward * (-1.5 + 3. * _z) + worldCamera.right * _x + worldCamera.up * _y + worldCamera.position), ret = true;
//	}
//	FreePointer(depthBuffer);
//	return ret;
//}

bool PathSelectionVisualization::select(int x, int  y, Eigen::Vector3f& out)
{
	GLfloat projectionMatrix[16];
	glGetFloatv(GL_PROJECTION_MATRIX, projectionMatrix);
	if (0) printf("GL Projection Matrix\n");
	if (0) for (int i = 0; i < 4; i++)printf("%f %f %f %f \n", projectionMatrix[i], projectionMatrix[4 + i], projectionMatrix[8 + i], projectionMatrix[12 + i]);
	float frustrumWidth = 2.f*nearPlane / projectionMatrix[0];
	float frustrumHeight = 2.f*nearPlane / projectionMatrix[5];
	if (0) printf("Frustrum Dimensions  %f x %f \n", frustrumWidth, frustrumHeight);

	bool ret = false;
	Pointer(float) depthBuffer = AllocPointer< float >(sizeof(float)* screenWidth * screenHeight);
	glReadPixels(0, 0, screenWidth, screenHeight, GL_DEPTH_COMPONENT, GL_FLOAT, depthBuffer);
	int x1 = (int)floor(x), y1 = (int)floor(y), x2 = x1 + 1, y2 = y1 + 1;
	float dx = x - x1, dy = y - y1;
	x1 = std::max< int >(0.f, std::min< int >(x1, screenWidth - 1));
	y1 = std::max< int >(0.f, std::min< int >(y1, screenHeight - 1));
	x2 = std::max< int >(0.f, std::min< int >(x2, screenWidth - 1));
	y2 = std::max< int >(0.f, std::min< int >(y2, screenHeight - 1));
	float z_depth =
		depthBuffer[(screenHeight - 1 - y1)*screenWidth + x1] * (1.f - dx) * (1.f - dy) +
		depthBuffer[(screenHeight - 1 - y1)*screenWidth + x2] * (dx)* (1.f - dy) +
		depthBuffer[(screenHeight - 1 - y2)*screenWidth + x1] * (1.f - dx) * (dy)+
		depthBuffer[(screenHeight - 1 - y2)*screenWidth + x2] * (dx)* (dy);
	if (z_depth < 1.f){
		float nz = z_depth*2.f - 1.f;
		float ez = -(2.f*(farPlane*nearPlane) / (farPlane - nearPlane)) / (-(farPlane + nearPlane) / (farPlane - nearPlane) + nz);
		float ex = (float(x) / float(screenWidth) - 0.5f)*frustrumWidth*(ez / nearPlane);
		float ey = -(float(y) / float(screenHeight) - 0.5f)*frustrumHeight*(ez / nearPlane);
		if (0) printf("Eye Coordinates %f %f %f \n", ex, ey, ez);
		Point3D<double> worldCoordinate = worldCamera.position + (worldCamera.right / Point3D<double>::Length(worldCamera.right)) *ex + (worldCamera.up / Point3D<double>::Length(worldCamera.up)) *ey + (worldCamera.forward / Point3D<double>::Length(worldCamera.forward)) *ez;
		if (0) printf("World Coordinates %f %f %f \n", worldCoordinate[0], worldCoordinate[1], worldCoordinate[2]);
		out = Eigen::Vector3f(worldCoordinate[0], worldCoordinate[1], worldCoordinate[2]);
		ret = true;
	}
	FreePointer(depthBuffer);
	return ret;
}
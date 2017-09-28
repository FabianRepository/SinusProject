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


#include <Src/IntersectableObject.h>
#include <Src/VectorIO.h>
#include <Util/Timer.h>
#include <Util/Visualization.h>
#include <Util/Image.h>
#include <Util/Camera.h>
#include <Eigen/Dense>
#include "ParametricCurve.h"

#ifndef M_PI
#define M_PI		3.14159265358979323846
#endif // M_PI

class ToolBuilderVisualization : public Visualization{
public:
	ToolBuilderVisualization();

	float scale;
	std::vector<Eigen::Vector3f> controlNodes;
	std::vector<double> controlRadius;
	std::vector<Eigen::Vector3f> curveSamples;
	std::vector<double> curveRadius;
	std::vector<double> sampleWeights;

	int numSegmentSamples;
	int numOutputNodes;


	int selectedNode;

	GLuint offscreen_depth_texture = 0;
	GLuint offscreen_color_texture = 0;
	GLuint offscreen_framebuffer_handle = 0;
	int offscreen_frame_width, offscreen_frame_height;
	void SetupOffScreenBuffer();
	void RenderOffScreenBuffer(Image<Point3D<float>> & image);

	double nearPlane;
	double farPlane;

	bool showTool;
	bool showControlNodes;
	bool showPlane;
	bool isNodeManipulation;

	Camera worldCamera;
	float zoom;
	int oldX, oldY, newX, newY;
	GLfloat lightAmbient[4], lightDiffuse[4], lightSpecular[4], lightPosition[4], shapeSpecular[4], shapeSpecularShininess;
	bool useLight;
	bool rotating, scaling, panning;

	void keyboardFunc(unsigned char key, int x, int y);
	void display(void);
	void mouseFunc(int button, int state, int x, int y);
	void motionFunc(int x, int y);

	static void ShowToolCallBack(Visualization* v, const char*){
		ToolBuilderVisualization* av = (ToolBuilderVisualization*)v;
		av->showTool = !av->showTool;
	}

	static void ShowControlNodesCallBack(Visualization* v, const char*){
		ToolBuilderVisualization* av = (ToolBuilderVisualization*)v;
		av->showControlNodes = !av->showControlNodes;
	}

	static void ShowPlaneCallBack(Visualization* v, const char*){
		ToolBuilderVisualization* av = (ToolBuilderVisualization*)v;
		av->showPlane = !av->showPlane;
	}

	static void ScreenshotCallBack(Visualization* v, const char* prompt);
	static void ExportCurveCallBack(Visualization* v, const char* prompt);
	bool select(int x, int  y, Eigen::Vector3f& out);
	void UpdateSelectedVertex(int x, int y);

	void UpdateCurve();
};


void ToolBuilderVisualization::SetupOffScreenBuffer(){
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

void ToolBuilderVisualization::RenderOffScreenBuffer(Image<Point3D<float>> & image){
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

void ToolBuilderVisualization::ScreenshotCallBack(Visualization* v, const char* prompt){
	Image<Point3D<float>> image;
	ToolBuilderVisualization* av = (ToolBuilderVisualization*)v;
	av->RenderOffScreenBuffer(image);
	image.write(prompt);
}

void UniformSampling(const std::vector<Eigen::Vector3f> & inputCurve, std::vector<Eigen::Vector3f> & outputCurve, const std::vector<double> & inputRadius, std::vector<double> & outputRadius, int numOutputNodes){
	int numInputNodes = inputCurve.size();
	std::vector<double> cumArcLength(numInputNodes, 0);
	double cumLength = 0;
	for (int i = 1; i < numInputNodes; i++){
		double segmentLength = (inputCurve[i] - inputCurve[i - 1]).norm();
		cumLength += segmentLength;
		cumArcLength[i] = cumLength;
	}
	outputCurve.resize(numOutputNodes);
	outputCurve[0] = inputCurve[0];
	outputCurve[numOutputNodes - 1] = inputCurve[numInputNodes - 1];

	outputRadius.resize(numOutputNodes);
	outputRadius[0] = inputRadius[0];
	outputRadius[numOutputNodes - 1] = inputRadius[numInputNodes - 1];

	int inputIndex = 0;
	for (int i = 1; i < numOutputNodes - 1; i++){
		double samplingPos = double(i) * cumLength / double(numOutputNodes - 1);
		while (samplingPos > cumArcLength[inputIndex]){
			inputIndex++;
		}
		double alpha = (samplingPos - cumArcLength[inputIndex - 1]) / (cumArcLength[inputIndex] - cumArcLength[inputIndex - 1]);
		if (alpha < 0 || alpha > 1){
			printf("Unexpected condition!\n");
			return;
		}
		outputCurve[i] = inputCurve[inputIndex - 1] * (1 - alpha) + inputCurve[inputIndex] * alpha;
		outputRadius[i] = inputRadius[inputIndex - 1] * (1 - alpha) + inputRadius[inputIndex] * alpha;
	}
}

void ToolBuilderVisualization::ExportCurveCallBack(Visualization* v, const char* prompt){
	ToolBuilderVisualization* av = (ToolBuilderVisualization*)v;
	av->UpdateCurve();
	std::vector<Eigen::Vector3f> outputSamples;
	std::vector<double> outputRadius;
	UniformSampling(av->curveSamples, outputSamples, av->curveRadius, outputRadius, av->numOutputNodes);
	
	for (int i = 0; i < outputSamples.size(); i++)outputSamples[i] *= av->scale;
	for (int i = 0; i <outputRadius.size(); i++) outputRadius[i] *= av->scale;

	char toolName[256];
	sprintf(toolName, "%s_positions.bin", prompt);
	WriteVector(outputSamples, toolName);
	sprintf(toolName, "%s_radius.bin", prompt);
	WriteVector(outputRadius, toolName);
}

void ToolBuilderVisualization::UpdateCurve(){
	int numNodes = controlNodes.size();
	int offset = 0;
	for (int i = 0; i < numNodes - 1; i++){
		Eigen::Vector3f prev = i > 0 ? controlNodes[i - 1] : 2.0 * controlNodes[0] - controlNodes[1];
		Eigen::Vector3f mid_prev = controlNodes[i];
		Eigen::Vector3f mid_next = controlNodes[i + 1];
		Eigen::Vector3f next = i < controlNodes.size() - 2 ? controlNodes[i + 2] : 2.0 * controlNodes[numNodes - 1] - controlNodes[numNodes - 2];
		
		double prev_radius = i > 0 ? controlRadius[i - 1] : 2.0 * controlRadius[0] - controlRadius[1];
		double mid_prev_radius = controlRadius[i];
		double mid_next_radius = controlRadius[i + 1];
		double next_radius = i < controlRadius.size() - 2 ? controlRadius[i + 2] : 2.0 * controlRadius[numNodes - 1] - controlRadius[numNodes - 2];

		
		for (int j = 0; j < numSegmentSamples + 1; j++){
			curveSamples[offset] = sampleWeights[4 * j] * prev + sampleWeights[4 * j + 1] * mid_prev + sampleWeights[4 * j + 2] * mid_next + sampleWeights[4 * j + 3] * next;
			curveRadius[offset] = sampleWeights[4 * j] * prev_radius + sampleWeights[4 * j + 1] * mid_prev_radius + sampleWeights[4 * j + 2] * mid_next_radius + sampleWeights[4 * j + 3] * next_radius;
			offset++;
		}
		if (i == numNodes - 2){
			int j = numSegmentSamples + 1;
			curveSamples[offset] = sampleWeights[4 * j] * prev + sampleWeights[4 * j + 1] * mid_prev + sampleWeights[4 * j + 2] * mid_next + sampleWeights[4 * j + 3] * next;
			curveRadius[offset] = sampleWeights[4 * j] * prev_radius + sampleWeights[4 * j + 1] * mid_prev_radius + sampleWeights[4 * j + 2] * mid_next_radius + sampleWeights[4 * j + 3] * next_radius;
		}
	}
}

void ToolBuilderVisualization::display(){
	
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

	if (showPlane){
		glBegin(GL_QUADS);
		glNormal3f(0, 0, 1);
		glColor3f(0.8, 0.8, 0.8);
		glVertex3f(-0.5, -0.5, 0);
		glVertex3f( 0.5, -0.5, 0);
		glVertex3f( 0.5,  0.5, 0);
		glVertex3f(-0.5,  0.5, 0);
		glEnd();
	}
	if (showTool){
		UpdateCurve();
		IntersectableCylindricalCurve curve(curveSamples, curveRadius);
		curve.Draw();
	}
	if (showControlNodes){
		glDisable(GL_DEPTH_TEST);
		GLUquadric* quad = gluNewQuadric();
		glColor3f(1.f, 0.f, 0.f);
		for (int i = 0; i < controlNodes.size(); i++){
			if (i != selectedNode)glColor3f(1.f, 0.f, 0.f);
			else glColor3f(0.f, 0.f, 1.f);
			glPushMatrix();
			glTranslatef(controlNodes[i][0], controlNodes[i][1], controlNodes[i][2]);
			gluSphere(quad, i != selectedNode ? 0.01f : 0.03f, 20, 20);
			glPopMatrix();
		}
		gluDeleteQuadric(quad);
		glEnable(GL_DEPTH_TEST);
	}
}


void ToolBuilderVisualization::keyboardFunc(unsigned char key, int x, int y)
{
	switch (key)
	{
		case '+': 
			if (selectedNode != -1) controlRadius[selectedNode] += 0.001;
			break;
		case '-':
			if (selectedNode != -1) controlRadius[selectedNode] = std::max<double>(0, controlRadius[selectedNode] - 0.001);
			break;
	}
}

void ToolBuilderVisualization::mouseFunc(int button, int state, int x, int y)
{
	newX = x; newY = y;

	rotating = scaling = panning = false;



	if (glutGetModifiers() & GLUT_ACTIVE_SHIFT){
		if (state == GLUT_DOWN) UpdateSelectedVertex(x, y);
		isNodeManipulation = true;
		if (button == GLUT_LEFT_BUTTON){
			panning = true;
		}
		else if (button == GLUT_RIGHT_BUTTON) scaling = true;
	}
	else{
		selectedNode = -1;
		isNodeManipulation = false;
		if (button == GLUT_LEFT_BUTTON){
			if (glutGetModifiers() & GLUT_ACTIVE_CTRL){
				 rotating = true;
			}
			else{
				panning = true;
			}
		}
		else if (button == GLUT_RIGHT_BUTTON) scaling = true;
	}
}

void ToolBuilderVisualization::motionFunc(int x, int y)
{
	oldX = newX, oldY = newY, newX = x, newY = y;
	int screenSize = std::min< int >(screenWidth, screenHeight);
	float rel_x = (newX - oldX) / (float)screenSize * 2;
	float rel_y = (newY - oldY) / (float)screenSize * 2;
	float pRight = -rel_x * zoom, pUp = rel_y * zoom;
	float pForward = rel_y * zoom;
	float rRight = rel_y, rUp = rel_x;

	if (isNodeManipulation){
		if (scaling && selectedNode != -1){
			Point3D<double> displacement = worldCamera.forward *pForward;
			controlNodes[selectedNode] += Eigen::Vector3f(displacement[0], displacement[1], displacement[2]);
		}
		else if (panning && selectedNode != -1){
			Point3D<double> displacement = -worldCamera.right * pRight - worldCamera.up * pUp;
			controlNodes[selectedNode] += Eigen::Vector3f(displacement[0], displacement[1], displacement[2]);
		}
	}
	else{
		if (rotating) worldCamera.rotateUp(rUp), worldCamera.rotateRight(rRight);
		else if (scaling) worldCamera.translate(worldCamera.forward *pForward);
		else if (panning) worldCamera.translate(worldCamera.right * pRight + worldCamera.up * pUp);
	}

	glutPostRedisplay();
}

ToolBuilderVisualization::ToolBuilderVisualization(){
	
	isNodeManipulation = false;
	showPlane = true;
	showTool = true;
	showControlNodes = true;
	useLight = true;
	selectedNode = -1;
	screenHeight = 600;
	screenWidth = 600;
	offscreen_frame_width = 1200;
	offscreen_frame_height = 1200;

	nearPlane = 0.01f;
	farPlane = 100.f;

	zoom = 1.f;
	rotating = scaling = panning = false;

	worldCamera.position = Point3D<float>(0.f, 0.f, -2.f);
	callBacks.push_back(KeyboardCallBack(this, 'e', "export curve", "File Name Prefix", ExportCurveCallBack));
	callBacks.push_back(KeyboardCallBack(this, 'K', "screenshot", "File Name", ScreenshotCallBack));
	//callBacks.push_back(KeyboardCallBack(this, 'c', "control nodes", ShowControlNodesCallBack));
	//callBacks.push_back(KeyboardCallBack(this, 't', " show tool", ShowToolCallBack));

	lightAmbient[0] = lightAmbient[1] = lightAmbient[2] = 0.25f, lightAmbient[3] = 1.f;
	lightDiffuse[0] = lightDiffuse[1] = lightDiffuse[2] = 0.70f, lightDiffuse[3] = 1.f;
	lightSpecular[0] = lightSpecular[1] = lightSpecular[2] = 1.00f, lightSpecular[3] = 1.f;
	lightPosition[3] = 0.f;
	shapeSpecular[0] = shapeSpecular[1] = shapeSpecular[2] = 0.0f, shapeSpecular[3] = 1.f;
	shapeSpecularShininess = 128;
}
//
////ORTHOGRAPHIC
////bool ToolBuilderVisualization::SelectPoint(int x, int  y, Point3D<float>& out){
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
//bool ToolBuilderVisualization::SelectPoint(int x, int  y, Eigen::Vector3f& out)
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
//void ToolBuilderVisualization::UpdateSelectedTriangle(int x, int y){
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



//ToolBuilderVisualization::ToolBuilderVisualization(void)
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
void ToolBuilderVisualization::UpdateSelectedVertex(int x, int y){
	Eigen::Vector3f point;
	if (select(x, y, point)){
		int vCount = controlNodes.size();
		int minIndex = -1;
		float minDistance = FLT_MAX;
		for (int i = 0; i < vCount; i++){
			float tDistance = (controlNodes[i] - point).norm();
			if (tDistance < minDistance){
				minDistance = tDistance;
				minIndex = i;
			}
		}
		selectedNode = minIndex;
		glutPostRedisplay();
	}
	else{
		selectedNode = -1;
	}
}

//bool ToolBuilderVisualization::select(int x, int  y, Eigen::Vector3f& out) //ORTHOGRAPHIC
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

bool ToolBuilderVisualization::select(int x, int  y, Eigen::Vector3f& out)
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

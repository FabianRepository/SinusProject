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

#include "Simulator.h"
#include <stdlib.h>
#include <GL/glew.h>
#include <GL/glut.h>
#include <Util/CmdLineParser.h>
#include <Util/Algebra.h>
#include <Util/Camera.h>
#include <Util/Visualization.h>

cmdLineParameter< char * > Input("input");
cmdLineParameter< char * > Fixed("fixed");
cmdLineParameter< char * > ToolPos("toolPos");
cmdLineParameter< char * > ToolRad("toolRad");
cmdLineParameter<  char* > Pose("pose");
cmdLineParameter<  char* > Path("path");
cmdLineParameter<  char* > Output("output");
cmdLineParameter< float > AxialScale("aScale", 1.0);
cmdLineParameter< float > RadialScale("rScale", 1.0);
cmdLineParameter<  char* > CameraConfig("camera");


cmdLineParameter< float > AdvanceStep("step", 0.33);
cmdLineParameter< float > CollisionSensitivity("sensitivity", 1.33);

cmdLineParameter<  int > CenterToolPeriod("centerP",1); 
cmdLineParameter<  int > EnforceAdvancePeriod("enforceP",3);
cmdLineParameter<  int > AvoidSurfacePeriod("avoidP",1);

cmdLineReadable Simulate("simulation");
cmdLineReadable* params[] =
{
	&Input, &Fixed, &ToolPos, &ToolRad, &Simulate, &CameraConfig, &Pose, &Path, &Output, &AxialScale, &RadialScale, &AdvanceStep, &CollisionSensitivity,
	&CenterToolPeriod, &EnforceAdvancePeriod, &AvoidSurfacePeriod,
	NULL
};

void ShowUsage(const char* ex)
{
	printf("Usage: %s [params]:\n", ex);

	printf("ROI and path parameters: \n");
	printf("\t[--%s <roi model>]\n", Input.name);
	printf("\t[--%s <fixed vertices>]\n", Fixed.name);
	printf("\t[--%s <path>]\n\n", Path.name);
	
	printf("Tool parameters: \n");
	printf("\t[--%s <tool node positions>]\n", ToolPos.name);
	printf("\t[--%s <tool node radius>]\n", ToolRad.name);
	printf("\t[--%s <tool initial pose>]\n", Pose.name);
	printf("\t[--%s <axial scale>][%f]\n", AxialScale.name, AxialScale.value);
	printf("\t[--%s <radial scale>][%f]\n\n", RadialScale.name, RadialScale.value);

	printf("Simulation parameters: \n");
	printf("\t[--%s <advance step size>][%f (mm)]\n", AdvanceStep.name, AdvanceStep.value);
	printf("\t[--%s <collision sensitivity>][%f (mm)]\n", CollisionSensitivity.name, CollisionSensitivity.value);
	printf("\t[--%s <center tool period>][%d]\n", CenterToolPeriod.name, CenterToolPeriod.value);
	printf("\t[--%s <enforce advance period>][%d]\n", EnforceAdvancePeriod.name, EnforceAdvancePeriod.value);
	printf("\t[--%s <avoid surface period>][%d]\n", AvoidSurfacePeriod.name, AvoidSurfacePeriod.value);
	printf("\t[--%s <run simulation>]\n\n", Simulate.name);

	printf("Auxiliar and output parameters: \n");
	printf("\t[--%s <camera  configuration>]\n", CameraConfig.name);
	printf("\t[--%s <output prefix>]\n", Output.name);

}

int main(int argc, char* argv[])
{
	cmdLineParse(argc - 1, argv + 1, params);
	if (!ToolPos.set || !Input.set || !Path.set){ ShowUsage(argv[0]); return EXIT_FAILURE; }

	if (!Simulator::Init(ToolPos.value, ToolRad.set ? ToolRad.value : NULL, Input.value, Path.value, Fixed.set? Fixed.value : NULL, Pose.set? Pose.value : NULL, AxialScale.value, RadialScale.value, AdvanceStep.value, CollisionSensitivity.value)){
		printf("Initialization failed! \n");
		return 0;
	}

	if (CameraConfig.set) Simulator::visualization.ReadSceneConfigurationCallBack((SimulatorVisualization*)&Simulator::visualization, CameraConfig.value);

	if (1){//Sort triangles
		std::set<std::pair<float, int>> triangleDistance;
		for (int i = 0; i < Simulator::visualization.triangles.size(); i++){
			Eigen::Vector3f triangleCentroid = (Simulator::visualization.vertices[Simulator::visualization.triangles[i][0]] +
				Simulator::visualization.vertices[Simulator::visualization.triangles[i][1]] +
				Simulator::visualization.vertices[Simulator::visualization.triangles[i][2]]) / 3.0;
			Point3D<double> cameraPos = Simulator::visualization.worldCamera.position;
			Eigen::Vector3f cameraRay = triangleCentroid - Eigen::Vector3f(cameraPos[0], cameraPos[1], cameraPos[2]);
			Point3D<double> cameraForward = Simulator::visualization.worldCamera.forward;
			float cameraProjection = cameraRay.dot(Eigen::Vector3f(cameraForward[0], cameraForward[1], cameraForward[2]));
			triangleDistance.insert(std::pair<float, int>(cameraProjection, i));
		}
		std::vector<TriangleIndex> sortedTriangles;
		for (auto iter = triangleDistance.begin(); iter != triangleDistance.end(); iter++){
			int tIndex = (*iter).second;
			sortedTriangles.push_back(Simulator::visualization.triangles[tIndex]);
		}
		if (sortedTriangles.size() != Simulator::visualization.triangles.size()){
			printf("ERROR : Sorted triangles does not match input triangles! \n");
			return 0;
		}
		Simulator::visualization.triangles = sortedTriangles;
		printf("Triangles sorted! \n");
	}

	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
	glutInitWindowSize(600, 600);
	glutInit(&argc, argv);
	char windowName[1024];
	sprintf(windowName, "Simulator");
	glutCreateWindow(windowName);

	if (glewInit() != GLEW_OK) fprintf(stderr, "[ERROR] glewInit failed\n"), exit(0);
	glutDisplayFunc(Simulator::Display);
	glutReshapeFunc(Simulator::Reshape);
	glutMouseFunc(Simulator::MouseFunc);
	glutMotionFunc(Simulator::MotionFunc);
	glutKeyboardFunc(Simulator::KeyboardFunc);

	if (Simulate.set){
		std::vector<Eigen::Vector3f> translations;
		std::vector<Eigen::Matrix3f> rotations;
		Simulator::visualization.useTransparency = true;
		Simulator::visualization.showFixedVertices = false;
		Simulator::visualization.showClosestPoints = false;
		Simulator::visualization.showToolCentroid = false;
		Simulator::visualization.showToolSamples = false;
		Simulator::visualization.showTargetDeformationPositions = false;
		Simulator::visualization.showCollision = false;
		Simulator::visualization.showMedialVectorField = false;
		Simulator::visualization.showClosestAndOppositePoints = false;
		Simulator::GeneratePath(Output.set? Output.value: NULL, CenterToolPeriod.value, EnforceAdvancePeriod.value, AvoidSurfacePeriod.value);
	}


	{
		Simulator::visualization.useTransparency = true;
		Simulator::visualization.showFixedVertices = false;
		Simulator::visualization.showMedialVectorField = true;
		Simulator::visualization.showClosestPoints = true;
		Simulator::visualization.showToolCentroid = false;
		Simulator::visualization.showToolSamples = false;
		Simulator::visualization.showTargetDeformationPositions = false;
		Simulator::visualization.showCollision = true;
		Simulator::visualization.showClosestAndOppositePoints = false;
	}

	glutMainLoop();

	return EXIT_SUCCESS;

}
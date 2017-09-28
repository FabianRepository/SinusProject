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

#include <Src/DeformableMesh.h>
#include "SimulatorVisualization.inl"
#include <Src/VectorIO.h>
#include <Src/CurveQuery.h>
#include <Src/ExtrinsicVectorField.h>

#define USE_NORMAL_DISTANCE 0

enum{
	VECTOR_FIELD_ADVANCE,
	RANDOM_ADVANCE,
	MIXED_ADVANCE,
	ADVANCE_MODE_COUNT
};


class PathInfo{
public:
	PathInfo(){
		axialRotationEnergy = 0;
		orthogonalRotationEnergy = 0;
		axialTranslationEnergy = 0;
		orthogonalTranslationEnergy = 0;
		cummulativeDistortionEnergy = 0;
		maxDistortionEnergy = 0;
	}
	float axialRotationEnergy;
	float orthogonalRotationEnergy;
	float axialTranslationEnergy;
	float orthogonalTranslationEnergy;
	float cummulativeDistortionEnergy;
	float maxDistortionEnergy;
	std::vector<Eigen::Vector3f> translations;
	std::vector<Eigen::Matrix3f> rotations;
};

int SavePathTransformations(const PathInfo & pathInfo, const char * prefix){
	char rotationsName[256];
	sprintf(rotationsName, "%s_Rotations.bin", prefix);
	WriteVector(pathInfo.rotations, rotationsName);
	char translationName[256];
	sprintf(translationName, "%s_Translations.bin", prefix);
	WriteVector(pathInfo.translations, translationName);
	return 1;
}

int SavePathStatistics(const PathInfo & pathInfo, const char * prefix){
	char statisticsName[256];
	sprintf(statisticsName, "%s_Statistics.txt", prefix);
	FILE * file = fopen(statisticsName, "w");
	if (!file){
		printf("Unable to open report file \n");
		return 0;
	}

	fprintf(file, "Axial rotation energy %f \n", pathInfo.axialRotationEnergy);
	fprintf(file, "Orthogonal rotation energy %f \n", pathInfo.orthogonalRotationEnergy);
	fprintf(file, "Axial translation energy %f \n", pathInfo.axialTranslationEnergy);
	fprintf(file, "Orthogonal translation energy %f \n", pathInfo.orthogonalTranslationEnergy);
	fprintf(file, "Cummulative distortion %f \n", pathInfo.cummulativeDistortionEnergy);
	fprintf(file, "Maximum distortion  %f\n", pathInfo.maxDistortionEnergy);

	fclose(file);
	return 1;
}

class Simulator
{
public:
	static DeformableMesh elasticObject;

	static IntersectableObject * toolObject;
	static std::vector<Eigen::Vector3f> toolNodes;
	static Eigen::Vector3f toolTip;
	static Eigen::Vector3f toolBottom;

	static std::vector<Eigen::Vector3f> pathNodes;
	static CurveQuery * pathQuery;
	static ExtrinsicVectorField * pathVectorField;

	static SimulatorVisualization visualization;
	static Eigen::Vector3f meshCentroid;
	static float meshScale;

	static std::vector<double> toolAxialMetric;

	static Eigen::Matrix3f axisAlignedFrame;
	static std::vector<Eigen::Vector3f> initialCollisionPosition;
	static std::vector<Eigen::Vector3f> targetCollisionPosition;
	static std::vector<Eigen::Vector3f> temporalTargetPosition;
	static std::vector<double> temporalFittingWeight;

	static float axialRotationWeight;
	static float orthogonalRotationWeight;
	static float axialTranslationWeight;
	static float orthogonalTranslationWeight;

	static float advanceStep;
	static float collisionSensitivity;

	static std::vector<bool> fixedVertices;

	static std::vector<double> toolRadiusVector;
	//static float toolHalfWidth;
	//static float toolLength;
	static float radialScale;
	static float axialScale;

	static int Init(const char * toolName, const char * toolRadiusName, const char* meshName, const char* pathName, const char * fixedVerticesName, const char * toolPoseName, const float _axialScale, const float _radialScale, const float _advanceStep, const float _collisionSensitity);
	static void Display(void){ visualization.Display(); }
	static void MouseFunc(int button, int state, int x, int y);
	static void MotionFunc(int x, int y);
	static void Reshape(int w, int h){ visualization.Reshape(w, h); }
	static void KeyboardFunc(unsigned char key, int x, int y){ visualization.KeyboardFunc(key, x, y); }

	static void ExportDeformedMesh(const char * outputName);
	static void ExportDeformedMeshCallBack(Visualization* v, const char* prompt);

	static void AvoidSurface(const int testingIterations, const double maxTranslationMagnitude, const double maxRotationMagnitude);
	static void	AvoidSurfaceCallBack(Visualization* v, const char*);

	static void CenterTool(const float centeringWeight = 0.5f);
	static void	CenterToolCallBack(Visualization* v, const char*);

	static void AdvanceVectorField( bool orthogonalCorrection = true);
	static void	AdvanceVectorFieldCallBack(Visualization* v, const char*);

	static void IdentifyCollisions(std::vector<int> & collisionIndices);
	static void IdentifyCollisionsCallBack(Visualization* v, const char*);

	static void ComputeDistortion(float & distortion);
	static void ComputeDistortionCallBack(Visualization* v, const char*);

	static void SolveCollision();
	static void SolveCollisionsCallBack(Visualization* v, const char*);

	static void ComputeClosestPoints(Eigen::Vector3f & toolPoint, Eigen::Vector3f & objectPoint);
	static void ComputeClosestPointsCallBack(Visualization* v, const char*);

	static int GeneratePath(const char* outputPrefix, const int centerToolPeriod, const int enforcedAdvancePeriod, const int avoidSurfacePeriod);
	//static double NormalizedDistanceToTarget();

	static double DistanceToTarget();
};

DeformableMesh								Simulator::elasticObject;
IntersectableObject	*						Simulator::toolObject;
SimulatorVisualization						Simulator::visualization;
float										Simulator::advanceStep = 0.33; //mm
float										Simulator::collisionSensitivity = 1.33; //mm
float										Simulator::axialScale = 1.00;
float										Simulator::radialScale = 1.00;

std::vector<double>							Simulator::toolRadiusVector;
CurveQuery *								Simulator::pathQuery;
std::vector<double>							Simulator::toolAxialMetric;
std::vector<Eigen::Vector3f>				Simulator::temporalTargetPosition;
std::vector<Eigen::Vector3f>				Simulator::initialCollisionPosition;
std::vector<Eigen::Vector3f>				Simulator::targetCollisionPosition;
std::vector<double>							Simulator::temporalFittingWeight;
ExtrinsicVectorField *						Simulator::pathVectorField;
Eigen::Vector3f								Simulator::toolTip;
Eigen::Vector3f								Simulator::toolBottom;
Eigen::Matrix3f								Simulator::axisAlignedFrame;
std::vector<Eigen::Vector3f>				Simulator::pathNodes;
std::vector<Eigen::Vector3f>				Simulator::toolNodes;
std::vector<bool>							Simulator::fixedVertices;
float										Simulator::axialRotationWeight = 1.0;
float										Simulator::orthogonalRotationWeight = 1.0;
float										Simulator::axialTranslationWeight = 1.0;
float										Simulator::orthogonalTranslationWeight = 1.0;
Eigen::Vector3f								Simulator::meshCentroid;
float										Simulator::meshScale;

int Simulator::Init(const char * toolName, const char * toolRadiusName, const char* meshName, const char* pathName, const char * fixedVerticesName, const char * toolPoseName, const float _axialScale, const float _radialScale, const float _advanceStep, const float _collisionSensitity){

	advanceStep = _advanceStep;
	collisionSensitivity = _collisionSensitity;
	axialScale = _axialScale;
	radialScale = _radialScale;

	//Setup Deformable mesh
	SimpleMesh mesh;
	if (!ReadSimpleMesh(mesh, meshName)){
		printf("Unable to read mesh!\n");
		return 0;
	}

	GetMeshCentroidAndScale(mesh, meshCentroid, meshScale);

	printf("Centroid %f %f %f \n", meshCentroid[0], meshCentroid[1], meshCentroid[2]);
	printf("Scale %f \n", meshScale);

	for (int i = 0; i < mesh.vertices.size(); i++) mesh.vertices[i] = (mesh.vertices[i] - meshCentroid) / meshScale;

	Eigen::MatrixXf	geodesicDescriptor;
	if (!ReadEigenMatrixf(geodesicDescriptor, "GeodesicDescriptor.vec")){
		printf("Constructing geodesic descriptors ...\n");

		int descriptorDimension = 20;

		std::vector<double> points(3 * mesh.vertices.size());
		for (int i = 0; i < mesh.vertices.size(); i++) for (int c = 0; c < 3; c++) points[3 * i + c] = (double)mesh.vertices[i][c];
		std::vector<unsigned> faces(3 * mesh.triangles.size());
		for (int i = 0; i < mesh.triangles.size(); i++) for (int c = 0; c < 3; c++) faces[3 * i + c] = (unsigned)mesh.triangles[i][c];

		geodesic::Mesh geoMesh;

		geoMesh.initialize_mesh_data(points, faces); //create internal mesh data structure including edges
		geodesic::GeodesicAlgorithmExact * geoAlgorithm = new geodesic::GeodesicAlgorithmExact(&geoMesh);

		GenerateGeodesicDescriptor(mesh, descriptorDimension, geodesicDescriptor, geoMesh, geoAlgorithm);
		printf("Save geodesic descriptors ...\n");
		SaveEigenMatrixf(geodesicDescriptor, "GeodesicDescriptor.vec");
	}
	else{
		printf("Geodesic descriptors succesfully read ...\n");
	}

	
	if (fixedVerticesName != NULL){
		std::vector<int> _fixedVertices;
		if (!ReadVector(_fixedVertices, fixedVerticesName)){
			printf("Unable to read fixed vertices! \n");
			return 0;
		}
		else{
			fixedVertices.resize(_fixedVertices.size());
			for (int i = 0; i < fixedVertices.size(); i++) fixedVertices[i] = bool(_fixedVertices[i]);
		}
	}
	else{
		fixedVertices.resize(mesh.vertices.size(), false);
	}
	visualization.isFixedVertex = fixedVertices;

	elasticObject = DeformableMesh(mesh.vertices, mesh.triangles, geodesicDescriptor, fixedVertices);

	printf("Hierrachy resolution: \n");
	for (int i = 0; i < elasticObject.hierarchy.numLevels; i++){
		printf("%d \n", elasticObject.hierarchy.hierarchyFineIndices[i].size());
	}

	visualization.vertices = mesh.vertices;
	visualization.triangles = mesh.triangles;
	visualization.normals = mesh.normals;
	visualization.colors.resize(visualization.vertices.size(), Eigen::Vector4f(.7, .7, .7,.8));
	visualization.centroid = meshCentroid;
	visualization.scale = meshScale;

	
	if (!ReadVector(pathNodes, pathName)){
		printf("Unable to read tool file! \n");
		return 0;
	}

	printf("Number of path nodes %d \n", pathNodes.size());

	float pathRealLength = (pathNodes[pathNodes.size() - 1] - pathNodes[0]).norm();
	for (int i = 0; i < pathNodes.size(); i++)pathNodes[i] = (pathNodes[i] - meshCentroid) / meshScale;
	
	pathQuery = new CurveQuery(pathNodes);
	pathVectorField = new CurveInducedField(pathNodes);

	visualization.pathNodes = pathNodes;

	//Setup rigid object

	if (!ReadVector(toolNodes, toolName)){
		printf("Unable to read tool file! \n");
		return 0;
	}

	printf("Number of tool nodes %d \n", toolNodes.size());

	Eigen::Vector3f inputToolCentroid = Eigen::Vector3f::Zero();
	for (int i = 0; i < toolNodes.size(); i++) inputToolCentroid += toolNodes[i];
	inputToolCentroid /= float(toolNodes.size());
	float toolScaleFactor = axialScale / meshScale;
	for (int i = 0; i < toolNodes.size(); i++)toolNodes[i] = (toolNodes[i] - inputToolCentroid) * toolScaleFactor;

	advanceStep /= meshScale;
	collisionSensitivity /= meshScale;
	toolRadiusVector.resize(toolNodes.size(), 1.0);
	if (toolRadiusName != NULL){
		if (!ReadVector(toolRadiusVector, toolRadiusName)){
			printf("Unable to read tool radius file! \n");
			return 0;
		}
		
	}

	float toolRadialFactor = radialScale / meshScale;
	for (int i = 0; i < toolRadiusVector.size(); i++) toolRadiusVector[i] *= toolRadialFactor;

	toolObject = new IntersectableCylindricalCurve(toolNodes, toolRadiusVector);
	toolObject->InitializeSamples(100);
	printf("Num samples %d \n", toolObject->samples.size());
	toolBottom = toolNodes[0];
	visualization.toolBottom = toolBottom;
	toolTip = toolNodes[toolNodes.size() - 1];
	visualization.toolTip = toolTip;
	Eigen::Vector3f toolPrincipalAxis = toolNodes[toolNodes.size() - 1] - toolNodes[0];

	toolPrincipalAxis /= toolPrincipalAxis.norm();
	Eigen::Vector3f orthBasis[2];
	orthBasis[0] = Eigen::Vector3f((float(rand()) / float(RAND_MAX)) - 0.5, (float(rand()) / float(RAND_MAX)) - 0.5, (float(rand()) / float(RAND_MAX)) - 0.5);
	orthBasis[0] = toolPrincipalAxis.cross(orthBasis[0]);
	orthBasis[0] /= orthBasis[0].norm();
	orthBasis[1] = toolPrincipalAxis.cross(orthBasis[0]);

	axisAlignedFrame.col(0) = toolPrincipalAxis;
	axisAlignedFrame.col(1) = orthBasis[0];
	axisAlignedFrame.col(2) = orthBasis[1];

	visualization.axisAlignedFrame = axisAlignedFrame;

	toolObject->UpdateTransformations();
	visualization.toolObject = toolObject;

	if (toolPoseName != NULL){
		visualization.ReadToolConfigurationCallBack((SimulatorVisualization*)&visualization, toolPoseName);
	}


	toolAxialMetric.resize(toolNodes.size() - 1);
	double cumMedialAxisLength = 0;
	for (int i = 0; i < toolNodes.size() - 1; i++){
		double segmentSquaredLength = (toolNodes[i + 1] - toolNodes[i]).squaredNorm();
		toolAxialMetric[i] = segmentSquaredLength;
		cumMedialAxisLength += sqrt(segmentSquaredLength);
	}
	for (int i = 0; i < toolNodes.size() - 1; i++){
		toolAxialMetric[i] /= (cumMedialAxisLength*cumMedialAxisLength);
	}
	double cumMass = 0;
	for (int i = 0; i < toolNodes.size() - 1; i++) cumMass += sqrt(toolAxialMetric[i]);
	if (abs(1.0 - cumMass) > 1e-10){
		printf("Precision error: Cum mass does not add up to 1! \n");
		return 0;
	}

	temporalTargetPosition.resize(mesh.vertices.size());
	temporalFittingWeight.resize(mesh.vertices.size(), 0);
	initialCollisionPosition.resize(mesh.vertices.size());
	targetCollisionPosition.resize(mesh.vertices.size());

	visualization.callBacks.push_back(Visualization::KeyboardCallBack(&visualization, 'v', "advance vfield", AdvanceVectorFieldCallBack));
	visualization.callBacks.push_back(Visualization::KeyboardCallBack(&visualization, 'g', "avoid surface", AvoidSurfaceCallBack));
	visualization.callBacks.push_back(Visualization::KeyboardCallBack(&visualization, 'k', "center tool", CenterToolCallBack));
	visualization.callBacks.push_back(Visualization::KeyboardCallBack(&visualization, 'E', "export mesh","File Name",ExportDeformedMeshCallBack));
	

	visualization.callBacks.push_back(Visualization::KeyboardCallBack(&visualization, 'c', "identify collision", IdentifyCollisionsCallBack));
	visualization.callBacks.push_back(Visualization::KeyboardCallBack(&visualization, 's', "solve collision", SolveCollisionsCallBack));
	visualization.callBacks.push_back(Visualization::KeyboardCallBack(&visualization, 'u', "closest points", ComputeClosestPointsCallBack));

	visualization.callBacks.push_back(Visualization::KeyboardCallBack(&visualization, 'd', "distortion", ComputeDistortionCallBack));

	return 1;
}

void Simulator::MouseFunc(int button, int state, int x, int y)
{
	visualization.panning = visualization.scaling = visualization.rotating = false;
	visualization.newX = x; visualization.newY = y;

	if (glutGetModifiers() & GLUT_ACTIVE_SHIFT && state == GLUT_DOWN){
		visualization.isToolActive = true;
	}
	else{
		visualization.isToolActive = false;
	}

	if (glutGetModifiers() & GLUT_ACTIVE_CTRL){
		visualization.scaling = true;
	}
	else{
		if (button == GLUT_LEFT_BUTTON){
			visualization.panning = true;
		}
		else if (button == GLUT_RIGHT_BUTTON){
			visualization.rotating = true;
		}
	}

	glutPostRedisplay();
}
void Simulator::MotionFunc(int x, int y){
	visualization.oldX = visualization.newX, visualization.oldY = visualization.newY, visualization.newX = x, visualization.newY = y;
	int screenSize = std::min< int >(visualization.screenWidth, visualization.screenHeight);
	float rel_x = (visualization.newX - visualization.oldX) / (float)screenSize * 2;
	float rel_y = (visualization.newY - visualization.oldY) / (float)screenSize * 2;
	float pRight = -rel_x * visualization.zoom, pUp = rel_y * visualization.zoom;
	float pForward = rel_y * visualization.zoom;
	float rRight = rel_y, rUp = rel_x;

	if (visualization.isToolActive){
		if (visualization.panning){
			Point3D<double> displacement = visualization.worldCamera.right * pRight + visualization.worldCamera.up * pUp;
			toolObject->objectToWorldTranslation += Eigen::Vector3f(displacement[0], displacement[1], displacement[2]);
		}
		else if (visualization.scaling){
			Point3D<double> displacement = visualization.worldCamera.forward *pForward;
			toolObject->objectToWorldTranslation += Eigen::Vector3f(displacement[0], displacement[1], displacement[2]);

		}
		else{
			Point3D<double> up = visualization.worldCamera.up;
			Eigen::Vector3f upDirection(up[0], up[1], up[2]);
			float upAngle = rel_x * 20 * M_PI / 180.0;

			Eigen::AngleAxisf upAxisRotation(upAngle, upDirection);
			Eigen::Matrix3f upMatrixRotation = upAxisRotation.matrix();

			Eigen::Vector3f centroid = toolObject->GetCenterOfMass();
			toolObject->objectToWorldRotation = upMatrixRotation * toolObject->objectToWorldRotation;
			toolObject->objectToWorldTranslation = upMatrixRotation * toolObject->objectToWorldTranslation + (Eigen::Matrix3f::Identity() - upMatrixRotation)*centroid;


			Point3D<double> right = visualization.worldCamera.right;
			Eigen::Vector3f rightDirection(right[0], right[1], right[2]);
			float rightAngle = rel_y * 20 * M_PI / 180.0;

			Eigen::AngleAxisf rigthAxisRotation(rightAngle, rightDirection);
			Eigen::Matrix3f rightMatrixRotation = rigthAxisRotation.matrix();

			centroid = toolObject->GetCenterOfMass();
			toolObject->objectToWorldRotation = rightMatrixRotation * toolObject->objectToWorldRotation;
			toolObject->objectToWorldTranslation = rightMatrixRotation * toolObject->objectToWorldTranslation + (Eigen::Matrix3f::Identity() - rightMatrixRotation)*centroid;
		}
		toolObject->UpdateTransformations();
	}
	else{
		if (visualization.panning) visualization.worldCamera.translate(visualization.worldCamera.right * pRight + visualization.worldCamera.up * pUp);
		else if (visualization.scaling)visualization.worldCamera.translate(visualization.worldCamera.forward *pForward);
		else visualization.worldCamera.rotateUp(rUp), visualization.worldCamera.rotateRight(rRight);
	}

	glutPostRedisplay();
}

void Simulator::ComputeClosestPoints(Eigen::Vector3f & toolPoint, Eigen::Vector3f & objectPoint){
	std::vector<Eigen::Vector3f> toolSamples;
	toolObject->GetSamples(toolSamples);

	double minSurfaceDistance = DBL_MAX;

	Eigen::VectorXf sqrD;
	Eigen::VectorXi I;
	Eigen::MatrixXf C;
	Eigen::MatrixXf P;
	P.resize(toolSamples.size(), 3);
	for (int s = 0; s < toolSamples.size(); s++) P.row(s) = toolSamples[s];

	//elasticObject.undeformedKdTree.squared_distance(elasticObject.undeformedVertexMatrix, elasticObject.triangleMatrix, P, sqrD, I, C);
	elasticObject.kdTree.squared_distance(elasticObject.vertexMatrix, elasticObject.triangleMatrix, P, sqrD, I, C);

	for (int s = 0; s < toolSamples.size(); s++){
		Eigen::Vector3f transformedSamplePos = P.row(s);
		Eigen::Vector3f closestSurfacePoint = C.row(s);
#if USE_NORMAL_DISTANCE
		Eigen::Vector3f surfaceNormal = elasticObject.undeformedTriangleNormals[I[s]];
		double normalDistance = surfaceNormal.dot(transformedSamplePos - closestSurfacePoint);
		if (normalDistance < minSurfaceDistance){
			minSurfaceDistance = normalDistance;
			toolPoint = transformedSamplePos;
			objectPoint = closestSurfacePoint;
		}
#else
		double distance = (transformedSamplePos - closestSurfacePoint).norm();
		if (distance < minSurfaceDistance){
			minSurfaceDistance = distance;
			toolPoint = transformedSamplePos;
			objectPoint = closestSurfacePoint;
		}
#endif
	}
}

void Simulator::ComputeClosestPointsCallBack(Visualization* v, const char*){
	Eigen::Vector3f toolPoint;
	Eigen::Vector3f objectPoint;
	ComputeClosestPoints(toolPoint, objectPoint);

	visualization.showClosestPoints = true;
	visualization.closestToolPoint = toolPoint;
	visualization.closestObjectPoint = objectPoint;

	glutPostRedisplay();
}

//Rotations along principal axis
void Simulator::AvoidSurface(const int testingIterations, const double maxTranslationMagnitude, const double maxRotationMagnitude){

	std::vector<Eigen::Vector3f> toolSamples;
	toolObject->GetSamples(toolSamples);

	double initialSurfaceDistance = DBL_MAX;

	{//Compute initial state

		Eigen::VectorXf sqrD;
		Eigen::VectorXi I;
		Eigen::MatrixXf C;
		Eigen::MatrixXf P;
		P.resize(toolSamples.size(), 3);
		for (int s = 0; s < toolSamples.size(); s++) P.row(s) = toolSamples[s];

		elasticObject.kdTree.squared_distance(elasticObject.vertexMatrix, elasticObject.triangleMatrix, P, sqrD, I, C);

		for (int s = 0; s < toolSamples.size(); s++){
			Eigen::Vector3f transformedSamplePos = P.row(s);
			Eigen::Vector3f closestSurfacePoint = C.row(s);
#if USE_NORMAL_DISTANCE
			Eigen::Vector3f surfaceNormal = elasticObject.triangleNormals[I[s]];
			double normalDistance = surfaceNormal.dot(transformedSamplePos - closestSurfacePoint);
			if (normalDistance < initialSurfaceDistance){
				initialSurfaceDistance = normalDistance;
			}
#else 
			double distance = (transformedSamplePos - closestSurfacePoint).norm();
			//double distance_sign = toolObject->isInterior(closestSurfacePoint)? -1.0 : 1.0;
			//distance *= distance_sign;
			if (distance < initialSurfaceDistance){
				initialSurfaceDistance = distance;
			}
#endif
			
		}
	}

	if (1) printf("Initial distance to surface %g \n", initialSurfaceDistance);

	Eigen::Vector3f toolCentroid = toolObject->GetCenterOfMass();
	Eigen::Matrix3f transformedAxisAlignedFrame = toolObject->objectToWorldRotation * axisAlignedFrame;

	clock_t startTimer = clock();

	Eigen::Vector3f optimalTranslation = Eigen::Vector3f::Zero();
	Eigen::Matrix3f optimalRotation = Eigen::Matrix3f::Identity();

	double optimalSurfaceDistance = initialSurfaceDistance;

	int threads = omp_get_num_procs();
#pragma omp parallel for num_threads( threads )
	for (int i = 0; i < testingIterations; i++){

		Eigen::Vector2f translationDirection((float(rand()) / float(RAND_MAX)) - 0.5, (float(rand()) / float(RAND_MAX)) - 0.5);
		translationDirection /= translationDirection.norm();
		float translationMagnitude = ((float(rand()) / float(RAND_MAX))) * maxTranslationMagnitude;

		Eigen::Vector3f translation = (transformedAxisAlignedFrame.col(1)*translationDirection[0] + transformedAxisAlignedFrame.col(2)*translationDirection[1])*translationMagnitude;

		float rotationAngle = (float(rand()) / float(RAND_MAX)) * maxRotationMagnitude *M_PI / 180.0;

		Eigen::AngleAxisf angleAxisRotation(rotationAngle, transformedAxisAlignedFrame.col(0));
		Eigen::Matrix3f rotation = angleAxisRotation.matrix();

		Eigen::VectorXf sqrD;
		Eigen::VectorXi I;
		Eigen::MatrixXf C;
		Eigen::MatrixXf P;
		P.resize(toolSamples.size(), 3);
		for (int s = 0; s < toolSamples.size(); s++) P.row(s) = rotation * (toolSamples[s] - toolCentroid) + toolCentroid + translation;

		//elasticObject.undeformedKdTree.squared_distance(elasticObject.undeformedVertexMatrix, elasticObject.triangleMatrix, P, sqrD, I, C);
		elasticObject.kdTree.squared_distance(elasticObject.vertexMatrix, elasticObject.triangleMatrix, P, sqrD, I, C);

		double surfaceDistance = DBL_MAX;

		for (int s = 0; s < toolSamples.size(); s++){
			Eigen::Vector3f transformedSamplePos = P.row(s);
			Eigen::Vector3f closestSurfacePoint = C.row(s);
#if USE_NORMAL_DISTANCE
			Eigen::Vector3f surfaceNormal = elasticObject.triangleNormals[I[s]];
			double normalDistance = surfaceNormal.dot(transformedSamplePos - closestSurfacePoint);
			if (normalDistance < surfaceDistance){
				surfaceDistance = normalDistance;
			}
#else
			double distance = (transformedSamplePos - closestSurfacePoint).norm();
			//double distance_sign = toolObject->isInterior(closestSurfacePoint) ? -1.0 : 1.0;
			//distance *= distance_sign;
			if (distance < surfaceDistance){
				surfaceDistance = distance;
			}
#endif
			//curveDistance = std::max<double>(curveDistance,pathQuery->GetTriangularDistance(transformedSamplePos));
		}
#pragma omp critical
		{
			if (surfaceDistance > optimalSurfaceDistance){
				optimalSurfaceDistance = surfaceDistance;
				optimalTranslation = translation;
				optimalRotation = rotation;
			}
		}
	}

	if (optimalSurfaceDistance > initialSurfaceDistance){
		printf("Final surface distance %g. Gain %g. \n", optimalSurfaceDistance, optimalSurfaceDistance - initialSurfaceDistance);

		if (0)printf("Advance time %.4f on %d iterations \n", double(clock() - startTimer) / CLOCKS_PER_SEC, testingIterations);
		if (0){
			printf("Translation: \n");
			printf("%f %f %f \n", optimalTranslation[0], optimalTranslation[1], optimalTranslation[2]);
			printf("Rotation: \n");
			for (int k = 0; k < 3; k++) printf("%f %f %f \n", optimalRotation(k, 0), optimalRotation(k, 1), optimalRotation(k, 2));
		}


		Eigen::Matrix3f initialToolRotation = toolObject->objectToWorldRotation;
		Eigen::Vector3f initialToolTranslation = toolObject->objectToWorldTranslation;

		toolObject->objectToWorldRotation = optimalRotation * toolObject->objectToWorldRotation;
		toolObject->objectToWorldTranslation = optimalRotation * (toolObject->objectToWorldTranslation - toolCentroid) + toolCentroid + optimalTranslation;
		toolObject->UpdateTransformations();

		std::vector<int> collisionIndices;
		IdentifyCollisions(collisionIndices);

		if (collisionIndices.size() > 0){
			printf("Motion discarded by collision!\n");
			toolObject->objectToWorldRotation = initialToolRotation;
			toolObject->objectToWorldTranslation = initialToolTranslation;
			toolObject->UpdateTransformations();
		}
	}
}

void Simulator::ExportDeformedMesh(const char * outputName){
	SimpleMesh outMesh = elasticObject.mesh;
	for (int i = 0; i < outMesh.vertices.size(); i++)outMesh.vertices[i] = outMesh.vertices[i] * meshScale + meshCentroid;
	WriteSimpleMesh(outMesh, outputName);
}

void Simulator::ExportDeformedMeshCallBack(Visualization* v, const char* prompt){
	ExportDeformedMesh(prompt);
}

void Simulator::CenterTool(const float centeringWeight){

	std::vector<Eigen::Vector3f> sourcePositions(toolNodes.begin(), toolNodes.end());
	for (int i = 0; i < sourcePositions.size(); i++) sourcePositions[i] = toolObject->objectToWorldRotation*sourcePositions[i] + toolObject->objectToWorldTranslation;
	std::vector<Eigen::Vector3f> targetPositions(sourcePositions.begin(), sourcePositions.end());
	std::vector<float> weights(sourcePositions.size(), (1.0 - centeringWeight) / double(sourcePositions.size()));

	Eigen::Vector3f worldToolTip = toolObject->GetWorldPosition(toolTip);
	Eigen::Vector3f toolTipTarget = pathQuery->GetClosestPoint(worldToolTip);

	sourcePositions.push_back(worldToolTip);
	targetPositions.push_back(toolTipTarget);
	weights.push_back(centeringWeight);

	Eigen::Vector3f optimalTranslation;
	Eigen::Matrix3f optimalRotation;
	RigidAlignment(sourcePositions, targetPositions, weights, optimalTranslation, optimalRotation);

	toolObject->objectToWorldRotation = optimalRotation * toolObject->objectToWorldRotation;
	toolObject->objectToWorldTranslation = optimalRotation * toolObject->objectToWorldTranslation + optimalTranslation;
	toolObject->UpdateTransformations();
}



#define USE_ELASTIC_SHAPE 1

void Simulator::AdvanceVectorField(bool orthogonalCorrection){

	std::vector<Eigen::Vector3f> sourcePositions(toolNodes.begin(), toolNodes.end());

	Eigen::VectorXf sqrD;
	Eigen::VectorXi I;
	Eigen::MatrixXf C;
	Eigen::MatrixXf P;
	P.resize(sourcePositions.size(), 3);
	for (int i = 0; i < sourcePositions.size(); i++){
		sourcePositions[i] = toolObject->objectToWorldRotation*sourcePositions[i] + toolObject->objectToWorldTranslation;
		for (int c = 0; c < 3; c++) P(i, c) = sourcePositions[i][c];
	}


#if USE_ELASTIC_SHAPE
	elasticObject.kdTree.squared_distance(elasticObject.vertexMatrix, elasticObject.triangleMatrix, P, sqrD, I, C);
#else
	elasticObject.undeformedKdTree.squared_distance(elasticObject.undeformedVertexMatrix, elasticObject.triangleMatrix, P, sqrD, I, C);
#endif
	
	visualization.closestPointsSurface.clear();
	visualization.oppositePointsSurface.clear();
	visualization.medialVectorField.resize(sourcePositions.size());
	visualization.medialPositions.resize(sourcePositions.size());

	double maxConstraintWeight = 0;
	double cummulativeConstraintWeight = 0;
	for (int i = 0; i < sourcePositions.size(); i++){

		visualization.medialPositions[i] = sourcePositions[i];
		Eigen::Vector3f medialAxisVectorField = pathVectorField->VectorField(sourcePositions[i]) * advanceStep;
		if (orthogonalCorrection){

			Eigen::Vector3f closestSurfacePoint(C(i, 0), C(i, 1), C(i, 2));
			visualization.closestPointsSurface.push_back(closestSurfacePoint);
			Eigen::Vector3f direction = sourcePositions[i] - closestSurfacePoint;
			double closesSurfacePointDistance = (direction).norm();
			Eigen::Vector3f orthogonalDirection = direction / closesSurfacePointDistance;

			double oppositeSurfacePointDistance = DBL_MAX;
			Eigen::Vector3f isect;

#if USE_ELASTIC_SHAPE
			bool intersects = elasticObject.rayIntersect(sourcePositions[i], orthogonalDirection, isect);
#else
			bool intersects = elasticObject.undeformedRayIntersect(sourcePositions[i], orthogonalDirection, isect);
#endif
			if (intersects){
				visualization.oppositePointsSurface.push_back(isect);
				if (closesSurfacePointDistance < toolRadiusVector[i]){
					printf("WARNING: Closest point within tool radius : %f  %f \n", closesSurfacePointDistance, toolRadiusVector[i]);
				}
				oppositeSurfacePointDistance = (sourcePositions[i] - isect).norm();
			}

			double half_difference = (oppositeSurfacePointDistance - closesSurfacePointDistance) / 2;
			if (half_difference < 0){
				printf("WARNING: half difference is not positive  %f \n", half_difference);
			}

			double orthogonalDisplacementFactor = std::min<double>(half_difference / (2.0*collisionSensitivity), 1.0);
			orthogonalDirection *= 2.0*collisionSensitivity* orthogonalDisplacementFactor;

			double orthogonalWeight;

			if (closesSurfacePointDistance > collisionSensitivity){
				orthogonalWeight = 0;
			}
			else{
				orthogonalWeight = 1.0;
			}
			visualization.medialVectorField[i] = medialAxisVectorField + orthogonalDirection * orthogonalWeight;
		}
		else{
			visualization.medialVectorField[i] = medialAxisVectorField;
		}
	}

	std::vector<double> initialLenghts(sourcePositions.size());
	double maxInitialLength = 0;
	for (int i = 0; i < sourcePositions.size(); i++){
		initialLenghts[i] = visualization.medialVectorField[i].norm();
		maxInitialLength = std::max<double>(initialLenghts[i], maxInitialLength);
	}

	if (0){ //Smooth vector field
		std::vector<Eigen::Triplet<double>> systemTriplets;
		Eigen::MatrixXd rhs;
		Eigen::MatrixXd x0;
		x0.resize(sourcePositions.size(), 3);
		rhs.resize(sourcePositions.size(), 3);
		for (int i = 0; i < sourcePositions.size(); i++)for (int c = 0; c < 3; c++)rhs(i, c) = 0;
		for (int i = 0; i < sourcePositions.size(); i++)for (int c = 0; c < 3; c++)x0(i, c) = visualization.medialVectorField[i][c];

		double stiffnessWeight = 1e-1;

		for (int i = 0; i < sourcePositions.size() - 1; i++){
			double metric = toolAxialMetric[i];
			double invMetric = 1.0 / metric;
			double mass = sqrt(metric);
			double stiffnessFactor = invMetric*mass*stiffnessWeight;
			double constraintWeight = 1.0;
			for (int k = 0; k < 2; k++){
				for (int l = 0; l < 2; l++){
					systemTriplets.push_back(Eigen::Triplet<double>(i + k, i + l, k == l ? stiffnessFactor : -stiffnessFactor));
					systemTriplets.push_back(Eigen::Triplet<double>(i + k, i + l, k == l ? mass*constraintWeight / 3.0 : mass*constraintWeight / 6.0));
					Eigen::Vector3f _rhs = visualization.medialVectorField[i + l] * mass*constraintWeight*(k == l ? 1.0 / 3.0 : 1.0 / 6.0);
					rhs(i + k, 0) += _rhs[0]; rhs(i + k, 1) += _rhs[1]; rhs(i + k, 2) += _rhs[2];
				}
			}
		}

		Eigen::SparseMatrix< double > systemMatrix;
		systemMatrix.resize(sourcePositions.size(), sourcePositions.size());
		systemMatrix.setFromTriplets(systemTriplets.begin(), systemTriplets.end());

		Eigen::MatrixXd Mx0 = systemMatrix * x0;

		Eigen::SimplicialLDLT< Eigen::SparseMatrix< double > > solver(systemMatrix);

		Eigen::MatrixXd solution = solver.solve(rhs);

		for (int i = 0; i < sourcePositions.size(); i++){
			visualization.medialVectorField[i] = Eigen::Vector3f(solution(i, 0), solution(i, 1), solution(i, 2));
		}

		std::vector<double> finalLenghts(sourcePositions.size());
		for (int i = 0; i < sourcePositions.size(); i++)finalLenghts[i] = visualization.medialVectorField[i].norm();

		for (int i = 0; i < sourcePositions.size(); i++) visualization.medialVectorField[i] *= (initialLenghts[i] / finalLenghts[i]);
	}

	for (int i = 0; i < sourcePositions.size(); i++) visualization.medialVectorField[i] *= (advanceStep / maxInitialLength);

	std::vector<float> weights(sourcePositions.size());
	std::vector<Eigen::Vector3f> targetPositions(sourcePositions.begin(), sourcePositions.end());
	for (int i = 0; i < targetPositions.size(); i++){
		targetPositions[i] += visualization.medialVectorField[i];
		weights[i] = visualization.medialVectorField[i].norm();
	}

	Eigen::Vector3f optimalTranslation;
	Eigen::Matrix3f optimalRotation;
	RigidAlignment(sourcePositions, targetPositions, weights, optimalTranslation, optimalRotation);

	toolObject->objectToWorldRotation = optimalRotation * toolObject->objectToWorldRotation;
	toolObject->objectToWorldTranslation = optimalRotation * toolObject->objectToWorldTranslation + optimalTranslation;
	toolObject->UpdateTransformations();
}


void Simulator::AvoidSurfaceCallBack(Visualization* v, const char*){
	AvoidSurface(100, collisionSensitivity / 16, 5);
	glutPostRedisplay();
}

void Simulator::AdvanceVectorFieldCallBack(Visualization* v, const char*){
	AdvanceVectorField();
	glutPostRedisplay();
}



void Simulator::CenterToolCallBack(Visualization* v, const char*){
	CenterTool();
	glutPostRedisplay();
}

void Simulator::IdentifyCollisions(std::vector<int> & collisionIndices){
	int collisionCount = 0;
	collisionIndices.clear();
	for (int i = 0; i < elasticObject.mesh.vertices.size(); i++){
		if (toolObject->isInterior(elasticObject.mesh.vertices[i])){
			collisionIndices.push_back(i);
			collisionCount++;
		}
	}
	if (1) printf("Total collisions %d \n", collisionCount);
}

void Simulator::IdentifyCollisionsCallBack(Visualization* v, const char*){
	std::vector<int> collisionIndices;
	IdentifyCollisions(collisionIndices);
	visualization.collisionIndices = collisionIndices;
	visualization.targetDeformationPositions.clear();
}

void Simulator::ComputeDistortion(float & distortion){
	std::vector<float> trianglePointwiseDistortion;
	distortion = elasticObject.ComputeDistortion(trianglePointwiseDistortion);
	printf("Total distortion %f \n", distortion);

	int vCount = elasticObject.mesh.vertices.size();
	int tCount = elasticObject.mesh.triangles.size();
	std::vector<float> vertexPointwiseDistortion(vCount, 0);
	std::vector<float> vertexCumWeight(vCount, 0);

	for (int t = 0; t < tCount; t++){
		float tArea = elasticObject.triangleAreas[t];
		float tDistortion = trianglePointwiseDistortion[t];
		for (int k = 0; k < 3; k++){
			vertexPointwiseDistortion[elasticObject.mesh.triangles[t][k]] += tArea*tDistortion;
			vertexCumWeight[elasticObject.mesh.triangles[t][k]] += tArea;
		}
	}

	for (int v = 0; v < vCount; v++){
		float vertexDistortion = vertexPointwiseDistortion[v] / vertexCumWeight[v];
		float distortionFactor;
		if (elasticObject.distortionMode == ISOMETRIC_DRICHLET){
			distortionFactor = std::min<float>(1.0, 10.0*vertexDistortion);
		}
		else if (elasticObject.distortionMode == ARAP ){
			distortionFactor = std::min<float>(1.0, vertexDistortion / 4.0);
		}
		else if (elasticObject.distortionMode == DISTANCE_REST_STATE){
			distortionFactor = std::min<float>(1.0, vertexDistortion*30.f); 
		}
		visualization.colors[v] = Eigen::Vector4f(1.0, 0, 0, .8)*distortionFactor + Eigen::Vector4f(.7, .7, .7, .8)*(1.0 - distortionFactor);
	}

	visualization.UpdateVertexBuffer();
}

void Simulator::ComputeDistortionCallBack(Visualization* v, const char*){
	float distortion;
	ComputeDistortion(distortion);
	glutPostRedisplay();
}

void Simulator::SolveCollision(){

	std::vector<double> fittingWeight;
	std::vector<int> vertexIndices;
	std::vector<Eigen::Vector3f> targetPositions;

	float temporalWeightDecrease = 5.0;
	float initialWeightFitting = 1.0;
	float padding = collisionSensitivity;
	for (int i = 0; i < elasticObject.mesh.vertices.size(); i++){
		if (toolObject->isInterior(elasticObject.mesh.vertices[i])){

			//Eigen::Vector3f normalIntersection;
			Eigen::Vector3f closesPathPoint = pathQuery->GetClosestPoint(elasticObject.mesh.vertices[i]);
			Eigen::Vector3f pathToSurfaceDirection = elasticObject.mesh.vertices[i] - closesPathPoint;
			Eigen::Vector3f pullingDirection = -elasticObject.smoothedNormals[i];
			if (pathToSurfaceDirection.dot(pullingDirection) < 0) pullingDirection *= -1;
			
			Eigen::Vector3f offset = elasticObject.mesh.vertices[i] + pullingDirection*padding;
			temporalFittingWeight[i] = initialWeightFitting;
			temporalTargetPosition[i] = offset;

			fittingWeight.push_back(initialWeightFitting);
			vertexIndices.push_back(i); 
			targetPositions.push_back(offset);
		}
		else{
			if (temporalFittingWeight[i] > 0){
				Eigen::Vector3f closestPoint = toolObject->closestSurfacePoint(elasticObject.mesh.vertices[i]);
				float distance = (closestPoint - elasticObject.mesh.vertices[i]).norm();
				if (distance < padding / 2.0){


					Eigen::Vector3f closesPathPoint = pathQuery->GetClosestPoint(elasticObject.mesh.vertices[i]);
					Eigen::Vector3f pathToSurfaceDirection = elasticObject.mesh.vertices[i] - closesPathPoint;
					Eigen::Vector3f pullingDirection = -elasticObject.smoothedNormals[i];
					if (pathToSurfaceDirection.dot(pullingDirection) < 0) pullingDirection *= -1;

					Eigen::Vector3f targetPos = elasticObject.mesh.vertices[i] + pullingDirection*padding;

					temporalFittingWeight[i] = initialWeightFitting;
					fittingWeight.push_back(temporalFittingWeight[i]);
					vertexIndices.push_back(i);
					targetPositions.push_back(targetPos);
				}
				else{
					Eigen::Vector3f restPoseVector = elasticObject.undeformedVertexPositions[i] - elasticObject.mesh.vertices[i];
					if (restPoseVector.norm() < (padding / 2.0)){
						temporalFittingWeight[i] = 0;
					}
					else{
						Eigen::Vector3f restPoseDirection = restPoseVector/restPoseVector.norm();
						Eigen::Vector3f targetPos = elasticObject.mesh.vertices[i] + restPoseDirection * (padding / 8.0);
						temporalFittingWeight[i] = initialWeightFitting;
						fittingWeight.push_back(temporalFittingWeight[i]);
						vertexIndices.push_back(i);
						targetPositions.push_back(targetPos);
					}
				}
			}

			//}
		}
	}




	if (1)printf("Selected vertices: %d \n", vertexIndices.size());
	visualization.targetDeformationPositions = targetPositions;

	if (vertexIndices.size() > 0){
		elasticObject.SolveDeformation(fittingWeight, vertexIndices, targetPositions);
		visualization.vertices = elasticObject.mesh.vertices;
		visualization.normals = elasticObject.mesh.normals;
		visualization.UpdateVertexBuffer();
	}
}

void Simulator::SolveCollisionsCallBack(Visualization* v, const char*){
	SolveCollision();
	glutPostRedisplay();
}


double Simulator::DistanceToTarget(){
	return pathQuery->GetTriangularDistance(toolObject->GetWorldPosition(toolTip));
}

int Simulator::GeneratePath(const char* outputPrefix, const int centerToolPeriod, const int enforcedAdvancePeriod, const int avoidSurfacePeriod){
	
	PathInfo pathInfo;

	std::vector<Eigen::Vector3f> & translations = pathInfo.translations;
	std::vector<Eigen::Matrix3f> & rotations = pathInfo.rotations;
	
	double distance = DistanceToTarget();

	translations.push_back(toolObject->objectToWorldTranslation);
	rotations.push_back(toolObject->objectToWorldRotation);

	printf("Initial Distance %f \n", distance);

	
	int _advanceIteration = 0; 
	int advanceItertation = 0;
	double targetDistance = 0.05;
	int maxAdvanceIterations = 1000;
	int maxCollisionSolveIterations = 10;

	float centeringDistanceStart = 1.0;
	float centeringDistanceEnd = 0.1;
	float centeringDistanceMaxWeight = 0.9; 

	bool excesiveCollision = false;

	Image<Point3D<float>> image;

	while (!excesiveCollision && distance > targetDistance && advanceItertation < maxAdvanceIterations){
		
		Eigen::Matrix3f currentToolRotation = toolObject->objectToWorldRotation;
		Eigen::Vector3f currentToolTranslation = toolObject->objectToWorldTranslation;
		Eigen::Vector3f currentToolCentroid = toolObject->GetCenterOfMass();
		Eigen::Matrix3f currentToolAxisAlignedFrame = currentToolRotation * axisAlignedFrame;
		Eigen::Vector3f alignedToolAxis = currentToolAxisAlignedFrame.col(0);

		AdvanceVectorField(!(_advanceIteration % enforcedAdvancePeriod == 0));
		if (_advanceIteration % avoidSurfacePeriod == 0) AvoidSurface(100, collisionSensitivity / 16, 5);
		if (_advanceIteration % centerToolPeriod == 0){
			float centerigWeight = distance > centeringDistanceStart ? 0 : std::min<float>(pow(((centeringDistanceStart - distance) / (centeringDistanceStart - centeringDistanceEnd)), 4.0), 1.0)*centeringDistanceMaxWeight;
			printf("Centering Weight %f \n", centerigWeight);
			CenterTool(centerigWeight);
		}

		SolveCollision();
		std::vector<int> collisionIndices;
		IdentifyCollisions(collisionIndices);
		printf("\t Collisions %d \n", collisionIndices.size());
		int collisionSolverIteration = 0;
		while (collisionIndices.size() && collisionSolverIteration < maxCollisionSolveIterations){
			SolveCollision();
			IdentifyCollisions(collisionIndices);
			printf("\t Collisions %d \n", collisionIndices.size());
			collisionSolverIteration++;
		}

		translations.push_back(toolObject->objectToWorldTranslation);
		rotations.push_back(toolObject->objectToWorldRotation);

		Eigen::Matrix3f newToolRotation = toolObject->objectToWorldRotation;
		Eigen::Vector3f newToolTranslation = toolObject->objectToWorldTranslation;

		Eigen::Matrix3f stepRotation = newToolRotation*currentToolRotation.inverse();
		Eigen::AngleAxisf angleAxisStepRotation(stepRotation);

		float stepRotationAngle = angleAxisStepRotation.angle();
		Eigen::Vector3f stepRotationAxis = angleAxisStepRotation.axis();
		Eigen::Vector3f axialRotationCommponent = alignedToolAxis * (stepRotationAxis.dot(alignedToolAxis));
		Eigen::Vector3f orthogonalRotationComponent = stepRotationAxis - axialRotationCommponent;

		Eigen::Vector3f stepTranslation = stepRotation*currentToolCentroid - currentToolCentroid + newToolTranslation - stepRotation*currentToolTranslation;
		Eigen::Vector3f axialTranslationComponent = alignedToolAxis * (stepTranslation.dot(alignedToolAxis));
		Eigen::Vector3f orthogonalTranslationComponent = stepTranslation - axialTranslationComponent;

		printf("Rotation Magnitude:\n Axial(%.6f) Orthogonal (%.6f) \n", axialRotationCommponent.squaredNorm() * stepRotationAngle, orthogonalRotationComponent.squaredNorm() *stepRotationAngle);
		printf("Translation Magnitude:\n Axial(%.6f) Orthogonal (%.6f) \n", axialTranslationComponent.norm(), orthogonalTranslationComponent.norm());
		 

		distance = DistanceToTarget();
		printf("Frame %06d. Current Distance %f \n", advanceItertation, distance);
		float distortion;
		ComputeDistortion(distortion);

		pathInfo.axialRotationEnergy += axialRotationCommponent.squaredNorm() * stepRotationAngle * axialRotationWeight;
		pathInfo.orthogonalRotationEnergy += orthogonalRotationComponent.squaredNorm() *stepRotationAngle * orthogonalRotationWeight;
		pathInfo.axialTranslationEnergy += axialTranslationComponent.norm() * axialTranslationWeight;
		pathInfo.orthogonalTranslationEnergy += orthogonalTranslationComponent.norm() * orthogonalTranslationWeight;

		pathInfo.cummulativeDistortionEnergy += distortion;
		pathInfo.maxDistortionEnergy = std::max<float>(pathInfo.maxDistortionEnergy, distortion);

		visualization.RenderOffScreenBuffer(image);
		char outputImage[256];
		sprintf(outputImage, "Frame-%06d.png", advanceItertation);
		image.write(outputImage);
		
		FILE * file;
		char outputPose[256];
		sprintf(outputPose, "Pose-%06d.bin", advanceItertation);
		file = fopen(outputPose, "wb");
		fwrite(&toolObject->objectToWorldTranslation, sizeof(Eigen::Vector3f), 1, file);
		fwrite(&toolObject->objectToWorldRotation, sizeof(Eigen::Matrix3f), 1, file);
		fclose(file);

		advanceItertation++;
		_advanceIteration++;

		if (collisionSolverIteration == maxCollisionSolveIterations){
			printf("WARNING: Collision non solved! \n");
			_advanceIteration = 0;
			if (collisionIndices.size() > 30) excesiveCollision = true;
		}
	}
	if (outputPrefix != NULL){
		SavePathStatistics(pathInfo, outputPrefix);
		SavePathTransformations(pathInfo, outputPrefix);
		char deformeshMeshName[256];
		sprintf(deformeshMeshName, "%s_DeformedROI.ply", outputPrefix);
		ExportDeformedMesh(deformeshMeshName);
	}
	return 1;
}
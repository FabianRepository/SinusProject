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

#include "FixedVerticesVisualization.inl"
#include <Src/Mesh.h>
#include <Src/VectorIO.h>
#include <igl/point_mesh_squared_distance.h>
#include <igl/AABB.h>
#include <igl/barycentric_coordinates.h>
#include <unordered_map>

class FixedVertices
{
public:
	static SimpleMesh roiMesh;
	static Eigen::Vector3f meshCentroid;
	static float meshScale;

	static SimpleMesh fixedMesh;
	static std::vector<double> closesPointDistance;
	static FixedVerticesVisualization visualization;

	static float fixedRadius;
	static bool isBoundaryFixed;

	static std::vector<bool> isBoundaryVertex;
	static std::vector<int> isFixedVertex;

	static void UpdateBrushAction( const int vertexIndex, const float radius, bool add = true);

	static void UpdateFixedVertices();

	static void UpdateFixedRadiusCallBack(Visualization* v, const char* prompt);
	static void ExportFixedIndicesCallBack(Visualization* v, const char* prompt);
	static void ToggleBoundaryVerticesCallBack(Visualization* v, const char* prompt);

	static int Init(const char* roiMeshName, const char* fixedMeshName);
	static void Display(void){ visualization.Display(); }
	static void MouseFunc(int button, int state, int x, int y);
	static void MotionFunc(int x, int y);
	static void Reshape(int w, int h){ visualization.Reshape(w, h); }
	static void KeyboardFunc(unsigned char key, int x, int y){ visualization.KeyboardFunc(key, x, y); }
};

FixedVerticesVisualization						FixedVertices::visualization;
SimpleMesh										FixedVertices::roiMesh;
SimpleMesh										FixedVertices::fixedMesh;
float											FixedVertices::fixedRadius = 0;
std::vector<int>								FixedVertices::isFixedVertex;
std::vector<double>								FixedVertices::closesPointDistance;
bool											FixedVertices::isBoundaryFixed = false;
std::vector<bool>								FixedVertices::isBoundaryVertex;
Eigen::Vector3f									FixedVertices::meshCentroid;
float											FixedVertices::meshScale;

int FixedVertices::Init(const char* roiMeshName, const char* fixedMeshName){
	if (!ReadSimpleMesh(roiMesh, roiMeshName)){
		printf("Unable to read mesh!\n");
		return 0;
	}

	GetMeshCentroidAndScale(roiMesh, meshCentroid, meshScale);
	for (int i = 0; i < roiMesh.vertices.size(); i++) roiMesh.vertices[i] = (roiMesh.vertices[i] - meshCentroid) / meshScale;

	printf("Centroid %f %f %f \n", meshCentroid[0], meshCentroid[1], meshCentroid[2]);
	printf("Scale %f \n", meshScale);

	visualization.vertices = roiMesh.vertices;
	visualization.triangles = roiMesh.triangles;
	visualization.normals = roiMesh.normals;
	visualization.colors.resize(roiMesh.vertices.size(), Eigen::Vector3f(.75, .75, .75));

	visualization.centroid = meshCentroid;
	visualization.scale = meshScale;

	visualization.callBacks.push_back(Visualization::KeyboardCallBack(&visualization, 's', "export fixed","File Name",ExportFixedIndicesCallBack));
	visualization.callBacks.push_back(Visualization::KeyboardCallBack(&visualization, 'b', "toggle boundary", ToggleBoundaryVerticesCallBack));
	visualization.callBacks.push_back(Visualization::KeyboardCallBack(&visualization, 'r', "fixed radius", "Radius(0.1 mm)", UpdateFixedRadiusCallBack));
	
	if (fixedMeshName != NULL){
		if (!ReadSimpleMesh(fixedMesh, fixedMeshName)){
			printf("Unable to read mesh!\n");
			return 0;
		}

		visualization.isFixedMeshInitialized = true;

		for (int i = 0; i < fixedMesh.vertices.size(); i++) fixedMesh.vertices[i] = (fixedMesh.vertices[i] - meshCentroid) / meshScale;
		visualization.fixedVertices = fixedMesh.vertices;
		visualization.fixedTriangles = fixedMesh.triangles;
		visualization.fixedNormals = fixedMesh.normals;
		visualization.fixedColors.resize(fixedMesh.vertices.size(), Eigen::Vector3f(.75, 1.0, .75));


		igl::AABB<Eigen::MatrixXd, 3> fixedMeshKDTree;
		Eigen::MatrixXd fixedMeshVertexMatrix;
		Eigen::MatrixXi fixedMeshTriangleMatrix;

		printf("Computing closest point distance...");
		fixedMeshVertexMatrix.resize(fixedMesh.vertices.size(), 3);
		for (int i = 0; i < fixedMesh.vertices.size(); i++)for (int c = 0; c < 3; c++) fixedMeshVertexMatrix(i, c) = fixedMesh.vertices[i][c];
		fixedMeshTriangleMatrix.resize(fixedMesh.triangles.size(), 3);
		for (int i = 0; i < fixedMesh.triangles.size(); i++)for (int c = 0; c < 3; c++) fixedMeshTriangleMatrix(i, c) = fixedMesh.triangles[i][c];
		fixedMeshKDTree.init(fixedMeshVertexMatrix, fixedMeshTriangleMatrix);

		Eigen::VectorXd sqrD;
		Eigen::VectorXi I;
		Eigen::MatrixXd C;
		Eigen::MatrixXd P;
		P.resize(roiMesh.vertices.size(), 3);
		for (int i = 0; i < roiMesh.vertices.size(); i++) for (int c = 0; c < 3; c++) P(i, c) = roiMesh.vertices[i][c];

		fixedMeshKDTree.squared_distance(fixedMeshVertexMatrix, fixedMeshTriangleMatrix, P, sqrD, I, C);
		closesPointDistance.resize(roiMesh.vertices.size());
		for (int i = 0; i < roiMesh.vertices.size(); i++) closesPointDistance[i] = sqrt(sqrD[i]);
	}
	else{
		visualization.isFixedMeshInitialized = false;
		closesPointDistance.resize(roiMesh.vertices.size(), DBL_MAX);
	}
	printf("Done! \n");

	printf("Computing boundary vertices ...");

	isBoundaryVertex.resize(roiMesh.vertices.size(), false);

	std::unordered_map<unsigned long long, int> edgeIndex;
	for (int i = 0; i < roiMesh.triangles.size(); i++){
		for (int k = 0; k < 3; k++){
			unsigned long long  edgeKey = SetMeshEdgeKey(roiMesh.triangles[i][k], roiMesh.triangles[i][(k + 1) % 3]);
			if (edgeIndex.find(edgeKey) == edgeIndex.end()){
				edgeIndex[edgeKey] = 3 * i + k;
			}
			else{
				printf("Non manifold mesh!! \n");
			}
		}
	}

	for (int i = 0; i < roiMesh.triangles.size(); i++){
		for (int k = 0; k < 3; k++){
			unsigned long long  edgeKey = SetMeshEdgeKey(roiMesh.triangles[i][(k + 1) % 3], roiMesh.triangles[i][k]);
			if (edgeIndex.find(edgeKey) == edgeIndex.end()){ //Boundary edge
				isBoundaryVertex[roiMesh.triangles[i][k]] = true;
				isBoundaryVertex[roiMesh.triangles[i][(k + 1) % 3]] = true;
			}
		}
	}
	printf("Done! \n");

	int boundaryCount = 0;
	for (int i = 0; i < isBoundaryVertex.size(); i++) if (isBoundaryVertex[i])boundaryCount++;
	printf("Boundary vertices %d of %d \n", boundaryCount, isBoundaryVertex.size());


	isFixedVertex.resize(roiMesh.vertices.size(), 0);
	UpdateFixedVertices();
	visualization.isFixedVertex = isFixedVertex;
}

void FixedVertices::UpdateFixedVertices(){
	for (int i = 0; i < roiMesh.vertices.size(); i++){
		if (closesPointDistance[i] < fixedRadius){
			isFixedVertex[i] = 1;
		}
		else if (isBoundaryVertex[i] && isBoundaryFixed){
			isFixedVertex[i] = 1;
		}
		else{
			isFixedVertex[i] = 0;
		}
	}
}

void FixedVertices::UpdateBrushAction(const int vertexIndex, const float radius, bool add){
	Eigen::Vector3f centerVertex = visualization.vertices[vertexIndex];
	float squaredRadius = radius * radius;
	for (int i = 0; i < visualization.vertices.size(); i++){
		if ((visualization.vertices[i] - centerVertex).squaredNorm() < squaredRadius){
			isFixedVertex[i] = add;
		}
	}
	visualization.isFixedVertex = isFixedVertex;
}

void FixedVertices::MouseFunc(int button, int state, int x, int y)
{
	if (state == GLUT_DOWN && glutGetModifiers() & GLUT_ACTIVE_SHIFT){

		visualization.isBrushActive = true;
		visualization.diskX = x;
		visualization.diskY = y;

		visualization.brushAdd = button == GLUT_LEFT_BUTTON;

		visualization.UpdateSelectedVertex(x, y);
		if (visualization.selectedVertex != -1) UpdateBrushAction(visualization.selectedVertex, visualization.brushRadius, visualization.brushAdd);
	}
	else{
		visualization.isBrushActive = false;

		visualization.panning = visualization.scaling = visualization.rotating = false;
		visualization.newX = x; visualization.newY = y;

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
	}

	glutPostRedisplay();
}
void FixedVertices::MotionFunc(int x, int y){
	
	if (visualization.isBrushActive){
		visualization.diskX = x;
		visualization.diskY = y;

		visualization.UpdateSelectedVertex(x, y);
		if (visualization.selectedVertex != -1) UpdateBrushAction(visualization.selectedVertex, visualization.brushRadius, visualization.brushAdd);
	}
	else{
		visualization.oldX = visualization.newX, visualization.oldY = visualization.newY, visualization.newX = x, visualization.newY = y;
		int screenSize = std::min< int >(visualization.screenWidth, visualization.screenHeight);
		float rel_x = (visualization.newX - visualization.oldX) / (float)screenSize * 2;
		float rel_y = (visualization.newY - visualization.oldY) / (float)screenSize * 2;
		float pRight = -rel_x * visualization.zoom, pUp = rel_y * visualization.zoom;
		float pForward = rel_y * visualization.zoom;
		float rRight = rel_y, rUp = rel_x;

		if (visualization.panning) visualization.worldCamera.translate(visualization.worldCamera.right * pRight + visualization.worldCamera.up * pUp);
		else if (visualization.scaling)visualization.worldCamera.translate(visualization.worldCamera.forward *pForward);
		else visualization.worldCamera.rotateUp(rUp), visualization.worldCamera.rotateRight(rRight);
	}

	glutPostRedisplay();
}

void FixedVertices::UpdateFixedRadiusCallBack(Visualization* v, const char* prompt){
	fixedRadius = atof(prompt)/meshScale;
	UpdateFixedVertices();
	visualization.isFixedVertex = isFixedVertex;
	glutPostRedisplay();
}

void FixedVertices::ExportFixedIndicesCallBack(Visualization* v, const char* prompt){
	WriteVector(isFixedVertex, prompt);
}

void FixedVertices::ToggleBoundaryVerticesCallBack(Visualization* v, const char* prompt){
	isBoundaryFixed = !isBoundaryFixed;
	UpdateFixedVertices();
	visualization.isFixedVertex = isFixedVertex;
	glutPostRedisplay();
}

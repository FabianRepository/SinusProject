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

#include <Src/IntersectableMesh.h>
#include "PathSelectionVisualization.inl"
#include <Src/Mesh.h>
#include <Src/VisibilityGrid.h>
#include <Src/VectorIO.h>
#include <Src/GeodesicDescriptors.h>
#include <igl/point_mesh_squared_distance.h>
#include <igl/AABB.h>
#include <igl/barycentric_coordinates.h>
#include <Util/KDtree.h>
#include <Geodesics/geodesic_algorithm_exact.h>

void UniformSampling(std::vector<Point3D<double>> & inputCurve, std::vector<Point3D<double>> & outputCurve, int numOutputNodes){
	int numInputNodes = inputCurve.size();
	std::vector<double> cumArcLength(numInputNodes, 0);
	double cumLength = 0;
	for (int i = 1; i < numInputNodes; i++){
		double segmentLength = Point3D<double>::Length(inputCurve[i] - inputCurve[i - 1]);
		cumLength += segmentLength;
		cumArcLength[i] = cumLength;
	}
	outputCurve.resize(numOutputNodes);
	outputCurve[0] = inputCurve[0];
	outputCurve[numOutputNodes - 1] = inputCurve[numInputNodes - 1];
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
	}
}

inline long long _HalfEdgeKey(int i1, int i2)
{
	return (((long long)i1) << 32) | ((long long)i2);
}
inline long long _EdgeKey(int i1, int i2)
{
	if (i1>i2) return _HalfEdgeKey(i1, i2);
	else		return _HalfEdgeKey(i2, i1);
}

bool SetTriangleAdjaceny(const std::vector< TriangleIndex> & triangles, int vCount, std::vector< TriangleIndex> & adjacency){
	std::unordered_map<long long, int> edgeToTriangle;
	for (int i = 0; i < triangles.size(); i++){
		for (int k = 0; k < 3; k++){
			long long edgeKey = _HalfEdgeKey(triangles[i][k], triangles[i][(k + 1) % 3]);
			if (edgeToTriangle.find(edgeKey) == edgeToTriangle.end()) edgeToTriangle[edgeKey] = i;
			else {
				printf("Non manifold mesh!! \n");
				return false;
			}
		}
	}

	adjacency.resize(triangles.size(), TriangleIndex(-1, -1, -1));
	for (int i = 0; i < triangles.size(); i++){
		for (int k = 0; k < 3; k++){
			long long edgeKey = _HalfEdgeKey(triangles[i][(k + 1) % 3], triangles[i][k]);
			if (edgeToTriangle.find(edgeKey) != edgeToTriangle.end()) adjacency[i][k] = edgeToTriangle[edgeKey];
		}
	}
	return true;
}

int AddComponent(std::vector<int> & vertexComponent, int vIndex, int currentComponent, const std::vector<std::vector<int>> & neighbours, int & componentSize){
	componentSize = 0;
	vertexComponent[vIndex] = currentComponent;
	std::queue<int> visitingQueue;
	visitingQueue.push(vIndex);
	componentSize++;

	while (!visitingQueue.empty()){
		int currentVertex = visitingQueue.front();
		visitingQueue.pop();
		const std::vector<int> & vertexNeighbours = neighbours[currentVertex];
		for (int i = 0; i < vertexNeighbours.size(); i++){
			if (vertexComponent[vertexNeighbours[i]] == -1) {
				vertexComponent[vertexNeighbours[i]] = currentComponent;
				visitingQueue.push(vertexNeighbours[i]);
				componentSize++;
			}
			else if (vertexComponent[vertexNeighbours[i]] == currentComponent){}
			else{
				printf("Unexpected Condition On A Connected Component. Expected %d. Obtained %d\n", currentComponent, vertexComponent[vertexNeighbours[i]]);
				return 0;
			}
		}
	}
	return 1;
}

void PreservedLargestComponent(const std::vector<TriangleIndex> & triangleAdjacency, std::vector<bool> & markedTriangle){


	int tCount = triangleAdjacency.size();

	std::vector<std::vector<int>> neighbours(tCount);
	for (int i = 0; i < tCount; i++){
		for (int k = 0; k < 3; k++){
			int neighbourIndex = triangleAdjacency[i][k];
			if (neighbourIndex != -1){
				if (markedTriangle[i] && markedTriangle[neighbourIndex]){
					neighbours[i].push_back(neighbourIndex);
					neighbours[neighbourIndex].push_back(i);
				}
			}
		}
	}

	std::vector<int> triangleComponent;
	triangleComponent.resize(tCount, -1);
	int currentComponent = -1;
	int largestComponentIndex = -1;
	int largestComponentSize = -1;
	for (int v = 0; v < tCount; v++){
		if (triangleComponent[v] == -1){
			currentComponent++;
			double currentComponentWeight;
			int currentComponentSize;
			AddComponent(triangleComponent, v, currentComponent, neighbours,currentComponentSize);
			if (currentComponentSize > largestComponentSize){
				largestComponentSize = currentComponentSize;
				largestComponentIndex = currentComponent;
			}
		}
	}

	int numComponents = currentComponent + 1;

	int removedTriangles = 0;
	int preservedTriangles = 0;
	for (int v = 0; v < tCount; v++){
		if (markedTriangle[v]){
			int vComponent = triangleComponent[v];
			if (vComponent != largestComponentIndex){
				markedTriangle[v] = false;
				removedTriangles++;
			}
			else{
				preservedTriangles++;
			}
		}
	}

	if (1) printf("Preserved %d. Removed  %d \n", preservedTriangles, removedTriangles);

}


class PathSelection
{
public:
	static SimpleMesh mesh;

	static Eigen::Vector3f meshCentroid;
	static float meshScale;
	static std::vector<TriangleIndex> triangleAdjacency;

	static int pathSource[3];
	static int pathTarget[3];
	static int croppedPathStartOffset;
	static int croppedPathEndOffset;

	static igl::AABB<Eigen::MatrixXd, 3> meshKDTree;
	static Eigen::MatrixXd vertexMatrix;
	static Eigen::MatrixXi triangleMatrix;

	static IntersectableMesh * intersectableMesh;
	static VisibilityGrid visibilityGrid;

	static PathSelectionVisualization visualization;

	static std::unordered_set<int> roiTriangles;
	static float visibilityRadius;
	static double geodesicPropagationDistance;

	static void UpdateBrushAction(const int selectedTriangle, const float radius, bool add);

	static void GeodesicGrowingCallBack(Visualization* v, const char* prompt);
	static void PreserveLargesComponentCallBack(Visualization* v, const char* prompt);
	static void ExtractVisibleMeshCallBack(Visualization* v, const char* prompt);
	static void UniformSamplingCallBack(Visualization* v, const char* prompt);
	static void	ExportPathCallBack(Visualization* v, const char* prompt);
	static void ExportRoiIndicesCallBack(Visualization* v, const char* prompt);
	static void ExportRoiMeshCallBack(Visualization* v, const char* prompt);
	static void	CropPathCallBack(Visualization* v, const char*);
	static void	AdvancePathStartCallBack(Visualization* v, const char*);
	static void	AdvancePathEndCallBack(Visualization* v, const char*);
	static void	ComputeShortestPathCallBack(Visualization* v, const char*);
	static void	SmoothPathCallBack(Visualization* v, const char*);

	static int Init(const char* meshName);
	static void Display(void){ visualization.Display(); }
	static void MouseFunc(int button, int state, int x, int y);
	static void MotionFunc(int x, int y);
	static void Reshape(int w, int h){ visualization.Reshape(w, h); }
	static void KeyboardFunc(unsigned char key, int x, int y){ visualization.KeyboardFunc(key, x, y); }
};

PathSelectionVisualization						PathSelection::visualization;
SimpleMesh										PathSelection::mesh;
int												PathSelection::pathSource[3];
int												PathSelection::pathTarget[3];
int												PathSelection::croppedPathStartOffset;
int												PathSelection::croppedPathEndOffset;
IntersectableMesh *								PathSelection::intersectableMesh;
VisibilityGrid									PathSelection::visibilityGrid;
igl::AABB<Eigen::MatrixXd, 3>					PathSelection::meshKDTree;
Eigen::MatrixXd									PathSelection::vertexMatrix;
Eigen::MatrixXi									PathSelection::triangleMatrix;
std::unordered_set<int>							PathSelection::roiTriangles;
float											PathSelection::visibilityRadius = 18;
double											PathSelection::geodesicPropagationDistance = 0.5;
Eigen::Vector3f									PathSelection::meshCentroid;
float											PathSelection::meshScale;
std::vector<TriangleIndex>						PathSelection::triangleAdjacency;

int PathSelection::Init(const char* meshName){
	if (!ReadSimpleMesh(mesh, meshName)){
		printf("Unable to read mesh!\n");
		return 0;
	}

	SetTriangleAdjaceny(mesh.triangles, mesh.vertices.size(), triangleAdjacency);

	GetMeshCentroidAndScale(mesh, meshCentroid, meshScale);
	for (int i = 0; i < mesh.vertices.size(); i++) mesh.vertices[i] = (mesh.vertices[i] - meshCentroid) / meshScale;

	printf("Centroid %f %f %f \n", meshCentroid[0], meshCentroid[1], meshCentroid[2]);
	printf("Scale %f \n", meshScale);

	//SaveEigenVector3f(centroid, "Centroid.bin");
	//SaveScalarf(scale, "Scale.bin");

	visualization.vertices = mesh.vertices;
	visualization.triangles = mesh.triangles;
	
	visualization.baricenters.resize(visualization.triangles.size());
	for (int i = 0; i < visualization.triangles.size(); i++) visualization.baricenters[i] = (visualization.vertices[visualization.triangles[i][0]] + visualization.vertices[visualization.triangles[i][1]] + visualization.vertices[visualization.triangles[i][2]]) / 3.0;

	visualization.normals = mesh.normals;
	visualization.colors.resize(visualization.vertices.size(), Eigen::Vector3f(.75, .75, .75));
	visualization.centroid = meshCentroid;
	visualization.scale = meshScale;

	intersectableMesh = new IntersectableMesh(visualization.vertices, visualization.triangles);
	int res = 128;
	visibilityGrid = VisibilityGrid(intersectableMesh, res);

	printf("Building IGL KD-Tree ...");
	vertexMatrix.resize(visualization.vertices.size(), 3);
	for (int i = 0; i < visualization.vertices.size(); i++)for (int c = 0; c < 3; c++) vertexMatrix(i, c) = visualization.vertices[i][c];
	triangleMatrix.resize(visualization.triangles.size(), 3);
	for (int i = 0; i < visualization.triangles.size(); i++)for (int c = 0; c < 3; c++) triangleMatrix(i, c) = visualization.triangles[i][c];
	meshKDTree.init(vertexMatrix, triangleMatrix);
	printf("Done! \n");

	visualization.callBacks.push_back(Visualization::KeyboardCallBack(&visualization, 'u', "compute path", ComputeShortestPathCallBack));
	visualization.callBacks.push_back(Visualization::KeyboardCallBack(&visualization, 'k', "uniform sampling", "Num Nodes(20)", UniformSamplingCallBack));
	visualization.callBacks.push_back(Visualization::KeyboardCallBack(&visualization, 'y', "smooth path", SmoothPathCallBack));
	visualization.callBacks.push_back(Visualization::KeyboardCallBack(&visualization, '+', "advance start", AdvancePathStartCallBack));
	visualization.callBacks.push_back(Visualization::KeyboardCallBack(&visualization, '-', "advance end", AdvancePathEndCallBack));
	visualization.callBacks.push_back(Visualization::KeyboardCallBack(&visualization, 'c', "crop path", CropPathCallBack));
	visualization.callBacks.push_back(Visualization::KeyboardCallBack(&visualization, 'p', "export path", "File Name", ExportPathCallBack));
	visualization.callBacks.push_back(Visualization::KeyboardCallBack(&visualization, 'i', "export roi indices", "File Name", ExportRoiIndicesCallBack));
	visualization.callBacks.push_back(Visualization::KeyboardCallBack(&visualization, 'o', "export roi mesh", "File Name", ExportRoiMeshCallBack));
	visualization.callBacks.push_back(Visualization::KeyboardCallBack(&visualization, 'v', "visible mesh", "Radius(18mm)", ExtractVisibleMeshCallBack));
	visualization.callBacks.push_back(Visualization::KeyboardCallBack(&visualization, 'g', "geodesic grow", "Radius(0.5mm)", GeodesicGrowingCallBack));
	visualization.callBacks.push_back(Visualization::KeyboardCallBack(&visualization, 'r', "preserve largest",PreserveLargesComponentCallBack));
	
}


void PathSelection::UpdateBrushAction(const int selectedTriangle, const float radius, bool add){
	
	std::queue<int> visitingQueue;
	visitingQueue.push(selectedTriangle);
	float squaredRadius = radius*radius;
	Eigen::Vector3f seedBaricenter = visualization.baricenters[selectedTriangle];
	std::unordered_set<int> alreadyVisited;
	alreadyVisited.insert(selectedTriangle);
	
	if (add) roiTriangles.insert(selectedTriangle);
	else roiTriangles.erase(selectedTriangle);

	while (!visitingQueue.empty()){
		int currentTriangle = visitingQueue.front();
		visitingQueue.pop();
		for (int n = 0; n < 3; n++){
			int neighbourTriangle = triangleAdjacency[currentTriangle][n];
			if (neighbourTriangle != -1){
				if (alreadyVisited.find(neighbourTriangle) == alreadyVisited.end()){
					alreadyVisited.insert(neighbourTriangle);
					Eigen::Vector3f neighbourBaricenter = visualization.baricenters[neighbourTriangle];
					if ((neighbourBaricenter - seedBaricenter).squaredNorm() < squaredRadius){
						visitingQueue.push(neighbourTriangle);
						if (add) roiTriangles.insert(neighbourTriangle);
						else roiTriangles.erase(neighbourTriangle);
					}
				}
			}
		}
	}

	visualization.roiTriangles = roiTriangles;
}





void PathSelection::MouseFunc(int button, int state, int x, int y)
{

	if (state == GLUT_DOWN && glutGetModifiers() & GLUT_ACTIVE_ALT){

		visualization.isBrushActive = true;
		visualization.diskX = x;
		visualization.diskY = y;

		visualization.brushAdd = button == GLUT_LEFT_BUTTON;

		visualization.UpdateSelectedTriangle(x, y);
		if (visualization.selectedTriangle != -1) UpdateBrushAction(visualization.selectedTriangle, visualization.brushRadius, visualization.brushAdd);
	}
	else{
		visualization.isBrushActive = false;
	}

	if (state == GLUT_DOWN && glutGetModifiers() & GLUT_ACTIVE_SHIFT){
		visualization.UpdateSelectedVertex(x, y);
		if (visualization.selectedVertex != -1){
			int i[3];
			visibilityGrid.bbox.GetClosestVoxel(visualization.vertices[visualization.selectedVertex], i);
			printf("i %d %d %d \n", i[0], i[1], i[2]);
			int ci[3];
			visibilityGrid.FindClosestValidVoxel(i, ci);
			printf("v %d %d %d \n", ci[0], ci[1], ci[2]);
			Eigen::Vector3f _vCenter = visibilityGrid.bbox.GetVoxelCenter(ci);
			Point3D<double> vCenter(_vCenter[0], _vCenter[1], _vCenter[2]);
			if (button == GLUT_LEFT_BUTTON){
				pathSource[0] = ci[0];
				pathSource[1] = ci[1];
				pathSource[2] = ci[2];

				visualization.pathStart = vCenter;
			}
			else if (button == GLUT_RIGHT_BUTTON){
				pathTarget[0] = ci[0];
				pathTarget[1] = ci[1];
				pathTarget[2] = ci[2];

				visualization.pathEnd = vCenter;
			}
		}
	}

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

	glutPostRedisplay();
}
void PathSelection::MotionFunc(int x, int y){

	if (visualization.isBrushActive){
		visualization.diskX = x;
		visualization.diskY = y;

		visualization.UpdateSelectedTriangle(x, y);
		if (visualization.selectedTriangle != -1) UpdateBrushAction(visualization.selectedTriangle, visualization.brushRadius, visualization.brushAdd);
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


void PathSelection::ComputeShortestPathCallBack(Visualization* v, const char*){
	if (pathSource[0] != -1 && pathSource[1] != -1 && pathSource[2] != -1 && pathTarget[0] != -1 && pathTarget[1] != -1 && pathTarget[2] != -1){
		std::vector<unsigned long> path;
		if (visibilityGrid.ComputeShortesPath(pathSource, pathTarget, path)){
			printf("Path found! Size %d. \n", path.size());
			std::vector<Point3D<double>> pathCenters;
			for (int s = 0; s < path.size(); s++){
				int ci[3];
				GetVoxelIndices(path[s], ci);
				if (0) printf("%d %d %d. \n", ci[0], ci[1], ci[2]);
				Eigen::Vector3f _pCenter = visibilityGrid.bbox.GetVoxelCenter(ci);
				Point3D<double> pCenter(_pCenter[0], _pCenter[1], _pCenter[2]);
				pathCenters.push_back(pCenter);
			}
			visualization.pathCenters = pathCenters;
		}
		else{
			printf("Path not found! \n");
		}
	}
	else{
		printf("Extrema not specified! \n");
	}
}

void PathSelection::SmoothPathCallBack(Visualization* v, const char*){
	//std::vector<Point3D<float>> closestPoints(visualization.pathCenters.size());
	if (visualization.pathCenters.size() > 0){
		Eigen::VectorXd sqrD;
		Eigen::VectorXi I;
		Eigen::MatrixXd C;
		Eigen::MatrixXd P;
		P.resize(visualization.pathCenters.size(), 3);
		for (int i = 0; i < visualization.pathCenters.size(); i++) for (int c = 0; c < 3; c++) P(i, c) = visualization.pathCenters[i][c];
		meshKDTree.squared_distance(vertexMatrix, triangleMatrix, P, sqrD, I, C);
		for (int i = 1; i < visualization.pathCenters.size() - 1; i++){
			Point3D<double> direction = Point3D<double>(C(i, 0), C(i, 1), C(i, 2)) - visualization.pathCenters[i];
			direction /= Point3D<double>::Length(direction);
			visualization.pathCenters[i] -= direction*0.001;
		}
		std::vector<Point3D<double>> vertexCopy = visualization.pathCenters;
		for (int i = 1; i < visualization.pathCenters.size() - 1; i++){
			visualization.pathCenters[i] = vertexCopy[i] * 0.5 + vertexCopy[i - 1] * 0.25 + vertexCopy[i + 1] * 0.25;
		}
		glutPostRedisplay();
	}
}

void PathSelection::CropPathCallBack(Visualization* v, const char*){
	if (visualization.pathCenters.size() > 0){
		std::vector<Point3D<double>> croppedPath(visualization.pathCenters.begin() + croppedPathStartOffset, visualization.pathCenters.end() - croppedPathEndOffset);
		croppedPathStartOffset = croppedPathEndOffset = 0;
		visualization.pathCenters = croppedPath;
		visualization.pathStart = visualization.pathCenters[0];
		visualization.pathEnd = visualization.pathCenters[visualization.pathCenters.size() - 1];
	}
	glutPostRedisplay();
}

void PathSelection::AdvancePathStartCallBack(Visualization* v, const char*){
	if (visualization.pathCenters.size() > 0){
		croppedPathStartOffset++;
		croppedPathStartOffset = std::min<int>(croppedPathStartOffset, visualization.pathCenters.size() - 1);
		visualization.pathStart = visualization.pathCenters[croppedPathStartOffset];
	}
	glutPostRedisplay();
}

void PathSelection::ExportPathCallBack(Visualization* v, const char* prompt){
	if (visualization.pathCenters.size() > 0){
		std::vector<Eigen::Vector3f> _pathCenters;
		_pathCenters.resize(visualization.pathCenters.size());
		for (int i = 0; i < visualization.pathCenters.size(); i++)_pathCenters[i] = Eigen::Vector3f(visualization.pathCenters[i][0], visualization.pathCenters[i][1], visualization.pathCenters[i][2]) * meshScale + meshCentroid;
		WriteVector(_pathCenters, prompt);
	}
	glutPostRedisplay();
}

void PathSelection::AdvancePathEndCallBack(Visualization* v, const char*){
	if (visualization.pathCenters.size() > 0){
		croppedPathEndOffset++;
		croppedPathEndOffset = std::min<int>(croppedPathEndOffset, visualization.pathCenters.size() - 1);
		visualization.pathEnd = visualization.pathCenters[visualization.pathCenters.size() - croppedPathEndOffset];
	}
	glutPostRedisplay();
}

void PathSelection::UniformSamplingCallBack(Visualization* v, const char* prompt){
	if (visualization.pathCenters.size() > 0){
		std::vector<Point3D<double>> newCenters;
		UniformSampling(visualization.pathCenters, newCenters, atoi(prompt));
		visualization.pathCenters = newCenters;
		glutPostRedisplay();
	}
}

void PathSelection::GeodesicGrowingCallBack(Visualization* v, const char* prompt){
	if (roiTriangles.size() > 0){
		geodesicPropagationDistance = atof(prompt)/meshScale;

		int vCount = mesh.vertices.size();
		std::vector<double> points(3 * vCount);
		for (int i = 0; i < vCount; i++) for (int c = 0; c < 3; c++) points[3 * i + c] = (double)mesh.vertices[i][c];
		int tCount = mesh.triangles.size();
		std::vector<unsigned> faces(3 * tCount);
		for (int i = 0; i < tCount; i++) for (int c = 0; c < 3; c++) faces[3 * i + c] = (unsigned)mesh.triangles[i][c];

		geodesic::Mesh geoMesh;
		geoMesh.initialize_mesh_data(points, faces); //create internal mesh data structure including edges
		geodesic::GeodesicAlgorithmExact * geoAlgorithm = new geodesic::GeodesicAlgorithmExact(&geoMesh);

		std::vector<geodesic::SurfacePoint> all_sources;
		std::unordered_set<int> visibleVertices;
		for (auto iter = roiTriangles.begin(); iter != roiTriangles.end(); iter++){
			for (int k = 0; k < 3; k++){
				int vIndex = mesh.triangles[*iter][k];
				if (visibleVertices.find(vIndex) == visibleVertices.end()){
					visibleVertices.insert(vIndex);
					geodesic::SurfacePoint sourcePoint(&geoMesh.vertices()[vIndex]);
					all_sources.push_back(sourcePoint);
				}
			}
		}

		printf("Source points %d \n", all_sources.size());

		int addedVertices = 0;
		geoAlgorithm->propagate(all_sources);	
		for (unsigned i = 0; i < geoMesh.vertices().size(); ++i)
		{
			geodesic::SurfacePoint p(&geoMesh.vertices()[i]);
			double distance;
			unsigned best_source = geoAlgorithm->best_source(p, distance);		
			if (distance < geodesicPropagationDistance && visibleVertices.find(i) == visibleVertices.end()){
				visibleVertices.insert(i);
				addedVertices++;
			}
		}
		int initialroiTriangles = roiTriangles.size();
		roiTriangles.clear();

		for (int i = 0; i < mesh.triangles.size(); i++){
			bool validTriangle = false;
			for (int k = 0; k < 3; k++){
				int vIndex = mesh.triangles[i][k];
				if (visibleVertices.find(vIndex) != visibleVertices.end()){
					validTriangle = true;
					break;
				}
			}
			if (validTriangle) roiTriangles.insert(i);
		}

		printf("Geodesic grow complete!. Added triangles  %d -> %d \n", initialroiTriangles,roiTriangles.size());
		visualization.roiTriangles = roiTriangles;
		glutPostRedisplay();
	}
}

void PathSelection::PreserveLargesComponentCallBack(Visualization* v, const char* prompt){
	std::vector<bool> isRoiTriangle(mesh.triangles.size(), false);
	for (auto iter = roiTriangles.begin(); iter != roiTriangles.end(); iter++)isRoiTriangle[*iter] = true;
	PreservedLargestComponent(triangleAdjacency, isRoiTriangle);
	roiTriangles.clear();
	for (int i = 0; i < isRoiTriangle.size(); i++) if (isRoiTriangle[i]) roiTriangles.insert(i);

	visualization.roiTriangles = roiTriangles;
	glutPostRedisplay();

}

void PathSelection::ExtractVisibleMeshCallBack(Visualization* v, const char* prompt){
	
	if (visualization.pathCenters.size() > 0){
		visibilityRadius = atof(prompt)/meshScale;

		roiTriangles.clear();

		KDtree<double, 3> querytree((double *)&visualization.pathCenters[0][0], visualization.pathCenters.size());

		int visibleTriangleCounter = 0;
		for (int t = 0; t < mesh.triangles.size(); t++){
			Eigen::Vector3f centroid = (mesh.vertices[mesh.triangles[t][0]] + mesh.vertices[mesh.triangles[t][1]] + mesh.vertices[mesh.triangles[t][2]]) / 3.0;
			Eigen::Vector3f tNormal = (mesh.vertices[mesh.triangles[t][1]] - mesh.vertices[mesh.triangles[t][0]]).cross(mesh.vertices[mesh.triangles[t][2]] - mesh.vertices[mesh.triangles[t][0]]);
			tNormal /= tNormal.norm();


			double queryPt[3] = { centroid[0], centroid[1], centroid[2] };

			const double * _closestPoint = querytree.closest_to_pt(&queryPt[0], DBL_MAX);
			Eigen::Vector3f closestPoint(_closestPoint[0], _closestPoint[1], _closestPoint[2]);

			if ((closestPoint - centroid).norm() < visibilityRadius){
				Eigen::Vector3f intersectionPoint;
				Eigen::Vector3f intersectionNormal;
				if (!intersectableMesh->edgeIntersect(centroid + tNormal*0.001, closestPoint, intersectionPoint, intersectionNormal)){
					roiTriangles.insert(t);
					visibleTriangleCounter++;
				}
				else{
					//Change the closest point to be in the normal direction
					float closestPointDistance = (closestPoint - centroid).norm();
					Eigen::Vector3f newQueryPt = centroid + tNormal * closestPointDistance;
					queryPt[0] = newQueryPt[0];
					queryPt[1] = newQueryPt[1];
					queryPt[2] = newQueryPt[2];
					
					_closestPoint = querytree.closest_to_pt(&queryPt[0], DBL_MAX);
					closestPoint[0] = _closestPoint[0];
					closestPoint[1] = _closestPoint[1];
					closestPoint[2] = _closestPoint[2];

					if ((closestPoint - centroid).norm() < visibilityRadius){
						if (!intersectableMesh->edgeIntersect(centroid + tNormal*0.001, closestPoint, intersectionPoint, intersectionNormal)){
							roiTriangles.insert(t);
							visibleTriangleCounter++;
						}
					}
				}
			}
		}

		//CPointNDList pathTreePoints(3, static_cast<int>(visualization.pathCenters.size()));
		//for (int i = 0; i < visualization.pathCenters.size(); i++)
		//{
		//	pathTreePoints.m_pPtND[i].m_fCoords[0] = visualization.pathCenters[i][0];
		//	pathTreePoints.m_pPtND[i].m_fCoords[1] = visualization.pathCenters[i][1];
		//	pathTreePoints.m_pPtND[i].m_fCoords[2] = visualization.pathCenters[i][2];
		//	pathTreePoints.m_pPtND[i].m_nIdx = i;
		//}
		//CKDTree pathKDTree;
		//pathKDTree.BuildKDTree(pathTreePoints, 100);

		//int visibleTriangleCounter = 0;
		//for (int t = 0; t < mesh.triangles.size(); t++){
		//	Eigen::Vector3f centroid = (mesh.vertices[mesh.triangles[t][0]] + mesh.vertices[mesh.triangles[t][1]] + mesh.vertices[mesh.triangles[t][2]]) / 3.0;
		//	Eigen::Vector3f tNormal = (mesh.vertices[mesh.triangles[t][1]] - mesh.vertices[mesh.triangles[t][0]]).cross(mesh.vertices[mesh.triangles[t][2]] - mesh.vertices[mesh.triangles[t][0]]);
		//	tNormal /= tNormal.norm();


		//	float queryPt[3] = { centroid[0], centroid[1], centroid[2] };
		//	PointND closestNode;
		//	pathKDTree.GetClosestPoint(queryPt, closestNode);
		//	int closestNodeIndex = closestNode.m_nIdx;

		//	Point3D<double> _closestPoint = visualization.pathCenters[closestNodeIndex];
		//	Eigen::Vector3f closestPoint(_closestPoint[0], _closestPoint[1], _closestPoint[2]);

		//	if ((closestPoint - centroid).norm() < visibilityRadius){
		//		Eigen::Vector3f intersectionPoint;
		//		Eigen::Vector3f intersectionNormal;
		//		if (!intersectableMesh->edgeIntersect(centroid + tNormal*0.001, closestPoint, intersectionPoint, intersectionNormal)){
		//			roiTriangles.insert(t);
		//			visibleTriangleCounter++;
		//		}
		//		else{
		//			//Change the closest point to be in the normal direction
		//			float closestPointDistance = (closestPoint - centroid).norm();
		//			Eigen::Vector3f newQueryPt = centroid + tNormal * closestPointDistance;
		//			queryPt[0] = newQueryPt[0];
		//			queryPt[1] = newQueryPt[1];
		//			queryPt[2] = newQueryPt[2];
		//			pathKDTree.GetClosestPoint(queryPt, closestNode);
		//			closestNodeIndex = closestNode.m_nIdx;

		//			_closestPoint = visualization.pathCenters[closestNodeIndex];
		//			closestPoint[0] = _closestPoint[0];
		//			closestPoint[1] = _closestPoint[1];
		//			closestPoint[2] = _closestPoint[2];

		//			if ((closestPoint - centroid).norm() < visibilityRadius){
		//				if (!intersectableMesh->edgeIntersect(centroid + tNormal*0.001, closestPoint, intersectionPoint, intersectionNormal)){
		//					roiTriangles.insert(t);
		//					visibleTriangleCounter++;
		//				}
		//			}
		//		}
		//	}
		//}

		if (1){//Preserve largest component
			std::vector<bool> isVisibleTriangle(mesh.triangles.size(), false);
			for (auto iter = roiTriangles.begin(); iter != roiTriangles.end(); iter++)isVisibleTriangle[*iter] = true;
			PreservedLargestComponent(triangleAdjacency, isVisibleTriangle);
			roiTriangles.clear();
			for (int i = 0; i < isVisibleTriangle.size(); i++) if (isVisibleTriangle[i]) roiTriangles.insert(i);
		}

		printf("Visible triangles %d of %d \n", visibleTriangleCounter, mesh.triangles.size());
		visualization.roiTriangles = roiTriangles;
		glutPostRedisplay();
	}
}

void PathSelection::ExportRoiIndicesCallBack(Visualization* v, const char* prompt){
	if (roiTriangles.size() > 0){
		std::vector<int> roiTriangleVector(roiTriangles.begin(), roiTriangles.end());
		WriteVector(roiTriangleVector, prompt);
	}
	glutPostRedisplay();
}

void PathSelection::ExportRoiMeshCallBack(Visualization* v, const char* prompt){
	if (roiTriangles.size() > 0){
		SimpleMesh roiMesh;
		std::vector<int> reducedVertexIndices(mesh.vertices.size(), -1);
		int lastReducedIndex = 0;
		for (auto iter = roiTriangles.begin(); iter != roiTriangles.end(); iter++){
			int corners[3];
			for (int k = 0; k < 3; k++){
				int vIndex = mesh.triangles[*iter][k];
				if (reducedVertexIndices[vIndex] == -1){
					reducedVertexIndices[vIndex] = lastReducedIndex;
					lastReducedIndex++;
					roiMesh.vertices.push_back(mesh.vertices[vIndex] * meshScale + meshCentroid);
				}
				corners[k] = reducedVertexIndices[vIndex];
			}
			roiMesh.triangles.push_back(TriangleIndex(corners[0], corners[1], corners[2]));
		}
		WriteSimpleMesh(roiMesh, prompt);
	}
	glutPostRedisplay();
}
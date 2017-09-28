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

#ifndef DEFORMABLE_MESH_INCLUDED
#define DEFORMABLE_MESH_INCLUDED

#include "DeformationSolver.h"
#include "IntersectableObject.h"
#include "RayTracer.h"
#include <igl/point_mesh_squared_distance.h>
#include <igl/AABB.h>
#include <igl/barycentric_coordinates.h>

Eigen::Matrix2f GetOrientedFrame(Eigen::Vector3f corners[3]){
	Eigen::Vector3f e0 = corners[1] - corners[0];
	Eigen::Vector3f e1 = corners[2] - corners[0];
	Eigen::Vector3f n = e0.cross(e1);
	n /= n.norm();
	Eigen::Vector3f d0 = e0 / e0.norm();
	Eigen::Vector3f d1 = n.cross(d0);
	d1 /= d1.norm();
	Eigen::Matrix2f frame;
	frame(0, 0) = e0.dot(d0);
	frame(1, 0) = e0.dot(d1);
	frame(0, 1) = e1.dot(d0);
	frame(1, 1) = e1.dot(d1);
	return frame;
}

enum{
	ISOMETRIC_DRICHLET,
	ARAP,
	DISTANCE_REST_STATE,
};

class DeformableMesh : virtual public IntersectableObject{
public:
	DeformableMesh(){}
	DeformableMesh(const std::vector<Eigen::Vector3f> & p_vertices, const std::vector<TriangleIndex> & p_triangles, const Eigen::MatrixXf & geodesicDescriptor, const std::vector<bool> & fixedVertices, const int numLevels = 5, const float p_fixedWeight = 1e3, const float p_rigidityWeight = 1e3, const int p_finestLevelSolved = 1){
		mesh.vertices = p_vertices;
		mesh.triangles = p_triangles;
		UpdateNormals(mesh);

		undeformedVertexPositions = mesh.vertices;
		undeformedVertexNormals = mesh.normals;

		rayQuery.Init(mesh.vertices, mesh.triangles);
		undeformedRayQuery.Init(mesh.vertices, mesh.triangles);

		//Initialize triangleNormals
		
		triangleNormals.resize(mesh.triangles.size());
		triangleAreas.resize(mesh.triangles.size());
		inverseReferenceFrames.resize(mesh.triangles.size());
		for (int t = 0; t < mesh.triangles.size(); t++){

			Eigen::Vector3f corners[3] = { mesh.vertices[mesh.triangles[t][0]], mesh.vertices[mesh.triangles[t][1]], mesh.vertices[mesh.triangles[t][2]] };
			Eigen::Matrix2f referenceFrame = GetOrientedFrame(corners);
			inverseReferenceFrames[t] = referenceFrame.inverse();

			Eigen::Vector3f d01 = mesh.vertices[mesh.triangles[t][1]] - mesh.vertices[mesh.triangles[t][0]];
			Eigen::Vector3f d02 = mesh.vertices[mesh.triangles[t][2]] - mesh.vertices[mesh.triangles[t][0]];
			Eigen::Vector3f n = d01.cross(d02);
			triangleAreas[t] = (n).norm() / 2.0;
			triangleNormals[t] = n / (n).norm();
		}
		undeformedTriangleNormals = triangleNormals;

		//Initialize radius
		Eigen::Vector3f  center(0, 0, 0);
		for (int v = 0; v < mesh.vertices.size(); v++) center += mesh.vertices[v];
		center /= (double)(mesh.vertices.size());
		radius = 0;
		for (int v = 0; v < mesh.vertices.size(); v++)radius = std::max<double>(radius, (center - mesh.vertices[v]).squaredNorm());
		radius = sqrt(radius);
		undeformedRadius = radius;
		//Initialize kdTree
		vertexMatrix.resize(mesh.vertices.size(), 3);
		for (int i = 0; i < mesh.vertices.size(); i++)for (int c = 0; c < 3; c++) vertexMatrix(i, c) = mesh.vertices[i][c];
		triangleMatrix.resize(mesh.triangles.size(), 3);
		for (int i = 0; i < mesh.triangles.size(); i++)for (int c = 0; c < 3; c++) triangleMatrix(i, c) = mesh.triangles[i][c];
		kdTree.init(vertexMatrix, triangleMatrix);

		undeformedVertexMatrix = vertexMatrix;
		undeformedKdTree = kdTree;

		//Initialize Hierarchy
		if (!hierarchy.ConstructHierachy(mesh, fixedVertices, numLevels, geodesicDescriptor)){
			printf("Unable to construct hierarchy! \n");
		}

		fixedWeight = p_fixedWeight;
		rigidityWeight = p_rigidityWeight;
		finestLevelSolved = p_finestLevelSolved;

		distortionMode = DISTANCE_REST_STATE;

		//Smoother
		GetSystemMatrices(mesh, massMatrix, stiffnessMatrix);
		float meshArea = GetMeshArea(mesh);
		float smoothingWeigth = 1e-3;
		smoother = new Eigen::SimplicialLDLT<Eigen::SparseMatrix<float>> (massMatrix + stiffnessMatrix*meshArea*smoothingWeigth);
		smoothedNormals.resize(mesh.normals.size());
		SmoothSignal(mesh.normals, smoothedNormals, true);
	}

	double undeformedRadius;
	double radius;
	MeshRayIntersection  rayQuery;
	MeshRayIntersection  undeformedRayQuery;
	SimpleMesh mesh;
	std::vector<Eigen::Vector3f> undeformedVertexPositions;
	std::vector<Eigen::Vector3f> undeformedVertexNormals;
	std::vector<Eigen::Vector3f> triangleNormals;
	std::vector<Eigen::Vector3f> undeformedTriangleNormals;
	igl::AABB<Eigen::MatrixXf, 3> kdTree;
	Eigen::MatrixXf vertexMatrix;
	Eigen::MatrixXi triangleMatrix;

	Eigen::MatrixXf undeformedVertexMatrix;
	igl::AABB<Eigen::MatrixXf, 3> undeformedKdTree;

	std::vector<Eigen::Matrix2f> inverseReferenceFrames;
	std::vector<float> triangleAreas;
	int distortionMode;

	//Deformable model
	HierarchicalModel hierarchy;
	float fixedWeight;
	float rigidityWeight;
	int finestLevelSolved;

	//Smoother
	Eigen::SparseMatrix<float> massMatrix;
	Eigen::SparseMatrix<float> stiffnessMatrix;
	Eigen::SimplicialLDLT<Eigen::SparseMatrix<float>> * smoother;
	std::vector<Eigen::Vector3f> smoothedNormals;

	void SmoothSignal(const std::vector<Eigen::Vector3f> & in, std::vector<Eigen::Vector3f> & out, bool normalize = false){
		Eigen::MatrixXf x0;
		x0.resize(in.size(), 3);
		for (int i = 0; i < in.size(); i++) x0.row(i) = in[i];
		Eigen::MatrixXf solution = smoother->solve(massMatrix * x0);
		for (int i = 0; i < out.size(); i++){
			out[i] = solution.row(i);
			if (normalize)out[i] / out[i].norm();
		}
	}

	void SolveDeformation(const std::vector<double> &fittingWeight, const std::vector<int> &targetVertexIndex, const std::vector<Eigen::Vector3f> &targetVertexPosition, bool verbose = false){
		std::vector<Eigen::Vector3f> referenceVertexPosition(targetVertexIndex.size());
		std::vector<Eigen::Vector3f> referenceVertexNormal(targetVertexIndex.size());

		for (int i = 0; i < targetVertexIndex.size(); i++){
			referenceVertexPosition[i] = undeformedVertexPositions[targetVertexIndex[i]];
			referenceVertexNormal[i] = undeformedVertexNormals[targetVertexIndex[i]];
		}

		HierarchicalDeformationSolver(hierarchy, fixedWeight, rigidityWeight, fittingWeight, mesh.vertices, targetVertexIndex, targetVertexPosition, referenceVertexPosition, referenceVertexNormal, finestLevelSolved,verbose);

		//Update normals
		UpdateNormals(mesh);
		SmoothSignal(mesh.normals, smoothedNormals, true);
		
		//Update ray tracer
		rayQuery.Terminate();
		rayQuery.Init(mesh.vertices, mesh.triangles);

		//Update triangleNormals
		for (int t = 0; t < mesh.triangles.size(); t++){
			Eigen::Vector3f d01 = mesh.vertices[mesh.triangles[t][1]] - mesh.vertices[mesh.triangles[t][0]];
			Eigen::Vector3f d02 = mesh.vertices[mesh.triangles[t][2]] - mesh.vertices[mesh.triangles[t][0]];
			Eigen::Vector3f n = d01.cross(d02);
			triangleNormals[t] = n / (n).norm();
		}

		//Update radius
		Eigen::Vector3f  center(0, 0, 0);
		for (int v = 0; v < mesh.vertices.size(); v++) center += mesh.vertices[v];
		center /= (double)(mesh.vertices.size());
		radius = 0;
		for (int v = 0; v < mesh.vertices.size(); v++)radius = std::max<double>(radius, (center - mesh.vertices[v]).norm());
		radius = sqrt(radius);

		//Update kdTree
		for (int i = 0; i < mesh.vertices.size(); i++)for (int c = 0; c < 3; c++) vertexMatrix(i, c) = mesh.vertices[i][c];
		kdTree.deinit();
		kdTree.init(vertexMatrix, triangleMatrix);
	}

	float ComputeDistortion(std::vector<float> & pointwiseDistortion, bool areaNormalize = true) const{

		pointwiseDistortion.resize(mesh.triangles.size());

		double integratedDistortion = 0;
		double cumArea = 0;
		for (int i = 0; i < mesh.triangles.size(); i++){
			Eigen::Vector3f corners[3] = { mesh.vertices[mesh.triangles[i][0]], mesh.vertices[mesh.triangles[i][1]], mesh.vertices[mesh.triangles[i][2]] };
			Eigen::Matrix2f deformedFrame = GetOrientedFrame(corners);
			Eigen::Matrix2f transformation = deformedFrame * inverseReferenceFrames[i];

			float distortion;
			if (distortionMode == ISOMETRIC_DRICHLET){
				distortion = (transformation.squaredNorm() + transformation.inverse().squaredNorm() - 4.0);
			}
			else if (distortionMode == ARAP){
				Eigen::JacobiSVD<Eigen::Matrix2f> mSVD(transformation, Eigen::ComputeFullU | Eigen::ComputeFullV);
				Eigen::Vector2f singularValues = mSVD.singularValues();
				distortion = (1.0 - singularValues[0])*(1.0 - singularValues[0]) + (1.0 - singularValues[1])*(1.0 - singularValues[1]);
			}
			else if (distortionMode == DISTANCE_REST_STATE){
				distortion = ((corners[0] + corners[1] + corners[2]) / 3.0 - (undeformedVertexPositions[mesh.triangles[i][0]] + undeformedVertexPositions[mesh.triangles[i][1]] + undeformedVertexPositions[mesh.triangles[i][2]]) / 3.0).norm() / undeformedRadius;
			}
			pointwiseDistortion[i] = distortion;
			integratedDistortion += distortion * triangleAreas[i];
			cumArea += triangleAreas[i];
		}
		if (areaNormalize) integratedDistortion /= cumArea;

		return integratedDistortion;
	}

	Eigen::Vector3f _GetCenterOfMass() const{
		Eigen::Vector3f cumCenter(0, 0, 0);
		double cumArea = 0;
		for (int t = 0; t < mesh.triangles.size(); t++){
			double area = ((mesh.vertices[mesh.triangles[t][1]] - mesh.vertices[mesh.triangles[t][0]]).cross(mesh.vertices[mesh.triangles[t][2]] - mesh.vertices[mesh.triangles[t][0]])).norm() / 2.0;
			Eigen::Vector3f center = (mesh.vertices[mesh.triangles[t][0]] + mesh.vertices[mesh.triangles[t][1]] + mesh.vertices[mesh.triangles[t][2]]) / 3.0;
			cumCenter += center*area;
			cumArea += area;
		}
		return (cumCenter / cumArea);
	}

	//Assume closed non self-intersection model
	bool _isInterior(const Eigen::Vector3f & p) const{
		Eigen::Vector3f direction(double(rand()) / double(RAND_MAX) - 0.5, double(rand()) / double(RAND_MAX) - 0.5, double(rand()) / double(RAND_MAX) - 0.5);
		direction /= (direction).norm();
		int tId;
		float u, v;
		bool intersect = rayQuery.IntersectRay(p, direction, radius, tId, u, v);
		if (intersect){
			return direction.dot(triangleNormals[tId]) > 0;
		}
		else{
			return false;
		}
	}
	Eigen::Vector3f _closestSurfacePoint(const Eigen::Vector3f & p) const{
		Eigen::VectorXf sqrD;
		Eigen::VectorXi I;
		Eigen::MatrixXf C;
		Eigen::MatrixXf P;
		P.resize(1, 3);
		P.row(0) = p;
		kdTree.squared_distance(vertexMatrix, triangleMatrix, P, sqrD, I, C);
		return C.row(0);
	}
	bool _edgeIntersect(const Eigen::Vector3f & p0, const Eigen::Vector3f & p1, Eigen::Vector3f & isect_p, Eigen::Vector3f & isect_n) const{
		Eigen::Vector3f direction = p1 - p0;
		double distance = (direction).norm();
		if (distance){
			direction /= distance;
			int tId;
			float u, v;
			return rayQuery.IntersectRay(p0, direction, distance, tId, u, v);
		}
		else{
			return false;
		}
	}
	/*void OpenGLInitialization(){
	UpdateVertexBuffer();
	UpdateFaceBuffer();
	}
	*/

	bool _undeformedRayIntersect(const Eigen::Vector3f & p, const Eigen::Vector3f & d, Eigen::Vector3f & isect_p) const{
		int tId;
		float u, v;
		bool intersects = undeformedRayQuery.IntersectRay(p, d, 2 * radius, tId, u, v);
		if (intersects){
			isect_p = undeformedVertexPositions[mesh.triangles[tId][0]] * (1.0 - u - v) + undeformedVertexPositions[mesh.triangles[tId][1]] * u + undeformedVertexPositions[mesh.triangles[tId][2]] * v;
			return true;
		}
		else{
			return false;
		}
	}


	bool undeformedRayIntersect(const Eigen::Vector3f & p, const Eigen::Vector3f & d, Eigen::Vector3f & isect_p){
		Eigen::Vector3f tp = worldToObjectRotation*p + worldToObjectTranslation;
		Eigen::Vector3f td = worldToObjectRotation*d;
		bool intersect = _undeformedRayIntersect(tp, td, isect_p);
		if (intersect){
			isect_p = objectToWorldRotation*isect_p + objectToWorldTranslation;
			return true;
		}
		else
			return false;
	}

	bool _rayIntersect(const Eigen::Vector3f & p, const Eigen::Vector3f & d, Eigen::Vector3f & isect_p) const{
		int tId;
		float u, v;
		bool intersects = rayQuery.IntersectRay(p, d, 2 * radius, tId, u, v);
		if (intersects){
			isect_p = mesh.vertices[mesh.triangles[tId][0]] * (1.0 - u - v) + mesh.vertices[mesh.triangles[tId][1]] * u + mesh.vertices[mesh.triangles[tId][2]] * v;
			return true;
		}
		else{
			return false;
		}
	}

	void _Draw(){
		if (!vertexBuffer || !normalBuffer) UpdateVertexBuffer();
		if (!faceBuffer) UpdateFaceBuffer();
		glBindBuffer(GL_ARRAY_BUFFER, vertexBuffer);
		glEnableClientState(GL_VERTEX_ARRAY);
		glVertexPointer(3, GL_FLOAT, 0, NULL);

		glBindBuffer(GL_ARRAY_BUFFER, normalBuffer);
		glEnableClientState(GL_NORMAL_ARRAY);
		glNormalPointer(GL_FLOAT, 0, NULL);

		glColor3f(0.75, 1.0, 0.75);

		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, faceBuffer);
		glDrawElements(GL_TRIANGLES, (GLsizei)(mesh.triangles.size() * 3), GL_UNSIGNED_INT, NULL);

		glDisableClientState(GL_VERTEX_ARRAY);
		glDisableClientState(GL_NORMAL_ARRAY);
		glBindBuffer(GL_ARRAY_BUFFER, 0);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

	}

	double _GetArea() const{
		printf("TODO\n");
		return  0;
	}
	Eigen::Vector3f _GetFirstMoment() const{
		printf("TODO\n");
		return Eigen::Vector3f(0, 0, 0);
	}
	Eigen::Matrix3f  _GetSecondMoment() const{
		printf("TODO\n");
		Eigen::Matrix3f secondMoment;
		return secondMoment;
	}

	Eigen::Vector3f _MinCorner() const{
		Eigen::Vector3f cumMin = mesh.vertices[0];
		for (int i = 0; i < mesh.vertices.size(); i++){
			for (int c = 0; c < 3; c++) cumMin[c] = std::min<double>(cumMin[c], mesh.vertices[i][c]);
		}
		return cumMin;
	}
	Eigen::Vector3f _MaxCorner() const{
		Eigen::Vector3f cumMax = mesh.vertices[0];
		for (int i = 0; i < mesh.vertices.size(); i++){
			for (int c = 0; c < 3; c++) cumMax[c] = std::max<double>(cumMax[c], mesh.vertices[i][c]);
		}
		return cumMax;
	}

	GLuint vertexBuffer = 0;
	GLuint normalBuffer = 0;
	GLuint faceBuffer = 0;
	void UpdateVertexBuffer();
	void UpdateFaceBuffer();

	void _InitializeSamples(const int numSamples){
		printf("WARNING: Not implemented! \n");
	}
};


void DeformableMesh::UpdateVertexBuffer(){
	if (!glIsBuffer(vertexBuffer)){
		glGenBuffers(1, &vertexBuffer);
	}
	if (!glIsBuffer(normalBuffer)){
		glGenBuffers(1, &normalBuffer);
	}
	int vCount = mesh.vertices.size();

	glBindBuffer(GL_ARRAY_BUFFER, vertexBuffer);
	glBufferData(GL_ARRAY_BUFFER, vCount * sizeof(Eigen::Vector3f), &mesh.vertices[0], GL_DYNAMIC_DRAW);
	glBindBuffer(GL_ARRAY_BUFFER, normalBuffer);
	glBufferData(GL_ARRAY_BUFFER, vCount * sizeof(Eigen::Vector3f), &mesh.normals[0], GL_DYNAMIC_DRAW);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
}

void DeformableMesh::UpdateFaceBuffer(){
	if (!glIsBuffer(faceBuffer)){
		glGenBuffers(1, &faceBuffer);
	}
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, faceBuffer);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, mesh.triangles.size() * sizeof(int)* 3, &mesh.triangles[0][0], GL_STATIC_DRAW);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
}

#endif //DEFORMABLE_MESH_INCLUDED
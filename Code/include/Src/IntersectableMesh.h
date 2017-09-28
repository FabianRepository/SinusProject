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

#ifndef INTERSECTABLE_MESH_INCLUDED
#define INTERSECTABLE_MESH_INCLUDED

#include "IntersectableObject.h"
#include "RayTracer.h"
#include <igl/point_mesh_squared_distance.h>
#include <igl/AABB.h>
#include <igl/barycentric_coordinates.h>

class IntersectableMesh : virtual public IntersectableObject{
public:
	IntersectableMesh(const std::vector<Eigen::Vector3f> & p_vertices, const std::vector<TriangleIndex> & p_triangles){
		mesh.vertices = p_vertices;
		mesh.triangles = p_triangles;
		UpdateNormals(mesh);
		rayQuery.Init(mesh.vertices, mesh.triangles);

		//Initialize tNormals
		tNormals.resize(mesh.triangles.size());
		for (int t = 0; t < mesh.triangles.size(); t++){
			Eigen::Vector3f d01 = mesh.vertices[mesh.triangles[t][1]] - mesh.vertices[mesh.triangles[t][0]];
			Eigen::Vector3f d02 = mesh.vertices[mesh.triangles[t][2]] - mesh.vertices[mesh.triangles[t][0]];
			Eigen::Vector3f n = d01.cross(d02);
			tNormals[t] = n / (n).norm();
		}

		//Initialize radius
		Eigen::Vector3f  center(0, 0, 0);
		for (int v = 0; v < mesh.vertices.size(); v++) center += mesh.vertices[v];
		center /= (double)(mesh.vertices.size());
		radius = 0;
		for (int v = 0; v < mesh.vertices.size(); v++)radius = std::max<double>(radius, (center - mesh.vertices[v]).squaredNorm());
		radius = sqrt(radius);

	}

	double radius;
	MeshRayIntersection  rayQuery;
	SimpleMesh mesh;
	std::vector<Eigen::Vector3f> tNormals;

	void Update(const std::vector<Eigen::Vector3f> & p_vertices){
		mesh.vertices = p_vertices;
		UpdateNormals(mesh);
		rayQuery.Terminate();
		rayQuery.Init(mesh.vertices, mesh.triangles);

		//Initialize tNormals
		for (int t = 0; t < mesh.triangles.size(); t++){
			Eigen::Vector3f d01 = mesh.vertices[mesh.triangles[t][1]] - mesh.vertices[mesh.triangles[t][0]];
			Eigen::Vector3f d02 = mesh.vertices[mesh.triangles[t][2]] - mesh.vertices[mesh.triangles[t][0]];
			Eigen::Vector3f n = d01.cross(d02);
			tNormals[t] = n / (n).norm();
		}

		//Initialize radius
		Eigen::Vector3f  center(0, 0, 0);
		for (int v = 0; v < mesh.vertices.size(); v++) center += mesh.vertices[v];
		center /= (double)(mesh.vertices.size());
		radius = 0;
		for (int v = 0; v < mesh.vertices.size(); v++)radius = std::max<double>(radius, (center - mesh.vertices[v]).norm());
		radius = sqrt(radius);
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
			return direction.dot(tNormals[tId]) > 0;
		}
		else{
			return false;
		}
	}
	Eigen::Vector3f _closestSurfacePoint(const Eigen::Vector3f & p) const{
		printf("Unimplemented! \n");
		return  Eigen::Vector3f();
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


void IntersectableMesh::UpdateVertexBuffer(){
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

void IntersectableMesh::UpdateFaceBuffer(){
	if (!glIsBuffer(faceBuffer)){
		glGenBuffers(1, &faceBuffer);
	}
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, faceBuffer);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, mesh.triangles.size() * sizeof(int)* 3, &mesh.triangles[0][0], GL_STATIC_DRAW);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
}

#endif// INTERSECTABLE_MESH_INCLUDED
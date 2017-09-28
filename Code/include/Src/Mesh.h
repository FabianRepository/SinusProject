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


#ifndef MESH_INCLUDED
#define MESH_INCLUDED

#include <Util/Ply.h>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>

unsigned long long SetMeshEdgeKey(const unsigned long i0, const unsigned long i1){
	return (((static_cast<unsigned long long>(i0) << 32) & 0xFFFFFFFF00000000) | (static_cast<unsigned long long>(i1) & 0x00000000FFFFFFFF));
}

void GetMeshEdgeIndices(unsigned long long key, unsigned long & i0, unsigned long & i1){
	i1 = static_cast<unsigned long>(key & 0x00000000FFFFFFFF);
	i0 = static_cast<unsigned long>((key >> 32) & 0x00000000FFFFFFFF);
}

class SimpleMesh{
public:
	std::vector<Eigen::Vector3f> vertices;
	std::vector<Eigen::Vector3f> normals;
	std::vector<TriangleIndex> triangles;
};

void GetMeshCentroidAndScale(const SimpleMesh & mesh, Eigen::Vector3f & centroid, float & scale){
	centroid = Eigen::Vector3f(0,0,0);
	for (int i = 0; i < mesh.vertices.size(); i++) {
		centroid += mesh.vertices[i];
	}
	centroid /= float(mesh.vertices.size());
	scale = 0;
	for (int i = 0; i < mesh.vertices.size(); i++) scale = std::max<float>(scale, (mesh.vertices[i] - centroid).norm());
}

void UpdateNormals(SimpleMesh & mesh){
	mesh.normals.resize(mesh.vertices.size(), Eigen::Vector3f(0.0, 0.0, 0.0));
	for (int t = 0; t < mesh.triangles.size(); t++){
		Eigen::Vector3f d01 = mesh.vertices[mesh.triangles[t][1]] - mesh.vertices[mesh.triangles[t][0]];
		Eigen::Vector3f d02 = mesh.vertices[mesh.triangles[t][2]] - mesh.vertices[mesh.triangles[t][0]];
		Eigen::Vector3f n = d01.cross(d02);
		for (int v = 0; v < 3; v++) mesh.normals[mesh.triangles[t][v]] += n;
	}

	for (int i = 0; i < mesh.normals.size(); i++){
		if (mesh.normals[i].norm()) mesh.normals[i] /= mesh.normals[i].norm();
	}
}

void UpdateNormals(const std::vector<Eigen::Vector3f> & vertices,const std::vector<TriangleIndex> & triangles, std::vector<Eigen::Vector3f> & normals){
	normals.clear();
	normals.resize(vertices.size(), Eigen::Vector3f(0.0, 0.0, 0.0));
	for (int t = 0; t < triangles.size(); t++){
		Eigen::Vector3f d01 = vertices[triangles[t][1]] - vertices[triangles[t][0]];
		Eigen::Vector3f d02 = vertices[triangles[t][2]] - vertices[triangles[t][0]];
		Eigen::Vector3f n = d01.cross(d02);
		for (int v = 0; v < 3; v++) normals[triangles[t][v]] += n;
	}

	for (int i = 0; i < normals.size(); i++){
		if (normals[i].norm()) normals[i] /= normals[i].norm();
	}
}


int ReadSimpleMesh(SimpleMesh & mesh, const char * fileName){
	mesh.vertices.clear();
	mesh.triangles.clear();
	int file_type;
	std::vector< PlyVertex< double > > ply_vertices;
	bool readFlags[PlyVertex< double >::ReadComponents];
	if (!PlyReadTriangles(fileName, ply_vertices, mesh.triangles, PlyVertex< double >::ReadProperties, readFlags, PlyVertex< double >::ReadComponents, file_type)) return 0;
	mesh.vertices.resize(ply_vertices.size());
	for (int i = 0; i < ply_vertices.size(); i++) mesh.vertices[i] = Eigen::Vector3f(ply_vertices[i].point[0], ply_vertices[i].point[1], ply_vertices[i].point[2]);
	UpdateNormals(mesh);
	return 1;
}

void WriteSimpleMesh(SimpleMesh & mesh, const char * fileName){
	std::vector< PlyVertex< float > > ply_vertices(mesh.vertices.size());
	for (int i = 0; i<mesh.vertices.size(); i++) ply_vertices[i].point = Point3D<float>(mesh.vertices[i][0], mesh.vertices[i][1], mesh.vertices[i][2]);
	PlyWriteTriangles(fileName, ply_vertices, mesh.triangles, PlyVertex< float >::WriteProperties, PlyVertex< float >::WriteComponents, PLY_BINARY_NATIVE);
}

Eigen::Matrix2f GetMetricFromEmbedding(const Eigen::Vector3f & corners_0, const Eigen::Vector3f & corners_1, const Eigen::Vector3f & corners_2){
	Eigen::Matrix2f g;
	g(0, 0) = (corners_1 - corners_0).dot(corners_1 - corners_0);
	g(1, 0) = g(0, 1) = (corners_1 - corners_0).dot(corners_2 - corners_0);
	g(1, 1) = (corners_2 - corners_0).dot(corners_2 - corners_0);

	return g;
}

Eigen::Matrix3f GetMassMatrix(const Eigen::Matrix2f & g, bool lump = false){
	double area = sqrt(g.determinant()) / 2.0;
	Eigen::Matrix3f mass;
	if (lump){
		mass(0, 0) = mass(1, 1) = mass(2, 2) = area*(1.0 / 3.0);
		mass(0, 1) = mass(0, 2) = mass(1, 0) = mass(1, 2) = mass(2, 0) = mass(2, 1) = 0.0;
	}
	else{
		mass(0, 0) = mass(1, 1) = mass(2, 2) = area*(1.0 / 6.0);
		mass(0, 1) = mass(0, 2) = mass(1, 0) = mass(1, 2) = mass(2, 0) = mass(2, 1) = area*(1.0 / 12.0);
	}
	return mass;
}

void GetVertexArea(const SimpleMesh & mesh, std::vector<double> & vertexArea){
	vertexArea.resize(mesh.vertices.size(), 0.0);
	for (int i = 0; i < mesh.triangles.size(); i++){
		double area = ((mesh.vertices[mesh.triangles[i][1]] - mesh.vertices[mesh.triangles[i][0]]).cross(mesh.vertices[mesh.triangles[i][2]] - mesh.vertices[mesh.triangles[i][0]])).norm() / 2.0;
		for (int j = 0; j < 3; j++) vertexArea[mesh.triangles[i][j]] += area / 3.0;
	}
}

double GetMeshArea(const SimpleMesh & mesh){
	double area = 0.0;
	for (int i = 0; i < mesh.triangles.size(); i++){
		area += ((mesh.vertices[mesh.triangles[i][1]] - mesh.vertices[mesh.triangles[i][0]]).cross(mesh.vertices[mesh.triangles[i][2]] - mesh.vertices[mesh.triangles[i][0]])).norm() / 2.0;
	}
	return area;
}


Eigen::Matrix3f GetStiffnessMatrix(const Eigen::Matrix2f & g){
	double area = sqrt(g.determinant()) / 2.0;
	Eigen::Matrix3f stiffness;
	Eigen::Vector2f d[3];
	d[0] = Eigen::Vector2f(-1.0, -1.0);
	d[1] = Eigen::Vector2f(1.0, 0.0);
	d[2] = Eigen::Vector2f(0.0, 1.0);
	Eigen::Matrix2f g_inv = g.inverse();
	for (int i = 0; i < 3; i++)for (int j = 0; j < 3; j++) stiffness(i, j) = d[i].dot(g_inv*d[j])*area;
	return stiffness;
}


void GetSystemMatrices(const SimpleMesh & mesh, Eigen::SparseMatrix<float> & massMatrix, Eigen::SparseMatrix<float> & stiffnessMatrix){
	std::vector<Eigen::Triplet<float>> massTriplets;
	std::vector<Eigen::Triplet<float>> stiffnessTriplets;
	
	for (int t = 0; t < mesh.triangles.size(); t++){
		Eigen::Matrix2f metric = GetMetricFromEmbedding(mesh.vertices[mesh.triangles[t][0]], mesh.vertices[mesh.triangles[t][1]], mesh.vertices[mesh.triangles[t][2]]);
		Eigen::Matrix3f mass = GetMassMatrix(metric);
		Eigen::Matrix3f stiffness = GetStiffnessMatrix(metric);

		int vIndex[3] = { mesh.triangles[t][0], mesh.triangles[t][1], mesh.triangles[t][2] };
		for (int k = 0; k < 3; k++)for (int l = 0; l < 3; l++){
			massTriplets.push_back(Eigen::Triplet<float>(vIndex[k], vIndex[l], mass(k, l)));
			stiffnessTriplets.push_back(Eigen::Triplet<float>(vIndex[k], vIndex[l], stiffness(k, l)));
		}
	}

	massMatrix.resize(mesh.vertices.size(), mesh.vertices.size());
	massMatrix.setFromTriplets(massTriplets.begin(), massTriplets.end());

	stiffnessMatrix.resize(mesh.vertices.size(), mesh.vertices.size());
	stiffnessMatrix.setFromTriplets(stiffnessTriplets.begin(), stiffnessTriplets.end());
}

#endif //MESH_INCLUDED
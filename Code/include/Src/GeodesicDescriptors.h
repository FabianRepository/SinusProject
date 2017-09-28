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

#ifndef GEODESIC_DESCRIPTORS_INCLUDED
#define GEODESIC_DESCRIPTORS_INCLUDED

#define GEODESIC_MAX_DISTANCE 100

#include "Mesh.h"
#include <queue>
#include <unordered_set>
#include <Geodesics/geodesic_algorithm_exact.h>

void GenerateGeodesicDescriptor(const SimpleMesh & mesh, const int dimension, Eigen::MatrixXf & descriptor, geodesic::Mesh & geoMesh, geodesic::GeodesicAlgorithmExact * geoAlgorithm){
	int vCount = mesh.vertices.size();
	descriptor.resize(vCount,dimension);

	Eigen::Vector3f average(0, 0, 0);
	for (int i = 0; i < mesh.vertices.size(); i++) average += mesh.vertices[i];
	average /= float(mesh.vertices.size());

	float farthestDistance = 0.f;
	int farthestIndex = -1;
	for (int i = 0; i < mesh.vertices.size(); i++){
		float distance = (mesh.vertices[i] - average).norm();
		if (distance >farthestDistance){
			farthestDistance = distance;
			farthestIndex = i;
		}
	}

	std::vector<int> referenceVertices(dimension,-1);
	referenceVertices[0] = farthestIndex;

	std::vector<double> minDistanceValue;
	minDistanceValue.resize(mesh.vertices.size(), DBL_MAX);

	for (int s = 0; s < dimension; s++){
		printf("%d of %d \r", s, dimension);
		geodesic::SurfacePoint source(&geoMesh.vertices()[referenceVertices[s]]);		//create source 
		std::vector<geodesic::SurfacePoint> all_sources(1, source);					//in general, there could be multiple sources, but now we have only one

		geoAlgorithm->propagate(all_sources);	//cover the whole mesh
		for (unsigned i = 0; i < geoMesh.vertices().size(); ++i)
		{
			geodesic::SurfacePoint p(&geoMesh.vertices()[i]);
			double distance;
			unsigned best_source = geoAlgorithm->best_source(p, distance);		//for a given surface point, find closets source and distance to this source
			distance = std::min<double>(distance, GEODESIC_MAX_DISTANCE);
			descriptor(i,s) = distance;
			minDistanceValue[i] = std::min<double>(distance, minDistanceValue[i]);
		}

		if (s < dimension - 1){
			float maxMinDistanceValue = 0.f;
			int maxMinDistanceIndex = -1;
			for (int i = 0; i < mesh.vertices.size(); i++) if (minDistanceValue[i] > maxMinDistanceValue){
				maxMinDistanceValue = minDistanceValue[i];
				maxMinDistanceIndex = i;
			}
			referenceVertices[s + 1] = maxMinDistanceIndex;
		}
	}
	//if (0){
	//	ColoredMesh descriptorVisualization;
	//	descriptorVisualization.vertices = mesh.vertices;
	//	descriptorVisualization.triangles = mesh.triangles;
	//	descriptorVisualization.colors.resize(mesh.vertices.size());
	//	for (int i = 0; i < descriptorVisualization.colors.size(); i++){
	//		descriptorVisualization.colors[i] = Point3D<double>(descriptor[i][0], descriptor[i][1], descriptor[i][2]);
	//	}

	//	for (int k = 0; k < 3; k++){
	//		double maxValue = 0;
	//		for (int v = 0; v < descriptorVisualization.colors.size(); v++) maxValue = std::max<double>(maxValue, descriptorVisualization.colors[v][k]);
	//		for (int v = 0; v < descriptorVisualization.colors.size(); v++) descriptorVisualization.colors[v][k] *= (255.0 / maxValue);
	//	}
	//	char outputName[256];
	//	sprintf(outputName, "geodesicDescriptor-%03d.ply", meshId);
	//	meshId++;
	//	WriteColoredMesh(descriptorVisualization, outputName);
	//}
	//if (0){
	//	ColoredMesh descriptorVisualization;
	//	descriptorVisualization.vertices = mesh.vertices;
	//	descriptorVisualization.triangles = mesh.triangles;
	//	descriptorVisualization.colors.resize(mesh.vertices.size());
	//	descriptorVisualization.colors[referenceVertices[0]] = Point3D<double>(255.0, 255.0, 255.0);
	//	descriptorVisualization.colors[referenceVertices[1]] = Point3D<double>(255.0, 0.0, 0.0);
	//	descriptorVisualization.colors[referenceVertices[2]] = Point3D<double>(0.0, 255.0, 0.0);
	//	char outputName[256];
	//	sprintf(outputName, "ReferenceVertices-%03d.ply", meshId);
	//	meshId++;
	//	WriteColoredMesh(descriptorVisualization, outputName);
	//}
}

#endif// GEODESIC_DESCRIPTORS_INCLUDED
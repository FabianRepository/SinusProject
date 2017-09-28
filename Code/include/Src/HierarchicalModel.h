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

#ifndef HIERARCHICAL_MODEL_INCLUDED
#define HIERARCHICAL_MODEL_INCLUDED
#include <Eigen/Sparse>
#include "Mesh.h"
#include "GeodesicDescriptors.h"
#include "RayTracer.h"
#include <unordered_set>

class RigidTransformation{
public:
	RigidTransformation(){
		rotation = Eigen::Matrix3f::Identity();
		translation = Eigen::Vector3f(0, 0, 0);
	}
	RigidTransformation(Eigen::Matrix3f p_rotation, Eigen::Vector3f p_translation){
		rotation = p_rotation;
		translation = p_translation;
	}
	Eigen::Matrix3f rotation;
	Eigen::Vector3f translation;
};


class IndexWeight{
public:
	IndexWeight(int p_index, float p_weight){
		index = p_index;
		weight = p_weight;
	}
	int index;
	float weight;
};


void TransformModel(const std::vector<Eigen::Vector3f> & restPosition, const std::vector<RigidTransformation> & transformations, const std::vector<std::vector<IndexWeight>> & parents, const std::vector<int> & fineIndices, std::vector<Eigen::Vector3f> & transformedPosition){
	std::vector<Eigen::Vector3f> nodePositions(fineIndices.size());
	for (int i = 0; i < fineIndices.size(); i++) nodePositions[i] = restPosition[fineIndices[i]];

	int vCount = restPosition.size();
	transformedPosition.resize(vCount);

	int threads = omp_get_num_procs();
#pragma omp parallel for num_threads( threads )
	for (int i = 0; i < vCount; i++){
		Eigen::Vector3f cumDeformation(0, 0, 0);
		for (int j = 0; j < parents[i].size(); j++){
			int parentIndex = parents[i][j].index;
			float parentWeight = parents[i][j].weight;
			cumDeformation += (transformations[parentIndex].rotation*(restPosition[i] - nodePositions[parentIndex]) + transformations[parentIndex].translation + nodePositions[parentIndex])*parentWeight;
		}
		transformedPosition[i] = cumDeformation;
	}
}

inline float ApproxGeodesicDistance(const int i, const int j, const Eigen::MatrixXf & geodesicDescriptor){
	return (geodesicDescriptor.row(i) - geodesicDescriptor.row(j)).lpNorm<Eigen::Infinity>();
}

struct IndexedAngle{
	float angle;
	int index;
};

bool CompareIndexedAngle(const IndexedAngle & a, const IndexedAngle & b){
	return a.angle < b.angle;
}

int ConstructProlongation(const SimpleMesh & mesh, const std::vector<bool> & fixedVertices, const Eigen::MatrixXf & geodesicDescriptor, const std::vector<std::vector<int>> & neighbours, const float radius, std::vector<int> & representatives,
	std::vector<std::vector<IndexWeight>> & nodeParents, std::vector<std::vector<IndexWeight>> & nodeSons,
	std::vector<std::vector<int>> & representativeNeighbours, std::vector<std::vector<TriangleIndex>> & representativeTriangles,
	std::vector<RigidTransformation> & frames, std::vector<RigidTransformation> & transformations, std::vector<bool> & fixedRepresentative){
	int vCount = neighbours.size();
	std::vector<bool> isDeepCovered(vCount, false);
	unsigned long long numRepresentatives = 0;

	nodeParents.resize(vCount);
	std::vector<float> cumParentWeight(vCount, 0);
	int averageNumParents = 0;

	std::vector<int> reducedIndex(vCount, -1);

	int fixedCount = 0;

	//Find representatives
	for (int i = 0; i < vCount; i++){
		if (!isDeepCovered[i]){

			std::vector<IndexWeight> currentSons;
			float cumSonsWeight = 0;
			reducedIndex[i] = numRepresentatives;
			representatives.push_back(i);

			bool isFixed = false;

			std::queue<int> visitingQueue;
			visitingQueue.push(i);
			std::unordered_set<int> alreadyAdded;
			alreadyAdded.insert(i);
			while (!visitingQueue.empty()){
				int currentVertex = visitingQueue.front();
				visitingQueue.pop();
				float distance = ApproxGeodesicDistance(i, currentVertex, geodesicDescriptor);

				if (distance <= 2.0*radius){
					
					float weight = radius > 0 ? 1.0 - (distance / (2.0*radius)) : 1.0;
					nodeParents[currentVertex].push_back(IndexWeight(numRepresentatives, weight));
					cumParentWeight[currentVertex] += weight;

					currentSons.push_back(IndexWeight(currentVertex, weight));
					cumSonsWeight += weight;

					if (fixedVertices[currentVertex]){
						isFixed = true;
					}

					averageNumParents++;

					if (distance <= radius){
						isDeepCovered[currentVertex] = true;
					}

					const std::vector<int> & vertexNeighbours = neighbours[currentVertex];
					for (int j = 0; j < vertexNeighbours.size(); j++){
						int neighbourVertex = vertexNeighbours[j];
						if (alreadyAdded.find(neighbourVertex) == alreadyAdded.end()){
							alreadyAdded.insert(neighbourVertex);
							visitingQueue.push(neighbourVertex);
						}
					}
				}
			}
			//Normalize sons weight
			for (int k = 0; k < currentSons.size(); k++){
				currentSons[k].weight /= cumSonsWeight;
			}

			if (currentSons.size() < 4 && !(radius == 0) ){
				printf("WARNING: Node with very few sons = %d! \n",currentSons.size());
			}

			nodeSons.push_back(currentSons);
			fixedRepresentative.push_back(isFixed);
			if (isFixed) fixedCount++;
			numRepresentatives++;
		}
	}

	if (radius == 0 && numRepresentatives != vCount){
		printf("ERROR: Num representatives unexpected! \n");
		return 0;
	}

	printf("Num representatives %d of %d \n", numRepresentatives, vCount);
	printf("Num fixed nodes %d of %d \n", fixedCount, numRepresentatives);
	printf("Averge num parents per node %d \n", averageNumParents / vCount);
	printf("Averge num sons per representative %d \n", averageNumParents / numRepresentatives);



	//Normalize parents weights
	for (int i = 0; i < vCount; i++){
		for (int j = 0; j < nodeParents[i].size(); j++){
			nodeParents[i][j].weight /= cumParentWeight[i];
		}
	}

	//Construct representative connectivity
	int avergeRepresentativeDegree = 0;
	if (radius != 0){
		representativeNeighbours.resize(numRepresentatives);
		for (int i = 0; i < numRepresentatives; i++){
			std::queue<int> visitingQueue;
			int root = representatives[i];
			visitingQueue.push(root);
			std::unordered_set<int> alreadyAdded;
			alreadyAdded.insert(root);
			while (!visitingQueue.empty()){
				int currentVertex = visitingQueue.front();
				visitingQueue.pop();
				float distance = ApproxGeodesicDistance(root, currentVertex, geodesicDescriptor);
				if (distance <= 2.0*radius){
					int j = reducedIndex[currentVertex];
					if (j != -1 && j != i && j < i){
						representativeNeighbours[i].push_back(j);
						representativeNeighbours[j].push_back(i);
						avergeRepresentativeDegree += 2;
					}
					const std::vector<int> & vertexNeighbours = neighbours[currentVertex];
					for (int k = 0; k < vertexNeighbours.size(); k++){
						int neighbourVertex = vertexNeighbours[k];
						if (alreadyAdded.find(neighbourVertex) == alreadyAdded.end()){
							alreadyAdded.insert(neighbourVertex);
							visitingQueue.push(neighbourVertex);
						}
					}
				}
			}
		}
	}
	else{
		representativeNeighbours = neighbours;
		for (int i = 0; i < representativeNeighbours.size(); i++)avergeRepresentativeDegree += representativeNeighbours[i].size();
	}

	avergeRepresentativeDegree /= numRepresentatives;
	printf("Averge representative degree %d \n", avergeRepresentativeDegree);

	frames.resize(numRepresentatives);
	transformations.resize(numRepresentatives);

	//Construct triangles //WARNING: This may modify the triangles of the finest level
	representativeTriangles.resize(numRepresentatives);
	for (int i = 0; i < numRepresentatives; i++){
		std::vector<int> currentNeighbours = representativeNeighbours[i];
		if (currentNeighbours.size() > 3){
			int currentFineIndex = representatives[i];
			Eigen::Vector3f currentNormal = mesh.normals[currentFineIndex];
			currentNormal /= currentNormal.norm();
			Eigen::Vector3f tangentialDirection(2.f * ((float(rand()) / FLT_MAX) - 0.5f), 2.f * ((float(rand()) / FLT_MAX) - 0.5f), 2.f * ((float(rand()) / FLT_MAX) - 0.5f));
			tangentialDirection = currentNormal.cross(tangentialDirection);
			tangentialDirection /= tangentialDirection.norm();
			Eigen::Vector3f biTangentialDirection = currentNormal.cross(tangentialDirection);
			biTangentialDirection /= biTangentialDirection.norm();

			Eigen::Matrix3f currentFrame;
			currentFrame.col(0) = currentNormal;
			currentFrame.col(1) = tangentialDirection;
			currentFrame.col(2) = biTangentialDirection;

			Eigen::Vector3f currentPos = mesh.vertices[currentFineIndex];

			frames[i] = RigidTransformation(currentFrame, currentPos);

			std::vector<IndexedAngle> angles(currentNeighbours.size());
			for (int n = 0; n < currentNeighbours.size(); n++){
				int j = currentNeighbours[n];
				int neighbourhFineIndex = representatives[j];
				Eigen::Vector3f neighbourPos = mesh.vertices[neighbourhFineIndex];
				Eigen::Vector3f direction = neighbourPos - currentPos;
				float projection[2] = { tangentialDirection.dot(direction), biTangentialDirection.dot(direction) };
				float angle = atan2(projection[0], projection[1]);
				angles[n].angle = angle;
				angles[n].index = j;
			}
			std::sort(angles.begin(), angles.end(),CompareIndexedAngle);
			for (int k = 0; k < angles.size(); k++) representativeTriangles[i].push_back(TriangleIndex(i, angles[k].index, angles[(k + 1) % currentNeighbours.size()].index));
		}
		else{
			if(0) printf("WARNING: Lonely vertex! \n");
		}
	}
	return 1;
}

class HierarchicalModel{
public:
	int ConstructHierachy(const SimpleMesh & mesh, const std::vector<bool> & fixedVertices, const int levels, const Eigen::MatrixXf & geodesicDescriptor);
	int numLevels;
	std::vector<std::vector<int>> hierarchyFineIndices;
	std::vector<std::vector<Eigen::Vector3f>> hierarchyReferencePosition;
	std::vector<std::vector<std::vector<IndexWeight>>> hierarchicalParents;
	std::vector<std::vector<std::vector<IndexWeight>>> hierarchicalSons;
	std::vector<std::vector<std::vector<int>>> hierarchyNeighbourhoods;
	std::vector<std::vector<std::vector<TriangleIndex>>> hierarchyTriangles;
	std::vector<std::vector<RigidTransformation>> hierarchicalTransformation;
	std::vector<std::vector<RigidTransformation>> hierarchicalFrames;
	std::vector<std::vector<bool>> hierarchicalFixedNodes;
	std::vector<float> hierarchicalSquaredEdgeLenght;
};

int HierarchicalModel::ConstructHierachy(const SimpleMesh & mesh, const std::vector<bool> & fixedVertices, const int levels, const Eigen::MatrixXf & geodesicDescriptor){

	numLevels = levels;

	int vCount = mesh.vertices.size();
	int tCount = mesh.triangles.size();

	std::vector<std::unordered_set<int>> _neighbours(vCount);
	double averageEdgeLenght = 0;
	for (int t = 0; t < tCount; t++){
		for (int k = 0; k < 3; k++){
			averageEdgeLenght += (mesh.vertices[mesh.triangles[t][(k + 1) % 3]] - mesh.vertices[mesh.triangles[t][k]]).norm();
			_neighbours[mesh.triangles[t][k]].insert(mesh.triangles[t][(k + 1) % 3]);
			_neighbours[mesh.triangles[t][k]].insert(mesh.triangles[t][(k + 2) % 3]);
		}
	}
	std::vector<std::vector<int>> neighbours(vCount);
	for (int i = 0; i < vCount; i++)neighbours[i] = std::vector<int>(_neighbours[i].begin(), _neighbours[i].end());

	averageEdgeLenght /= float(3 * mesh.triangles.size());
	
	printf("Average edge length %g\n", averageEdgeLenght);

	hierarchyFineIndices.resize(numLevels);
	hierarchyNeighbourhoods.resize(numLevels);
	hierarchicalParents.resize(numLevels);
	hierarchicalSons.resize(numLevels);
	hierarchyTriangles.resize(numLevels);
	hierarchicalTransformation.resize(numLevels);
	hierarchicalFrames.resize(numLevels);
	hierarchicalFixedNodes.resize(numLevels);

	for (int l = 0; l < numLevels; l++){
		printf("Level %d \n", l);
		float radius = l > 0 ? averageEdgeLenght * float(1 << l) : 0;
		if (!ConstructProlongation(mesh, fixedVertices, geodesicDescriptor, neighbours, radius, hierarchyFineIndices[l], hierarchicalParents[l], hierarchicalSons[l], hierarchyNeighbourhoods[l], hierarchyTriangles[l], hierarchicalFrames[l], hierarchicalTransformation[l], hierarchicalFixedNodes[l])){
			printf("Failed construction of level %d \n",l);
			return 0;
		}
	}

	hierarchyReferencePosition.resize(numLevels);
	for (int l = 0; l < numLevels; l++){
		hierarchyReferencePosition[l].resize(hierarchyFineIndices[l].size());
		for (int i = 0; i < hierarchyFineIndices[l].size(); i++){
			hierarchyReferencePosition[l][i] = mesh.vertices[hierarchyFineIndices[l][i]];
		}
	}

	hierarchicalSquaredEdgeLenght.resize(numLevels);
	printf("Squared edge lengths: \n");
	for (int l = 0; l < numLevels; l++){
		float cumSquaredEdgeLenght = 0;
		const std::vector<std::vector<int>> & neighbourhoods = hierarchyNeighbourhoods[l];
		const std::vector<Eigen::Vector3f> & nodeReferencePosition = hierarchyReferencePosition[l];
		for (int j = 0; j < neighbourhoods.size(); j++){
			Eigen::Vector3f refPosCurrent = nodeReferencePosition[j];
			for (int n = 0; n < neighbourhoods[j].size(); n++){
				int k = neighbourhoods[j][n];
				Eigen::Vector3f refPosNeighbour = nodeReferencePosition[k];
				cumSquaredEdgeLenght += (refPosCurrent - refPosNeighbour).norm();
			}
		}
		hierarchicalSquaredEdgeLenght[l] = cumSquaredEdgeLenght;
		printf("Level %d = %g \n", l, cumSquaredEdgeLenght);
	}

	return 1;
}

#endif //HIERARCHICAL_MODEL_INCLUDED
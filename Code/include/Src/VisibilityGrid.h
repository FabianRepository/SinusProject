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


#ifndef VISIBILITY_GRID_INCLUDED
#define VISIBILITY_GRID_INCLUDED
#include "IntersectableObject.h"

template<class T>
class Grid3D{
public:
	Grid3D(){}
	Grid3D(int p_x_dim, int p_y_dim, int p_z_dim){
		resize(p_x_dim, p_y_dim, p_z_dim);
	}
	void resize(int p_x_dim, int p_y_dim, int p_z_dim){
		x_dim = p_x_dim;
		y_dim = p_y_dim;
		z_dim = p_z_dim;
		xy_dim = x_dim * y_dim;
		xyz_dim = x_dim * y_dim *z_dim;
		data.resize(xyz_dim);
	}
	void clear(){
		data.clear();
	}
	T& operator () (int i[3]){ return data[xy_dim*i[2] + x_dim*i[1] + i[0]]; }
	T& operator () (int x, int y, int z){ return data[xy_dim*z + x_dim*y + x]; }
	const T& operator() (int i[3]) const { return data[xy_dim*i[2] + x_dim*i[1] + i[0]]; }
	const T& operator() (int x, int y, int z) const { return data[xy_dim*z + x_dim*y + x]; }
	std::vector<T> data;
	int x_dim, y_dim, z_dim, xy_dim, xyz_dim;
};


class BBox{
public:
	BBox(){}
	BBox(IntersectableObject * object, int p_res){
		res = p_res;
		Eigen::Vector3f minCorner = object->MinCorner();
		Eigen::Vector3f maxCorner = object->MaxCorner();
		Eigen::Vector3f diagonal = maxCorner - minCorner;
		double scale = std::max<double>(std::max<double>(diagonal[0], diagonal[1]), diagonal[2]);
		step = scale / double(res);
		Eigen::Vector3f center = (maxCorner + minCorner) / 2.0;
		bboxMin = center - Eigen::Vector3f(scale, scale, scale) / 2.0;
	}
	Eigen::Vector3f GetVoxelCenter(int i[3]);
	void GetClosestVoxel(const Eigen::Vector3f & p, int i[3]);
	Eigen::Vector3f bboxMin;
	int res;
	double step;
};

Eigen::Vector3f BBox::GetVoxelCenter(int i[3]){
	return bboxMin + (Eigen::Vector3f(i[0], i[1], i[2]) + Eigen::Vector3f(0.5, 0.5, 0.5)) * step;
}
void BBox::GetClosestVoxel(const Eigen::Vector3f & p, int i[3]){
	Eigen::Vector3f diff = ((p - bboxMin) / step) - Eigen::Vector3f(0.5, 0.5, 0.5);
	i[0] = std::min<int>(std::max<int>(0, int(round(diff[0]))), res - 1);
	i[1] = std::min<int>(std::max<int>(0, int(round(diff[1]))), res - 1);
	i[2] = std::min<int>(std::max<int>(0, int(round(diff[2]))), res - 1);
}

class VisibilityGrid{
public:
	VisibilityGrid(){}
	VisibilityGrid(IntersectableObject * p_object, int p_res){
		object = p_object;
		res = p_res;
		SetExtertiorVoxel();
		DiscardExternalVoxel(4);
		SetVisibilityEdges();
	}

	void SetExtertiorVoxel();
	void DiscardExternalVoxel(const int visibleCorners);
	void SetVisibilityEdges();

	bool ComputeShortesPath(int i0[3], int i1[3], std::vector<unsigned long> & path);
	bool FindClosestValidVoxel(int i0[3], int i1[3]);

	int res;
	IntersectableObject *object;
	Grid3D<int> isExterior;
	std::unordered_set<unsigned long long> occludedEdges;
	BBox bbox;
};

void VisibilityGrid::SetVisibilityEdges(){
	for (int l = 0; l < res; l++)for (int j = 0; j < res; j++)for (int k = 0; k < res; k++){
		int i[3] = { l, j, k };
		if (isExterior(i)){
			Eigen::Vector3f currentVoxelCenter = bbox.GetVoxelCenter(i);
			unsigned long currentVoxelKey = SetVoxelKey(i);
			for (int c = 0; c < 3; c++)for (int di = 0; di < 2; di++){
				int pi[3] = { i[0], i[1], i[2] };
				pi[c] = std::min<int>(std::max<int>(0, i[c] + 2 * di - 1), res - 1);
				if (isExterior(pi)){
					Eigen::Vector3f neighbourVoxelCenter = bbox.GetVoxelCenter(pi);
					Eigen::Vector3f psect;
					Eigen::Vector3f nsect;
					if (object->edgeIntersect(currentVoxelCenter, neighbourVoxelCenter, psect, nsect)){
						unsigned long neighbourVoxelKey = SetVoxelKey(pi);
						unsigned long long ocluddedEdgeKey = SetEdgeKey(currentVoxelKey, neighbourVoxelKey);
						occludedEdges.insert(ocluddedEdgeKey);
						ocluddedEdgeKey = SetEdgeKey(neighbourVoxelKey, currentVoxelKey);
						occludedEdges.insert(ocluddedEdgeKey);
					}
				}
			}
		}
	}
	printf("Number of occluded edges %d \n", occludedEdges.size());
}

bool VisibilityGrid::FindClosestValidVoxel(int i0[3], int i1[3]){
	if (isExterior(i0)){
		i1[0] = i0[0];
		i1[1] = i0[1];
		i1[2] = i0[2];
		return true;
	}

	unsigned long sourceKey = SetVoxelKey(i0);

	std::unordered_set<unsigned long> addedVoxels;
	addedVoxels.insert(sourceKey);

	std::queue<unsigned long> voxelQueue;
	voxelQueue.push(sourceKey);

	int maxIterValue = 11 * 11 * 11;
	int iter = 0;
	while (iter < maxIterValue){
		iter++;
		unsigned long voxelKey = voxelQueue.front();
		voxelQueue.pop();
		int i[3];
		GetVoxelIndices(voxelKey, i);
		for (int c = 0; c < 3; c++)for (int di = 0; di < 2; di++){
			int pi[3] = { i[0], i[1], i[2] };
			pi[c] = std::min<int>(std::max<int>(0, i[c] + 2 * di - 1), res - 1);
			if (isExterior(pi)){
				i1[0] = pi[0];
				i1[1] = pi[1];
				i1[2] = pi[2];
				return true;
			}
			else{
				unsigned long neighbourKey = SetVoxelKey(pi);
				if (addedVoxels.find(neighbourKey) == addedVoxels.end()){
					addedVoxels.insert(neighbourKey);
					voxelQueue.push(neighbourKey);
				}
			}
		}
	}
	return false;
}

void VisibilityGrid::SetExtertiorVoxel(){
	isExterior.clear();
	isExterior.resize(res, res, res);
	bbox = BBox(object, res);

	for (int i = 0; i < res; i++)for (int j = 0; j < res; j++)for (int k = 0; k < res; k++){
		int i3[3] = { i, j, k };
		Eigen::Vector3f voxelCenter = bbox.GetVoxelCenter(i3);
		isExterior(i3) = object->isInterior(voxelCenter) ? 0 : 1;
	}

	int exteriorCount = 0;
	for (int i = 0; i < res; i++)for (int j = 0; j < res; j++)for (int k = 0; k < res; k++) if (isExterior(i, j, k)) exteriorCount++;
	printf("Initial exterior voxels %d \n", exteriorCount);
}

void VisibilityGrid::DiscardExternalVoxel(const int visibleCorners){
	Eigen::Vector3f bboxCorners[8];
	for (int i = 0; i < 2; i++)for (int j = 0; j < 2; j++)for (int k = 0; k < 2; k++){
		bboxCorners[4 * i + 2 * j + k] = bbox.bboxMin + Eigen::Vector3f(bbox.step*i, bbox.step*j, bbox.step*k);
	}

	for (int i = 0; i < res; i++)for (int j = 0; j < res; j++)for (int k = 0; k < res; k++){
		int i3[3] = { i, j, k };
		if (isExterior(i3)){
			Eigen::Vector3f voxelCenter = bbox.GetVoxelCenter(i3);
			int count = 0;
			for (int s = 0; s < 8; s++){
				Eigen::Vector3f psect;
				Eigen::Vector3f nsect;
				if (!object->edgeIntersect(voxelCenter, bboxCorners[s], psect, nsect)){
					count++;
				}
				if (count >= visibleCorners){
					isExterior(i3) = 0;
					break;
				}
			}
		}
	}

	int exteriorCount = 0;
	for (int i = 0; i < res; i++)for (int j = 0; j < res; j++)for (int k = 0; k < res; k++) if (isExterior(i, j, k)) exteriorCount++;
	printf("Exterior voxels after discarding by corner visibility %d \n", exteriorCount);
}

class VoxelPropagationData{
public:
	unsigned long parentKey;
	bool fromSource;
	int depth;
};


bool VisibilityGrid::ComputeShortesPath(int i0[3], int i1[3], std::vector<unsigned long> & path){
	if (!isExterior(i0) || !isExterior(i1)){
		printf("Voxels are not exterior! \n");
		return false;
	}

	path.clear();

	unsigned long sourceKey = SetVoxelKey(i0);
	VoxelPropagationData source;
	source.fromSource = true;
	source.parentKey = -1;
	source.depth = 0;

	unsigned long targetKey = SetVoxelKey(i1);
	VoxelPropagationData target;
	target.fromSource = false;
	target.parentKey = -1;
	target.depth = 0;

	std::unordered_map<unsigned long, VoxelPropagationData> addedVoxels;
	addedVoxels[sourceKey] = source;
	addedVoxels[targetKey] = target;


	std::queue<unsigned long> voxelQueue;
	voxelQueue.push(sourceKey);
	voxelQueue.push(targetKey);

	int currentDepth = 0;
	int diameter = 3 * res;
	while (currentDepth < diameter){
		unsigned long voxelKey = voxelQueue.front();
		VoxelPropagationData & currentVoxel = addedVoxels[voxelKey];
		currentDepth = currentVoxel.depth;
		voxelQueue.pop();
		int i[3];
		GetVoxelIndices(voxelKey, i);
		for (int c = 0; c < 3; c++) for (int di = 0; di < 2; di++){
			int pi[3] = { i[0], i[1], i[2] };
			pi[c] = std::min<int>(std::max<int>(0, i[c] + 2 * di - 1), res - 1);
			if (isExterior(pi)){
				unsigned long neighbourKey = SetVoxelKey(pi);
				unsigned long long edgeKey = SetEdgeKey(voxelKey, neighbourKey);
				if (occludedEdges.find(edgeKey) == occludedEdges.end()){
					if (addedVoxels.find(neighbourKey) == addedVoxels.end()){
						VoxelPropagationData neighbour;
						neighbour.fromSource = currentVoxel.fromSource;
						neighbour.parentKey = voxelKey;
						neighbour.depth = currentDepth + 1;
						addedVoxels[neighbourKey] = neighbour;
						voxelQueue.push(neighbourKey);
					}
					else{
						VoxelPropagationData neighbour = addedVoxels[neighbourKey];
						if (neighbour.fromSource != currentVoxel.fromSource){
							std::vector<unsigned long> inverseCurrentPath;
							unsigned long pathKey = voxelKey;
							while (pathKey != -1){
								inverseCurrentPath.push_back(pathKey);
								pathKey = addedVoxels[pathKey].parentKey;
							}

							std::vector<unsigned long> inverseNeighbourPath;
							pathKey = neighbourKey;
							while (pathKey != -1){
								inverseNeighbourPath.push_back(pathKey);
								pathKey = addedVoxels[pathKey].parentKey;
							}

							const std::vector<unsigned long> & inverseFromSource = currentVoxel.fromSource ? inverseCurrentPath : inverseNeighbourPath;
							const std::vector<unsigned long> & inverseFromTarget = currentVoxel.fromSource ? inverseNeighbourPath : inverseCurrentPath;

							int pathSize = inverseFromSource.size() + inverseFromTarget.size();
							for (int l = inverseFromSource.size() - 1; l >= 0; l--) path.push_back(inverseFromSource[l]);
							for (int l = 0; l < inverseFromTarget.size(); l++) path.push_back(inverseFromTarget[l]);

							return true;
						}
					}
				}
			}

		}
	}

	return false;
}

#endif //VISIBILITY_GRID_INCLUDED
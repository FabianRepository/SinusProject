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

#include <Util/Geometry.h>
#include <Util/KDtree.h>
#include <Eigen/Dense>

class ExtrinsicVectorField{
public:
	virtual Eigen::Vector3f VectorField(const Eigen::Vector3f & p) = 0;
	virtual Eigen::Vector3f GetClosestPoint(const Eigen::Vector3f & p) = 0;
	virtual void GetClosestPoint(const Eigen::Vector3f & p, int & segementIndex, double & t) =0;
};

class CurveInducedField : public ExtrinsicVectorField{
public:
	CurveInducedField(std::vector<Eigen::Vector3f> & p_nodes);
	std::vector<Eigen::Vector3f> nodes;
	std::vector<double> normalizedLength;
	std::vector<Eigen::Vector3f> normalizedDirection;
	double curveLenght;
	//CKDTree queryKDTree;
	KDtree<float, 3> * queryKDTree;
	Eigen::Vector3f GetClosestPoint(const Eigen::Vector3f & p);
	void GetClosestPoint(const Eigen::Vector3f & p, int & segementIndex, double & t);
	double GetClosestTime(const Eigen::Vector3f & p);
	double GetTriangularDistance(const Eigen::Vector3f & p, bool toTarget = true);
	Eigen::Vector3f VectorField(const Eigen::Vector3f & p);
};

CurveInducedField::CurveInducedField(std::vector<Eigen::Vector3f> & p_nodes){
	nodes = p_nodes;

	normalizedLength.resize(nodes.size(), 0);
	for (int i = 1; i < nodes.size(); i++){
		double segmentLegnth = (nodes[i] - nodes[i - 1]).norm();
		normalizedLength[i] = normalizedLength[i-1] + segmentLegnth;
	}
	curveLenght = normalizedLength[nodes.size() - 1];
	for (int i = 0; i < nodes.size(); i++) normalizedLength[i] /= curveLenght;

	normalizedDirection.resize(nodes.size());
	normalizedDirection[0] = (nodes[1] - nodes[0]) / (nodes[1] - nodes[0]).norm();
	for (int i = 1; i < normalizedDirection.size() - 1; i++) normalizedDirection[i] = (nodes[i + 1] - nodes[i - 1]) / (nodes[i + 1] - nodes[i - 1]).norm();
	normalizedDirection[nodes.size() - 1] = (nodes[nodes.size() - 1] - nodes[nodes.size() - 2]) / ((nodes[nodes.size() - 1] - nodes[nodes.size() - 2])).norm();

	//CPointNDList kdTreePoints(3, static_cast<int>(nodes.size()));
	//for (int i = 0; i < nodes.size(); i++)
	//{
	//	kdTreePoints.m_pPtND[i].m_fCoords[0] = nodes[i][0];
	//	kdTreePoints.m_pPtND[i].m_fCoords[1] = nodes[i][1];
	//	kdTreePoints.m_pPtND[i].m_fCoords[2] = nodes[i][2];
	//	kdTreePoints.m_pPtND[i].m_nIdx = i;
	//}
	//queryKDTree.BuildKDTree(kdTreePoints, 100);

	queryKDTree = new KDtree<float, 3>(&nodes[0][0], nodes.size());
}

Eigen::Vector3f CurveInducedField::GetClosestPoint(const Eigen::Vector3f & p){
	int segmentIndex;
	double t;
	GetClosestPoint(p, segmentIndex, t);
	//printf("num nodes = %d, i = %d ,t = %f \n", nodes.size(), segmentIndex, t);
	Eigen::Vector3f closestPoint = nodes[segmentIndex] *(1.0 -  t) + nodes[segmentIndex + 1] * t;
	return closestPoint;
}

void CurveInducedField::GetClosestPoint(const Eigen::Vector3f & p, int & segementIndex, double & t){
	float queryPt[3] = { p[0], p[1], p[2] };
	//PointND closestNode;
	//queryKDTree.GetClosestPoint(queryPt, closestNode);
	//int closestNodeIndex = closestNode.m_nIdx;

	//{//For Debuggin!!
	//	float testPt[3] = { nodes[10][0], nodes[10][1], nodes[10][2] };
	//	const float * _closestPoint = queryKDTree->closest_to_pt(&testPt[0], FLT_MAX);
	//	uintptr_t closest_ptr = reinterpret_cast<std::uintptr_t>(_closestPoint);
	//	uintptr_t start_ptr = reinterpret_cast<std::uintptr_t>(&nodes[0][0]);
	//	unsigned int closestNodeIndex = (closest_ptr - start_ptr) / (3 * sizeof(float));
	//	printf("Here \n");
	//	printf("%d \n", closestNodeIndex);
	//	printf("%llu \n", closest_ptr);
	//	printf("%llu \n", start_ptr);
	//	printf("%f %f %f \n", nodes[10][0], nodes[10][1], nodes[10][2]);
	//	printf("%f %f %f \n", nodes[closestNodeIndex][0], nodes[closestNodeIndex][1], nodes[closestNodeIndex][2]);
	//}

	const float * _closestPoint = queryKDTree->closest_to_pt(&queryPt[0], FLT_MAX);
	uintptr_t closest_ptr = reinterpret_cast<std::uintptr_t>(_closestPoint);
	uintptr_t start_ptr = reinterpret_cast<std::uintptr_t>(&nodes[0][0]);
	unsigned int closestNodeIndex = (closest_ptr - start_ptr) / (3 * sizeof(float));

	Eigen::Vector3f ray = p - nodes[closestNodeIndex];
	double nodeDistance2 =(ray).squaredNorm();

	int _segmentIndex = closestNodeIndex < nodes.size() - 1? closestNodeIndex : nodes.size() - 2;
	double _t = closestNodeIndex < nodes.size() - 1 ? 0 : 1;

	double minDistance2 = nodeDistance2;

	if (closestNodeIndex > 0){//Check previous segment
		Eigen::Vector3f segmentDirection = nodes[closestNodeIndex - 1] - nodes[closestNodeIndex];
		double segmentLenght = (segmentDirection).norm();
		Eigen::Vector3f segmenNormal = segmentDirection / segmentLenght;
		double absoluteProjection = ray.dot(segmenNormal);
		double normalizedProjection = absoluteProjection / segmentLenght;

		if (normalizedProjection > 0 && normalizedProjection < 1.0){
			double orthogonalDistance2 = nodeDistance2 - absoluteProjection*absoluteProjection;
			if (orthogonalDistance2 < minDistance2){
				minDistance2 = orthogonalDistance2;
				_segmentIndex = closestNodeIndex - 1;
				_t = 1.0 - normalizedProjection;
			}
		}
	}

	if (closestNodeIndex < nodes.size() - 1){//Check next segment
		Eigen::Vector3f segmentDirection = nodes[closestNodeIndex + 1] - nodes[closestNodeIndex];
		double segmentLenght = (segmentDirection).norm();
		Eigen::Vector3f segmenNormal = segmentDirection / segmentLenght;
		double absoluteProjection = ray.dot(segmenNormal);
		double normalizedProjection = absoluteProjection / segmentLenght;

		if (normalizedProjection > 0 && normalizedProjection < 1.0){
			double orthogonalDistance2 = nodeDistance2 - absoluteProjection*absoluteProjection;
			if (orthogonalDistance2 < minDistance2){
				minDistance2 = orthogonalDistance2;
				_segmentIndex = closestNodeIndex;
				_t = normalizedProjection;
			}
		}
	}

	segementIndex = _segmentIndex;
	t = _t;
}

Eigen::Vector3f CurveInducedField::VectorField(const Eigen::Vector3f & p){
	//int segmentIndex;
	//double t;
	//GetClosestPoint(p, segmentIndex, t); 
	//Eigen::Vector3f direction = segmentIndex < nodes.size() - 1 ? normalizedDirection[segmentIndex] * (1.0 - t) + normalizedDirection[segmentIndex + 1] * t : normalizedDirection[nodes.size() - 1];
	//direction /= (direction).norm();
	//return direction;

	//Eigen::Vector3f closestPoint = nodes[segmentIndex] * (1.0 - t) + nodes[segmentIndex + 1] * t;
	//Eigen::Vector3f ray = closestPoint - p;
	//double rayLenght = ray.norm();
	//Eigen::Vector3f normalizedRay = ray / rayLenght;

	//double factor = std::min<double>(1.0, rayLenght*20.0 / curveLenght);
	//Eigen::Vector3f vector = normalizedRay*factor + normalizedDirection[segmentIndex] * (1.0 - factor);
	////printf("Ray lenght %f . Curve length %f \n", rayLenght, curveLenght);
	////printf("Axial weight %f : Directional weight %f \n", factor, 1.0 - factor);

	//return (vector / Eigen::Vector3f::Length(vector));

	int segmentIndex;
	double t;
	GetClosestPoint(p, segmentIndex, t);
	Eigen::Vector3f closestPoint = nodes[segmentIndex] * (1.0 - t) + nodes[segmentIndex + 1] * t;
	Eigen::Vector3f ray = closestPoint - p;
	double rayLenght = ray.norm();
	Eigen::Vector3f normalizedRay = ray / rayLenght;
	
	double factor = std::min<double>(1.0, rayLenght*5.0 / curveLenght);
	Eigen::Vector3f curveDirection = normalizedDirection[segmentIndex] * (1.0 - t) + normalizedDirection[segmentIndex + 1] * t;
	Eigen::Vector3f vector = normalizedRay*factor + curveDirection*(1.0 - factor);

	return (vector /vector.norm());
}

double CurveInducedField::GetClosestTime(const Eigen::Vector3f & p){
	int segmentIndex;
	double t;
	GetClosestPoint(p, segmentIndex, t);
	double closestTime = normalizedLength[segmentIndex] * (1.0 - t) + normalizedLength[segmentIndex + 1] * t;
	return closestTime;
}

double CurveInducedField::GetTriangularDistance(const Eigen::Vector3f & p, bool toTarget){
	int segmentIndex;
	double t;
	GetClosestPoint(p, segmentIndex, t);
	double closestTime = normalizedLength[segmentIndex] * (1.0 - t) + normalizedLength[segmentIndex + 1] * t;
	Eigen::Vector3f closestPoint = nodes[segmentIndex] * (1.0 - t) + nodes[segmentIndex + 1] * t;
	if (toTarget){
		return ((1.0 - closestTime) * curveLenght + (closestPoint - p).norm());
	}
	else{
		return (closestTime * curveLenght + (closestPoint - p).norm());
	}
}
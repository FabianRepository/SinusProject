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

#ifndef DEFORMATION_SOLVER_INCLUDED
#define DEFORMATION_SOLVER_INCLUDED


#include <ceres/ceres.h>
#include <ceres/rotation.h>
#include"HierarchicalModel.h"
#include <omp.h>
#include <ctime>

const int numVariablesPerNode = 6;
//
// Regularity Costs
//
struct RigidityCostFunctor
{
	RigidityCostFunctor(double tX, double tY, double tZ)
	: m_tX(tX), m_tY(tY), m_tZ(tZ)
	{}

	template <typename T>
	bool operator() (const T* const R1, const T* const T1, const T* const T2, T* residual) const {
		T point[3] = { T(m_tX), T(m_tY), T(m_tZ) };
		T rotatedPoint[3];
		ceres::AngleAxisRotatePoint(R1, point, rotatedPoint);
		residual[0] = rotatedPoint[0] + T1[0] - T2[0] - T(m_tX);
		residual[1] = rotatedPoint[1] + T1[1] - T2[1] - T(m_tY);
		residual[2] = rotatedPoint[2] + T1[2] - T2[2] - T(m_tZ);
		return true;
	}

private:
	const double m_tX;
	const double m_tY;
	const double m_tZ;
};

void SetRigidityResiduals(const std::vector<std::vector<int>> & neighbourhoods, const std::vector<Eigen::Vector3f> & nodeReferencePosition, double rigidityWeight, Eigen::VectorXd &vars, ceres::Problem &problem)
{
	for (int j = 0; j < neighbourhoods.size(); j++)
	{
		Eigen::Vector3f refPosCurrent = nodeReferencePosition[j];
		for (int n = 0; n < neighbourhoods[j].size(); n++)
		{
			int k = neighbourhoods[j][n];
			int jVariableOffset = int(j*numVariablesPerNode);
			int kVariableOffset = int(k*numVariablesPerNode);
			Eigen::Vector3f diff = nodeReferencePosition[k] - refPosCurrent;
			ceres::CostFunction* cost_function = new ceres::AutoDiffCostFunction<RigidityCostFunctor, 3 /*#residuals*/, 3 /*rotate*/,3 /*jTranslate*/, 3/*kTranslate*/>(new RigidityCostFunctor(diff[0], diff[1], diff[2]));
			problem.AddResidualBlock(cost_function, new ceres::ScaledLoss(NULL, rigidityWeight, ceres::TAKE_OWNERSHIP), &(vars[jVariableOffset]) /*rotate*/, &(vars[jVariableOffset + 3]) /*jTranslate*/, &(vars[kVariableOffset + 3]) /*kTranslate*/);
		}
	}
}

struct FixedNodeRotationCostFunctor
{
	FixedNodeRotationCostFunctor(){}

	template <typename T>
	bool operator() (const T* const R, T* residual) const {
		residual[0] = R[0];
		residual[1] = R[1];
		residual[2] = R[2];
		return true;
	}
};

void SetFixedNodeRotationResiduals(const std::vector<bool> & fixedNodes, double fixedWeight, Eigen::VectorXd &vars, ceres::Problem &problem)
{
	for (int j = 0; j < fixedNodes.size(); j++)
	{
		if (fixedNodes[j]){
			int jVariableOffset = int(j*numVariablesPerNode);
			ceres::CostFunction* cost_function = new ceres::AutoDiffCostFunction<FixedNodeRotationCostFunctor, 3 /*#residuals*/, 3 /*rotate*/>(new FixedNodeRotationCostFunctor());
			problem.AddResidualBlock(cost_function, new ceres::ScaledLoss(NULL, fixedWeight, ceres::TAKE_OWNERSHIP), &(vars[jVariableOffset]) /*rotate*/);
		}
	}
}


struct FixedNodeTranslationCostFunctor
{
	FixedNodeTranslationCostFunctor(){}

	template <typename T>
	bool operator() (const T* const R, T* residual) const {
		residual[0] = R[0];
		residual[1] = R[1];
		residual[2] = R[2];
		return true;
	}
};

void SetFixedNodeTranslationResiduals(const std::vector<bool> & fixedNodes, double fixedWeight, Eigen::VectorXd &vars, ceres::Problem &problem)
{
	for (int j = 0; j < fixedNodes.size(); j++)
	{
		if (fixedNodes[j]){
			int jVariableOffset = int(j*numVariablesPerNode);
			ceres::CostFunction* cost_function = new ceres::AutoDiffCostFunction<FixedNodeTranslationCostFunctor, 3 /*#residuals*/, 3 /*translation*/>(new FixedNodeTranslationCostFunctor());
			problem.AddResidualBlock(cost_function, new ceres::ScaledLoss(NULL, fixedWeight, ceres::TAKE_OWNERSHIP), &(vars[jVariableOffset + 3]) /*translation*/);
		}
	}
}

struct FitCostFunctor
{
	typedef
	ceres::DynamicAutoDiffCostFunction < FitCostFunctor, 12 /*batch size of derivative evaluation*/ >
	CostFunction;

	FitCostFunctor(
		const int numNodes,
		const std::vector<Eigen::Vector3f> &diff,
		const std::vector<Eigen::Vector3f> &nodePos,
		const std::vector<double> &nodeWeight,
		const double fittingWeight,
		const Eigen::Vector3f &cpPos,
		const Eigen::Vector3f &cpNormal)
		: m_numNodes(numNodes)
		, m_diff(diff)
		, m_nodePos(nodePos)
		, m_nodeWeight(nodeWeight)
		, m_fittinWeight(fittingWeight)
		, m_cpPos(cpPos)
		, m_cpNormal(cpNormal)
	{}

	template <typename T>
	bool operator() (T const* const* RT, T* residual) const
	{
		T pX = -T(m_cpPos[0]);
		T pY = -T(m_cpPos[1]);
		T pZ = -T(m_cpPos[2]);
		for (int i = 0; i < m_numNodes; i++)
		{
			int ri = i * 2;     // start index of rotation
			int ti = i * 2 + 1; // start index of translation
			T point[3] = { T(m_diff[i][0]), T(m_diff[i][1]), T(m_diff[i][2]) };
			T rotatedPoint[3];
			ceres::AngleAxisRotatePoint(RT[ri], point, rotatedPoint);
			pX += ((rotatedPoint[0]) + RT[ti][0] + T(m_nodePos[i][0])) * T(m_nodeWeight[i]);
			pY += ((rotatedPoint[1]) + RT[ti][1] + T(m_nodePos[i][1])) * T(m_nodeWeight[i]);
			pZ += ((rotatedPoint[2]) + RT[ti][2] + T(m_nodePos[i][2])) * T(m_nodeWeight[i]);
		}
		// Point Constraints
		residual[0] = pX * 0.316 * m_fittinWeight; // Following [Li et al. 2009], we weight the point-to-point constraint 10% as important as the point-to-normal constraint
		residual[1] = pY * 0.316 * m_fittinWeight; // i.e., sqrt(0.1) ~= 0.316. (taking square root here since Ceres wants the "signed cost", not "energy")
		residual[2] = pZ * 0.316 * m_fittinWeight;
		// Normal Constraint
		residual[3] = (pX * T(m_cpNormal[0])) + (pY * T(m_cpNormal[1])) + (pZ * T(m_cpNormal[2]))*m_fittinWeight;
		return true;
	}

	static const int c_numResiduals = 4; // i.e., 3 for point constraints and 1 for normal 

private:
	const int m_numNodes;                   // Number of neighboring nodes (usually between 2~5)
	const std::vector<Eigen::Vector3f> m_diff;     // Difference vectors between the current point and its neighboring nodes
	const std::vector<Eigen::Vector3f> m_nodePos;  // Positions of the neighboring nodes
	const std::vector<double> m_nodeWeight; // Weights of the neighboring nodes
	const double m_fittinWeight;			//Fitting weight
	const Eigen::Vector3f m_cpPos;                 // Position of the closest point
	const Eigen::Vector3f m_cpNormal;              // Normal of the closest point
};

void SetFitResiduals(
	const std::vector<Eigen::Vector3f>          &nodeReferencePosition,
	const std::vector<std::vector<IndexWeight>> &vertexParents,
	const std::vector<double>					& fittingWeight,
	const std::vector<int>						&targetVertexIndex,
	const std::vector<Eigen::Vector3f>			&targetVertexPosition,
	const std::vector<Eigen::Vector3f>			&referenceVertexPosition,
	const std::vector<Eigen::Vector3f>			&targetVertexNormal,
	Eigen::VectorXd								&vars,
	ceres::Problem								&problem)
{
	for (int k = 0; k < targetVertexIndex.size(); k++)
	{
		const int fineIndex = targetVertexIndex[k];
		const std::vector<IndexWeight> & parents = vertexParents[fineIndex];
		const Eigen::Vector3f  referencePosition = referenceVertexPosition[k];

		int numNodes = int(parents.size());
		std::vector<double*> parameter_blocks(numNodes * 2); // rotation and translation
		std::vector<Eigen::Vector3f> diff(numNodes);
		std::vector<Eigen::Vector3f> nPos(numNodes);
		std::vector<double> nWeight(numNodes);

		for (int i = 0; i < numNodes; i++)
		{
			int   nodeIdx = parents[i].index;
			double nodeWeight = parents[i].weight;
			diff[i] = referencePosition - nodeReferencePosition[nodeIdx];
			int baseIdx = nodeIdx*numVariablesPerNode;
			parameter_blocks[i * 2] = &(vars[baseIdx]);         // rotation
			parameter_blocks[i * 2 + 1] = &(vars[baseIdx + 3]); // translation
			nPos[i] = nodeReferencePosition[nodeIdx];
			nWeight[i] = nodeWeight;
		}

		FitCostFunctor* fitCostFunctor = new FitCostFunctor(numNodes, diff, nPos, nWeight, fittingWeight[k], targetVertexPosition[k], targetVertexNormal[k]);
		FitCostFunctor::CostFunction* cost_function = new FitCostFunctor::CostFunction(fitCostFunctor);

		for (int i = 0; i < numNodes; i++)
		{
			cost_function->AddParameterBlock(3); //for rotations
			cost_function->AddParameterBlock(3); //for translations
		}

		cost_function->SetNumResiduals(FitCostFunctor::c_numResiduals);

		problem.AddResidualBlock(cost_function, NULL, parameter_blocks);
	}
}

double ComputeRigidityError(const std::vector<std::vector<int>> & neighbourhoods, const std::vector<Eigen::Vector3f> & nodeReferencePosition, double rigidityWeight,const Eigen::VectorXd &vars){

	double cumRigiditySquaredError = 0;

	for (int j = 0; j < neighbourhoods.size(); j++)
	{
		Eigen::Vector3f refPosCurrent = nodeReferencePosition[j];
		for (int n = 0; n < neighbourhoods[j].size(); n++)
		{
			int k = neighbourhoods[j][n];
			Eigen::Vector3f refPosNeighbour = nodeReferencePosition[k];

			int jVariableOffset = int(j*numVariablesPerNode);
			int kVariableOffset = int(k*numVariablesPerNode);
			Eigen::Vector3f angleAxisj(vars[jVariableOffset], vars[jVariableOffset + 1], vars[jVariableOffset + 2]);
			Eigen::Matrix3f rotationj;
			ceres::AngleAxisToRotationMatrix((float *)&angleAxisj, (float *)&rotationj);
			Eigen::Vector3f translationj(vars[jVariableOffset + 3], vars[jVariableOffset + 4], vars[jVariableOffset + 5]);
			Eigen::Vector3f translationk(vars[kVariableOffset + 3], vars[kVariableOffset + 4], vars[kVariableOffset + 5]);

			Eigen::Vector3f residual = rotationj*(refPosNeighbour - refPosCurrent) + translationj + refPosCurrent - (translationk + refPosNeighbour);
			cumRigiditySquaredError += residual.squaredNorm();
		}
	}

	return 0.5*cumRigiditySquaredError*rigidityWeight;
}

double ComputeFitError(const std::vector<Eigen::Vector3f> &nodeReferencePosition, const std::vector<std::vector<IndexWeight>> &vertexParents, const std::vector<double> & fittingWeight, const std::vector<int> &targetVertexIndex,
	const std::vector<Eigen::Vector3f> &targetVertexPosition, const std::vector<Eigen::Vector3f> &referenceVertexPosition, const std::vector<Eigen::Vector3f> &targetVertexNormal, const Eigen::VectorXd &vars)
{
	double cumDistanceSquaredError = 0;
	double cumNormalSquaredError = 0;
	for (int k = 0; k < targetVertexIndex.size(); k++)
	{
		const int fineIndex = targetVertexIndex[k];
		const std::vector<IndexWeight> & parents = vertexParents[fineIndex];
		const Eigen::Vector3f  referencePosition = referenceVertexPosition[k];

		Eigen::Vector3f predictedPosition = Eigen::Vector3f::Zero();

		for (int i = 0; i < parents.size(); i++)
		{
			int variableOffset = parents[i].index*numVariablesPerNode;
			Eigen::Vector3f angleAxis(vars[variableOffset], vars[variableOffset + 1], vars[variableOffset + 2]);
			Eigen::Matrix3f rotation;
			ceres::AngleAxisToRotationMatrix((float *)&angleAxis, (float *)&rotation);
			Eigen::Vector3f translation(vars[variableOffset + 3], vars[variableOffset + 4], vars[variableOffset + 5]);

			Eigen::Vector3f  refNodePosition = nodeReferencePosition[parents[i].index];
			predictedPosition += (rotation*(referencePosition - refNodePosition) + translation + refNodePosition)*parents[i].weight;
		}
		Eigen::Vector3f residual = predictedPosition - targetVertexPosition[k];
		cumDistanceSquaredError += residual.squaredNorm() * 0.316 * 0.316;
		double normalResidual = residual.dot(targetVertexNormal[k]);
		cumNormalSquaredError += normalResidual*normalResidual;
		cumDistanceSquaredError *= (fittingWeight[k] * fittingWeight[k]);
		cumNormalSquaredError *= (fittingWeight[k] * fittingWeight[k]);
	}
	return 0.5*(cumDistanceSquaredError + cumNormalSquaredError);
}


double ComputeFixedNodeError(const std::vector<bool> & fixedNodes, double fixedWeight, const Eigen::VectorXd &vars){
	double cumSquaredError = 0;
	for (int j = 0; j < fixedNodes.size(); j++)
	{
		if (fixedNodes[j]){
			int variableOffset = int(j*numVariablesPerNode);
			for (int k = 0; k < numVariablesPerNode; k++) cumSquaredError += (vars[variableOffset + k] * vars[variableOffset + k]);
		}
	}
	return 0.5*cumSquaredError*fixedWeight;
}

void SolveDeformation(const std::vector<bool> & fixedNodes, const std::vector<Eigen::Vector3f> & nodeReferencePosition, const std::vector<std::vector<int>> & neighbourhoods, const std::vector<std::vector<IndexWeight>> &vertexParents, std::vector<RigidTransformation> & transformations,
	const double fixedWeight, const double rigidityWeight, const std::vector<double> & fittingWeight, const std::vector<int> &targetVertexIndex, const std::vector<Eigen::Vector3f> &targetVertexPosition, const std::vector<Eigen::Vector3f>	&referenceVertexPosition, const std::vector<Eigen::Vector3f> &targetVertexNormal, bool verbose = false){


	int numNodes = nodeReferencePosition.size();
	if (verbose)printf("Num nodes %d \n", numNodes);

	ceres::Problem problem;

	Eigen::VectorXd vars;
	vars.resize(numVariablesPerNode*numNodes);

	for (int i = 0; i < transformations.size(); i++){
		Eigen::Vector3f angleAxis;
		Eigen::Matrix3f rotation = transformations[i].rotation; 
		Eigen::Vector3f translation = transformations[i].translation;
		ceres::RotationMatrixToAngleAxis((float *)&rotation, (float *)&angleAxis);
		
		vars[i*numVariablesPerNode] = angleAxis[0];
		vars[i*numVariablesPerNode + 1] = angleAxis[1];
		vars[i*numVariablesPerNode + 2] = angleAxis[2];

		vars[i*numVariablesPerNode + 3] = translation[0];
		vars[i*numVariablesPerNode + 4] = translation[1];
		vars[i*numVariablesPerNode + 5] = translation[2];
	}
	
	if (verbose){
		double rigidityError = ComputeRigidityError(neighbourhoods, nodeReferencePosition, rigidityWeight, vars);
		double fittingError = ComputeFitError(nodeReferencePosition, vertexParents, fittingWeight, targetVertexIndex, targetVertexPosition, referenceVertexPosition, targetVertexNormal, vars);
		double fixedError = ComputeFixedNodeError(fixedNodes, fixedWeight, vars);
		printf("Initial Error %g : Rigid %g  - Fit %g - Fixed %g \n", rigidityError + fittingError + fixedError, rigidityError, fittingError, fixedError);
	}

	if (fixedNodes.size() != transformations.size()){
		printf("Unexpected node count! \n");
	}

	SetRigidityResiduals(neighbourhoods, nodeReferencePosition, rigidityWeight, vars, problem);
	SetFixedNodeRotationResiduals(fixedNodes, fixedWeight, vars, problem);
	SetFixedNodeTranslationResiduals(fixedNodes, fixedWeight, vars, problem);
	SetFitResiduals(nodeReferencePosition, vertexParents, fittingWeight, targetVertexIndex, targetVertexPosition, referenceVertexPosition, targetVertexNormal, vars, problem);
	

	ceres::Solver::Options options;
	options.linear_solver_type = ceres::SPARSE_NORMAL_CHOLESKY;
	options.num_threads = omp_get_num_procs(); // Use OpenMP to parallelize
	options.num_linear_solver_threads = omp_get_num_procs(); // Use OpenMP to parallelize
	options.minimizer_progress_to_stdout = false;
	options.function_tolerance = 5e-3;

	ceres::Solver::Summary summary;
	ceres::Solve(options, &problem, &summary);

	std::string report = summary.BriefReport();
	const char * _report = report.c_str();
	if (verbose)printf("Ceres Summary : %s \n", _report);

	if (verbose){
		double rigidityError = ComputeRigidityError(neighbourhoods, nodeReferencePosition, rigidityWeight, vars);
		double fittingError = ComputeFitError(nodeReferencePosition, vertexParents, fittingWeight, targetVertexIndex, targetVertexPosition, referenceVertexPosition, targetVertexNormal, vars);
		double fixedError = ComputeFixedNodeError(fixedNodes, fixedWeight, vars);
		printf("Initial Error %g : Rigid %g  - Fit %g - Fixed %g \n", rigidityError + fittingError + fixedError, rigidityError, fittingError, fixedError);
	}

	for (int i = 0; i < transformations.size(); i++){
		Eigen::Vector3f angleAxis(vars[i*numVariablesPerNode], vars[i*numVariablesPerNode + 1], vars[i*numVariablesPerNode + 2]);
		Eigen::Matrix3f rotation;
		Eigen::Vector3f translation(vars[i*numVariablesPerNode + 3], vars[i*numVariablesPerNode + 4], vars[i*numVariablesPerNode + 5]);
		ceres::AngleAxisToRotationMatrix((float *)&angleAxis,(float *)&rotation);
		transformations[i].rotation = rotation;
		transformations[i].translation = translation;
	}
}

void RigidAlignment(const std::vector<Eigen::Vector3f> & sourcePts, const std::vector<Eigen::Vector3f> & targetPts, const std::vector<float> & weights, Eigen::Vector3f & optimalTranslation, Eigen::Matrix3f & optimalRotation)
{
	if (sourcePts.size() != targetPts.size()) printf("Target points does not match source points! \n");

	int vCount = (int)sourcePts.size();
	Eigen::Vector3f sCentroid(0.f, 0.f, 0.f);
	Eigen::Vector3f tCentroid(0.f, 0.f, 0.f);
	float cumWeight = 0.f;
	for (int i = 0; i < vCount; i++){
		float weight = weights[i];
		sCentroid += sourcePts[i] * weight;
		tCentroid += targetPts[i] * weight;
		cumWeight += weight;
	}

	sCentroid /= cumWeight;
	tCentroid /= cumWeight;

	Eigen::Matrix3f varianceMatrix = Eigen::Matrix3f::Zero();
	for (int i = 0; i < vCount; i++){
		Eigen::Vector3f sourceVector = sourcePts[i] - sCentroid;
		Eigen::Vector3f targetVector = targetPts[i] - tCentroid;
		for (int k = 0; k < 3; k++) for (int l = 0; l < 3; l++) varianceMatrix(k,l) += weights[i] * sourceVector[k] * targetVector[l];
	}

	Eigen::JacobiSVD<Eigen::Matrix3f> _svd(varianceMatrix, Eigen::ComputeFullU | Eigen::ComputeFullV);
	Eigen::Matrix3f U = _svd.matrixU();
	Eigen::Matrix3f V = _svd.matrixV();

	if (U.determinant()*V.determinant() < 0.f) for (int k = 0; k < 3; k++) U(k, 2) *= -1.f;
	Eigen::Matrix3f U_transpose = U.transpose();
	optimalRotation = V*U_transpose;
	optimalTranslation = tCentroid - optimalRotation*sCentroid;
}

//TODO: Make sure each node has at least 4 sons to ensure rigid transformation is well defined?

void ProjectDeformation(const std::vector<Eigen::Vector3f> & vertexReferencePositions, const std::vector<Eigen::Vector3f> & vertexDeformedPositions, const std::vector<Eigen::Vector3f> & nodeReferencePositions, const std::vector<std::vector<IndexWeight>> &nodeSons, std::vector<RigidTransformation> & transformations){
	
	int threads = omp_get_num_procs();
#pragma omp parallel for num_threads( threads )
	for (int i = 0; i < nodeSons.size(); i++){
		const Eigen::Vector3f nodePosition = nodeReferencePositions[i];
		const std::vector<IndexWeight> & currentSons = nodeSons[i];
		std::vector<Eigen::Vector3f> sonReferencePositions(currentSons.size());
		std::vector<Eigen::Vector3f> sonTargetPositions(currentSons.size());
		std::vector<float> sonsWeights(currentSons.size());
		for (int j = 0; j < currentSons.size(); j++){
			sonReferencePositions[j] = vertexReferencePositions[currentSons[j].index] - nodePosition;
			sonTargetPositions[j] = vertexDeformedPositions[currentSons[j].index] - nodePosition;
			sonsWeights[j] = currentSons[j].weight;
		}
		RigidAlignment(sonReferencePositions, sonTargetPositions, sonsWeights, transformations[i].translation, transformations[i].rotation);
	}
}

void HierarchicalDeformationSolver(HierarchicalModel & hierarchy, const double fixedWeight, const double rigidityWeight, const std::vector<double> & fittingWeight, std::vector<Eigen::Vector3f> &vertexPosition, const std::vector<int> &targetVertexIndex, const std::vector<Eigen::Vector3f> &targetVertexPosition, const std::vector<Eigen::Vector3f> &referenceVertexPosition, const std::vector<Eigen::Vector3f> &targetVertexNormal, const int finestLevelSoved, bool verbose = false){

	std::vector<Eigen::Vector3f> transformedPositions = vertexPosition;
	double cummulativeTime = 0;
	for (int l = hierarchy.numLevels - 1; l >= finestLevelSoved; l--){
		if (verbose)printf("Level %d ... \n", l);
		clock_t startTimer = clock();
		ProjectDeformation(hierarchy.hierarchyReferencePosition[0], transformedPositions, hierarchy.hierarchyReferencePosition[l], hierarchy.hierarchicalSons[l], hierarchy.hierarchicalTransformation[l]);
		SolveDeformation(hierarchy.hierarchicalFixedNodes[l], hierarchy.hierarchyReferencePosition[l], hierarchy.hierarchyNeighbourhoods[l], hierarchy.hierarchicalParents[l], hierarchy.hierarchicalTransformation[l], fixedWeight, rigidityWeight / hierarchy.hierarchicalSquaredEdgeLenght[l], fittingWeight, targetVertexIndex, targetVertexPosition, referenceVertexPosition, targetVertexNormal, verbose);
		TransformModel(hierarchy.hierarchyReferencePosition[0], hierarchy.hierarchicalTransformation[l], hierarchy.hierarchicalParents[l], hierarchy.hierarchyFineIndices[l], transformedPositions);
		double partialTime = clock() - startTimer;
		cummulativeTime += partialTime;
		if (verbose)printf("Time %.4f \n", double(partialTime) / CLOCKS_PER_SEC);
	}
	vertexPosition = transformedPositions;
	if (verbose)printf("Total Time %.4f \n", double(cummulativeTime) / CLOCKS_PER_SEC);
}

#endif //DEFORMATION_SOLVER_INCLUDED
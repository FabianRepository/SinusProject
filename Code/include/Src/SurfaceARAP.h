#ifndef SURFACE_ARAP_INCLUDED
#define SURFACE_ARAP_INCLUDED
#include "Mesh.h"
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/StdVector>
#include <Eigen/Sparse>
#include <Eigen/SVD>
#include <omp.h>


Eigen::Matrix2f TriMetricTensor(const Eigen::Vector3f v0, const  Eigen::Vector3f v1, const  Eigen::Vector3f v2){
	Eigen::Vector3f d[2];
	d[0] = v1 - v0;
	d[1] = v2 - v0;
	Eigen::Matrix2f g;
	for (int i = 0; i < 2; i++)for (int j = 0; j < 2; j++){
		g(i, j) = d[i].dot(d[j]);
	}
	return g;
}

Eigen::Matrix3f TriStiffnesMatrix(const Eigen::Matrix2f & g){
	Eigen::Vector2f d[3];
	d[0] = Eigen::Vector2f(-1.0, -1.0);
	d[1] = Eigen::Vector2f(1.0, 0.0);
	d[2] = Eigen::Vector2f(0.0, 1.0);
	Eigen::Matrix3f s;
	Eigen::Matrix2f  g_inv = g.inverse();
	double triArea = sqrt(g.determinant());
	for (int i = 0; i < 3; i++)for (int j = 0; j < 3; j++) s(i, j) = -(d[i].dot(g_inv*d[j]))*triArea;
	return s;
}

void TriangleCotangentWeights(const std::vector<Eigen::Vector3f> & vertices, const std::vector<std::vector<TriangleIndex>> & triangles, std::vector<std::vector<Eigen::Vector3f>> & cotangentWeights){
	cotangentWeights.resize(triangles.size());
	for (int t = 0; t < triangles.size(); t++){
		cotangentWeights[t].resize(triangles[t].size());
		for (int n = 0; n < triangles[t].size(); n++){
			Eigen::Matrix2f g = TriMetricTensor(vertices[triangles[t][n][0]], vertices[triangles[t][n][1]], vertices[triangles[t][n][2]]);
			if (g.determinant() <= 0.0){
				printf("Non Positive Area Triangle!. %f\n", g.determinant());
				printf("Mass matrix assigned to 0.000001*Id\n");
				g(0, 0) = g(1, 1) = 0.000001;
				g(0, 1) = g(1, 0) = 0.0;
			}
			Eigen::Matrix3f s = TriStiffnesMatrix(g);
			for (int i = 0; i < 3; i++) cotangentWeights[t][n][i] = s(i, (i + 1) % 3);
		}
	}
}

class ARAPModel{
public:
	void Initialize(const std::vector<Eigen::Vector3f> & p_referenceVertices, const std::vector<std::vector<TriangleIndex>> & p_triangles, const std::vector<int> & p_fixedIndices, const float & p_softScale);
	void UpdteSoftWeight(const float & p_softScale);
	void Solve(const std::vector<Eigen::Vector3f> & softConstraints, std::vector<Eigen::Vector3f> & currentVertices);
	double Energy(const std::vector<Eigen::Vector3f> & currentVertices, std::vector<double> & energy);

	std::vector<Eigen::Vector3f> referenceVertices;
	std::vector<std::vector<TriangleIndex>> triangles;
	std::vector<std::vector<Eigen::Vector3f>> cotangentWeights;
	std::vector<float> softWeights;
	std::vector<int> variableIndex;
	std::vector<int> fixedIndices;

	Eigen::SimplicialLDLT<Eigen::SparseMatrix<float>> ARAPCholesky;
	Eigen::SparseMatrix<float> stiffnessMatrix;
	int freeVarCount;
};

double ARAPModel::Energy(const std::vector<Eigen::Vector3f> & currentVertices, std::vector<double> & energy){
	energy.resize(triangles.size());
	double cumEnergy = 0;
	for (int t = 0; t < triangles.size(); t++){
		Eigen::Matrix3f scatter = Eigen::Matrix3f::Zero();
		for (int n = 0; n < triangles[t].size(); n++){
			for (int i = 0; i < 3; i++){
				Eigen::Vector3f restEdge = referenceVertices[triangles[t][n][(i + 1) % 3]] - referenceVertices[triangles[t][n][i]];
				Eigen::Vector3f newEdge = currentVertices[triangles[t][n][(i + 1) % 3]] - currentVertices[triangles[t][n][i]];
				scatter += cotangentWeights[t][n][i] * restEdge*newEdge.transpose();
			}
		}
		Eigen::JacobiSVD<Eigen::Matrix3f> mSVD(scatter, Eigen::ComputeFullU | Eigen::ComputeFullV);
		Eigen::Matrix3f U = mSVD.matrixU();
		Eigen::Matrix3f V = mSVD.matrixV();
		if (U.determinant()*V.determinant() < 0.0) for (int k = 0; k < 3; k++) U(k, 2) *= -1.0;
		Eigen::Matrix3f rotation = V*U.transpose();

		double cumNeighbourhoodEnergy = 0;

		for (int n = 0; n < triangles[t].size(); n++){
			for (int j = 0; j < 3; j++){
				Eigen::Vector3f newEdge = currentVertices[triangles[t][n][j]] - currentVertices[triangles[t][n][(j + 1) % 3]];
				Eigen::Vector3f restEdge = referenceVertices[triangles[t][n][j]] - referenceVertices[triangles[t][n][(j + 1) % 3]];
				Eigen::Vector3f rotatedEdge = rotation * restEdge;
				cumNeighbourhoodEnergy += (newEdge - rotatedEdge).squaredNorm() * cotangentWeights[t][n][j];
			}
		}

		if (cumNeighbourhoodEnergy <0) printf("Negative neighbourhood Energy!!\n");

		energy[t] = cumNeighbourhoodEnergy;
		cumEnergy += cumNeighbourhoodEnergy;
	}
	return cumEnergy;
}

void ARAPModel::Solve(const std::vector<Eigen::Vector3f> & softConstraints, std::vector<Eigen::Vector3f> & currentVertices){

		Eigen::MatrixXf rhs = Eigen::MatrixXf::Zero(freeVarCount, 3);
		for (int t = 0; t < triangles.size(); t++){
			Eigen::Matrix3f scatter = Eigen::Matrix3f::Zero();
			for (int n = 0; n < triangles[t].size(); n++){
				for (int i = 0; i < 3; i++){
					Eigen::Vector3f restEdge = referenceVertices[triangles[t][n][(i + 1) % 3]] - referenceVertices[triangles[t][n][i]];
					Eigen::Vector3f newEdge = currentVertices[triangles[t][n][(i + 1) % 3]] - currentVertices[triangles[t][n][i]];
					scatter += cotangentWeights[t][n][i]*restEdge*newEdge.transpose();
				}
			}
			Eigen::JacobiSVD<Eigen::Matrix3f> mSVD(scatter, Eigen::ComputeFullU | Eigen::ComputeFullV);
			Eigen::Matrix3f U = mSVD.matrixU();
			Eigen::Matrix3f V = mSVD.matrixV();
			if (U.determinant()*V.determinant() < 0.0) for (int k = 0; k < 3; k++) U(k, 2) *= -1.0;
			Eigen::Matrix3f rotation = V*U.transpose();

			for (int n = 0; n < triangles[t].size(); n++){
				for (int i = 0; i < 3; i++){
					for (int j = 0; j < 3; j++){
						int vP = variableIndex[triangles[t][n][j]];
						int vN = variableIndex[triangles[t][n][(j + 1) % 3]];
						Eigen::Vector3f restEdge = referenceVertices[triangles[t][n][j]] - referenceVertices[triangles[t][n][(j + 1) % 3]];
						Eigen::Vector3f rotatedEdge = rotation * restEdge * cotangentWeights[t][n][j];
						if (vP != -1){
							rhs.row(vP) += rotatedEdge;
							if (vN == -1) rhs.row(vP) += currentVertices[triangles[t][n][(j + 1) % 3]] * cotangentWeights[t][n][j];
						}
						if (vN != -1){
							rhs.row(vN) -= rotatedEdge;
							if (vP == -1) rhs.row(vN) += currentVertices[triangles[t][n][j]] * cotangentWeights[t][n][j];
						}
					}
				}
			}
		}

		int threads = omp_get_num_procs();

#pragma omp parallel for num_threads( threads )
		for (int i = 0; i < referenceVertices.size(); i++){
			int vi = variableIndex[i];
			if(vi != -1) rhs.row(vi) += softConstraints[i] * softWeights[i];
		}

		Eigen::MatrixXf solution = ARAPCholesky.solve(rhs);

#pragma omp parallel for num_threads( threads )
		for (int i = 0; i < referenceVertices.size(); i++){
			int vi = variableIndex[i];
			if (vi != -1) for (int k = 0; k < 3; k++) currentVertices[i][k] = solution(vi, k);
		}
}

void ARAPModel::UpdteSoftWeight(const float & p_softScale){
	int vCount = referenceVertices.size();
	softWeights.resize(vCount, p_softScale / float(vCount));

	std::vector<Eigen::Triplet<double>> softWeighTriplets;
	softWeighTriplets.reserve(freeVarCount);
	for (int i = 0; i < vCount; i++){
		int vi = variableIndex[i];
		if (vi != -1) softWeighTriplets.push_back(Eigen::Triplet<double>(vi, vi, softWeights[i]));
	}

	Eigen::SparseMatrix<float> softWeightsMatrix;
	softWeightsMatrix.resize(freeVarCount, freeVarCount);
	softWeightsMatrix.setFromTriplets(softWeighTriplets.begin(), softWeighTriplets.end());

	ARAPCholesky.factorize(stiffnessMatrix + softWeightsMatrix);
}


void ARAPModel::Initialize(const std::vector<Eigen::Vector3f> & p_referenceVertices, const std::vector<std::vector<TriangleIndex>> & p_triangles, const std::vector<int> & p_fixedIndices, const float & p_softScale){
	fixedIndices = p_fixedIndices;
	referenceVertices = p_referenceVertices;
	triangles = p_triangles;

	TriangleCotangentWeights(referenceVertices, triangles, cotangentWeights);

	int vCount = referenceVertices.size();
	softWeights.resize(vCount, p_softScale / float(vCount));

	{//Set free variables
		std::vector<bool> isFixed(vCount, false);
		for (int i = 0; i < fixedIndices.size(); i++)isFixed[fixedIndices[i]] = true;

		variableIndex.resize(vCount, -1);
		int varCounter = 0;
		for (int i = 0; i < vCount; i++) if (!isFixed[i]){
			variableIndex[i] = varCounter;
			varCounter++;
		}

		freeVarCount = varCounter;

		if (freeVarCount != (vCount - fixedIndices.size()))printf("Variable counters unexpected! \n");

	}

	{//Set matrices
		std::vector<Eigen::Triplet<double>> stiffnessTriplets;

		for (int t = 0; t < triangles.size(); t++){
			for (int n = 0; n < triangles[t].size(); n++){
				for (int i = 0; i < 3; i++){
					for (int j = 0; j < 3; j++){
						int vP = variableIndex[triangles[t][n][j]];
						int vN = variableIndex[triangles[t][n][(j + 1) % 3]];
						if (vP != -1){
							stiffnessTriplets.push_back(Eigen::Triplet<double>(vP, vP, cotangentWeights[t][n][j]));
							if (vN != -1) stiffnessTriplets.push_back(Eigen::Triplet<double>(vP, vN, -cotangentWeights[t][n][j]));
						}
						if (vN != -1){
							stiffnessTriplets.push_back(Eigen::Triplet<double>(vN, vN, cotangentWeights[t][n][j]));
							if (vP != -1) stiffnessTriplets.push_back(Eigen::Triplet<double>(vN, vP, -cotangentWeights[t][n][j]));
						}
					}
				}
			}
		}


		stiffnessMatrix.resize(freeVarCount, freeVarCount);
		stiffnessMatrix.setFromTriplets(stiffnessTriplets.begin(), stiffnessTriplets.end());


		std::vector<Eigen::Triplet<double>> softWeighTriplets;
		softWeighTriplets.reserve(freeVarCount);
		for (int i = 0; i < vCount; i++){
			int vi = variableIndex[i];
			if (vi != -1) softWeighTriplets.push_back(Eigen::Triplet<double>(vi, vi, softWeights[i]));
		}

		Eigen::SparseMatrix<float> softWeightsMatrix;
		softWeightsMatrix.resize(freeVarCount, freeVarCount);
		softWeightsMatrix.setFromTriplets(softWeighTriplets.begin(), softWeighTriplets.end());

		ARAPCholesky.analyzePattern(stiffnessMatrix + softWeightsMatrix);
		ARAPCholesky.factorize(stiffnessMatrix + softWeightsMatrix);
	}
}

#endif //SURFACE_ARAP_INCLUDED
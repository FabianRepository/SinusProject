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

#ifndef VECTOR_IO_INCLUDED
#define VECTOR_IO_INCLUDED
#include<cstdlib>
#include<cstdio>
#include<vector>

template<typename T>
int ReadVector(std::vector<T> & vec, const char * fileName){

	vec.clear();
	FILE * file;
	file = fopen(fileName, "rb");
	if (!file) { printf("Unable to read %s \n", fileName); return 0; }
	int vecSize;
	fread(&vecSize, sizeof(int), 1, file);
	vec.resize(vecSize);
	fread(&vec[0], sizeof(T), vecSize, file);
	fclose(file);
	return 1;
}

template<typename T>
void WriteVector(const std::vector<T> & vec, const char * fileName){

	FILE * file;
	file = fopen(fileName, "wb");
	int vecSize = (int)vec.size();
	fwrite(&vecSize, sizeof(int), 1, file);
	fwrite(&vec[0], sizeof(T), vecSize, file);
	fclose(file);
}

void SaveEigenMatrixf(const Eigen::MatrixXf & matrix, const char * fileName){
	FILE * file;
	file = fopen(fileName, "wb");
	int numRows = (int)matrix.rows();
	int numCols = (int)matrix.cols();
	fwrite(&numRows, sizeof(int), 1, file);
	fwrite(&numCols, sizeof(int), 1, file);
	fwrite(&matrix(0, 0), sizeof(float), numRows*numCols, file);
	fclose(file);
}

int ReadEigenMatrixf(Eigen::MatrixXf & matrix, const char * fileName){
	FILE * file;
	file = fopen(fileName, "rb");
	if (!file) { printf("Unable to read %s \n", fileName); return 0; }
	int numRows;
	int numCols;
	fread(&numRows, sizeof(int), 1, file);
	fread(&numCols, sizeof(int), 1, file);
	matrix.resize(numRows, numCols);
	fread(&matrix(0, 0), sizeof(float), numRows*numCols, file);
	fclose(file);
	return 1;
}

void SaveEigenVectorf(const Eigen::VectorXf & vector, const char * fileName){
	FILE * file;
	file = fopen(fileName, "wb");
	int size = (int)vector.size();
	fwrite(&size, sizeof(int), 1, file);
	fwrite(&vector[0], sizeof(float), size, file);
	fclose(file);
}

int ReadEigenVectorf(Eigen::VectorXf & vector, const char * fileName){
	FILE * file;
	file = fopen(fileName, "rb");
	if (!file) { printf("Unable to read %s \n", fileName); return 0; }
	int size;
	fread(&size, sizeof(int), 1, file);
	vector.resize(size);
	fread(&vector[0], sizeof(float), size, file);
	fclose(file);
	return 1;
}

void SaveEigenVector3f(const Eigen::Vector3f & vector, const char * fileName){
	FILE * file;
	file = fopen(fileName, "wb");
	fwrite(&vector[0], sizeof(float), 3, file);
	fclose(file);
}

int ReadEigenVector3f(Eigen::Vector3f & vector, const char * fileName){
	FILE * file;
	file = fopen(fileName, "rb");
	if (!file) { printf("Unable to read %s \n", fileName); return 0; }
	fread(&vector[0], sizeof(float), 3, file);
	fclose(file);
	return 1;
}


void SaveScalarf(const float & value, const char * fileName){
	FILE * file;
	file = fopen(fileName, "wb");
	fwrite(&value, sizeof(float), 1, file);
	fclose(file);
}

int ReadScalarf(float & value, const char * fileName){
	FILE * file;
	file = fopen(fileName, "rb");
	if (!file) { printf("Unable to read %s \n", fileName); return 0; }
	fread(&value, sizeof(float), 1, file);
	fclose(file);
	return 1;
}

#endif //VECTOR_IO_INCLUDED
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

#include "ToolBuilderVisualization.inl"
#include "ParametricCurve.h"

class ToolBuilder
{
public:

	static ToolBuilderVisualization visualization;
	static int Init(const int numNodes, const int numSegmentSamples, const int numOutputSamples, const float scale);
	static void Display(void){ visualization.Display(); }
	static void MouseFunc(int button, int state, int x, int y){ visualization.mouseFunc(button, state, x, y); }
	static void MotionFunc(int x, int y){ visualization.motionFunc(x, y); }
	static void Reshape(int w, int h){ visualization.Reshape(w, h); }
	static void KeyboardFunc(unsigned char key, int x, int y){ visualization.KeyboardFunc(key, x, y); }
};

ToolBuilderVisualization						ToolBuilder::visualization;


int ToolBuilder::Init(const int numNodes, const int numSegmentSamples, const int numOutputSamples, const float scale){
	visualization.scale = scale;
	visualization.controlNodes.resize(numNodes);
	visualization.controlRadius.resize(numNodes,0.02);
	visualization.numSegmentSamples = numSegmentSamples;
	int totalCurveSamples = numNodes + (numNodes - 1) * numSegmentSamples;
	visualization.curveSamples.resize(totalCurveSamples);
	visualization.curveRadius.resize(totalCurveSamples);
	visualization.sampleWeights.resize(4 * (numSegmentSamples+2));
	for (int i = 0; i < numSegmentSamples + 2; i++){
		double t = double(i) / double(numSegmentSamples + 1);
		double weights[4];
		//CardinalCubicBspline(weights, t);
		UniformCubicBspline(weights, t);
		for (int k = 0; k < 4; k++){
			visualization.sampleWeights[4 * i + k] = weights[k];
		}
	}

	for (int i = 0; i < numNodes; i++) visualization.controlNodes[i] = Eigen::Vector3f(0, double(i) / double(numNodes - 1), 0) - Eigen::Vector3f(0, 0.5, 0);	
	visualization.UpdateCurve();
	visualization.numOutputNodes = numOutputSamples;
	return 1;
}
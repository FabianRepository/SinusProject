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

#include "PathSelection.h"
#include <stdlib.h>
#include <GL/glew.h>
#include <GL/glut.h>
#include <Util/CmdLineParser.h>
#include <Util/Algebra.h>
#include <Util/Camera.h>
#include <Util/Visualization.h>


cmdLineParameter< char * > Input("input");
cmdLineReadable* params[] =
{
	&Input,
	NULL
};

void ShowUsage(const char* ex)
{
	printf("Usage: %s [params]:\n", ex);
	printf("\t[--%s <input model>]\n", Input.name);
}

int main(int argc, char* argv[])
{

	cmdLineParse(argc - 1, argv + 1, params);
	if (!Input.set){ ShowUsage(argv[0]); return EXIT_FAILURE; }

	PathSelection::Init(Input.value);

	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
	glutInitWindowSize(600, 600);
	glutInit(&argc, argv);
	char windowName[1024];
	sprintf(windowName, "Select path and ROI");
	glutCreateWindow(windowName);

	if (glewInit() != GLEW_OK) fprintf(stderr, "[ERROR] glewInit failed\n"), exit(0);
	glutDisplayFunc(PathSelection::Display);
	glutReshapeFunc(PathSelection::Reshape);
	glutMouseFunc(PathSelection::MouseFunc);
	glutMotionFunc(PathSelection::MotionFunc);
	glutKeyboardFunc(PathSelection::KeyboardFunc);
	glutMainLoop();

	return EXIT_SUCCESS;

}
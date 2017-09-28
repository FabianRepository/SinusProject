#ifndef RASTERIZE_INCLUDED
#define RASTERIZE_INCLUDED

#include "Geometry.h"
#include "CubeGrid.h"

/////////////////Rasterizer

void Rasterize(Point3D<float> v1, CubeGrid< char >& grid, float scale = float(1.));

void Rasterize(Point3D<float> v1, Point3D<float> v2, CubeGrid< char >& grid, float scale = float(1.));

void Rasterize(Point3D<float> v1, Point3D<float> v2, Point3D<float> v3, CubeGrid< char >& grid, float scale = float(1.));


void Rasterize(const std::vector< Point3D<float> >& vertices, const std::vector< TriangleIndex >& triangles, CubeGrid< char >& grid, float scale = float(1.), int threads = 1);

void RasterizeEdges(const std::vector< Point3D<float> >& vertices, const std::vector< TriangleIndex >& triangles, CubeGrid< char >& grid, float scale = float(1.), int threads = 1);

#include "Rasterize.inl"
#endif //RASTERIZE_INCLUDED
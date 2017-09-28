/* -*- C++ -*-
Copyright (c) 2006, Michael Kazhdan and Matthew Bolitho
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
#ifndef GEOMETRY_INCLUDED
#define GEOMETRY_INCLUDED
#include <cmath>
#include <cassert>
#include <complex>
#include <vector>
#include "Algebra.h"
#include <cfloat>
template< class Real > Real Random( void );

template< class Real , int Dim >
class Point : public InnerProductSpace< Real , Point< Real , Dim > >
{
public:
	/////////////////////////////////
	// Inner product space methods //
	void Add            ( const Point& p );
	void Scale          ( Real s );
	Real InnerProduct   ( const Point< Real , Dim >& p ) const;
	/////////////////////////////////

	Real coords[Dim];
	Point ( void ) { memset( coords , 0 , sizeof(Real)*Dim ); }
	template<class Real2>
	operator Point< Real2, Dim > ( void ) const
	{
		Point< Real2, Dim > p;
		for( int d=0 ; d<Dim ; d++ ) p.coords[d] = Real2( coords[d] ); 
		return p;
	}
	Real& operator [] (int idx) { return coords[idx]; }
	const Real& operator [] (int idx) const { return coords[idx]; }

	////////////////////////////////////////////////////////////////////////////
	/*! Writes an ASCII representation of the point to an output stream.
	//  @param[in]  os   output stream
	//  @param[in]  p    the vector to be output
	//  @return     reference to the output stream for output operator chaining
	*///////////////////////////////////////////////////////////////////////////
	friend std::ostream & operator<<(std::ostream &os, const Point &p)
	{
		os << "[";
		for (int i = 0; i < Dim; ++i)   {
			if (i)
				os << ", ";
			os << p[i];
		}

		os << "]";

		return os;
	}

};

////////////////////////////////////////////////////////////////////////////////
/*! Convenience function to create Point<Real, 2>s
//  @param[in]  x   the first component
//  @param[in]  x   the second component
//  @return     Point<Real, 2> initialized with the two components
*///////////////////////////////////////////////////////////////////////////////
template<typename Real>
Point<Real, 2> MakePoint2(Real x, Real y)
{
	Point<Real, 2> result;
	result[0] = x;
	result[1] = y;
	return result;
}

////////////////////////////////////////////////////////////////////////////////
/*! Complex-real vector dot product
//  @param[in]  a   complex vector lhs
//  @param[in]  b   real vector rhs
//  @return     (complex result) a dot b
*///////////////////////////////////////////////////////////////////////////////
template <class Real, int Dim>
inline std::complex<Real> operator*(const Point<std::complex<Real>, Dim> &a
	, const Point<Real              , Dim> &b)
{
	std::complex<Real> dot = 0;
	for (int i = 0; i < Dim; ++i)   {
		dot += a[i] * b[i];
	}
	return dot;
}


template< class Real >
struct Point2D
{
	Point2D( void ){ coords[0] = coords[1] = (Real)0.; }
	Point2D( Real x , Real y ){ coords[0] = x , coords[1] = y; }
	template< class Real2 > Point2D( const Point2D< Real2 >& p ){ coords[0] = (Real)p.coords[0] , coords[1] = (Real)p.coords[1]; }
	template< class Real2 > Point2D( const Point< Real2 , 2 >& p ){ coords[0] = (Real)p.coords[0] , coords[1] = (Real)p.coords[1]; }
	template< class Real2 > Point2D& operator = ( const Point2D< Real2 >& p ){ coords[0] = (Real)p.coords[0] , coords[1] = (Real)p.coords[1] ; return *this; }
	template< class Real2 > Point2D& operator = ( const Point< Real2 , 2 >& p ){ coords[0] = (Real)p.coords[0] , coords[1] = (Real)p.coords[1] ; return *this; }
	Real& operator[] ( int i ) { return coords[i]; }
	const Real& operator[] ( int i ) const { return coords[i]; }

	Point2D& operator += ( const Point2D& p ){ coords[0] += p.coords[0] , coords[1] += p.coords[1] ; return *this; }
	Point2D& operator -= ( const Point2D& p ){ coords[0] -= p.coords[0] , coords[1] -= p.coords[1] ; return *this; }
	Point2D& operator *= ( Real s ){ coords[0] *= s , coords[1] *= s ; return *this; }
	Point2D& operator /= ( Real s ){ coords[0] /= s , coords[1] /= s ; return *this; }

	Point2D operator - ( void ) const { return Point2D( -coords[0] , -coords[1] ); }

	Point2D operator + ( const Point2D& p ) const { return Point2D( coords[0] + p.coords[0] , coords[1] + p.coords[1] ); }
	Point2D operator - ( const Point2D& p ) const { return Point2D( coords[0] - p.coords[0] , coords[1] - p.coords[1] ); }
	Point2D operator * ( Real s ) const { return Point2D( coords[0] * s , coords[1] * s ); }
	Point2D operator / ( Real s ) const { return Point2D( coords[0] / s , coords[1] / s ); }

	Real squareNorm( void ) const { return coords[0]*coords[0] + coords[1]*coords[1]; }
	static Real Dot( const Point2D& p1 , const Point2D& p2 ){ return p1.coords[0]*p2.coords[0] + p1.coords[1]*p2.coords[1]; }
	static Real SquareNorm( const Point2D& p ){ return p.coords[0]*p.coords[0] + p.coords[1]*p.coords[1]; }
	static Real Length( const Point2D& p ){ return (Real)sqrt( SquareNorm(p) ); }
	Real coords[2];
};

template< class Real >
struct Point3D
{
	Point3D( void ){ coords[0] = coords[1] = coords[2] = (Real)0.; }
	Point3D( Real x , Real y , Real z ){ coords[0] = x , coords[1] = y , coords[2] = z; }
	template< class Real2 > Point3D( const Point3D< Real2 >& p ){ coords[0] = (Real)p.coords[0] , coords[1] = (Real)p.coords[1] , coords[2] = (Real)p.coords[2]; }
	template< class Real2 > Point3D( const Point< Real2 , 3 >& p ){ coords[0] = (Real)p.coords[0] , coords[1] = (Real)p.coords[1] , coords[2] = (Real)p.coords[2]; }
	template< class Real2 > Point3D& operator = ( const Point3D< Real2 >& p ){ coords[0] = (Real)p.coords[0] , coords[1] = (Real)p.coords[1] , coords[2] = (Real)p.coords[2] ; return *this; }
	template< class Real2 > Point3D& operator = ( const Point< Real2 , 3 >& p ){ coords[0] = (Real)p.coords[0] , coords[1] = (Real)p.coords[1] , coords[2] = (Real)p.coords[2] ; return *this; }
	Real& operator[] ( int i ) { return coords[i]; }
	const Real& operator[] ( int i ) const { return coords[i]; }

	Point3D& operator += ( const Point3D& p ){ coords[0] += p.coords[0] , coords[1] += p.coords[1] , coords[2] += p.coords[2] ; return *this; }
	Point3D& operator -= ( const Point3D& p ){ coords[0] -= p.coords[0] , coords[1] -= p.coords[1] , coords[2] -= p.coords[2] ; return *this; }
	Point3D& operator *= ( Real s ){ coords[0] *= s , coords[1] *= s , coords[2] *= s ; return *this; }
	Point3D& operator /= ( Real s ){ coords[0] /= s , coords[1] /= s , coords[2] /= s ; return *this; }

	Point3D operator - ( void ) const { return Point3D( -coords[0] , -coords[1] , -coords[2] ); }

	Point3D operator + ( const Point3D& p ) const { return Point3D( coords[0] + p.coords[0] , coords[1] + p.coords[1] , coords[2] + p.coords[2] ); }
	Point3D operator - ( const Point3D& p ) const { return Point3D( coords[0] - p.coords[0] , coords[1] - p.coords[1] , coords[2] - p.coords[2] ); }
	Point3D operator * ( Real s ) const { return Point3D( coords[0] * s , coords[1] * s , coords[2] * s ); }
	Point3D operator / ( Real s ) const { return Point3D( coords[0] / s , coords[1] / s , coords[2] / s ); }

	Real squareNorm( void ) const { return coords[0]*coords[0] + coords[1]*coords[1] + coords[2]*coords[2]; }
	static Real Dot( const Point3D& p1 , const Point3D& p2 ){ return p1.coords[0]*p2.coords[0] + p1.coords[1]*p2.coords[1] + p1.coords[2]*p2.coords[2]; }
	static Real SquareNorm( const Point3D& p ){ return p.coords[0]*p.coords[0] + p.coords[1]*p.coords[1] + p.coords[2]*p.coords[2]; }
	static Real Length( const Point3D& p ){ return (Real)sqrt( SquareNorm(p) ); }
	static Point3D CrossProduct( const Point3D& p1 , const Point3D& p2 )
	{
		return Point3D( p1.coords[1]*p2.coords[2] - p1.coords[2]*p2.coords[1], - p1.coords[0]*p2.coords[2] + p1.coords[2]*p2.coords[0] , p1.coords[0]*p2.coords[1] - p1.coords[1]*p2.coords[0] );
	}
	Real coords[3];
};

class TriangleIndex
{
protected:
	unsigned int v[3];
public:
	TriangleIndex() { v[0] = v[1] = v[2] = 0; }
	TriangleIndex( unsigned int v0 , unsigned int v1 , unsigned int v2 ){ v[0] = v0; v[1] = v1; v[2] = v2; }
	unsigned int &operator[]( unsigned int idx ) { return v[idx]; }
	unsigned int  operator[]( unsigned int idx ) const { return v[idx]; }
};

//////////////////////////////
// MinimalAreaTriangulation //
//////////////////////////////
template <class Real>
class MinimalAreaTriangulation
{
	double* bestTriangulation;
	int* midPoint;
	double GetArea( const int& i , const int& j , const std::vector<Point3D<Real> >& vertices );
	void GetTriangulation( const int& i , const int& j , const std::vector<Point3D<Real> >& vertices,std::vector<TriangleIndex>& triangles , int& idx);
public:
	MinimalAreaTriangulation(void);
	~MinimalAreaTriangulation(void);
	double GetArea(const std::vector<Point3D<Real> >& vertices);
	void GetTriangulation( const std::vector<Point3D<Real> >& vertices , std::vector<TriangleIndex>& triangles );
};

template <class Real>
MinimalAreaTriangulation<Real>::MinimalAreaTriangulation(void)
{
	bestTriangulation=NULL;
	midPoint=NULL;
}
template <class Real>
MinimalAreaTriangulation<Real>::~MinimalAreaTriangulation(void)
{
	if(bestTriangulation)
		delete[] bestTriangulation;
	bestTriangulation=NULL;
	if(midPoint)
		delete[] midPoint;
	midPoint=NULL;
}
template <class Real>
void MinimalAreaTriangulation<Real>::GetTriangulation( const std::vector<Point3D<Real> >& vertices , std::vector<TriangleIndex>& triangles )
{
	triangles.resize( vertices.size() - 2 );
	if( vertices.size()==3 )
	{
		triangles[0][0]=0;
		triangles[0][1]=1;
		triangles[0][2]=2;
		return;
	}
	else if( vertices.size()==4 )
	{
		TriangleIndex tIndex[2][2];
		Real area[2];

		area[0]=area[1]=0;

		tIndex[0][0][0]=0;
		tIndex[0][0][1]=1;
		tIndex[0][0][2]=2;
		tIndex[0][1][0]=2;
		tIndex[0][1][1]=3;
		tIndex[0][1][2]=0;

		tIndex[1][0][0]=0;
		tIndex[1][0][1]=1;
		tIndex[1][0][2]=3;
		tIndex[1][1][0]=3;
		tIndex[1][1][1]=1;
		tIndex[1][1][2]=2;

		Point3D<Real> n,p1,p2;
		for(int i=0;i<2;i++)
			for(int j=0;j<2;j++)
				for(int k=0;k<3;k++)
				{
					p1.coords[k]=vertices[tIndex[i][j][1]].coords[k]-vertices[tIndex[i][j][0]].coords[k];
					p2.coords[k]=vertices[tIndex[i][j][2]].coords[k]-vertices[tIndex[i][j][0]].coords[k];
					n = Point3D<Real>::CrossProduct(p1,p2);
					area[i] += Real( Point3D<Real>::Length(n) );
				}
		if(area[0]>area[1])
		{
			triangles[0]=tIndex[1][0];
			triangles[1]=tIndex[1][1];
		}
		else
		{
			triangles[0]=tIndex[0][0];
			triangles[1]=tIndex[0][1];
		}
		return;
	}

	if(bestTriangulation) delete[] bestTriangulation;
	if(midPoint) delete[] midPoint;
	bestTriangulation=NULL;
	midPoint=NULL;
	size_t eCount=vertices.size();
	bestTriangulation=new double[eCount*eCount];
	midPoint=new int[eCount*eCount];
	for (unsigned int i = 0; i < eCount * eCount; i++)
        bestTriangulation[i] = -1;
	memset(midPoint,-1,sizeof(int)*eCount*eCount);
	GetArea(0,1,vertices);
//	triangles.clear();
	int idx = 0;
//	GetTriangulation(0,1,vertices,triangles);
	GetTriangulation( 0 , 1 , vertices , triangles , idx );
}
template <class Real>
double MinimalAreaTriangulation<Real>::GetArea(const std::vector<Point3D<Real> >& vertices)
{
	if(bestTriangulation)
		delete[] bestTriangulation;
	if(midPoint)
		delete[] midPoint;
	bestTriangulation=NULL;
	midPoint=NULL;
	size_t eCount=vertices.size();
	bestTriangulation=new double[eCount*eCount];
	midPoint=new int[eCount*eCount];
	for(int i=0;i<eCount*eCount;i++)
		bestTriangulation[i]=-1;
	memset(midPoint,-1,sizeof(int)*eCount*eCount);
	return GetArea(0,1,vertices);
}
template<class Real>
void MinimalAreaTriangulation<Real>::GetTriangulation( const int& i , const int& j , const std::vector<Point3D<Real> >& vertices , std::vector<TriangleIndex>& triangles , int& idx )
{
	TriangleIndex tIndex;
	size_t eCount=vertices.size();
	int ii=i;
	if( i<j ) ii+=(int)eCount;
	if( j+1>=ii ) return;
	ii=midPoint[i*eCount+j];
	if(ii>=0)
	{
		tIndex[0]=i;
		tIndex[1]=j;
		tIndex[2]=ii;
		triangles[idx++] = tIndex;
		GetTriangulation( i , ii , vertices , triangles , idx );
		GetTriangulation( ii , j , vertices , triangles , idx );
	}
}

template<class Real>
double MinimalAreaTriangulation<Real>::GetArea(const int& i,const int& j,const std::vector<Point3D<Real> >& vertices)
{
	double a=FLT_MAX,temp;
	size_t eCount=vertices.size();
	size_t idx=i*eCount+j;
	size_t ii=i;
	if(i<j) ii+=eCount;
	if(j+1>=(int) ii)
	{
		bestTriangulation[idx]=0;
		return 0;
	}
	int mid=-1;
	for(unsigned int r=j+1;r<ii;r++)
	{
		int rr=r%eCount;
		size_t idx1=i*eCount+rr,idx2=rr*eCount+j;
		Point3D<Real> p,p1,p2;
		for(int k=0;k<3;k++)
		{
			p1.coords[k]=vertices[i].coords[k]-vertices[rr].coords[k];
			p2.coords[k]=vertices[j].coords[k]-vertices[rr].coords[k];
		}
		//CrossProduct(p1,p2,p);
		p = Point3D<Real>::CrossProduct(p1,p2);
		temp=Point3D<Real>::Length(p);

		if(bestTriangulation[idx1]>0)
		{
			temp+=bestTriangulation[idx1];
			if(temp>a)
				continue;
			if(bestTriangulation[idx2]>0)
				temp+=bestTriangulation[idx2];
			else
				temp+=GetArea(rr,j,vertices);
		}
		else
		{
			if(bestTriangulation[idx2]>0)
				temp+=bestTriangulation[idx2];
			else
				temp+=GetArea(rr,j,vertices);
			if(temp>a)
				continue;
			temp+=GetArea(i,rr,vertices);
		}

		if(temp<a)
		{
			a=temp;
			mid=rr;
		}
	}
	bestTriangulation[idx]=a;
	midPoint[idx]=mid;

	return a;
}

#endif // GEOMETRY_INCLUDED

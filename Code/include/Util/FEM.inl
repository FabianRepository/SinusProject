#include <math.h>
#include <vector>
#include <hash_map>
#include <algorithm>
#include "lineqn.h"
#include "Hash.h"

/////////
// FEM //
/////////
template< class Real >
inline SquareMatrix< Real , 2 > FEM::MakeConformal( const SquareMatrix< Real , 2 >& tensor ){ return tensor * sqrt( tensor.determinant() / tensor.determinant() ); }
template< class Real >
inline SquareMatrix< Real , 2 > FEM::MakeAuthalic( const SquareMatrix< Real , 2 >& tensor ){ return SquareMatrix< Real , 2 >::Identity() * sqrt( tensor.determinant() / tensor.determinant() ); }
template< class Real >
inline Point2D< Real > FEM::Rotate90( const SquareMatrix< Real , 2 >& tensor , Point2D< Real > v )
{
	Point2D< Real > w = tensor.inverse() * Point2D< Real >( -v[1] , v[0] );
	Real vNorm2 = Point2D< Real >::Dot( tensor * v , v ) , wNorm2 = Point2D< Real >::Dot( tensor * w , w );
	if( wNorm2 ) return w * (Real)sqrt( vNorm2 / wNorm2 );
	else         return w;
}
template< class Real >
inline void FEM::SetMassMatrix( const SquareMatrix< Real , 2 >& tensor , SquareMatrix< Real , 3 >& massMatrix , bool lump )
{
	// F_0(x,y) = (1-x-y) ; grad( F_0 ) = ( -1 , -1 )
	// F_1(x,y) = x       ; grad( F_1 ) = (  1 ,  0 )
	// F_2(x,y) = y       ; grad( F_2 ) = (  0 ,  1 )

	// < F_0 , F_0 > = \int_0^1 \int_0^{1-y}( 1 + x^2 + y^2 - 2*x - 2*y + 2*x*y ) dx dy
	//               = \int_0^1 [ (1-y) + 1/3*(1-y)^3 + (1-y)*y^2 - (1-y)^2 - 2*(1-y)*y + (1-y)^2*y ] dy
	//               = \int_0^1 [ (1-y) + 1/3*(1-3*y+3*y^2-y^3) + (y^2-y^3) - (1-2*y+y^2) - 2*(y-y^2) + (y-2*y^2+y^3) ] dy
	//               = [ 1 - 1/2  +  1/3 - 1/2 + 1/3 - 1/12  +  1/3 - 1/4  -  1 + 1 - 1/3  -  1 + 2/3  +  1/2 - 2/3 + 1/4 ]
	//               = [ 1/3 - 1/2 + 1/3 - 1/12 ]
	//               = 1/12
	// < F_0 , F_1 > = \int_0^1 \int_0^{1-y}( x - x^2 - x*y ) dx dy
	//               = \int_0^1 [ 1/2*(1-y)^2 - 1/3*(1-y)^3 - 1/2*(1-y)^2*y ] dy
	//               = \int_0^1 [ 1/2*(1-2*y+y^2) - 1/3*(1-3*y+3*y^2-y*3) - 1/2*(y-2*y^2+y^3) ] dy
	//               = [ 1/2 - 1/2 + 1/6 - 1/3 + 1/2 - 1/3 + 1/12 - 1/4 + 1/3 - 1/8 ] dy
	//               = [ 1/6 + 1/4 - 1/3 + 1/12 - 1/8 ]
	//               = 1/24
	// < F_1 , F_1 > = \int_0^1 \int_0^{1-y}( x^2 ) dx dy
	//               = \int_0^1 [ 1/3*(1-y)^3 ] dy
	//               = \int_0^1 [ 1/3*(1-3*y+3*y^2-y^3) ] dy
	//               = [ 1/3 - 1/2 + 1/3 - 1/12 ]
	//               = 1/12
	// < F_1 , F_2 > = \int_0^1 \int_0^{1-y}( x*y ) dx dy
	//               = \int_0^1 [ 1/2*(1-y)^2*y ] dy
	//               = \int_0^1 [ 1/2*(y-2*y^2+y^3) ] dy
	//               = [ 1/4 - 1/3 + 1/8 ]
	//               = 1/24

	if( !tensor.determinant() )
	{
		fprintf( stderr , "[WARNING] Vanishing metric tensor determinant\n" );
		massMatrix *= 0;
		return;
	}
	massMatrix(0,0) = massMatrix(1,1) = massMatrix(2,2) = ( lump ? 1./6 : 1./12 );
	massMatrix(1,0) = massMatrix(0,1) = massMatrix(1,2) = massMatrix(2,1) = massMatrix(0,2) = massMatrix(2,0) = ( lump ? 0. : 1./24 );

	double sqrDet = sqrt( tensor.determinant() );
	massMatrix *= sqrDet;
}
template< class Real >
inline void FEM::SetMassAndStiffnessMatrices( const SquareMatrix< Real , 2 >& tensor , SquareMatrix< Real , 3 >& massMatrix , SquareMatrix< Real , 3 >& stiffnessMatrix , bool lump )
{
	// F_0(x,y) = (1-x-y) ; grad( F_0 ) = ( -1 , -1 )
	// F_1(x,y) = x       ; grad( F_1 ) = (  1 ,  0 )
	// F_2(x,y) = y       ; grad( F_2 ) = (  0 ,  1 )

	// < F_0 , F_0 > = \int_0^1 \int_0^{1-y}( 1 + x^2 + y^2 - 2*x - 2*y + 2*x*y ) dx dy
	//               = \int_0^1 [ (1-y) + 1/3*(1-y)^3 + (1-y)*y^2 - (1-y)^2 - 2*(1-y)*y + (1-y)^2*y ] dy
	//               = \int_0^1 [ (1-y) + 1/3*(1-3*y+3*y^2-y^3) + (y^2-y^3) - (1-2*y+y^2) - 2*(y-y^2) + (y-2*y^2+y^3) ] dy
	//               = [ 1 - 1/2  +  1/3 - 1/2 + 1/3 - 1/12  +  1/3 - 1/4  -  1 + 1 - 1/3  -  1 + 2/3  +  1/2 - 2/3 + 1/4 ]
	//               = [ 1/3 - 1/2 + 1/3 - 1/12 ]
	//               = 1/12
	// < F_0 , F_1 > = \int_0^1 \int_0^{1-y}( x - x^2 - x*y ) dx dy
	//               = \int_0^1 [ 1/2*(1-y)^2 - 1/3*(1-y)^3 - 1/2*(1-y)^2*y ] dy
	//               = \int_0^1 [ 1/2*(1-2*y+y^2) - 1/3*(1-3*y+3*y^2-y*3) - 1/2*(y-2*y^2+y^3) ] dy
	//               = [ 1/2 - 1/2 + 1/6 - 1/3 + 1/2 - 1/3 + 1/12 - 1/4 + 1/3 - 1/8 ] dy
	//               = [ 1/6 + 1/4 - 1/3 + 1/12 - 1/8 ]
	//               = 1/24
	// < F_1 , F_1 > = \int_0^1 \int_0^{1-y}( x^2 ) dx dy
	//               = \int_0^1 [ 1/3*(1-y)^3 ] dy
	//               = \int_0^1 [ 1/3*(1-3*y+3*y^2-y^3) ] dy
	//               = [ 1/3 - 1/2 + 1/3 - 1/12 ]
	//               = 1/12
	// < F_1 , F_2 > = \int_0^1 \int_0^{1-y}( x*y ) dx dy
	//               = \int_0^1 [ 1/2*(1-y)^2*y ] dy
	//               = \int_0^1 [ 1/2*(y-2*y^2+y^3) ] dy
	//               = [ 1/4 - 1/3 + 1/8 ]
	//               = 1/24

	if( !tensor.determinant() )
	{
		fprintf( stderr , "[WARNING] Vanishing metric tensor determinant\n" );
		massMatrix *= 0 , stiffnessMatrix *= 0;
		return;
	}
	massMatrix(0,0) = massMatrix(1,1) = massMatrix(2,2) = ( lump ? 1./6 : 1./12 );
	massMatrix(1,0) = massMatrix(0,1) = massMatrix(1,2) = massMatrix(2,1) = massMatrix(0,2) = massMatrix(2,0) = ( lump ? 0. : 1./24 );

	SquareMatrix< Real , 2 > iTensor = tensor.inverse();
	Point2D< Real > grad[3];
	grad[0][0] = -1. , grad[0][1] = -1.;
	grad[1][0] =  1. , grad[1][1] =  0.;
	grad[2][0] =  0. , grad[2][1] =  1.;
	for( int i=0 ; i<3 ; i++ ) for( int j=0 ; j<3 ; j++ ) stiffnessMatrix(i,j) = Point2D< Real >::Dot( grad[i] , iTensor * grad[j] ) / 2.;

	double sqrDet = sqrt( tensor.determinant() );
	massMatrix *= sqrDet;
	stiffnessMatrix *= sqrDet;
}
template< class Real >
inline Point2D< Real > FEM::Gradient( const SquareMatrix< Real , 2 >& tensor , const Real values[3] ){ return tensor.inverse() * Point2D< Real >( values[1]-values[0] , values[2]-values[0] ); }

///////////////
// FEM::Mesh //
///////////////
inline int FEM::Mesh::vCount( void ) const
{
	int count = 0;
	for( unsigned int i=0 ; i<tCount ; i++ ) for( int j=0 ; j<3 ; j++ ) count = std::max< int >( count , triangles[i][j] );
	return count+1;
}
template< class Data >
inline int FEM::Mesh::Subdivide( const FEM::Mesh& inMesh , const Data* inVertexData , FEM::Mesh& outMesh , Data** outVertexData )
{
#define EDGE_KEY( i1 , i2 ) ( (i1)>(i2) ? ( ( (long long) (i1) )<<32 ) | ( (long long) (i2) ) : ( ( (long long) (i2) )<<32 ) | ( (long long) (i1) ) )

	std::vector< Data > vertexData;
	int vCount = inMesh.vCount();
	vertexData.resize( vCount );
	for( int i=0 ; i<vCount ; i++ ) vertexData[i] = inVertexData[i];

	hash_map< long long , int > vMap;
	std::vector< Data > _vertexData = vertexData;
	std::vector< TriangleIndex > _triangles;
	for( int i=0 ; i<inMesh.tCount ; i++ )
	{
		long long keys[] = { EDGE_KEY( inMesh.triangles[i][0] , inMesh.triangles[i][1] ) , EDGE_KEY( inMesh.triangles[i][1] , inMesh.triangles[i][2] ) , EDGE_KEY( inMesh.triangles[i][2] , inMesh.triangles[i][0] ) };
		int eIndex[3];
		for( int j=0 ; j<3 ; j++ )
			if( vMap.find( keys[j] )==vMap.end() ) vMap[ keys[j] ] = eIndex[j] = (int)_vertexData.size() , _vertexData.push_back( (vertexData[inMesh.triangles[i][j]] + vertexData[ inMesh.triangles[i][(j+1)%3] ] )/2.f );
			else eIndex[j] = vMap[ keys[j] ];
			_triangles.push_back( TriangleIndex( eIndex[0] , eIndex[1] , eIndex[2] ) );
			_triangles.push_back( TriangleIndex( inMesh.triangles[i][0] , eIndex[0] , eIndex[2] ) );
			_triangles.push_back( TriangleIndex( inMesh.triangles[i][1] , eIndex[1] , eIndex[0] ) );
			_triangles.push_back( TriangleIndex( inMesh.triangles[i][2] , eIndex[2] , eIndex[1] ) );
	}

	*outVertexData = new Data[ _vertexData.size() ];
	for( int i=0 ; i<_vertexData.size() ; i++ ) (*outVertexData)[i] = _vertexData[i];
	outMesh.triangles = new TriangleIndex[ _triangles.size() ];
	for( int i=0 ; i<_triangles.size() ; i++ ) outMesh.triangles[i] = _triangles[i];
	outMesh.tCount = _triangles.size();
	return (int)_vertexData.size();

#undef EDGE_KEY
}
template< class Data >
inline int FEM::Mesh::LoopSubdivide( const FEM::Mesh& inMesh , const Data* inVertexData , FEM::Mesh& outMesh , Data** outVertexData )
{
#define EDGE_KEY( i1 , i2 ) ( (i1)>(i2) ? ( ( (long long) (i1) )<<32 ) | ( (long long) (i2) ) : ( ( (long long) (i2) )<<32 ) | ( (long long) (i1) ) )

	int tCount = (int)inMesh.tCount , vCount = (int)inMesh.vCount();
	std::vector< Data > vertexData( vCount );
	std::vector< int > valence( vCount , 0 );
	for( int i=0 ; i<vCount ; i++ ) vertexData[i] = inVertexData[i];
	for( int i=0 ; i<tCount ; i++ ) for( int j=0 ; j<3 ; j++ ) valence[ inMesh.triangles[i][j] ]++;

	hash_map< long long , int > vMap;
	std::vector< Data > _vertexData( vCount );
	std::vector< TriangleIndex > _triangles;
	for( int i=0 ; i<vCount ; i++ )
	{
		double beta = valence[i]==3 ? 3./16 : 3. / ( valence[i]*8 );
		_vertexData[i] = vertexData[i] * ( 1. - valence[i] * beta );
	}
	for( int i=0 ; i<tCount ; i++ ) for( int j=0 ; j<3 ; j++ )
	{
		double beta = valence[ inMesh.triangles[i][j] ]==3 ? 3./16 : 3. / ( valence[ inMesh.triangles[i][j] ]*8 );
		_vertexData[ inMesh.triangles[i][j] ] += vertexData[ inMesh.triangles[i][(j+1)%3] ] * beta;
	}
	for( int i=0 ; i<tCount ; i++ )
	{
		long long keys[] = { EDGE_KEY( inMesh.triangles[i][0] , inMesh.triangles[i][1] ) , EDGE_KEY( inMesh.triangles[i][1] , inMesh.triangles[i][2] ) , EDGE_KEY( inMesh.triangles[i][2] , inMesh.triangles[i][0] ) };
		int eIndex[3];
		for( int j=0 ; j<3 ; j++ )
		{
			if( vMap.find( keys[j] )==vMap.end() )  vMap[ keys[j] ] = eIndex[j] = (int)_vertexData.size() , _vertexData.push_back( Data() );
			else eIndex[j] = vMap[ keys[j] ];
			_vertexData[ eIndex[j] ] += vertexData[inMesh.triangles[i][j]] * (3./8) + vertexData[ inMesh.triangles[i][(j+2)%3] ] * (1./8);
		}
		_triangles.push_back( TriangleIndex( eIndex[0] , eIndex[1] , eIndex[2] ) );
		_triangles.push_back( TriangleIndex( inMesh.triangles[i][0] , eIndex[0] , eIndex[2] ) );
		_triangles.push_back( TriangleIndex( inMesh.triangles[i][1] , eIndex[1] , eIndex[0] ) );
		_triangles.push_back( TriangleIndex( inMesh.triangles[i][2] , eIndex[2] , eIndex[1] ) );
	}

	*outVertexData = new Data[ _vertexData.size() ];
	for( int i=0 ; i<_vertexData.size() ; i++ ) (*outVertexData)[i] = _vertexData[i];
	outMesh.triangles = new TriangleIndex[ _triangles.size() ];
	for( int i=0 ; i<_triangles.size() ; i++ ) outMesh.triangles[i] = _triangles[i];
	outMesh.tCount = _triangles.size();
	return (int)_vertexData.size();

#undef EDGE_KEY
}
template< class Vertex , class Real >
inline void FEM::Mesh::centerAndArea( const Vertex* vertices , Point3D< Real >& center , Real& area ) const
{
	center[0] = center[1] = center[2] = area = 0;
	for( int i=0 ; i<tCount ; i++ )
	{
		Point3D< Real > v[] = { Point3D< Real >( vertices[ triangles[i][0] ] ) , Point3D< Real >( vertices[ triangles[i][1] ] ) , Point3D< Real >( vertices[ triangles[i][2] ] ) };
		Point3D< Real > N = Point3D< Real >::CrossProduct( v[1]-v[0] , v[2]-v[0] );
		Real a = (Real)Length( N ) / 2.;
		Point3D< Real > c = ( v[0] + v[1] + v[2] ) / 3.;
		center += c*a;
		area += a;
	}
	center /= area;
}
template< class Vertex , class Real >
inline SquareMatrix< Real , 2 >* FEM::Mesh::getFirstFundamentalForm( const Vertex* vertices ) const
{
	SquareMatrix< Real , 2 >* I = new SquareMatrix< Real , 2 >[ tCount ];
	setFirstFundamentalForm( vertices , I );
	return I;
}
template< class Vertex , class Real >
inline void FEM::Mesh::setFirstFundamentalForm( const Vertex* vertices , SquareMatrix< Real , 2 >* I ) const
{
#pragma omp parallel for
	for( int i=0 ; i<tCount ; i++ )
	{
		Point3D< Real > e[] = { Point3D< Real >( vertices[ triangles[i][1] ] ) - Point3D< Real >( vertices[ triangles[i][0] ] ) , Point3D< Real >( vertices[ triangles[i][2] ] ) - Point3D< Real >( vertices[ triangles[i][0] ] ) };
		for( int j=0 ; j<2 ; j++ ) for( int k=0 ; k<2 ; k++ ) I[i](j,k) = Point3D< Real >::Dot( e[j] , e[k] );
	}
}

template< class Vertex , class Real >
inline SquareMatrix< Real , 2 >* FEM::Mesh::getSecondFundamentalForm( const Vertex* vertices , int nSmooth ) const
{
	SquareMatrix< Real , 2 >* II = new SquareMatrix< Real , 2 >[ tCount ];
	setSecondFundamentalForm( vertices , II , nSmooth );
	return II;
}

template< class Vertex , class Real >
inline void FEM::Mesh::setSecondFundamentalForm( const Vertex* vertices , SquareMatrix< Real , 2 >* II , int nSmooth ) const
{
	// Code borrowed from TriMesh2: http://www.cs.princeton.edu/gfx/proj/trimesh2/
	// Author: Szymon Rusinkiewicz

	// Get the vertex count
	int vCount = this->vCount();

	// Compute vertex and triangle unit normals
	std::vector< Point3D< Real > > vNormals( vCount ) , tNormals( tCount );
	for( int i=0 ; i<tCount ; i++ )
	{
		Point3D< Real > v[] = { Point3D< Real >( vertices[ triangles[i][0] ] ) , Point3D< Real >( vertices[ triangles[i][1] ] ) , Point3D< Real >( vertices[ triangles[i][2] ] ) };
		Point3D< Real > normal = Point3D< Real >::CrossProduct( v[1]-v[0] , v[2]-v[0] );
		tNormals[i] = normal;
		for( int j=0 ; j<3 ; j++ ) vNormals[ triangles[i][j] ] += normal;
	}
	// Smooth out the vertex normal field
	for( int i=0 ; i<nSmooth ; i++ )
	{
		std::vector< Point3D< Real > > _vNormals( vCount );
		for( int j=0 ; j<tCount ; j++ ) for( int k=0 ; k<3 ; k++ ) for( int l=0 ; l<3 ; l++ ) _vNormals[ triangles[j][k] ] += vNormals[ triangles[j][l] ];
		vNormals = _vNormals;
	}
#pragma omp parallel for
	for( int i=0 ; i<vNormals.size() ; i++ )
	{
		double l = Point3D< Real >::Length( vNormals[i] );
		if( l ) vNormals[i] /= l;
	}
	for( int i=0 ; i<tNormals.size() ; i++ )
	{
		double l = Point3D< Real >::Length( tNormals[i] );
		if( l ) tNormals[i] /= l;
	}

#pragma omp parallel for
	for( int i=0 ; i<tCount ; i++ )
	{
		int idx[] = { triangles[i][0] , triangles[i][1] , triangles[i][2] };
		Point3D< Real > v[] = { Point3D< Real >( vertices[ idx[0] ] ) , Point3D< Real >( vertices[ idx[1] ] ) , Point3D< Real >( vertices[ idx[2] ] ) };
		// build the n-t-b frame
		Point3D< Real > e[3] = { v[2] - v[1] , v[0] - v[2] , v[1] - v[0] };
		Point3D< Real > t = e[0];
		t /= Point3D< Real >::Length( t );
		Point3D< Real > n = tNormals[i];
		n /= Point3D< Real >::Length( n );
		Point3D< Real > b = Point3D< Real >::CrossProduct( n , t );
		b /= Point3D< Real >::Length( b );

		// build the linear system to solve
		double m[3] = { 0, 0, 0 };
		double w[3][3] = { {0,0,0}, {0,0,0}, {0,0,0} };
		for( int j=0 ; j<3 ; j++ )
		{
			double u = Point3D< Real >::Dot( e[j] , t );
			double v = Point3D< Real >::Dot( e[j] , b );
			w[0][0] += u*u;
			w[0][1] += u*v;
			w[2][2] += v*v;
			Point3D< Real > dn = vNormals[ idx[(j+2)%3] ] - vNormals[ idx[(j+1)%3] ];
			double dnu = Point3D< Real >::Dot( dn , t );
			double dnv = Point3D< Real >::Dot( dn , b );
			m[0] += dnu*u;
			m[1] += dnu*v + dnv*u;
			m[2] += dnv*v;
		}
		w[1][1] = w[0][0] + w[2][2];
		w[1][2] = w[0][1];

		// Solve by Least-Square 
		double diag[3];
		if( !ldltdc< double , 3 >( w , diag ) ) continue;
		ldltsl< double , 3 >( w , diag , m , m );

		double disc = 4*m[1]*m[1] + m[0]*m[0] - 2*m[0]*m[2] + m[2]*m[2];
		if( disc<0 ) printf( "Negative discriminant: %g\n" , disc ) , disc = 0;
		else disc = sqrt( disc );
		double kV[2];
		Point3D< Real > kD[2];
		kV[0] = ( 0.5 ) * ( m[0] + m[2] - disc );
		kV[1] = ( 0.5 ) * ( m[0] + m[2] + disc );
		
		if( kV[0]==kV[1] )  kD[0] = t , kD[1] = b;
		else
		{
			double KminusM2, x1, x2;
			KminusM2 = (kV[0] - m[2]);
			//x1 = KminusM2	/ sqrt( m[1]*m[1] + KminusM2*KminusM2 );
			//x2 = m[1]		/ sqrt( m[1]*m[1] + KminusM2*KminusM2 );
			x1 = m[1]-m[2]+kV[1];
			x2 = m[1]-m[0]+kV[1];
			kD[1]  = ( t * x1 ) + ( b * x2 );
			kD[1] /= Point3D< Real >::Length( kD[1] );

			kD[0] = Point3D< Real >::CrossProduct( n , kD[1] );
			kD[0] /= Point3D< Real >::Length( kD[0] );
		}
		{
			Point3D< Real > e[] = { v[1]-v[0] , v[2]-v[0] };
			for( unsigned int j=0 ; j<2 ; j++ ) for( unsigned int k=0 ; k<2 ; k++ )
			{
				II[i](j,k) = 0.;
				for( unsigned int l=0 ; l<2 ; l++ ) II[i](j,k) += Point3D< Real >::Dot( e[j] , kD[l] ) * Point3D< Real >::Dot( e[k] , kD[l] ) * kV[l];
			}
		}
	}
}
template< class Vertex , class Real >
inline SquareMatrix< Real , 2 >* FEM::Mesh::getSecondFundamentalForm( const Vertex* vertices , const Point3D< Real >* vNormals ) const
{
	SquareMatrix< Real , 2 >* II = new SquareMatrix< Real , 2 >[ tCount ];
	setSecondFundamentalForm( vertices , vNormals , II );
	return II;
}
template< class Vertex , class Real >
inline void FEM::Mesh::setSecondFundamentalForm( const Vertex* vertices , const Point3D< Real >* vNormals , SquareMatrix< Real , 2 >* II ) const
{
	// Code borrowed from TriMesh2: http://www.cs.princeton.edu/gfx/proj/trimesh2/
	// Author: Szymon Rusinkiewicz

#pragma omp parallel for
	for( int i=0 ; i<tCount ; i++ )
	{
		int idx[] = { triangles[i][0] , triangles[i][1] , triangles[i][2] };
		Point3D< Real > v[] = { Point3D< Real >( vertices[ idx[0] ] ) , Point3D< Real >( vertices[ idx[1] ] ) , Point3D< Real >( vertices[ idx[2] ] ) };
		Point3D< Real > e[3] = { v[2] - v[1] , v[0] - v[2] , v[1] - v[0] };
		Point3D< Real > dN[] = { vNormals[ idx[2] ] - vNormals[ idx[1] ] , vNormals[ idx[0] ] - vNormals[ idx[2] ] , vNormals[ idx[1] ] - vNormals[ idx[0] ] };

		// build the n-t-b frame
		Point3D< Real > n = Point3D< Real >::CrossProduct( e[1] , e[2] );
		Point3D< Real > t = e[0];
		Point3D< Real > b = Point3D< Real >::CrossProduct( n , t );
		t /= Point3D< Real >::Length( t ) , n /= Point3D< Real >::Length( n ) , b /= Point3D< Real >::Length( b );

		// Solve for the matrix A minimizing:
		// E(A) = \sum_{e,f \in T} [e^t A f - < dN(e) , f >]^2
		//      = \sum_{e,f \in T} [e^t A f]^2 - 2 e^t A f < dN(e) , f > + constant
		// Taking the derivatives w.r.t. A_{ij} gives:
		// 0 = \sum_{e,f \in T} 2[e^t A f] * e[i] * f[j] - 2 e[i] * f[j] < dN(e) , f >
		// <=>
		// \sum_{e,f \in T} e[i] * f[j] * < dN(e) , f > = \sum_{e,f \in T} [e^t A f] * e[i] * f[j]

		SquareMatrix< double , 2 > A;
		{
			int index[3][2] = { {0,0} , {1,1} , {0,1} };
			Point< double , 3 > _x , _b;
			SquareMatrix< double , 3 > _L;
			for( int j=0 ; j<3 ; j++ )
			{
				Point2D< double > _e( Point3D< Real >::Dot( e[j] , t ) , Point3D< Real >::Dot( e[j] , b ) );
				for( int k=0 ; k<3 ; k++ )
				{
					Point2D< double > _f( Point3D< Real >::Dot( e[k] , t ) , Point3D< Real >::Dot( e[k] , b ) );
					for( int l=0 ; l<3 ; l++ )
					{
						_b[l] += _e[ index[l][0] ] * _f[ index[l][1] ] * ( Point3D< Real >::Dot( dN[j] , e[k] ) + Point3D< Real >::Dot( dN[k] , e[j] ) ) / 2.;
						for( int m=0 ; m<3 ; m++ ) _L(l,m) += _e[ index[m][0] ] * _f[ index[m][1] ] * _e[ index[l][0] ] * _f[ index[l][1] ];
					}
				}
			}
			_x = _L.inverse() * _b;
			for( int i=0 ; i<3 ; i++ ) A( index[i][0] , index[i][1] ) = _x[i];
			A(1,0) = A(0,1);
		}

		// Characteristic polynomial: P(x) = x^2 - x * A.tr + A.det
		// x = ( - A.tr +/- [A.tr^2 - 4*A.det]^{1/2} ) / 2
		double tr = A.trace() , det = A.determinant();
		double disc = tr * tr - 4. * det;

		if( disc<0 ) printf( "Negative discriminant: %g\n" , disc ) , disc = 0;
		else disc = sqrt( disc );
		double kValues[] = { ( 0.5 ) * ( tr - disc ) , ( 0.5 ) * ( tr + disc ) };	// principal curvature values
		Point3D< Real > kDirections[2];												// principal curvature vectors
		
		if( kValues[0]==kValues[1] ) kDirections[0] = t , kDirections[1] = b;
		else
		{
			// Solve for v = (1,s=y/x)^t such that (A-k)v = 0
			// Setting B = A-k
			// => E(s) = || B(1,s)^t ||^2  = (1,s) B^t B (1,s)^t
			//         = B^2(0,0) + [ B^2(1,0) + B^2(0,1) ] * s + s^2 * B^2(1,1)
			// =>     0 = B^2(0,1) + B^2(1,0) + 2*s * B^2(1,1)
			//        s = -B^2(0,1) + B^2(1,0) / ( 2 * B^2(1,1) )
			// => (x,y) = ( -2* B^2(1,1) , B^2(0,1) + B^2(1,0) )
			// Choose the larger eigenvalue to avoid choosing a vanishing direction
			int which = fabs( kValues[1] ) > fabs( kValues[0] ) ? 1 : 0;
			{
				SquareMatrix< double , 2 > B = A;
				B(0,0) -= kValues[which] , B(1,1) -= kValues[which];
				B = SquareMatrix< double , 2 >( B.transpose() ) * B;
				Point2D< Real > v( (Real)( -2. * B(1,1) ) , (Real)( B(1,0) + B(0,1) ) );
				kDirections[ which     ] = t*v[0] +  b*v[1];
				kDirections[(which+1)%2] = Point3D< Real >::CrossProduct( n , kDirections[which] );
			}
			kDirections[0] /= Point3D< Real >::Length( kDirections[0] ) , kDirections[1] /= Point3D< Real >::Length( kDirections[1] );
		}
		Point3D< Real > basis[] = { v[1]-v[0] , v[2]-v[0] };
		for( unsigned int j=0 ; j<2 ; j++ ) for( unsigned int k=0 ; k<2 ; k++ )
		{
			II[i](j,k) = 0.;
			for( unsigned int l=0 ; l<2 ; l++ ) II[i](j,k) += Point3D< Real >::Dot( basis[j] , kDirections[l] ) * Point3D< Real >::Dot( basis[k] , kDirections[l] ) * kValues[l];
		}
	}
}
template< class CType >
inline void FEM::Mesh::setBasisToElementMatrix( SparseMatrix< CType , int >& P ) const
{
	P.resize( (int)tCount*3 );
	for( int i=0 ; i<tCount ; i++ ) for( int j=0 ; j<3 ; j++ )
	{
		int idx = i*3+j;
		P.SetRowSize( idx , 1 );
		P[idx][0] = MatrixEntry< CType , int >( triangles[i][j] , 1. );
	}
}

template< class Vertex , class Real >
FEM::Mesh::Edge< Real >* FEM::Mesh::getEdgeXForms( const Vertex* vertices ) const
{
	FEM::Mesh::Edge< Real >* edgeXForms = new FEM::Mesh::Edge< Real >[ tCount*3 ];
	setEdgeXForms( vertices , edgeXForms );
	return edgeXForms;
}
template< class Vertex , class Real >
void FEM::Mesh::setEdgeXForms( const Vertex* vertices , FEM::Mesh::Edge< Real >* edgeXForms ) const
{
	auto unfold = [] ( int idx , Point3D< Real > tangents[2][3] , FEM::Mesh::Edge< Real >* edgeXForms , const TriangleIndex* triangles , const Vertex* vertices )
	{
		if( edgeXForms[idx].opposite==-1 ) fprintf( stderr , "[ERROR] Boundary edge (TriangleMesh::unfold)\n" ) , exit( 0 );

		int tIdx[] = { idx/3 , edgeXForms[idx].opposite/3 };
		int v[] = { (int)triangles[idx/3][(idx+1)%3] , (int)triangles[idx/3][(idx+2)%3] };

		// Compute the 3D positions of the triangles
		for( int i=0 ; i<2 ; i++ ) for( int j=0 ; j<3 ; j++ ) tangents[i][j] = Point3D< Real >( vertices[ triangles[tIdx[i]][j] ] );

		int ii[] = { -1 , -1 };
		// Get the index of vertex v[0]
		{
			for( int i=0 ; i<2 ; i++ ) for( int j=0 ; j<3 ; j++ )  if( triangles[tIdx[i]][j]==v[0] ) ii[i] = j;
			if( ii[0]==-1 || ii[1]==-1 ) fprintf( stderr , "[ERROR] Failed to find edge vertex in triangle: %d %d -- %d\n" , ii[0] , ii[1] , v[0] ) , exit( 0 );
		}

		// Translate both triangles so that vertex v[0] is at the origin;
		for( int i=0 ; i<2 ; i++ ) for( int j=0 ; j<3 ; j++ ) tangents[i][j] -= Point3D< Real >( vertices[ triangles[tIdx[i]][ii[i]] ] );
		// Rotate both triangles so that they lie in the same plane
		for( int i=0 ; i<2 ; i++ )
		{
			SquareMatrix< Real , 3 > R;
			Point3D< Real > frame[3];
			frame[0] = Point3D< Real >( vertices[ v[1] ] ) - Point3D< Real >( vertices[ v[0] ] );
			frame[1] = Point3D< Real >::CrossProduct( tangents[i][1]-tangents[i][0] , tangents[i][2]-tangents[i][0] );
			frame[0] /= Real( sqrt( Point3D< Real >::SquareNorm( frame[0] ) ) );
			frame[1] /= Real( sqrt( Point3D< Real >::SquareNorm( frame[1] ) ) );
			frame[2] = Point3D< Real >::CrossProduct( frame[0] , frame[1] );
			for( int j=0 ; j<3 ; j++ ) for( int k=0 ; k<3 ; k++ ) R(j,k) = frame[k][j];
			for( int j=0 ; j<3 ; j++ ) tangents[i][j] = R * tangents[i][j];
		}
	};
	hash_map< long long , int > edgeMap;
	for( int t=0 ; t<tCount ; t++ ) for( int v=0 ; v<3 ; v++ )
	{
		int idx = t*3 + ( (v+2)%3 );
		long long key = HalfEdgeKey( triangles[t][v] , triangles[t][(v+1)%3] );
		if( edgeMap.find(key)!=edgeMap.end() ) fprintf( stderr , "[ERROR] Edge is occupied\n" ) , exit( 0 );
		edgeMap[ key ] = idx;
	}
#pragma omp parallel for
	for( int t=0 ; t<tCount ; t++ ) for( int v=0 ; v<3 ; v++ )
	{
		int idx = t*3 + ( (v+2)%3 );
		hash_map< long long , int >::const_iterator iter = edgeMap.find( HalfEdgeKey( triangles[t][(v+1)%3] , triangles[t][v] ) );
		if( iter==edgeMap.end() ) edgeXForms[idx].opposite = -1;
		else
		{
			edgeXForms[idx].opposite = iter->second;

			Point3D< Real > p[2][3];
			int tt = iter->second/3 , vv = iter->second%3;

			unfold( idx , p , edgeXForms , triangles , vertices );
			// b -> c = p[0][0] + ( p[0][1]-p[0][0] ) * b[0] + ( p[0][2]-p[0][0] ) * b[1] -> d = m^{-1}( c - p[1][0] )
			// b -> p[0][0] + M1 * b -> M2 * ( p[0][0] - p[1][0] + M1 * b ) -> M2 * ( p[0][0]-p[1][0] ) + M2 * M1 * b
			// xForm(b) = M(b) + t
			// xForm^{-1}(b) = M^{-1}(b-t)
			Point3D< Real > d[2][3];
			SquareMatrix< Real , 3 > M[2];
			for( int k=0 ; k<2 ; k++ )
			{
				d[k][0] = p[k][1]-p[k][0] , d[k][1] = p[k][2]-p[k][0] , d[k][2] = Point3D< Real >::CrossProduct( d[k][0] , d[k][1] );
				for( int i=0 ; i<3 ; i++ ) for( int j=0 ; j<3 ; j++ ) M[k](i,j) = d[k][i][j];
			}
			Point3D< Real > translate = M[1].inverse() * ( p[0][0] - p[1][0] );
			SquareMatrix< Real , 3 > xForm = M[1].inverse() * M[0];
			for( int i=0 ; i<2 ; i++ )
			{
				edgeXForms[idx].xForm.constant[i] = translate[i];
				for( int j=0 ; j<2 ; j++ ) edgeXForms[idx].xForm.linear(i,j) = xForm(i,j);
			}
		}
	}
}
template< class Real >
FEM::Mesh::EdgeCoordinateXForm< Real > FEM::Mesh::exp( const FEM::Mesh::Edge< Real >* edgeXForms , FEM::Mesh::HermiteSamplePoint< Real >& p , Real eps ) const
{
	EdgeCoordinateXForm< Real > xForm = EdgeCoordinateXForm< Real >::Identity();
	if( !Point2D< Real >::SquareNorm( p.v ) ) return xForm;
	const int MAX_ITERS = 10000;
	int count = 0;
	int inEdge = -1;
	while( count<MAX_ITERS )
	{
		// Intersect the ray p + s*v with each of the three edges
		// Bottom edge:   p[1] + s * v[1] = 0                       => s = -p[1]/v[1]
		// Left edge:     p[0] + s * v[0] = 0                       => s = -p[0]/v[0]
		// Diagonal edge: p[1] + s * v[1] = 1 - ( p[0] + s * v[0] ) => s = ( 1 - p[0]  - p[1] ) / ( v[1] + v[0] )
		Real maxS = 0;
		int idx = -1;
		{
			Real s[] = { -p.p[1] / p.v[1]  , -p.p[0] / p.v[0] , ( Real(1.) - p.p[0]  - p.p[1] ) / ( p.v[1] + p.v[0] ) };
			if( inEdge!=2 && s[0]>0 ){ Real foo = p.p[0] + p.v[0] * s[0] ; if( foo>=-eps && foo<=1+eps ) if( s[0]>maxS ) idx = 2 , maxS = s[0]; }
			if( inEdge!=1 && s[1]>0 ){ Real foo = p.p[1] + p.v[1] * s[1] ; if( foo>=-eps && foo<=1+eps ) if( s[1]>maxS ) idx = 1 , maxS = s[1]; }
			if( inEdge!=0 && s[2]>0 ){ Real foo = p.p[0] + p.v[0] * s[2] ; if( foo>=-eps && foo<=1+eps ) if( s[2]>maxS ) idx = 0 , maxS = s[2]; }
		}
		if( idx==-1 ) fprintf( stderr , "[ERROR] Ray does not intersect triangle[%d]: (%f %f) (%g %g) [%g/%g]\n" , count , p.p[0] , p.p[1] , p.v[0] , p.v[1] , Point2D< Real >::SquareNorm(p.v) , eps*eps ) , exit( 0 );
		if( maxS>1 ) // The end-point is within the triangle
		{
			p.p += p.v , p.v -= p.v;
			return xForm;
		}
		else // The end-point is outside the triangle
		{
			const Edge< Real >& edge = edgeXForms[ p.tIdx*3 + idx ];

			p.p += p.v*maxS , p.v -= p.v*maxS , p.tIdx = edge.opposite/3;
			p.p = edge.xForm( p.p ) , p.v = edge.xForm.linear * p.v;
			inEdge = edge.opposite%3;
			xForm = edge.xForm * xForm;
		}
		count++;
	}
	fprintf( stderr , "[WARNING] Failed to converge exp after %d iterations\n" , MAX_ITERS );
	return xForm;
}

template< class Real >
FEM::Mesh::EdgeCoordinateXForm< Real > FEM::Mesh::flow( const FEM::Mesh::Edge< Real >* edgeXForms , const Point2D< Real >* vf , Real delta , FEM::Mesh::SamplePoint< Real >& p , Real eps ) const
{
	EdgeCoordinateXForm< Real > xForm = EdgeCoordinateXForm< Real >::Identity();
	const int MAX_ITERS = 1000000;
	int count = 0;
	int inEdge = -1;
	Point2D< Real > v = vf[ p.tIdx ];
	Real direction = (delta<0) ? (Real)-1. : (Real)1.;
	delta *= direction;
	while( count<MAX_ITERS )
	{
		v *= direction;
		if( !Point2D< Real >::SquareNorm( v ) ) return xForm;
		// Intersect the ray p + s*v with each of the three edges
		// Bottom edge:   p[1] + s * v[1] = 0                       => s = -p[1]/v[1]
		// Left edge:     p[0] + s * v[0] = 0                       => s = -p[0]/v[0]
		// Diagonal edge: p[1] + s * v[1] = 1 - ( p[0] + s * v[0] ) => s = ( 1 - p[0]  - p[1] ) / ( v[1] + v[0] )
		Real maxD = 0;
		int idx = -1;
		{
			Real s[] = { -p.p[1] / v[1]  , -p.p[0] / v[0] , ( Real(1.) - p.p[0]  - p.p[1] ) / ( v[1] + v[0] ) };
			if( inEdge!=2 && s[0]>0 ){ Real foo = p.p[0] + v[0] * s[0] ; if( foo>=-eps && foo<=1+eps ) if( s[0]>maxD ) idx = 2 , maxD = s[0]; }
			if( inEdge!=1 && s[1]>0 ){ Real foo = p.p[1] + v[1] * s[1] ; if( foo>=-eps && foo<=1+eps ) if( s[1]>maxD ) idx = 1 , maxD = s[1]; }
			if( inEdge!=0 && s[2]>0 ){ Real foo = p.p[0] + v[0] * s[2] ; if( foo>=-eps && foo<=1+eps ) if( s[2]>maxD ) idx = 0 , maxD = s[2]; }
		}
		if( idx==-1 ) return xForm;
		if( maxD<1e-6 ) return xForm;
		if( maxD>delta ) // The end-point is within the triangle
		{
			p.p += v;
			return xForm;
		}
		else // The end-point is outside the triangle
		{
			const Edge< Real >& edge = edgeXForms[ p.tIdx*3 + idx ];
			p.p += v*maxD , p.tIdx = edge.opposite/3 , delta -= maxD;
			p.p = edge.xForm( p.p ); // v = edge.xForm.linear * v;
			v = vf[ p.tIdx ];
			inEdge = edge.opposite%3;


			xForm = edge.xForm * xForm;
		}
		count++;
	}
	fprintf( stderr , "[WARNING] Failed to converge flow after %d iterations\n" , MAX_ITERS );
	return xForm;
}
template< class Real >
FEM::Mesh::EdgeCoordinateXForm< Real > FEM::RiemannianMesh< Real >::flow_( const FEM::Mesh::Edge< Real >* edgeXForms , const Point2D< Real >* vf , Real flowTime , FEM::Mesh::SamplePoint< Real >& p , Real maxStepSize , Real eps , std::vector< SamplePoint< Real > >* path ) const
{
	EdgeCoordinateXForm< Real > xForm = EdgeCoordinateXForm< Real >::Identity();
	int MAX_ITERS = 1000000;
	int count = 0;
	int inEdge = -1;
	Real direction = (flowTime<0) ? (Real)-1. : (Real)1.;
	Real stepSizeLeft = maxStepSize;
	Point2D< Real > v = vf[ p.tIdx ] * direction;
	flowTime *= direction;
	if( path ) path->push_back( p );
	while( count<MAX_ITERS )
	{
		if( !Point2D< Real >::SquareNorm( v ) ) return xForm;
		// Intersect the ray p + s * v with each of the three edges
		// Bottom edge:   p[1] + s * v[1] = 0                       => s = -p[1]/v[1]
		// Left edge:     p[0] + s * v[0] = 0                       => s = -p[0]/v[0]
		// Diagonal edge: p[1] + s * v[1] = 1 - ( p[0] + s * v[0] ) => s = ( 1 - p[0]  - p[1] ) / ( v[1] + v[0] )
		Real s = 0;
		int idx = -1;
		{
			Real _s[] = { -p.p[1] / v[1]  , -p.p[0] / v[0] , ( Real(1.) - p.p[0]  - p.p[1] ) / ( v[1] + v[0] ) };
			if( inEdge!=2 && _s[0]>0 ){ Real foo = p.p[0] + v[0] * _s[0] ; if( foo>=-eps && foo<=1+eps ) if( _s[0]>s ) idx = 2 , s = _s[0]; }
			if( inEdge!=1 && _s[1]>0 ){ Real foo = p.p[1] + v[1] * _s[1] ; if( foo>=-eps && foo<=1+eps ) if( _s[1]>s ) idx = 1 , s = _s[1]; }
			if( inEdge!=0 && _s[2]>0 ){ Real foo = p.p[0] + v[0] * _s[2] ; if( foo>=-eps && foo<=1+eps ) if( _s[2]>s ) idx = 0 , s = _s[2]; }
		}
		if( idx==-1 ) return xForm;
		Real squareStepSize = Point2D< Real >::Dot( v , g[p.tIdx] * v ) * s * s;
		bool updateVector = false;
		if( maxStepSize>0 && squareStepSize>stepSizeLeft*stepSizeLeft )
		{
			s = stepSizeLeft / (Real)sqrt( Point2D< Real >::Dot( v , g[p.tIdx] * v ) );
			updateVector = true;
		}

		// If we can finish the flow
		if( flowTime<s )
		{
			p.p += v * flowTime;
			if( path ) path->push_back( p );
			return xForm;
		}
		// If we do not cross a boundary, change direction
		else if( updateVector )
		{
			p.p += v * s , flowTime -= s;

			// If the the vectors are oppositely oriented, terminate the flow
			if( Point2D< Real >::Dot( v , g[p.tIdx] * vf[p.tIdx] )*direction < 0 ) return xForm;

			v = vf[ p.tIdx ] * direction;
			stepSizeLeft = maxStepSize;
			if( path ) path->push_back( p );
			inEdge = -1;
		}
		// If we cross the boundary, transport the direction of the flow
		else // The end-point is outside the triangle
		{
			// Advance along the flow until you hit the edge
			p.p += v*s , flowTime -= s;

			const Edge< Real >& edge = edgeXForms[ p.tIdx*3 + idx ];
			// Switch into the next triangle
			p.tIdx = edge.opposite/3;
			p.p = edge.xForm( p.p );
			v = edge.xForm.linear * v;

			// Mark the edge we came in on
			inEdge = edge.opposite%3;

			// Accumulate the transformations
			xForm = edge.xForm * xForm;

			stepSizeLeft -= (Real)sqrt( squareStepSize );
			if( path ) path->push_back( p );
		}
		count++;
	}
	fprintf( stderr , "[WARNING] Failed to converge flow after %d iterations\n" , MAX_ITERS );
	return xForm;
}
template< class Real >
FEM::Mesh::EdgeCoordinateXForm< Real > FEM::RiemannianMesh< Real >::flow( const FEM::Mesh::Edge< Real >* edgeXForms , const Point2D< Real >* vf , Real delta , FEM::Mesh::SamplePoint< Real >& p , Real eps , std::vector< SamplePoint< Real > >* path ) const
{
	EdgeCoordinateXForm< Real > xForm = EdgeCoordinateXForm< Real >::Identity();
	int MAX_ITERS = 1000000;
	int count = 0;
	int inEdge = -1;
	Real direction = (delta<0) ? (Real)-1. : (Real)1.;
	Point2D< Real > v = vf[ p.tIdx ] * direction;
	delta *= direction;
	if( path ) path->push_back( p );
	while( count<MAX_ITERS )
	{
		if( !Point2D< Real >::SquareNorm( v ) ) return xForm;
		// Intersect the ray p + s*`v with each of the three edges
		// Bottom edge:   p[1] + s * v[1] = 0                       => s = -p[1]/v[1]
		// Left edge:     p[0] + s * v[0] = 0                       => s = -p[0]/v[0]
		// Diagonal edge: p[1] + s * v[1] = 1 - ( p[0] + s * v[0] ) => s = ( 1 - p[0]  - p[1] ) / ( v[1] + v[0] )
		Real maxD = 0;
		int idx = -1;
		{
			Real s[] = { -p.p[1] / v[1]  , -p.p[0] / v[0] , ( Real(1.) - p.p[0]  - p.p[1] ) / ( v[1] + v[0] ) };
			if( inEdge!=2 && s[0]>0 ){ Real foo = p.p[0] + v[0] * s[0] ; if( foo>=-eps && foo<=1+eps ) if( s[0]>maxD ) idx = 2 , maxD = s[0]; }
			if( inEdge!=1 && s[1]>0 ){ Real foo = p.p[1] + v[1] * s[1] ; if( foo>=-eps && foo<=1+eps ) if( s[1]>maxD ) idx = 1 , maxD = s[1]; }
			if( inEdge!=0 && s[2]>0 ){ Real foo = p.p[0] + v[0] * s[2] ; if( foo>=-eps && foo<=1+eps ) if( s[2]>maxD ) idx = 0 , maxD = s[2]; }
		}
		if( idx==-1 ) return xForm;
		else if( maxD>delta ) // The end-point is within the triangle
		{
			p.p += v*delta;
			if( path ) path->push_back( p );
			return xForm;
		}
		else // The end-point is outside the triangle
		{
			// Advance along the flow until you hit the edge
			p.p += v*maxD , delta -= maxD;
			const Edge< Real >& edge = edgeXForms[ p.tIdx*3 + idx ];

			// If the the vectors on the two sides of the edge are oppositely oriented, terminate the flow
//			if( Point2D< Real >::Dot( edge.xForm.linear * v , g[edge.opposite/3] * vf[edge.opposite/3] )*direction < 0 ) return xForm;

			// Switch into the next triangle
			p.tIdx = edge.opposite/3;
			p.p = edge.xForm( p.p );
			v = vf[ p.tIdx ] * direction;

			// Mark the edge we came in on
			inEdge = edge.opposite%3;

			// Accumulate the transformations
			xForm = edge.xForm * xForm;
			if( path ) path->push_back( p );
		}
		count++;
	}
	fprintf( stderr , "[WARNING] Failed to converge flow after %d iterations\n" , MAX_ITERS );
	return xForm;
}
template< class Real >
Real FEM::RiemannianMesh< Real >::flow( const FEM::Mesh::Edge< Real >* edgeXForms , const Point2D< Real >* vf , Real delta , FEM::Mesh::SamplePoint< Real >& p , FEM::Mesh::EdgeCoordinateXForm< Real >& xForm , Real eps ) const
{
	Real distance = (Real)0;
	xForm = EdgeCoordinateXForm< Real >::Identity();
	const int MAX_ITERS = 1000000;
	int count = 0;
	int inEdge = -1;
	Point2D< Real > v = vf[ p.tIdx ];
	Real direction = (delta<0) ? (Real)-1. : (Real)1.;
	delta *= direction;
	while( count<MAX_ITERS )
	{
		v *= direction;
		if( !Point2D< Real >::SquareNorm( v ) ) return distance;
		// Intersect the ray p + s*v with each of the three edges
		// Bottom edge:   p[1] + s * v[1] = 0                       => s = -p[1]/v[1]
		// Left edge:     p[0] + s * v[0] = 0                       => s = -p[0]/v[0]
		// Diagonal edge: p[1] + s * v[1] = 1 - ( p[0] + s * v[0] ) => s = ( 1 - p[0]  - p[1] ) / ( v[1] + v[0] )
		Real maxD = 0;
		int idx = -1;
		{
			Real s[] = { -p.p[1] / v[1]  , -p.p[0] / v[0] , ( Real(1.) - p.p[0]  - p.p[1] ) / ( v[1] + v[0] ) };
			if( inEdge!=2 && s[0]>0 ){ Real foo = p.p[0] + v[0] * s[0] ; if( foo>=-eps && foo<=1+eps ) if( s[0]>maxD ) idx = 2 , maxD = s[0]; }
			if( inEdge!=1 && s[1]>0 ){ Real foo = p.p[1] + v[1] * s[1] ; if( foo>=-eps && foo<=1+eps ) if( s[1]>maxD ) idx = 1 , maxD = s[1]; }
			if( inEdge!=0 && s[2]>0 ){ Real foo = p.p[0] + v[0] * s[2] ; if( foo>=-eps && foo<=1+eps ) if( s[2]>maxD ) idx = 0 , maxD = s[2]; }
		}
		Real vLength = (Real)sqrt( Point2D< Real >::Dot( v , g[p.tIdx] * v ) );
		if( idx==-1 ) return distance;
		if( maxD>delta ) // The end-point is within the triangle
		{
			distance += vLength * delta;
			p.p += v*delta;
			return distance;
		}
		else // The end-point is outside the triangle
		{
			const Edge< Real >& edge = edgeXForms[ p.tIdx*3 + idx ];

			// If the the vectors on the two sides of the edge are oppositely oriented, terminate the flow
			if( Point2D< Real >::Dot( edge.xForm.linear * v , g[edge.opposite/3] * vf[edge.opposite/3] )*direction < 0 ) return distance;
			distance += vLength * maxD;
			p.p += v*maxD , p.tIdx = edge.opposite/3 , delta -= maxD;
			p.p = edge.xForm( p.p );
			v = vf[ p.tIdx ];
			inEdge = edge.opposite%3;

			xForm = edge.xForm * xForm;
		}
		count++;
	}
	fprintf( stderr , "[WARNING] Failed to converge flow after %d iterations\n" , MAX_ITERS );
	return distance;
}

/////////////////////////
// FEM::RiemannianMesh //
/////////////////////////
template< class Real >
inline void FEM::RiemannianMesh< Real >::makeUnitArea( void )
{
	double scale = 0;
#pragma omp parallel for reduction( + : scale )
	for( int i=0 ; i<tCount ; i++ ) scale += sqrt( g[i].determinant() );
	scale = 2. / scale;
#pragma omp parallel for
	for( int i=0 ; i<tCount ; i++ ) g[i] *= (Real)scale;
}
template< class Real >
inline Real FEM::RiemannianMesh< Real >::area( void ) const
{
	Real area = 0;
#pragma omp parallel for reduction( + : area )
	for( int i=0 ; i<tCount ; i++ ) area += (Real)sqrt( g[i].determinant() );
	return area / (Real)2.;
}
template< class Real >
inline Real FEM::RiemannianMesh< Real >::area( int idx ) const { return (Real)sqrt( g[idx].determinant() ) / (Real)2.; }

template< class Real >
template< class Vertex >
void FEM::RiemannianMesh< Real >::setInnerProduct( const Vertex* vertices )
{
	if( g ) delete[] g;
	g = new SquareMatrix< Real , 2 >[ tCount ];
	Mesh::setFirstFundamentalForm( vertices , g );
}
template< class Real >
inline Point2D< Real >* FEM::RiemannianMesh< Real >::getGradients( const Real* vertexValues ) const
{
	Point2D< Real >* gradients = new Point2D< Real >[ tCount ];
	setGradients( vertexValues , gradients );
	return gradients;
}
template< class Real >
inline void FEM::RiemannianMesh< Real >::setGradients( const Real* vertexValues , Point2D< Real >* gradients ) const
{
#pragma omp parallel for
	for( int i=0 ; i<tCount ; i++ )
	{
		Real values[3];
		for( int j=0 ; j<3 ; j++ ) values[j] = vertexValues[ triangles[i][j] ];
		gradients[i] = FEM::Gradient( g[i] , values );
	}
}
template< class Real >
inline Real* FEM::RiemannianMesh< Real >::getDivergence( const Point2D< Real >* vectorField ) const
{
	Real* divergence = new Real[ vCount() ];
	setDivergence( vectorField , divergence );
	return divergence;
}
template< class Real >
inline void FEM::RiemannianMesh< Real >::setDivergence( const Point2D< Real >* vectorField , Real* divergence ) const
{
	int vCount = this->vCount();
	Point2D< Real > grads[3];
	grads[0][0] = -1. , grads[0][1] = -1.;
	grads[1][0] =  1. , grads[1][1] =  0.;
	grads[2][0] =  0. , grads[2][1] =  1.;
#pragma omp parallel for
	for( int i=0 ; i<vCount ; i++ ) divergence[i] = 0;
#pragma omp parallel for
	for( int i=0 ; i<tCount ; i++ )
	{
		double area = sqrt( g[i].determinant() ) / 2.;
		for( int j=0 ; j<3 ; j++ )
		{
#pragma omp atomic
			divergence[ triangles[i][j] ] += (Real)( Point2D< Real >::Dot( grads[j] , vectorField[i] ) * area );
		}
	}
}
template< class Real >
template< class Data >
inline Data* FEM::RiemannianMesh< Real >::getProlongation( const Data* faceData ) const
{
	Data* vertexData = new Data[ vCount() ];
	setProlongation( faceData , vertexData );
	return vertexData;
}
template< class Real >
template< class Data >
inline void FEM::RiemannianMesh< Real >::setProlongation( const Data* faceData , Data* vertexData ) const
{
	int vCount = this->vCount();
	double* areas = new double[ vCount ];
#pragma omp parallel for
	for( int i=0 ; i<vCount ; i++ ) areas[i] = 0 , vertexData[i] *= (Real)0;
	for( int i=0 ; i<tCount ; i++ )
	{
		double area = sqrt( g[i].determinant() ) / 2.;
		for( int j=0 ; j<3 ; j++ )
		{
			areas[ triangles[i][j] ] += area;
			vertexData[ triangles[i][j] ] += faceData[i] * (Real)area;
		}
	}
#pragma omp parallel for
	for( int i=0 ; i<vCount ; i++ ) vertexData[i] /= (Real)areas[i];
	delete[] areas;
}
template< class Real , int Dim >
SquareMatrix< Real , Dim > __SPDSquareRoot( const SquareMatrix< Real , Dim >& A )
{
	Real _A[Dim][Dim] , _d[Dim];
	for( int i=0 ; i<Dim ; i++ ) for( int j=0 ; j<Dim ; j++ ) _A[i][j] = A(i,j);
	eigdc< Real , Dim >( _A , _d );
	SquareMatrix< Real , Dim > U , D;
	for( int i=0 ; i<Dim ; i++ ) for( int j=0 ; j<Dim ; j++ ) U(i,j) = _A[i][j];
	for( int i=0 ; i<Dim ; i++ ) D(i,i) = (Real)sqrt( _d[i] );
	return SquareMatrix< Real , Dim >( U.transpose() ) * D * U;
}
template< class Real >
template< class Vertex >
inline void FEM::RiemannianMesh< Real >::convertGradients( const Vertex* vertices , const Point3D< Real >* inGradients , Point2D< Real >* outGradients ) const
{
	// Solve for A, giving the coefficients of the 2D gradient in the dual basis to {e1,e2}, such that:
	// A^t g A = D^{-1} D D^{-1}
	// => A = g^{-1/2} * D^{-1/2}
	// where the square-root is symmetric
#pragma omp parallel for
	for( int i=0 ; i<tCount ; i++ )
	{
		SquareMatrix< Real , 2 > dot , A;
		Point3D< Real > e1 = Point3D< Real >( vertices[ triangles[i][1] ] ) - Point3D< Real >( vertices[ triangles[i][0] ] );
		Point3D< Real > e2 = Point3D< Real >( vertices[ triangles[i][2] ] ) - Point3D< Real >( vertices[ triangles[i][0] ] );
		dot(0,0) = Point3D< Real >::Dot(e1,e1);
		dot(0,1) = dot(1,0) = Point3D< Real >::Dot(e1,e2);
		dot(1,1) = Point3D< Real >::Dot(e2,e2);
		A = __SPDSquareRoot( g[i].inverse() ) * __SPDSquareRoot( dot.inverse() );
		outGradients[i] = A * Point2D< Real >( Point3D< Real >::Dot( e1 , inGradients[i] ) , Point3D< Real >::Dot( e2 , inGradients[i] ) );
	}
}
template< class Real >
template< class Vertex >
inline void FEM::RiemannianMesh< Real >::convertGradients( const Vertex* vertices , const Point2D< Real >* inGradients , Point3D< Real >* outGradients ) const
{
	// Solve for A, giving the coefficients of the 3D gradient in the basis {e1,e2}, such that:
	// A^t D A = g
	// => A = D^{-1/2} * g^{1/2}
	// where the square-root is symmetric
#pragma omp parallel for
	for( int i=0 ; i<tCount ; i++ )
	{
		SquareMatrix< Real , 2 > dot , A;
		Point3D< Real > e1 = Point3D< Real >( vertices[ triangles[i][1] ] ) - Point3D< Real >( vertices[ triangles[i][0] ] );
		Point3D< Real > e2 = Point3D< Real >( vertices[ triangles[i][2] ] ) - Point3D< Real >( vertices[ triangles[i][0] ] );
		dot(0,0) = Point3D< Real >::Dot(e1,e1);
		dot(0,1) = dot(1,0) = Point3D< Real >::Dot(e1,e2);
		dot(1,1) = Point3D< Real >::Dot(e2,e2);
		A = __SPDSquareRoot( dot.inverse() ) * __SPDSquareRoot( g[i] );
		Point2D< Real > v = A * inGradients[i];
		outGradients[i] = e1*v[0] + e2*v[1];
	}
}
template< class Real >
inline Real* FEM::RiemannianMesh< Real >::getConstraints( const Real* vertexValues , const Point2D< Real >* triangleGradients , bool lump , double vScale , double gScale ) const
{
	int vCount = this->vCount();
	Real* constraints = new Real[ vCount ];
	setConstraints( vertexValues , triangleGradients , lump , constraints , vScale , gScale );
	return constraints;
}
template< class Real >
inline void FEM::RiemannianMesh< Real >::setConstraints( const Real* vertexValues , const Point2D< Real >* triangleGradients , bool lump , Real* constraints , double vScale , double gScale ) const
{
	int vCount = this->vCount();
	SquareMatrix< Real , 3 > massMatrix;
	Point2D< Real > grads[3];
	grads[0][0] = -1. , grads[0][1] = -1.;
	grads[1][0] =  1. , grads[1][1] =  0.;
	grads[2][0] =  0. , grads[2][1] =  1.;
#pragma omp parallel for
	for( int i=0 ; i<vCount ; i++ ) constraints[i] = 0;
	for( int i=0 ; i<tCount ; i++ )
	{
		if( vertexValues )
		{
			FEM::SetMassMatrix( g[i] , massMatrix , lump );
			for( int j=0 ; j<3 ; j++ ) for( int k=0 ; k<3 ; k++ ) constraints[ triangles[i][k] ] += (Real)( massMatrix(j,k) * vertexValues[ triangles[i][j] ] * vScale );
		}
		if( triangleGradients )
		{
			double area = sqrt( g[i].determinant() ) / 2.;
			for( int j=0 ; j<3 ; j++ ) constraints[ triangles[i][j] ] += (Real)( Point2D< Real >::Dot( grads[j] , triangleGradients[i] ) * area * gScale );
		}
	}
}

template< class Real >
inline void FEM::RiemannianMesh< Real >::setElementSystemMatrices( SparseMatrix< Real , int >* mass , SparseMatrix< Real , int >* stiffness , bool lump ) const
{
	if( mass ) mass->resize( tCount*3 );
	if( stiffness ) stiffness->resize( tCount*3 );
#pragma omp parallel for
	for( int i=0 ; i<tCount ; i++ )
	{
		SquareMatrix< Real , 3 > massMatrix , stiffnessMatrix;
		SetMassAndStiffnessMatrices( g[i] , massMatrix , stiffnessMatrix , lump );

		if( mass )
			for( int j=0 ; j<3 ; j++ )
			{
				int idx = i*3+j;
				mass->SetRowSize( idx , 3 );
				for( int k=0 ; k<3 ; k++ ) (*mass)[idx][k] = MatrixEntry< Real , int >( i*3+k , massMatrix(j,k) );
			}
		if( stiffness )
			for( int j=0 ; j<3 ; j++ )
			{
				int idx = i*3+j;
				stiffness->SetRowSize( idx , 3 );
				for( int k=0 ; k<3 ; k++ ) (*stiffness)[idx][k] = MatrixEntry< Real , int >( i*3+k , stiffnessMatrix(j,k) );
			}
		}
}
template< class Real >
inline void FEM::RiemannianMesh< Real >::setElementSystemMatrix( SparseMatrix< Real , int >& M , double mWeight , double sWeight , bool lump ) const
{
	M.resize( (int)tCount*3 );
	for( int i=0 ; i<tCount ; i++ ) for( int j=0 ; j<3 ; j++ ) M.SetRowSize( i*3+j , 3 );
#pragma omp parallel for
	for( int i=0 ; i<tCount ; i++ )
	{
		SquareMatrix< Real , 3 > massMatrix , stiffnessMatrix;
		SetMassAndStiffnessMatrices( g[i] , massMatrix , stiffnessMatrix , lump );
		for( int j=0 ; j<3 ; j++ )
		{
			int idx = i*3+j;
			for( int k=0 ; k<3 ; k++ ) M[idx][k] = MatrixEntry< Real , int >( i*3+k , massMatrix(j,k)*mWeight + stiffnessMatrix(j,k)*sWeight );
		}
	}
}
template< class Real >
inline void FEM::RiemannianMesh< Real >::setBasisSystemMatrices( SparseMatrix< Real , int >* mass , SparseMatrix< Real , int >* stiffness , bool lump ) const
{
	int vCount = this->vCount();
	SparseMatrix< Real , int > P , R;
	Mesh::setBasisToElementMatrix( P );
	Transpose( P , R , vCount );
	setBasisSystemMatrices( P , R , mass , stiffness , lump );
}
template< class Real >
inline void FEM::RiemannianMesh< Real >::setBasisSystemMatrix( SparseMatrix< Real , int >& M , double mWeight , double sWeight , bool lump ) const
{
	SparseMatrix< Real , int > P , R;
	Mesh::setBasisToElementMatrix( P );
	Transpose( P , R , vCount() );
	setBasisSystemMatrix( P , R , M , mWeight , sWeight , lump );
}
template< class Real >
inline void FEM::RiemannianMesh< Real >::setBasisSystemMatrices( const SparseMatrix< Real , int >& P , const SparseMatrix< Real , int >& R , SparseMatrix< Real , int >* mass , SparseMatrix< Real , int >* stiffness , bool lump ) const
{
	SparseMatrix< Real , int > elementM , elementS;
	setElementSystemMatrices( (mass!=NULL) ? &elementM : NULL , (stiffness!=NULL) ? &elementS : NULL , lump );

	if( mass )
	{
		SparseMatrix< Real , int > temp;
		Multiply( elementM , P , temp , omp_get_max_threads() ) , Multiply( R , temp , *mass , omp_get_max_threads() );
	}
	if( stiffness )
	{
		SparseMatrix< Real , int > temp;
		Multiply( elementS , P , temp , omp_get_max_threads() ) , Multiply( R , temp , *stiffness , omp_get_max_threads() );
	}
}
template< class Real >
inline void FEM::RiemannianMesh< Real >::setBasisSystemMatrix( const SparseMatrix< Real , int >& P , const SparseMatrix< Real , int >& R , SparseMatrix< Real , int >& M , double mWeight , double sWeight , bool lump ) const
{
	SparseMatrix< Real , int > temp;
	setElementSystemMatrix( M , mWeight , sWeight , lump );
	Multiply( M , P , temp , omp_get_max_threads() ) , Multiply( R , temp , M , omp_get_max_threads() );
}
template< class Real >
inline Real FEM::RiemannianMesh< Real >::getIntegral( const Real* coefficients ) const
{
	Real integral = (Real)0;
#pragma omp parallel for
	for( int i=0 ; i<tCount ; i++ )
	{
		SquareMatrix< Real , 3 > mass;
		FEM::SetMassMatrix( g[i] , mass , false );
		for( int j=0 ; j<3 ; j++ )
		{
			Real sum = (Real)0;
			for( int k=0 ; k<3 ; k++ ) sum += mass(j,k);
#pragma omp atomic
			integral += coefficients[ triangles[i][j] ] * sum;
		}
	}
	return integral;
}
template< class Real >
inline Real FEM::RiemannianMesh< Real >::getDotProduct( const Real* coefficients1 , const Real* coefficients2 , bool lump ) const
{
	Real dotProduct = (Real)0;
#pragma omp parallel for reduction( + : dotProduct )
	for( int i=0 ; i<tCount ; i++ )
	{
		SquareMatrix< Real , 3 > mass;
		FEM::SetMassMatrix( g[i] , mass , lump );
		for( int j=0 ; j<3 ; j++ ) for( int k=0 ; k<3 ; k++ ) dotProduct += mass(j,k) * coefficients1[ triangles[i][j] ] * coefficients2[ triangles[i][k] ];
	}
	return dotProduct;
}
template< class Real >
inline void FEM::RiemannianMesh< Real >::rotate90( Point2D< Real >* vf ) const
{
#pragma omp parallel for
	for( int i=0 ; i<tCount ; i++ ) vf[i] = FEM::Rotate90( g[i] , vf[i] );
}
////////////////////////
// BasisSystemBuilder //
////////////////////////
template< class Real >
inline FEM::BasisSystemBuilder< Real >::BasisSystemBuilder( const FEM::Mesh& mesh , bool lump , SparseMatrix< Real , int >& M )
{
	_lump = lump;
	_elementIndices = new _TriangleMatrixIndex[ mesh.tCount ];
	_elementMatrices = new SquareMatrix< Real , 3 >[ mesh.tCount ];

	// Generate the basis matrix with the appropriate connectivity
	{
		// First construct the element matrix
		SparseMatrix< Real , int > element;
		element.resize( mesh.tCount*3 );
		for( int i=0 ; i<mesh.tCount ; i++ ) for( int j=0 ; j<3 ; j++ )
		{
			int idx = i*3+j;
			element.SetRowSize( idx , 3 );
			for( int k=0 ; k<3 ; k++ ) element[idx][k] = MatrixEntry< Real , int >( i*3+k , 1. );
		}

		int vCount = mesh.vCount();
		SparseMatrix< Real , int > P , R , temp;
		mesh.setBasisToElementMatrix( P );
		Transpose( P , R , vCount );
		Multiply( element , P , temp , omp_get_max_threads() ) , Multiply( R , temp , M , omp_get_max_threads() );
		for( int i=0 ; i<M.rows ; i++ ) for( int j=0 ; j<M.rowSizes[i] ; j++ ) M[i][j].Value = (Real)0.;
	}

	hash_map< long long , std::pair< int , int > > edgeMap;

#define HALF_EDGE_KEY( v1 , v2 ) ( ( ( (long long) v1)<<32 ) | ( (long long) v2) )

	// Construct the map from half-edges to matrix indices
	for( int i=0 ; i<M.rows ; i++ ) for( int j=0 ; j<M.rowSizes[i] ; j++ )
		edgeMap[ HALF_EDGE_KEY( i , M[i][j].N ) ] = std::pair< int , int >( i , j );

#pragma omp parallel for
	for( int i=0 ; i<mesh.tCount ; i++ ) for( int j=0 ; j<3 ; j++ ) for( int k=0 ; k<3 ; k++ )
		_elementIndices[i].idx[j][k] = edgeMap[ HALF_EDGE_KEY( mesh.triangles[i][j] , mesh.triangles[i][k] ) ];

#undef HALF_EDGE_KEY
}

template< class Real >
inline FEM::BasisSystemBuilder< Real >::~BasisSystemBuilder( void )
{
	if( _elementIndices ) delete[] _elementIndices , _elementIndices = NULL;
	if( _elementMatrices ) delete[] _elementMatrices , _elementMatrices = NULL;
}
template< class Real >
inline void FEM::BasisSystemBuilder< Real >::update( const RiemannianMesh< Real >& mesh , double mWeight , double sWeight , SparseMatrix< Real , int >& M )
{
#pragma omp parallel for
	for( int i=0 ; i<mesh.tCount ; i++ )
	{
		SquareMatrix< Real , 3 > massMatrix , stiffnessMatrix;
		SetMassAndStiffnessMatrices( mesh.g[i] , massMatrix , stiffnessMatrix , _lump );
		_elementMatrices[i] = massMatrix*mWeight + stiffnessMatrix*sWeight;
	}
#pragma omp parallel for
	for( int i=0 ; i<M.rows ; i++ ) for( int j=0 ; j<M.rowSizes[i] ; j++ ) M[i][j].Value = 0;

	for( int i=0 ; i<mesh.tCount ; i++ )
	{
		const _TriangleMatrixIndex& idx = _elementIndices[i];
		for( int j=0 ; j<3 ; j++ ) for( int k=0 ; k<3 ; k++ ) M[ idx.idx[j][k].first ][ idx.idx[j][k].second ].Value += _elementMatrices[i](j,k);
	}
}

#ifndef FEM_INCLUDED
#define FEM_INCLUDED

#include <string.h>
#include "SparseMatrix.h"
#include "Geometry.h"

namespace FEM
{
	// Compute the 3x3 element matrix for a right triangle with the given metric tensor
	template< class Real > void SetMassMatrix( const SquareMatrix< Real , 2 >& tensor , SquareMatrix< Real , 3 >& massMatrix , bool lump );
	template< class Real > void SetMassAndStiffnessMatrices( const SquareMatrix< Real , 2 >& tensor , SquareMatrix< Real , 3 >& massMatrix , SquareMatrix< Real , 3 >& stiffnessMatrix , bool lump );
	template< class Real > SquareMatrix< Real , 2 > MakeConformal( const SquareMatrix< Real , 2 >& tensor );
	template< class Real > SquareMatrix< Real , 2 > MakeAuthalic ( const SquareMatrix< Real , 2 >& tensor );
	template< class Real > Point2D< Real > Rotate90( const SquareMatrix< Real , 2 >& tensor , Point2D< Real > v );
	template< class Real > Point2D< Real > Gradient( const SquareMatrix< Real , 2 >& tensor , const Real values[3] );

	struct Mesh
	{
		TriangleIndex* triangles;
		size_t tCount;

		Mesh( void ) { triangles = NULL , tCount = 0; }

		inline int vCount( void ) const;

		// Compute the first fundamental form using the embedding
		template< class Vertex , class Real > SquareMatrix< Real , 2 >* getFirstFundamentalForm( const Vertex* vertices                               ) const;
		template< class Vertex , class Real > void                      setFirstFundamentalForm( const Vertex* vertices , SquareMatrix< Real , 2 >* I ) const;

		// Estimate the second fundamental form using the normal variation
		template< class Vertex , class Real > SquareMatrix< Real , 2 >* getSecondFundamentalForm( const Vertex* vertices ,                                int nSmooth=0 ) const;
		template< class Vertex , class Real > void                      setSecondFundamentalForm( const Vertex* vertices , SquareMatrix< Real , 2 >* II , int nSmooth=0 ) const;
		template< class Vertex , class Real > SquareMatrix< Real , 2 >* getSecondFundamentalForm( const Vertex* vertices , const Point3D< Real > *vNormals                                ) const;
		template< class Vertex , class Real > void                      setSecondFundamentalForm( const Vertex* vertices , const Point3D< Real > *vNormals , SquareMatrix< Real , 2 >* II ) const;

		// Compute the prolongation matrix expressing basis coefficients in terms of element coefficients
		template< class CType > void setBasisToElementMatrix( SparseMatrix< CType , int >& P ) const;

		// Perform 1-to-4 subdivision
		template< class Data > static int Subdivide( const Mesh& inMesh , const Data* inVertexData , Mesh& outMesh , Data** outVertexData );
		// Perform Loop-subdivision (assumes mesh is oriented and water-tight)
		template< class Data > static int LoopSubdivide( const Mesh& inMesh , const Data* inVertexData , Mesh& outMesh , Data** outVertexData );

		// Compute the center of mass and area of the triangle mesh
		template< class Vertex , class Real > void centerAndArea( const Vertex* vertices , Point3D< Real >& center , Real& area ) const;

		// Structures for calculating geodesics
		template< class Real >
		struct SamplePoint
		{
			int tIdx;
			Point2D< Real > p;
			SamplePoint( void ){ ; }
			SamplePoint( int tIdx , const Point2D< Real > p ){ this->tIdx = tIdx , this->p = p; }
		};
		template< class Real >
		struct HermiteSamplePoint : public SamplePoint< Real >
		{
			using SamplePoint< Real >::tIdx;
			using SamplePoint< Real >::p;
			Point2D< Real > v;
			HermiteSamplePoint( void ){ ; }
			HermiteSamplePoint( int tIdx , const Point2D< Real >& p , const Point2D< Real > v=Point2D< Real >() ){ this->tIdx = tIdx , this->p = p , this->v = v; }
			HermiteSamplePoint( const typename FEM::Mesh::SamplePoint< Real >& p , const Point2D< Real > v=Point2D< Real >() ){ tIdx = p.tIdx , this->p = p.p , this->v = v; }
		};
		template< class Real >
		struct EdgeCoordinateXForm : public Group< EdgeCoordinateXForm< Real > >
		{
			SquareMatrix< Real , 2 > linear;
			Point2D< Real > constant;
			// (A,s) * p = A*p + s
			Point2D< Real > operator() ( const Point2D< Real >& p ) const { return linear*p + constant; }
			Point2D< Real > operator * ( const Point2D< Real >& p ) const { return linear*p + constant; }

			// (A,s) * (B,t) * p = (A,s) * (B*p + t) = A*B*p + (A*t + s)
			void SetIdentity( void ){ constant = Point2D< Real >() , linear = SquareMatrix< Real , 2 >::Identity(); }
			void Multiply( const EdgeCoordinateXForm& xForm ){ constant += linear * xForm.constant , linear *= xForm.linear; }
			void Invert( void ){ linear = linear.inverse() , constant = - linear * constant; }
		};
		template< class Real >
		struct Edge
		{
			int opposite;
			EdgeCoordinateXForm< Real > xForm;
			Edge( void ) { opposite = -1; }
		};

		template< class Vertex , class Real > Edge< Real >* getEdgeXForms( const Vertex* vertices                            ) const;
		template< class Vertex , class Real > void          setEdgeXForms( const Vertex* vertices , Edge< Real >* edgeXForms ) const;
		template< class Real > EdgeCoordinateXForm< Real > exp( const Edge< Real >* edgeXForms , HermiteSamplePoint< Real >& p , Real eps=(Real)0 ) const;
		template< class Real > EdgeCoordinateXForm< Real > flow( const Edge< Real >* edgeXForms , const Point2D< Real >* vf , Real delta , SamplePoint< Real >& p , Real eps=(Real)0 ) const;
	};
	// This structure represents a Riemmanian mesh, with the triangles giving the connectivity and the square (symmetric) matrices giving the metric
	template< class Real >
	struct RiemannianMesh : public Mesh
	{
		SquareMatrix< Real , 2 >* g;

		RiemannianMesh( void ){ g = NULL; }

		template< class Vertex > void setInnerProduct( const Vertex* vertices );
		void makeUnitArea( void );
		Real area( void ) const;
		Real area( int idx ) const;


		// Compute the gradients of a function
		Point2D< Real >* getGradients( const Real* vertexValues                              ) const;
		void             setGradients( const Real* vertexValues , Point2D< Real >* gradients ) const;

		// Compute the divergence of a vector field
		Real* getDivergence( const Point2D< Real >* vectorField                    ) const;
		void  setDivergence( const Point2D< Real >* vectorField , Real* divergence ) const;

		// Average per-triangle values into the vertices
		template< class Data > Data* getProlongation( const Data* faceData                    ) const;
		template< class Data > void  setProlongation( const Data* faceData , Data* vertexData ) const;

		// Transform 3D gradients (w.r.t. the standard inner product) into 2D gradients (w.r.t. to the inner product on the mesh) and vice-versa
		template< class Vertex > void convertGradients( const Vertex* vertices , const Point3D< Real >* inGradients , Point2D< Real >* outGradients ) const;
		template< class Vertex > void convertGradients( const Vertex* vertices , const Point2D< Real >* inGradients , Point3D< Real >* outGradients ) const;

		// Computes the sum of the value and gradient constraints. If either is NULL, it doesn't contribute.
		Real* getConstraints( const Real* vertexValues , const Point2D< Real >* triangleGradients , bool lump ,                     double vScale=1. , double gScale=1. ) const;
		void  setConstraints( const Real* vertexValues , const Point2D< Real >* triangleGradients , bool lump , Real* constraints , double vScale=1. , double gScale=1. ) const;

		// Compute the system matrix, expressed in terms of element coefficients
		void setElementSystemMatrices( SparseMatrix< Real , int >* mass , SparseMatrix< Real , int >* stiffness , bool lump ) const;
		void setElementSystemMatrix  ( SparseMatrix< Real , int >& M , double mWeight , double sWeight ,          bool lump ) const;

		// Compute the system matrix, expressed in terms of basis coefficients
		void setBasisSystemMatrices(                                                                             SparseMatrix< Real , int >* mass , SparseMatrix< Real , int >* stiffness , bool lump ) const;
		void setBasisSystemMatrix  (                                                                             SparseMatrix< Real , int >& M , double mWeight , double sWeight ,          bool lump ) const;
		void setBasisSystemMatrices( const SparseMatrix< Real , int >& P , const SparseMatrix< Real , int >& R , SparseMatrix< Real , int >* mass , SparseMatrix< Real , int >* stiffness , bool lump ) const;
		void setBasisSystemMatrix  ( const SparseMatrix< Real , int >& P , const SparseMatrix< Real , int >& R , SparseMatrix< Real , int >& M , double mWeight , double sWeight ,          bool lump ) const;

		// Integrate the piecewise linear function over the mesh
		Real getIntegral( const Real* coefficients ) const;
		Real getDotProduct( const Real* c1 , const Real* c2 , bool lump ) const;

		void rotate90( Point2D< Real >* vf ) const;

		EdgeCoordinateXForm< Real > flow_( const Edge< Real >* edgeXForms , const Point2D< Real >* vf , Real delta , SamplePoint< Real >& p , Real maxStepSize , Real eps=(Real)0 , std::vector< SamplePoint< Real > >* path=NULL ) const;
		EdgeCoordinateXForm< Real > flow( const Edge< Real >* edgeXForms , const Point2D< Real >* vf , Real delta , SamplePoint< Real >& p , Real eps=(Real)0 , std::vector< SamplePoint< Real > >* path=NULL ) const;
		Real flow( const Edge< Real >* edgeXForms , const Point2D< Real >* vf , Real delta , SamplePoint< Real >& p , EdgeCoordinateXForm< Real >& xForm , Real eps=(Real)0 ) const;
	};


	template< class Real >
	struct BasisSystemBuilder
	{
		BasisSystemBuilder( const Mesh& mesh , bool lump , SparseMatrix< Real , int >& M );
		~BasisSystemBuilder( void );
		void update( const RiemannianMesh< Real >& mesh , double mWeight , double sWeight , SparseMatrix< Real , int >& M );
	private:
		bool _lump;
		struct _TriangleMatrixIndex{ std::pair< int , int > idx[3][3]; };
		_TriangleMatrixIndex* _elementIndices;
		SquareMatrix< Real , 3 >* _elementMatrices;
	};
}
#include "FEM.inl"

#endif // FEM_INCLUDED
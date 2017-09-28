#ifndef SPARSE_MATRIX_INTERFACE_INCLUDED
#define SPARSE_MATRIX_INTERFACE_INCLUDED

#define USE_SIZE_T_INDEX 0
#define FORCE_TWO_BYTE_ALIGNMENT 1
#include "Array.h"
#include <vector>


#if FORCE_TWO_BYTE_ALIGNMENT
#pragma pack(push)
#pragma pack(2)
#endif // FORCE_TWO_BYTE_ALIGNMENT
template< class T , class IndexType >
struct MatrixEntry
{
	MatrixEntry( void )       { N =-1 , Value = 0; }
	MatrixEntry( int i )      { N = i , Value = 0; }
	MatrixEntry( int n , T v ){ N = n , Value = v; }
#if USE_SIZE_T_INDEX
	size_t N;
#else // !USE_SIZE_T_INDEX
	IndexType N;
#endif // USE_SIZE_T_INDEX
	T Value;
};
#if FORCE_TWO_BYTE_ALIGNMENT
#pragma pack(pop)
#endif // FORCE_TWO_BYTE_ALIGNMENT


#define MULTIPLY_ADD 1
#define MULTIPLY_NEGATE 2
template< class T , class const_iterator > class SparseMatrixInterface
{
public:
	virtual const_iterator begin( int row ) const = 0;
	virtual const_iterator end  ( int row ) const = 0;
	virtual size_t Rows   ( void )          const = 0;
	virtual size_t RowSize( size_t idx )    const = 0;

	size_t Entries( void ) const;

	double SquareNorm( void ) const;
	double SquareASymmetricNorm( void ) const;
	double SquareASymmetricNorm( int& idx1 , int& idx2 ) const;

	template< class T2 > void Multiply           (  ConstPointer( T2 ) In , Pointer( T2 ) Out ) const;
	template< class T2 > void Multiply           (       Pointer( T2 ) In , Pointer( T2 ) Out ) const { return Multiply( ( ConstPointer( T2 ) ) In , Out ); }
	template< class T2 > void MultiplyParallel   ( ConstPointer( T2 )  In , Pointer( T2 ) Out , int threads , int multiplyFlag=MULTIPLY_ADD ) const;
	template< class T2 > void MultiplyParallel   (      Pointer( T2 )  In , Pointer( T2 ) Out , int threads , int multiplyFlag=MULTIPLY_ADD ) const { MultiplyParallel( ( ConstPointer(T2) )( In ) , Out , threads , multiplyFlag ); }
	template< class T2 > void MultiplyScaledParallel   ( T scale , ConstPointer( T2 )  In , Pointer( T2 ) Out , int threads , int multiplyFlag=MULTIPLY_ADD ) const;
	template< class T2 > void MultiplyScaledParallel   ( T scale ,      Pointer( T2 )  In , Pointer( T2 ) Out , int threads , int multiplyFlag=MULTIPLY_ADD ) const { MultiplyScaledParallel( scale , ( ConstPointer(T2) )( In ) , Out , threads , multiplyFlag ); }
	template< class T2 > void MultiplyTranspose  ( ConstPointer( T2 )  In , Pointer( T2 ) Out , int outDim , int multiplyFlag , T (*TransposeFunction)( const T& )=NULL ) const;
	template< class T2 > void MultiplyTranspose  (      Pointer( T2 )  In , Pointer( T2 ) Out , int outDim , int multiplyFlag , T (*TransposeFunction)( const T& )=NULL ) const { MultiplyTranspose( ( ConstPointer( T2 ) ) In , Out , outDim , multiplyFlag , TransposeFunction ); }
	template< class T2 > void BMinusMX           ( ConstPointer( T2 ) B , ConstPointer( T2 ) X , Pointer( T2 ) out ) const;
	template< class T2 > void BMinusMX           (      Pointer( T2 ) B , ConstPointer( T2 ) X , Pointer( T2 ) out ) const { return BMinusMX( ( ConstPointer( T2 ) ) B ,                        X , out ); }
	template< class T2 > void BMinusMX           ( ConstPointer( T2 ) B ,      Pointer( T2 ) X , Pointer( T2 ) out ) const { return BMinusMX(                        B , ( ConstPointer( T2 ) ) X , out ); }
	template< class T2 > void BMinusMX           (      Pointer( T2 ) B ,      Pointer( T2 ) X , Pointer( T2 ) out ) const { return BMinusMX( ( ConstPointer( T2 ) ) B , ( ConstPointer( T2 ) ) X , out ); }

	template< class T2 , class T3 , bool unitDiagonal , bool StrippedDiagonal , bool UseSOR > void SolveGaussSeidel           ( ConstPointer( T3 ) Diagonal , ConstPointer( T2 ) b , int iters , Pointer( T2 ) Solution , T sor ) const;
	template< class T2 , class T3 , bool unitDiagonal , bool StrippedDiagonal , bool UseSOR > void SolveGaussSeidelAndResidual( ConstPointer( T3 ) Diagonal , ConstPointer( T2 ) b , int iters , Pointer( T2 ) Solution , Pointer( T2 ) Residual , T sor ) const;

	template< class T2 , class T3 , bool UseSOR > void SolveGaussSeidel( ConstPointer( T3 ) Diagonal , ConstPointer( T2 ) b , int iters , Pointer( T2 ) Solution , T sor , bool unitDiagonal , bool strippedDiagonal ) const
	{
		if( strippedDiagonal )
			if( unitDiagonal ) SolveGaussSeidel< T2 , T3 , true  , true  , UseSOR >( Diagonal , b , iters , Solution , sor );
			else               SolveGaussSeidel< T2 , T3 , false , true  , UseSOR >( Diagonal , b , iters , Solution , sor );
		else
			if( unitDiagonal ) SolveGaussSeidel< T2 , T3 , true  , false , UseSOR >( Diagonal , b , iters , Solution , sor );
			else               SolveGaussSeidel< T2 , T3 , false , false , UseSOR >( Diagonal , b , iters , Solution , sor );
	}
	template< class T2 , class T3 , bool UseSOR > void SolveGaussSeidelAndResidual( ConstPointer( T3 ) Diagonal , ConstPointer( T2 ) b , int iters , Pointer( T2 ) Solution , Pointer( T2 ) Residual , T sor , bool unitDiagonal , bool strippedDiagonal ) const
	{
		if( strippedDiagonal )
			if( unitDiagonal ) SolveGaussSeidelAndResidual< T2 , T3 , true  , true  , UseSOR >( Diagonal , b , iters , Solution , Residual , sor );
			else               SolveGaussSeidelAndResidual< T2 , T3 , false , true  , UseSOR >( Diagonal , b , iters , Solution , Residual , sor );
		else
			if( unitDiagonal ) SolveGaussSeidelAndResidual< T2 , T3 , true  , false , UseSOR >( Diagonal , b , iters , Solution , Residual , sor );
			else               SolveGaussSeidelAndResidual< T2 , T3 , false , false , UseSOR >( Diagonal , b , iters , Solution , Residual , sor );
	}
	template< class T2 > void MultiplyTransposeParallel( ConstPointer( T2 ) in ,                       Pointer( T2 ) out , const std::vector< int >& entriesPerLine , int linesPerSlice , const std::vector< int >& slicesPerThread , int sliceDependence ) const;
	template< class T2 > void MultiplyParallel         ( ConstPointer( T2 ) in ,                       Pointer( T2 ) out , const std::vector< int >& entriesPerLine , int linesPerSlice , const std::vector< int >& slicesPerThread ) const;
	template< class T2 > void BMinusMXParallel         ( ConstPointer( T2 ) B , ConstPointer( T2 ) X , Pointer( T2 ) out , const std::vector< int >& entriesPerLine , int linesPerSlice , const std::vector< int >& slicesPerThread ) const;

	template< class T2 > void MultiplyTransposeParallel(      Pointer( T2 ) in ,                       Pointer( T2 ) out , const std::vector< int >& entriesPerLine , int linesPerSlice , const std::vector< int >& slicesPerThread , int sliceDependence ) const { MultiplyTransposeParallel( ( ConstPointer( T2 ) ) in ,                              out , entriesPerLine , linesPerSlice , slicesPerThread , sliceDependence ); }
	template< class T2 > void MultiplyParallel         (      Pointer( T2 ) in ,                       Pointer( T2 ) out , const std::vector< int >& entriesPerLine , int linesPerSlice , const std::vector< int >& slicesPerThread )                       const { MultiplyParallel         ( ( ConstPointer( T2 ) ) in ,                              out , entriesPerLine , linesPerSlice , slicesPerThread ); }
	template< class T2 > void BMinusMXParallel         (      Pointer( T2 ) B , ConstPointer( T2 ) X , Pointer( T2 ) out , const std::vector< int >& entriesPerLine , int linesPerSlice , const std::vector< int >& slicesPerThread )                       const { BMinusMXParallel         ( ( ConstPointer( T2 ) ) B ,                        X ,    out , entriesPerLine , linesPerSlice , slicesPerThread ); }
	template< class T2 > void BMinusMXParallel         ( ConstPointer( T2 ) B ,      Pointer( T2 ) X , Pointer( T2 ) out , const std::vector< int >& entriesPerLine , int linesPerSlice , const std::vector< int >& slicesPerThread )                       const { BMinusMXParallel         (                        B , ( ConstPointer( T2 ) ) X ,    out , entriesPerLine , linesPerSlice , slicesPerThread ); }
	template< class T2 > void BMinusMXParallel         (      Pointer( T2 ) B ,      Pointer( T2 ) X , Pointer( T2 ) out , const std::vector< int >& entriesPerLine , int linesPerSlice , const std::vector< int >& slicesPerThread )                       const { BMinusMXParallel         ( ( ConstPointer( T2 ) ) B , ( ConstPointer( T2 ) ) X ,    out , entriesPerLine , linesPerSlice , slicesPerThread ); }

	template< class T2 > void MultiplyTransposeParallel( ConstPointer( T2 ) in ,                       Pointer( T2 ) out , const std::vector< int >& entriesPerSlice , const std::vector< int >& slicesPerThread , int sliceDependence ) const;
	template< class T2 > void MultiplyParallel         ( ConstPointer( T2 ) in ,                       Pointer( T2 ) out , const std::vector< int >& entriesPerSlice , const std::vector< int >& slicesPerThread ) const;
	template< class T2 > void BMinusMXParallel         ( ConstPointer( T2 ) B , ConstPointer( T2 ) X , Pointer( T2 ) out , const std::vector< int >& entriesPerSlice , const std::vector< int >& slicesPerThread ) const;

	template< class T2 > void MultiplyTransposeParallel(      Pointer( T2 ) in ,                       Pointer( T2 ) out , const std::vector< int >& entriesPerSlice , const std::vector< int >& slicesPerThread , int sliceDependence ) const { MultiplyTransposeParallel( ( ConstPointer( T2 ) ) in , out , entriesPerSlice , slicesPerThread , sliceDependence ); }
	template< class T2 > void MultiplyParallel         (      Pointer( T2 ) in ,                       Pointer( T2 ) out , const std::vector< int >& entriesPerSlice , const std::vector< int >& slicesPerThread )                       const { MultiplyParallel         ( ( ConstPointer( T2 ) ) in , out , entriesPerSlice , slicesPerThread ); }
	template< class T2 > void BMinusMXParallel         (      Pointer( T2 ) B , ConstPointer( T2 ) X , Pointer( T2 ) out , const std::vector< int >& entriesPerSlice , const std::vector< int >& slicesPerThread )                       const { BMinusMXParallel         ( ( ConstPointer( T2 ) ) B ,                        X ,    out , entriesPerSlice , slicesPerThread ); }
	template< class T2 > void BMinusMXParallel         ( ConstPointer( T2 ) B ,      Pointer( T2 ) X , Pointer( T2 ) out , const std::vector< int >& entriesPerSlice , const std::vector< int >& slicesPerThread )                       const { BMinusMXParallel         (                        B , ( ConstPointer( T2 ) ) X ,    out , entriesPerSlice , slicesPerThread ); }
	template< class T2 > void BMinusMXParallel         (      Pointer( T2 ) B ,      Pointer( T2 ) X , Pointer( T2 ) out , const std::vector< int >& entriesPerSlice , const std::vector< int >& slicesPerThread )                       const { BMinusMXParallel         ( ( ConstPointer( T2 ) ) B , ( ConstPointer( T2 ) ) X ,    out , entriesPerSlice , slicesPerThread ); }

	template< class T2 , class T3 , bool UnitDiagonal , bool StrippedDiagonal , bool UseSOR >
	void SolveGaussSeidelParallel
		(
		ConstPointer( T3 ) Diagonal , ConstPointer( T2 ) b , int iters , Pointer( T2 ) Solution , T sor ,
		const std::vector< int >& entriesPerSlice , const std::vector< int >& slicesPerThread , int sliceDependence
		) const;
	template< class T2 , class T3 , bool UnitDiagonal , bool StrippedDiagonal , bool UseSOR >
	void SolveGaussSeidelAndResidualParallel
		(
		ConstPointer( T3 ) Diagonal , ConstPointer( T2 ) b , int iters , Pointer( T2 ) Solution , Pointer( T2 ) Residual , T sor ,
		const std::vector< int >& entriesPerSlice , const std::vector< int >& slicesPerThread , int sliceDependence
		) const;

	template< class T2 , class T3 , bool UseSOR >
	void SolveGaussSeidelParallel
		(
		ConstPointer( T3 ) Diagonal , ConstPointer( T2 ) b , int iters , Pointer( T2 ) Solution , T sor ,
		const std::vector< int >& entriesPerSlice , const std::vector< int >& slicesPerThread , int sliceDependence , bool unitDiagonal , bool strippedDiagonal
		) const
	{
		if( strippedDiagonal )
			if( unitDiagonal ) SolveGaussSeidelParallel< T2 , T3 , true  , true  , UseSOR >( Diagonal , b , iters , Solution , sor , entriesPerSlice , slicesPerThread , sliceDependence );
			else               SolveGaussSeidelParallel< T2 , T3 , false , true  , UseSOR >( Diagonal , b , iters , Solution , sor , entriesPerSlice , slicesPerThread , sliceDependence );
		else
			if( unitDiagonal ) SolveGaussSeidelParallel< T2 , T3 , true  , false , UseSOR >( Diagonal , b , iters , Solution , sor , entriesPerSlice , slicesPerThread , sliceDependence );
			else               SolveGaussSeidelParallel< T2 , T3 , false , false , UseSOR >( Diagonal , b , iters , Solution , sor , entriesPerSlice , slicesPerThread , sliceDependence );
	}
	template< class T2 , class T3 , bool UseSOR >
	void SolveGaussSeidelAndResidualParallel
		(
		ConstPointer( T3 ) Diagonal , ConstPointer( T2 ) b , int iters , Pointer( T2 ) Solution , Pointer( T2 ) Residual ,  T sor ,
		const std::vector< int >& entriesPerSlice , const std::vector< int >& slicesPerThread , int sliceDependence , bool unitDiagonal , bool strippedDiagonal
		) const
	{
		if( strippedDiagonal )
			if( unitDiagonal ) SolveGaussSeidelAndResidualParallel< T2 , T3 , true  , true  , UseSOR >( Diagonal , b , iters , Solution , Residual , sor , entriesPerSlice , slicesPerThread , sliceDependence );
			else               SolveGaussSeidelAndResidualParallel< T2 , T3 , false , true  , UseSOR >( Diagonal , b , iters , Solution , Residual , sor , entriesPerSlice , slicesPerThread , sliceDependence );
		else
			if( unitDiagonal ) SolveGaussSeidelAndResidualParallel< T2 , T3 , true  , false , UseSOR >( Diagonal , b , iters , Solution , Residual , sor , entriesPerSlice , slicesPerThread , sliceDependence );
			else               SolveGaussSeidelAndResidualParallel< T2 , T3 , false , false , UseSOR >( Diagonal , b , iters , Solution , Residual , sor , entriesPerSlice , slicesPerThread , sliceDependence );
	}
};

#include "SparseMatrixInterface.inl"
#endif // SPARSE_MATRIX_INTERFACE_INCLUDED

/*
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

#ifndef __SPARSEMATRIX_HPP
#define __SPARSEMATRIX_HPP

#include "SparseMatrixInterface.h"
#include "Array.h"

template <class T>
struct MatrixEntry2
{
	MatrixEntry2( void )	{ inN = outN = -1; value = 0; }
	int inN , outN;
	T value;
};
template< class T , class IndexType > class SparseMatrix : public SparseMatrixInterface< T , ConstPointer( MatrixEntry< T ,IndexType > ) >
{
	template< class T2 , class IndexType2 > friend class SparseMatrix;
	bool _contiguousMemory;
	int _maxColumnsPerRow;
public:
	typedef SparseMatrixInterface< T , ConstPointer( MatrixEntry< T , IndexType > ) > Interface;

	int rows;
	Pointer( int ) rowSizes;
	Pointer( Pointer( MatrixEntry< T , IndexType > ) ) m_ppElements;

	SparseMatrix( );
	SparseMatrix( const SparseMatrix& M );
	template< class T2 , class IndexType2 >
	SparseMatrix( const SparseMatrix< T2 , IndexType2 >& M );
	~SparseMatrix();
	SparseMatrix< T , IndexType >& operator = ( const SparseMatrix< T , IndexType >& M );
	template< class T2 , class IndexType2 >
	SparseMatrix< T , IndexType >& operator = ( const SparseMatrix< T2 , IndexType2 >& M );

	template< class T2 > void operator()( const T2* in , T2* out ) const;

	template< class T2 , class IndexType2 >
	SparseMatrix< T , IndexType >& copy( const SparseMatrix< T2 , IndexType2 >& M );

	inline ConstPointer( MatrixEntry< T , IndexType > ) begin( int row )    const { return m_ppElements[row]; }
	inline ConstPointer( MatrixEntry< T , IndexType > ) end  ( int row )    const { return m_ppElements[row]+rowSizes[row]; }
	inline size_t Rows                              ( void )       const { return rows; }
	inline size_t RowSize                           ( size_t idx ) const { return rowSizes[idx]; }

	SparseMatrix( int rows );
	void resize	( int rows );
	SparseMatrix( int rows , int maxColumnsPerRow );
	void resize	( int rows , int maxColumnsPerRow );
	void SetRowSize( int row , int count );
	inline      Pointer( MatrixEntry< T , IndexType > ) operator[] ( int idx )       { return m_ppElements[idx]; }
	inline ConstPointer( MatrixEntry< T , IndexType > ) operator[] ( int idx ) const { return m_ppElements[idx]; }

	bool isContiguous( void ) const;
	void MakeContiguous( void );

	double SquareNorm(void) const;
	double ASymmetricSquareNorm( void ) const;
	double AHermitianSquareNorm( void ) const;

	/** Sets the column index of all allocated entries to -1 so that they are
	*  treated as non-existent. This is needed because SetRowSize() uses
	*  malloc instead of new and MatrixEntry's constructor is never called. */
	void invalidateEntries();

	/** Adds a scalar value to an element in the Matrix, using a new element if
	*  necessary. If no pre-allocated space for a new element exists, false is
	*  returned.
	*  WARNING: no check is done to remove entries that become zero after
	*  addition */
	bool addScalarToEntry( T s , IndexType i , IndexType j );
};
template< class T1 , class T2 , class T3 , class IndexType >
bool TransposeMultiply( const SparseMatrix< T1 , IndexType >& At , const SparseMatrix< T2 , IndexType >& B , SparseMatrix< T3 , IndexType >& out , T1 (*TransposeFunction)( const T1& )=NULL );
template< class T1 , class T2 , class T3 , class IndexType >
bool TransposeMultiply( const SparseMatrix< T1 , IndexType >& At , const SparseMatrix< T2 , IndexType >& B , SparseMatrix< T3 , IndexType >& out , int outRows , T1 (*TransposeFunction)( const T1& )=NULL );
template< class A_T , class A_const_iterator , class B_T , class B_const_iterator , class Out_T , class Out_IndexType >
bool Multiply( const SparseMatrixInterface< A_T , A_const_iterator >& A , const SparseMatrixInterface< B_T , B_const_iterator >& B , SparseMatrix< Out_T , Out_IndexType >& out , int threads = 1 );
template< class T , class In_const_iterator , class Out_IndexType >
bool Transpose( const SparseMatrixInterface< T , In_const_iterator >& At , SparseMatrix< T , Out_IndexType >& A ,               T (*TransposeFunction)( const T& )=NULL );
template< class T , class In_const_iterator , class Out_IndexType >
bool Transpose( const SparseMatrixInterface< T , In_const_iterator >& At , SparseMatrix< T , Out_IndexType >& A , int outRows , T (*TransposeFunction)( const T& )=NULL );

template< class T , unsigned int Radius > class BandedMatrix
{
	template< class T2 , unsigned int Radius2 > friend class BandedMatrix;
	size_t _rows;
	Pointer( T ) _entries;
public:
	BandedMatrix( void );
	BandedMatrix( size_t rows );
	BandedMatrix( size_t rows , const T& clearValue );
	template< class T2 > BandedMatrix( const BandedMatrix< T2 , Radius >& M );
	template< class T2 > BandedMatrix& operator = ( const BandedMatrix< T2 , Radius >& M );
	~BandedMatrix();

	template< class T2 > void multiply ( ConstPointer( T2 ) in , Pointer( T2 ) out , int threads=1 ) const;
	template< class T2 > void multiply2( ConstPointer( T2 ) in , Pointer( T2 ) out , int threads=1 ) const;

	inline size_t rows( void ) const { return _rows; }
	inline size_t entries( void ) const { return _rows * ( 2 * Radius + 1 ); }

	void resize( size_t rows );
	void resize( size_t rows , const T& clearValue );
	inline      Pointer( T ) operator[] ( unsigned int idx )       { return _entries + idx * ( 2 * Radius + 1 ); }
	inline ConstPointer( T ) operator[] ( unsigned int idx ) const { return _entries + idx * ( 2 * Radius + 1 ); }
	double squareNorm( void ) const;
};
#include "SparseMatrix.inl"
#endif /* __SPARSEMATRIX_HPP */

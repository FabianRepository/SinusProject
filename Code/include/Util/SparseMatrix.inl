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

#include <float.h>
#include <complex>
#include "Hash.h"

///////////////////
//  SparseMatrix //
///////////////////
///////////////////////////////////////
// SparseMatrix Methods and Memebers //
///////////////////////////////////////

template< class T , class IndexType >
SparseMatrix< T , IndexType >::SparseMatrix( void )
{
	_maxColumnsPerRow = 0;
	_contiguousMemory = false;
	rowSizes = NullPointer< int >( );
	rows = 0;
	m_ppElements= NullPointer< Pointer( MatrixEntry< T , IndexType > ) >( );
}
////////////////////////////////////////////////////////////////////////////////
/*! Sets the column index of all allocated entries to -1 so that they are
//  treated as non-existent. This is needed because SetRowSize() uses malloc
//  instead of new and MatrixEntry's constructor is never called.
*///////////////////////////////////////////////////////////////////////////////
template <class T, class IndexType>
void SparseMatrix<T, IndexType>::invalidateEntries()
{
    IndexType numRows = Rows();
    for (IndexType r = 0; r < numRows; ++r) {
        IndexType numColumns = RowSize(r);
        MatrixEntry<T, IndexType> *row = m_ppElements[r];
        for (IndexType c = 0; c < numColumns; ++c) {
            row[c].N = -1;
        }
    }
}

////////////////////////////////////////////////////////////////////////////////
/*! Adds a scalar value to an element in the Matrix, using a new element if
//  necesary. If no pre-allocated space for a new element exists, false is
//  returned.
//  WARNING: no check is done to remove entries that become zero after addition
//  @param[in]  s   The scalar to add to the destination element in the matrix
//  @param[in]  i   The destination element's row
//  @param[in]  j   The destination element's column
//  @return     true if successful, false if no pre-allocated space exists for
//              the creation of a new nonzero element
*///////////////////////////////////////////////////////////////////////////////
template <class T, class IndexType>
bool SparseMatrix<T, IndexType>::addScalarToEntry(T s, IndexType i,
          IndexType j)
{
    // Don't add unset entries (possibly change this to use a tolerance)
    if ((s.real() == 0) && (s.imag() == 0))
        return true;

    MatrixEntry<T, IndexType> *row = m_ppElements[i];
    IndexType rSize = RowSize(i);

    bool success = false;
    int availableIdx = -1;
    for (IndexType k = 0; !success && (k < rSize); ++k) {
        if (row[k].N == j) {
            row[k].Value += s;
            success = true;
        }
        if ((availableIdx == -1) && (row[k].N == (unsigned int) -1))
            availableIdx = k;
    }

    if (!success && (availableIdx != -1))   {
        row[availableIdx].Value = s;
        row[availableIdx].N = j;
        success = true;
    }

    return success;
}

template< class T , class IndexType >
SparseMatrix< T , IndexType >::SparseMatrix( int rows )
{
	_maxColumnsPerRow = 0;
	_contiguousMemory = false;
	rows = 0;
	rowSizes = NullPointer< int >( );
	m_ppElements= NullPointer< Pointer( MatrixEntry< T , IndexType > ) >( );
	resize( rows );
}
template< class T , class IndexType >
SparseMatrix< T , IndexType >::SparseMatrix( int rows , int maxColumnsPerRow )
{
	_maxColumnsPerRow = maxColumnsPerRow;
	_contiguousMemory = false;
	rows = 0;
	rowSizes = NullPointer< int >( );
	m_ppElements= NullPointer< Pointer( MatrixEntry< T , IndexType > ) >( );
	resize( rows , maxColumnsPerRow );
}
template< class T , class IndexType >
SparseMatrix< T , IndexType >::SparseMatrix( const SparseMatrix& M )
{
	_maxColumnsPerRow = 0;
	_contiguousMemory = false;
	rowSizes = NullPointer< int >( );
	rows = 0;
	m_ppElements = NullPointer< Pointer( MatrixEntry< T , IndexType > ) >( );
	resize( M.rows );
	for( int i=0 ; i<rows ; i++ )
	{
		SetRowSize( i , M.rowSizes[i] );
		for( int j=0 ; j<rowSizes[i] ; j++ ) m_ppElements[i][j] = M.m_ppElements[i][j];
	}
}
template< class T , class IndexType >
template< class T2 , class IndexType2 >
SparseMatrix< T , IndexType >::SparseMatrix( const SparseMatrix< T2 , IndexType2 >& M )
{
	_maxColumnsPerRow = 0;
	_contiguousMemory = false;
	rowSizes = NULL;
	rows = 0;
	m_ppElements = NULL;
	resize( M.rows );
	for( int i=0 ; i<rows ; i++ )
	{
		SetRowSize( i , M.rowSizes[i] );
		for( int j=0 ; j<rowSizes[i] ; j++ ) m_ppElements[i][j] = MatrixEntry< T , IndexType >( M.m_ppElements[i][j].N , T( M.m_ppElements[i][j].Value ) );
	}
}

template< class T , class IndexType >
double SparseMatrix< T , IndexType >::SquareNorm(void) const
{
	double l2=0;
	for(int i=0;i<rows;i++)
		for(int j=0;j<rowSizes[i];j++)
			l2+=m_ppElements[i][j].Value*m_ppElements[i][j].Value;
	return l2;
}
template< class T , class IndexType >
double SparseMatrix< T , IndexType >::ASymmetricSquareNorm(void) const
{
	double l2 = 0;
	for( int i=0 ; i<rows ; i++ ) for( int j=0 ; j<rowSizes[i] ; j++ )
	{
		double t1=0 , t2=0;
		int N = m_ppElements[i][j].N;
		if( N==i ) continue;
		// [WARNING] multi-counting 
		for( int k=0 ; k<rowSizes[i] ; k++ ) if( m_ppElements[i][k].N==N ) t1 += m_ppElements[i][k].Value;
		for( int k=0 ; k<rowSizes[N] ; k++ ) if( m_ppElements[N][k].N==i ) t2 += m_ppElements[N][k].Value;
		l2 += (t1-t2)*(t1-t2);
	}
	return l2;
}
template< class T , class IndexType >
double SparseMatrix< T , IndexType >::AHermitianSquareNorm(void) const
{
	double l2=0;
	for(int i=0;i<rows;i++)
		for(int j=0;j<rowSizes[i];j++)
		{
            std::complex< double > t , t1 = 0. , t2 = 0.;
			int N=m_ppElements[i][j].N;
			if( N==i ) continue;
			t1 = m_ppElements[i][j].Value;
			for(int k=0;k<rowSizes[N];k++)
				if( m_ppElements[N][k].N==i )
				{
					t2 = m_ppElements[N][k].Value;
					t2 = std::complex< double >( t2.real() , -t2.imag() );
					break;
				}
			t = t1-t2;
			l2 += t.real()*t.real() + t.imag()*t.imag();
		}
	return l2;
}
template< class T , class IndexType >
template< class T2 , class IndexType2 >
SparseMatrix< T , IndexType >& SparseMatrix< T , IndexType >::copy( const SparseMatrix< T2 , IndexType2 >& M  )
{
	resize( M.rows );
	for ( int i=0 ; i<rows ; i++)
	{
		SetRowSize( i , M.rowSizes[i] );
		for( int j=0 ; j<rowSizes[i] ; j++ )
		{
			int idx = M.m_ppElements[i][j].N;
			m_ppElements[i][j] = MatrixEntry< T , IndexType >( idx , T( M[i][j].Value ) );
		}
	}
	return *this;
}

template< class T , class IndexType >
SparseMatrix< T , IndexType >& SparseMatrix< T , IndexType >::operator = ( const SparseMatrix< T , IndexType >& M )
{
	resize( M.rows );
	for (int i=0; i<rows; i++)
	{
		SetRowSize( i , M.rowSizes[i] );
		for( int j=0 ; j<rowSizes[i] ; j++ ) m_ppElements[i][j]=M.m_ppElements[i][j];
	}
	return *this;
}
template< class T , class IndexType >
template< class T2 , class IndexType2 >
SparseMatrix< T , IndexType >& SparseMatrix< T , IndexType >::operator = (const SparseMatrix< T2 , IndexType2 >& M)
{
	resize( M.rows );
	for ( int i=0 ; i<rows ; i++)
	{
		SetRowSize( i , M.rowSizes[i] );
		for( int j=0 ; j<rowSizes[i] ; j++ ) m_ppElements[i][j] = MatrixEntry< T , IndexType >( M.m_ppElements[i][j].N , T( M.m_ppElements[i][j].Value ) );
	}
	return *this;
}
template< class T , class IndexType >
template< class T2 >
void SparseMatrix< T , IndexType >::operator() ( const T2* in , T2* out ) const { Interface::MultiplyParallel( in , out , omp_get_num_procs() , 0 ); }


template< class T , class IndexType >
SparseMatrix< T , IndexType >::~SparseMatrix()
{
	resize( 0 );
}
template< class T , class IndexType >
bool SparseMatrix< T , IndexType >::isContiguous( void ) const { return _contiguousMemory; }
template< class T , class IndexType >
void SparseMatrix< T , IndexType >::MakeContiguous( void )
{
	if( _contiguousMemory ) return;
	Pointer( MatrixEntry< T , IndexType > ) entries = AllocPointer< MatrixEntry< T , IndexType > >( Interface::Entries() , __FILE__ , __LINE__ );

	for( int i=0 ; i<rows ; i++ )
	{
		memcpy( entries , m_ppElements[i] , sizeof( MatrixEntry< T , IndexType > ) * rowSizes[i] );
		FreePointer( m_ppElements[i] );
		m_ppElements[i] = entries;
		entries += rowSizes[i];
	}
	_contiguousMemory = true;
}

template< class T , class IndexType >
void SparseMatrix< T , IndexType >::resize( int r )
{
	if( rows>0 )
	{
		if( !_contiguousMemory ){ for( int i=0 ; i<rows ; i++ ) FreePointer( m_ppElements[i] ); }
		else                                                    FreePointer( m_ppElements[0] );
		FreePointer( m_ppElements );
		FreePointer( rowSizes );
	}
	_maxColumnsPerRow = 0;
	_contiguousMemory = false;
	rows = r;
	if( r )
	{
		rowSizes = AllocPointer< int >( r );
		memset( rowSizes , 0 , sizeof(int)*r );
		m_ppElements = AllocPointer< Pointer( MatrixEntry< T , IndexType > ) >( r );
		for( int i=0 ; i<r ; i++ ) m_ppElements[i] = NullPointer< MatrixEntry< T , IndexType > >( );
	}
}
template< class T , class IndexType >
void SparseMatrix< T , IndexType >::resize( int r , int maxColumnsPerRow )
{
	if( r==rows && maxColumnsPerRow<=_maxColumnsPerRow )
	{
		memset( rowSizes , 0 , sizeof(int)*r );
		return;
	}
	if( rows>0 )
	{
		if( !_contiguousMemory ){ for( int i=0 ; i<rows ; i++ ) FreePointer( m_ppElements[i] ); }
		else                                                    FreePointer( m_ppElements[0] );
		FreePointer( m_ppElements );
		FreePointer( rowSizes );
	}
	_contiguousMemory = true;
	rows = r;
	if( r )
	{
		rowSizes = AllocPointer< int >( r ) , memset( rowSizes , 0 , sizeof(int)*r );
		m_ppElements = AllocPointer< Pointer( MatrixEntry< T , IndexType > ) >( r );
		m_ppElements[0] = AllocPointer< MatrixEntry< T , IndexType > >( rows * maxColumnsPerRow );
		for( int i=1 ; i<r ; i++ ) m_ppElements[i] = m_ppElements[i-1] + maxColumnsPerRow;
	}
}

template< class T , class IndexType >
void SparseMatrix< T , IndexType >::SetRowSize( int row , int count )
{
	if( _contiguousMemory )
	{
		fprintf( stderr , "Cannot set row/column size in contiguous-memory mode\n" );
		return;
	}
	if( row>=0 && row<rows )
	{
		FreePointer( m_ppElements[row] );
		if( count>0 )
		{
			m_ppElements[ row ] = AllocPointer< MatrixEntry< T , IndexType > >( count );
			memset( m_ppElements[ row ] , 0 , sizeof(MatrixEntry< T , IndexType >)*count );
		}
		rowSizes[row ]=count;
	}
}

// Given matrices At and B, compute A * B.
template< class T1 , class T2 , class T3 , class IndexType >
bool TransposeMultiply( const SparseMatrix< T1 , IndexType >& At , const SparseMatrix< T2 , IndexType >& B , SparseMatrix< T3 , IndexType >& out , T1 (*TransposeFunction)( const T1& ) )
{
	int aRows = 0;
	for( int i=0 ; i<At.rows ; i++ ) for( int j=0 ; j<At.rowSizes[i] ; j++ ) if( aRows<=At[i][j].N ) aRows = At[i][j].N+1;
	return TransposeMultiply( At , B , out , aRows , TransposeFunction );
}
template< class T1 , class T2 , class T3 , class IndexType >
bool TransposeMultiply( const SparseMatrix< T1 , IndexType >& At , const SparseMatrix< T2 , IndexType >& B , SparseMatrix< T3 , IndexType >& out , int aRows , T1 (*TransposeFunction)( const T1& ) )
{
	int _aRows = 0     , aCols = At.rows;
	int bRows = B.rows , bCols = 0;
	for( int i=0 ; i<At.rows ; i++ ) for( int j=0 ; j<At.rowSizes[i] ; j++ ) if( _aRows<=At[i][j].N ) _aRows = At[i][j].N+1;
	for( int i=0 ; i< B.rows ; i++ ) for( int j=0 ; j< B.rowSizes[i] ; j++ ) if(  bCols<= B[i][j].N )  bCols =  B[i][j].N+1;
	if( _aRows>aRows )
	{
		fprintf( stderr , "[Error] prescribed output dimension too low: %d < %d\n" , aRows , _aRows );
		return false;
	}
	if( bCols>At.rows )
	{
		fprintf( stderr , "[Error] Matrix sizes do not support multiplication %d x %d * %d x %d\n" , aRows , aCols , bRows , bCols );
		return false;
	}

	std::vector< hash_map< IndexType , T3 > > rows;
	rows.resize( _aRows );
	for( int i=0 ; i<bRows ; i++ )
		for( int j=0 ; j<B.rowSizes[i] ; j++ )
		{
			IndexType idx1 = B[i][j].N;
			T2 BValue = B[i][j].Value;
			if( idx1<0 ) continue;
			for( int k=0 ; k<At.rowSizes[i] ; k++ )
			{
				T3 temp;
				IndexType idx2 = At[i][k].N;
				T1 AValue = At[i][k].Value;
				if( TransposeFunction ) temp = T3( TransposeFunction( AValue ) * BValue );
				else                    temp = T3(                    AValue   * BValue ); // temp = At( idx2 , idx1 ) * B( i , idx1 ) = A( idx1 , idx2 ) * B( i , idx1 )
				typename hash_map< IndexType , T3 >::iterator iter = rows[idx2].find(idx1);
				if( iter==rows[idx2].end() ) rows[idx2][idx1] = temp;
				else iter->second += temp;
			}
		}

	out.resize( aRows );
#pragma omp parallel for
	for( int i=0 ; i<rows.size() ; i++ )
	{
		out.SetRowSize( i , rows[i].size() );
		out.rowSizes[i] = 0;
		for( typename hash_map< IndexType , T3 >::iterator iter=rows[i].begin() ; iter!=rows[i].end() ; iter++ )
			out[i][ out.rowSizes[i]++ ] = MatrixEntry< T3 , IndexType >( iter->first , iter->second );
	}
	return true;
}
template< class A_T , class A_const_iterator , class B_T , class B_const_iterator , class Out_T , class Out_IndexType >
bool Multiply( const SparseMatrixInterface< A_T , A_const_iterator >& A , const SparseMatrixInterface< B_T , B_const_iterator >& B , SparseMatrix< Out_T , Out_IndexType >& out , int threads )
{
	size_t aCols = 0 , aRows = A.Rows();
	size_t bCols = 0 , bRows = B.Rows();
	for( int i=0 ; i<A.Rows() ; i++ ) for( A_const_iterator iter=A.begin(i) ; iter!=A.end(i) ; iter++ ) if( aCols<iter->N ) aCols = iter->N+1;
	for( int i=0 ; i<B.Rows() ; i++ ) for( B_const_iterator iter=B.begin(i) ; iter!=B.end(i) ; iter++ ) if( bCols<iter->N ) bCols = iter->N+1;
	if( bRows<aCols )
	{
		fprintf( stderr , "[Error] Matrix sizes do not support multiplication %lld x %lld * %lld x %lld\n" , (unsigned long long)aRows , (unsigned long long)aCols , (unsigned long long)bRows , (unsigned long long)bCols );
		return false;
	}

	out.resize( (int)aRows );
#pragma omp parallel for num_threads( threads )
	for( int i=0 ; i<aRows ; i++ )
	{
		hash_map< Out_IndexType , Out_T > row;
		for( A_const_iterator iterA=A.begin(i) ; iterA!=A.end(i) ; iterA++ )
		{
			Out_IndexType idx1 = iterA->N;
if( idx1==-1 ) continue;
			A_T AValue = iterA->Value;
			if( idx1<0 ) continue;
			for( B_const_iterator iterB=B.begin(idx1) ; iterB!=B.end(idx1) ; iterB++ )
			{
				Out_IndexType idx2 = iterB->N;
if( idx2==-1 ) continue;
				B_T BValue = iterB->Value;
				Out_T temp = Out_T( BValue * AValue ); // temp = A( i , idx1 ) * B( idx1 , idx2 )
				typename hash_map< Out_IndexType , Out_T >::iterator iter = row.find(idx2);
				if( iter==row.end() ) row[idx2] = temp;
				else iter->second += temp;
			}
		}
		out.SetRowSize( i , (int)row.size() );
		out.rowSizes[i] = 0;
		for( typename hash_map< Out_IndexType , Out_T >::iterator iter=row.begin() ; iter!=row.end() ; iter++ )
			out[i][ out.rowSizes[i]++ ] = MatrixEntry< Out_T , Out_IndexType >( iter->first , iter->second );
	}
	return true;
}
template< class T , class In_const_iterator , class Out_IndexType >
bool Transpose( const SparseMatrixInterface< T , In_const_iterator >& At , SparseMatrix< T , Out_IndexType >& A , T (*TransposeFunction)( const T& ) )
{
	int aRows = 0 , aCols = At.Rows();
	for( int i=0 ; i<At.Rows() ; i++ ) for( In_const_iterator iter=At.begin(i) ; iter!=At.end(i) ; iter++ ) if( aRows<=iter->N ) aRows = iter->N+1;

	A.resize( aRows );
	for( int i=0 ; i<aRows ; i++ ) A.rowSizes[i] = 0;
	for( int i=0 ; i<At.Rows() ; i++ ) for( In_const_iterator iter=At.begin(i) ; iter!=At.end(i) ; iter++ ) A.rowSizes[ iter->N ]++;
	for( int i=0 ; i<A.rows ; i++ )
	{
		int t = A.rowSizes[i];
		A.rowSizes[i] = 0;
		A.SetRowSize( i , t );
		A.rowSizes[i] = 0;
	}
	if( TransposeFunction )
		for( int i=0 ; i<At.Rows() ; i++ ) for( In_const_iterator iter=At.begin(i) ; iter!=At.end(i) ; iter++ )
		{
			int ii = iter->N;
			A[ii][ A.rowSizes[ii]++ ] = MatrixEntry< T , Out_IndexType >( i , TransposeFunction( iter->Value ) );
		}
	else
		for( int i=0 ; i<At.Rows() ; i++ ) for( In_const_iterator iter=At.begin(i) ; iter!=At.end(i) ; iter++ )
		{
			int ii = iter->N;
			A[ii][ A.rowSizes[ii]++ ] = MatrixEntry< T , Out_IndexType >( i , iter->Value );
		}
	return true;
}
template< class T , class In_const_iterator , class Out_IndexType >
bool Transpose( const SparseMatrixInterface< T , In_const_iterator >& At , SparseMatrix< T , Out_IndexType >& A , int aRows , T (*TransposeFunction)( const T& ) )
{
	size_t _aRows = 0 , aCols = At.Rows();
	for( int i=0 ; i<At.Rows() ; i++ ) for( In_const_iterator iter=At.begin(i) ; iter!=At.end(i) ; iter++ ) if( aCols<=iter->N ) _aRows = iter->N+1;
	if( _aRows>aRows )
	{
		fprintf( stderr , "[Error] prescribed output dimension too low: %d < %lld\n" , aRows , (unsigned long long)_aRows );
		return false;
	}

	A.resize( aRows );
	for( int i=0 ; i<aRows ; i++ ) A.rowSizes[i] = 0;
	for( int i=0 ; i<At.Rows() ; i++ ) for( In_const_iterator iter=At.begin(i) ; iter!=At.end(i) ; iter++ ) A.rowSizes[ iter->N ]++;
	for( int i=0 ; i<A.rows ; i++ )
	{
		int t = A.rowSizes[i];
		A.rowSizes[i] = 0;
		A.SetRowSize( i , t );
		A.rowSizes[i] = 0;
	}
	if( TransposeFunction )
		for( int i=0 ; i<At.Rows() ; i++ ) for( In_const_iterator iter=At.begin(i) ; iter!=At.end(i) ; iter++ )
		{
			int ii = iter->N;
			A[ii][ A.rowSizes[ii]++ ] = MatrixEntry< T , Out_IndexType >( i , TransposeFunction( iter->Value ) );
		}
	else
		for( int i=0 ; i<At.Rows() ; i++ ) for( In_const_iterator iter=At.begin(i) ; iter!=At.end(i) ; iter++ )
		{
			int ii = iter->N;
			A[ii][ A.rowSizes[ii]++ ] = MatrixEntry< T , Out_IndexType >( i , iter->Value );
		}
	return true;
}

inline std::complex< double > operator + ( std::complex< double > c1 , std::complex< float > c2 )
{
	return std::complex< double >( c1.real()+c2.real() , c1.imag()+c2.imag() );
}
inline std::complex< double > operator + ( std::complex< float > c1 , std::complex< double > c2 )
{
	return std::complex< double >( c1.real()+c2.real() , c1.imag()+c2.imag() );
}
inline std::complex< double > operator * ( std::complex< double > c1 , std::complex< float > c2 )
{
	return std::complex< double >( c1.real()*c2.real() - c1.imag()*c2.imag() , c1.real()*c2.imag() + c1.imag()*c2.real() );
}
inline std::complex< double > operator * ( std::complex< float > c1 , std::complex< double > c2 )
{
	return std::complex< double >( c1.real()*c2.real() - c1.imag()*c2.imag() , c1.real()*c2.imag() + c1.imag()*c2.real() );
}
inline std::complex< double > operator * ( std::complex< double > c1 , float c2 )
{
	return std::complex< double >( c1.real()*c2 , c1.imag()*c2 );
}
inline std::complex< double > operator * ( std::complex< float > c1 , double c2 )
{
	return std::complex< double >( c1.real()*c2 , c1.imag()*c2 );
}
template< class T , unsigned int Radius > BandedMatrix< T , Radius >::~BandedMatrix( void ){ _rows = 0 ; FreePointer( _entries ); }
template< class T , unsigned int Radius > BandedMatrix< T , Radius >::BandedMatrix( void ){ _rows = 0 , _entries = NullPointer< T >(); }
template< class T , unsigned int Radius > BandedMatrix< T , Radius >::BandedMatrix( size_t rows ){ _rows = 0 , _entries = NullPointer< T >() ; resize( rows ); }
template< class T , unsigned int Radius > BandedMatrix< T , Radius >::BandedMatrix( size_t rows , const T& clearValue ){ _rows = 0 , _entries = NullPointer< T >() ; resize( rows , clearValue ); }
template< class T , unsigned int Radius >
template< class T2 >
BandedMatrix< T , Radius >::BandedMatrix( const BandedMatrix< T2 , Radius >& M )
{
	_rows = 0 ; _entries = NullPointer< T >();
	resize( M._rows );
	for( size_t i=0 ; i<entries() ; i++ ) _entries[i] = (T)( M._entries[i] );
}
template< class T , unsigned int Radius >
template< class T2 >
BandedMatrix< T , Radius >& BandedMatrix< T , Radius >::operator = ( const BandedMatrix< T2 , Radius >& M )
{
	resize( M._rows );
	for( size_t i=0 ; i<entries() ; i++ ) _entries[i] = (T)( M._entries[i] );
	return *this;
}

template< class T , unsigned int Radius > void BandedMatrix< T , Radius >::resize( size_t rows )
{
	if( rows==_rows ) return;
	FreePointer( _entries );
	_rows = 0;
	if( rows )
	{
		_rows = rows;
		_entries = AllocPointer< T >( rows * ( 2 * Radius + 1 ) );
		if( !_entries ) fprintf( stderr , "[ERROR] Failed to allocate BandedMatrix::_entries ( %d x %d )\n" , rows , 2*Radius+1 ) , exit( 0 );
	}
}
template< class T , unsigned int Radius > void BandedMatrix< T , Radius >::resize( size_t rows , const T& clearValue )
{
	resize( rows );
	for( size_t i=0 ; i<entries() ; i++ ) _entries[i] = clearValue;
}
template< class T , unsigned int Radius >
template< class T2 >
void BandedMatrix< T , Radius >::multiply( ConstPointer( T2 ) in , Pointer( T2 ) out , int threads ) const
{
	for( int i=0 ; i<Radius && i<_rows-Radius ; i++ )
	{
		T2 sum(0);
		const T* __entries = _entries + i * ( 2 * Radius + 1 );
		size_t ii = i + _rows - Radius;
		for( int j=0 ; j<=2*Radius ; j++ ) sum += (T2)( in[(ii+j)%_rows] * __entries[j] );
		out[i] = sum;
	}
	if( Radius==1 )
	{
#pragma omp parallel for num_threads( threads )
		for( int i=1 ; i<_rows-1 ; i++ )
		{
			ConstPointer( T ) __entries = _entries + i * 3;
			ConstPointer( T2 ) _in = in + i - 1;
			out[i] = (T2)( _in[0] * __entries[0] + _in[1] * __entries[1] + _in[2] * __entries[2] );
		}
	}
	else
	{
#pragma omp parallel for num_threads( threads )
		for( int i=Radius ; i<_rows-Radius ; i++ )
		{
			T2 sum(0);
			ConstPointer( T ) __entries = _entries + i * ( 2 * Radius + 1 );
			ConstPointer( T2 ) _in = in + i - Radius;
			for( int j=0 ; j<=2*Radius ; j++ ) sum += (T2)( _in[j] * __entries[j] );
			out[i] = sum;
		}
	}
	for( int i=(int)_rows-Radius ; i<_rows ; i++ )
	{
		T2 sum(0);
		const T* __entries = _entries + i * ( 2 * Radius + 1 );
		int ii = (int)( i + _rows - Radius );
		for( int j=0 ; j<=2*Radius ; j++ ) sum += (T2)( in[(ii+j)%_rows] * (T2)__entries[j] );
		out[i] = sum;
	}
}
template< class T , unsigned int Radius >
template< class T2 >
void BandedMatrix< T , Radius >::multiply2( ConstPointer( T2 ) in , Pointer( T2 ) out , int threads ) const
{
	for( int i=0 ; i<Radius && i<_rows-Radius ; i++ )
	{
		T2 sum0(0) , sum1(0);
		const T* __entries = _entries + i * ( 2 * Radius + 1 );
		size_t ii = i + _rows - Radius;
		for( int j=0 ; j<=2*Radius ; j++ )
		{
			int iii = (int)( (ii + j)%_rows );
			sum0 += (T2)( in[ iii<<1   ] * __entries[j] );
			sum1 += (T2)( in[(iii<<1)|1] * __entries[j] );
		}
		out[ i<<1   ] = sum0;
		out[(i<<1)|1] = sum1;
	}
	if( Radius==1 )
	{
#pragma omp parallel for num_threads( threads )
		for( int i=1 ; i<_rows-1 ; i++ )
		{
			ConstPointer( T ) __entries = _entries + i * 3;
			ConstPointer( T2 ) _in = in + (i-1)*2;
			out[ i<<1   ] = (T2)( _in[0] * __entries[0] + _in[2] * __entries[1] + _in[4] * __entries[2] );
			out[(i<<1)|1] = (T2)( _in[1] * __entries[0] + _in[3] * __entries[1] + _in[5] * __entries[2] );
		}
	}
	else
	{
#pragma omp parallel for num_threads( threads )
		for( int i=Radius ; i<_rows-Radius ; i++ )
		{
			T2 sum0(0) , sum1(0);
			ConstPointer( T ) __entries = _entries + i * ( 2 * Radius + 1 );
			ConstPointer( T2 ) _in = in + (i-Radius)*2;
			for( int j=0 ; j<=2*Radius ; j++ ) sum0 += (T2)( _in[j<<1] * __entries[j] ) ,  sum1 += (T2)( _in[(j<<1)|1] * __entries[j] );
			out[ i<<1   ] = sum0;
			out[(i<<1)|1] = sum1;
		}
	}
	for( int i=(int)_rows-Radius ; i<_rows ; i++ )
	{
		T2 sum0(0) , sum1(0);
		const T* __entries = _entries + i * ( 2 * Radius + 1 );
		int ii = i + (int)(_rows-Radius);
		for( int j=0 ; j<=2*Radius ; j++ )
		{
			int iii = (ii+j)%_rows;
			sum0 += (T2)( in[ iii<<1   ] * __entries[j] );
			sum1 += (T2)( in[(iii<<1)|1] * __entries[j] );
		}
		out[ i<<1   ] = sum0;
		out[(i<<1)|1] = sum1;
	}
}
template< class T , unsigned int Radius >
double BandedMatrix< T , Radius >::squareNorm( void ) const
{
	double n2 = 0;
	for( int i=0 ; i<entries() ; i++ ) n2 += _entries[i] * _entries[i];
	return n2;
}

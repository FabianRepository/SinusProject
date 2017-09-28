
template< class T , class const_iterator > size_t SparseMatrixInterface< T , const_iterator >::Entries( void ) const
{
	size_t entries = 0;
	for( size_t i=0 ; i<Rows() ; i++ ) entries += RowSize( i );
	return entries;
}
template< class T , class const_iterator > double SparseMatrixInterface< T , const_iterator >::SquareNorm( void ) const
{
	double n=0;
	for( size_t i=0 ; i<Rows() ; i++ )
	{
		const_iterator e = end( i );
		for( const_iterator iter = begin( i ) ; iter!=e ; iter++ ) n += iter->Value * iter->Value;
	}
	return n;

}
template< class T , class const_iterator > double SparseMatrixInterface< T , const_iterator >::SquareASymmetricNorm( void ) const
{
	double n=0;
	for( size_t i=0 ; i<Rows() ; i++ )
	{
		const_iterator e = end( i );
		for( const_iterator iter1 = begin( i ) ; iter1!=e ; iter1++ )
		{
			int j = iter1->N;
			const_iterator e = end( j );
			double value = 0;
			for( const_iterator iter2 = begin( j ) ; iter2!=e ; iter2++ )
			{
				int k = iter2->N;
				if( k==i ) value += iter2->Value;
			}
			n += (iter1->Value-value) * (iter1->Value-value);
		}
	}
	return n;
}
template< class T , class const_iterator > double SparseMatrixInterface< T , const_iterator >::SquareASymmetricNorm( int& idx1 , int& idx2 ) const
{
	double n=0;
	double max=0;
	for( size_t i=0 ; i<Rows() ; i++ )
	{
		const_iterator e = end( i );
		for( const_iterator iter = begin( i ) ; iter!=e ; iter++ )
		{
			int j = iter->N;
			const_iterator e = end( j );
			double value = 0;
			for( const_iterator iter2 = begin( j ) ; iter2!=e ; iter2++ )
			{
				int k = iter2->N;
				if( k==i ) value += iter2->Value;
			}
			double temp = (iter->Value-value) * (iter->Value-value);
			n += temp;
			if( temp>=max ) idx1 = i , idx2 = j , max=temp;
		}
	}
	return n;
}
template< class T , class const_iterator >
template< class T2>
void SparseMatrixInterface< T , const_iterator >::Multiply( ConstPointer( T2 ) In , Pointer( T2 ) Out ) const
{
	ConstPointer( T2 ) in = In;
	for( size_t i=0 ; i<Rows() ; i++ )
	{
		T2 temp = Out[i];
		ConstPointer( T2 ) _in = in;
		const_iterator e = end( i );
		for( const_iterator iter = begin( i ) ; iter!=e ; iter++ ) temp += _in[ iter->N ] * iter->Value;
		Out[i] = temp;
	}
}
template< class T , class const_iterator >
template< class T2 >
void SparseMatrixInterface< T , const_iterator >::MultiplyParallel( ConstPointer( T2 ) In , Pointer( T2 ) Out , int threads , int multiplyFlag ) const
{
	ConstPointer( T2 ) in = In;
#pragma omp parallel for num_threads( threads )
	for( int i=0 ; i<Rows() ; i++ )
	{
		T2 temp;
		memset( &temp , 0 , sizeof(T2) );
		ConstPointer( T2 ) _in = in;
		const_iterator e = end( i );
		for( const_iterator iter = begin( i ) ; iter!=e ; iter++ ) temp += (T2)( _in[ iter->N ] * iter->Value );
		if( multiplyFlag & MULTIPLY_NEGATE ) temp = -temp;
		if( multiplyFlag & MULTIPLY_ADD ) Out[i] += temp;
		else                              Out[i]  = temp;
	}
}
template< class T , class const_iterator >
template< class T2 >
void SparseMatrixInterface< T , const_iterator >::MultiplyScaledParallel( T scale , ConstPointer( T2 ) In , Pointer( T2 ) Out , int threads , int multiplyFlag ) const
{
	ConstPointer( T2 ) in = In;
#pragma omp parallel for num_threads( threads )
	for( int i=0 ; i<Rows() ; i++ )
	{
		T2 temp;
		memset( &temp , 0 , sizeof(T2) );
		ConstPointer( T2 ) _in = in;
		const_iterator e = end( i );
		for( const_iterator iter = begin( i ) ; iter!=e ; iter++ ) temp += _in[ iter->N ] * iter->Value;
		temp *= scale;
		if( multiplyFlag & MULTIPLY_NEGATE ) temp = -temp;
		if( multiplyFlag & MULTIPLY_ADD ) Out[i] += temp;
		else                              Out[i]  = temp;
	}
}
template< class T , class const_iterator >
template< class T2 >
void SparseMatrixInterface< T , const_iterator >::MultiplyTranspose( ConstPointer( T2 ) In , Pointer( T2 ) Out , int outDim , int multiplyFlag , T (*TransposeFunction)( const T& ) ) const
{
	T2* out = &Out[0];
	if( !(multiplyFlag & MULTIPLY_ADD) ) for( int i=0 ; i<outDim ; i++ ) out[i] *= 0;
	if( multiplyFlag & MULTIPLY_NEGATE ) for( int i=0 ; i<outDim ; i++ ) out[i] = -out[i];
	for( size_t i=0 ; i<Rows() ; i++ )
	{
		T2* _out = out;
		const_iterator e = end( i );
		if( TransposeFunction ) for( const_iterator iter = begin( i ) ; iter!=e ; iter++ ) _out[ iter->N ] += In[i] * TransposeFunction( iter->Value );
		else                    for( const_iterator iter = begin( i ) ; iter!=e ; iter++ ) _out[ iter->N ] += In[i] * iter->Value;
	}
	if( multiplyFlag&MULTIPLY_NEGATE ) for( int i=0 ; i<outDim ; i++ ) out[i] = -out[i];
}
template< class T , class const_iterator >
template< class T2 >
void SparseMatrixInterface< T , const_iterator >::BMinusMX( ConstPointer( T2 ) B , ConstPointer( T2 ) X , Pointer( T2 ) Out ) const
{
	ConstPointer( T2 ) x = X;
	for( size_t i=0 ; i<Rows() ; i++ )
	{
		T2 temp;
		temp *= 0;
		ConstPointer( T2 ) _x = x;
		const_iterator e = end( i );
		for( const_iterator iter = begin( i ) ; iter!=e ; iter++ ) temp += _x[ iter->N ] * iter->Value;
		Out[i] = B[i] - temp;
	}

}

template< class T , class const_iterator >
template< class T2 , class T3 , bool UnitDiagonal , bool StrippedDiagonal , bool UseSOR >
void SparseMatrixInterface< T , const_iterator >::SolveGaussSeidel( ConstPointer( T3 ) diagonal , ConstPointer( T2 ) b , int iters , Pointer( T2 ) solution , T sor ) const
{
	for( int i=0 ; i<iters ; i++ )
		if( UnitDiagonal )
			for( size_t j=0 ; j<Rows() ; j++ )
			{
				T2 temp;
				temp *= 0;
				Pointer( T2 ) _solution = solution;
				const_iterator e = end( j );
				for( const_iterator iter = begin( j ) ; iter!=e ; iter++ ) temp += _solution[ iter->N ] * iter->Value;
				if( UseSOR )
					if(StrippedDiagonal) solution[j]  = solution[j]*( T(1.0)-sor ) + (b[j]-temp) * sor;
					else                 solution[j] +=                              (b[j]-temp) * sor;
				else
					if(StrippedDiagonal) solution[j]  = (b[j]-temp);
					else                 solution[j] += (b[j]-temp);
			}
		else
			for( size_t j=0 ; j<Rows() ; j++ )
			{
				T2 temp;
				temp *= 0;
				Pointer( T2 ) _solution = solution;
				T dValue = T(1.) / diagonal[j];
				const_iterator e = end( j );
				for( const_iterator iter = begin( j ) ; iter!=e ; iter++ ) temp += _solution[ iter->N ] * iter->Value;
				if( UseSOR )
					if(StrippedDiagonal) solution[j]  = solution[j]*( T(1.0)-sor ) + (b[j]-temp) * dValue * sor;
					else                 solution[j] +=                              (b[j]-temp) * dValue * sor;
				else
					if(StrippedDiagonal) solution[j]  = (b[j]-temp) * dValue;
					else                 solution[j] += (b[j]-temp) * dValue;
			}
}

template< class T , class const_iterator >
template< class T2 , class T3 , bool UnitDiagonal , bool StrippedDiagonal , bool UseSOR >
void SparseMatrixInterface< T , const_iterator >::SolveGaussSeidelAndResidual( ConstPointer( T3 ) diagonal , ConstPointer( T2 ) b , int iters , Pointer( T2 ) solution , Pointer( T2 ) residual , T sor ) const
{
	for( int i=0 ; i<iters ; i++ )
		if( UnitDiagonal )
			for( size_t j=0 ; j<Rows() ; j++ )
			{
				T2 temp;
				temp *= 0;
				Pointer( T2 ) _solution = solution;
				const_iterator e = end( j );
				for( const_iterator iter = begin( j ) ; iter!=e ; iter++ ) temp += _solution[ iter->N ] * iter->Value;
				if( UseSOR )
					if( StrippedDiagonal ) solution[j]  = solution[j]*T(1.0-sor) + (b[j]-temp) * sor;
					else                   solution[j] +=                          (b[j]-temp) * sor;
				else
					if( StrippedDiagonal ) solution[j]  = (b[j]-temp);
					else                   solution[j] += (b[j]-temp);
			}
		else
			for( size_t j=0 ; j<Rows() ; j++ )
			{
				T2 temp;
				temp *= 0;
				Pointer( T2 ) _solution = solution;
				T dValue = T(1.) / diagonal[j];
				const_iterator e = end( j );
				for( const_iterator iter = begin( j ) ; iter!=e ; iter++ ) temp += _solution[ iter->N ] * iter->Value;
				if( UseSOR )
					if( StrippedDiagonal ) solution[j]  = solution[j]*T(1.0-sor) + (b[j]-temp) * dValue * sor;
					else                   solution[j] +=                          (b[j]-temp) * dValue * sor;
				else
					if( StrippedDiagonal ) solution[j]  = (b[j]-temp) * dValue;
					else                   solution[j] += (b[j]-temp) * dValue;
			}
	if( UnitDiagonal )
		for( size_t j=0 ; j<Rows() ; j++ )
		{
			T2 temp;
			temp *= 0;
			Pointer( T2 ) _solution = solution;
			const_iterator e = end( j );
			for( const_iterator iter = begin( j ) ; iter!=e ; iter++ ) temp += _solution[ iter->N ] * iter->Value;
			if( StrippedDiagonal ) residual[j] = b[j] - temp - solution[j];
			else                   residual[j] = b[j] - temp;
		}
	else
		for( size_t j=0 ; j<Rows() ; j++ )
		{
			T2 temp;
			temp *= 0;
			Pointer( T2 ) _solution = solution;
			T dValue = diagonal[j];
			const_iterator e = end( j );
			for( const_iterator iter = begin( j ) ; iter!=e ; iter++ ) temp += _solution[ iter->N ] * iter->Value;
			if( StrippedDiagonal ) residual[j] = b[j] - temp - ( solution[j] * dValue );
			else                   residual[j] = b[j] - temp;
		}
}

template< class T , class const_iterator >
template< class T2 >
void SparseMatrixInterface< T , const_iterator >::MultiplyParallel( 
	ConstPointer( T2 ) in , Pointer( T2 ) out ,
	const std::vector< int >& entriesPerSlice , const std::vector< int >& slicesPerThread ) const
{
	ConstPointer( T2 ) _in = &in[0];
	int threads = slicesPerThread.size();
	if( threads<1 ) return Multiply( in , out );
	std::vector< int > startEntry;

	startEntry.resize( entriesPerSlice.size() );
	startEntry[0] = 0;
	for( int i=1 ; i<entriesPerSlice.size() ; i++ ) startEntry[i] = startEntry[i-1] + entriesPerSlice[i-1];
#pragma omp parallel for num_threads( threads )
	for( int s=0 ; s<threads ; s++ )
	{
		int startSlice = 0 , endSlice;
		for( int i=0 ; i<s ; i++ ) startSlice += slicesPerThread[i];
		endSlice = startSlice + slicesPerThread[s];
		for( int i=startSlice ; i<endSlice ; i++ )
		{
			int start = startEntry[i];
			int stop  = startEntry[i] + entriesPerSlice[i];
			for( size_t j=start ; j<stop ; j++ )
			{
				T2 temp;
				temp *= 0;
				ConstPointer( T2 ) in = _in;
				const_iterator e = end( j );
				for( const_iterator iter = begin( j ) ; iter!=e ; iter++ ) temp += in[ iter->N ] * iter->Value;
				out[j] += temp;
			}
		}
	}
}
template< class T , class const_iterator >
template< class T2 >
void SparseMatrixInterface< T , const_iterator >::MultiplyParallel( 
	ConstPointer( T2 ) in , Pointer( T2 ) out ,
	const std::vector< int >& entriesPerLine , int linesPerSlice , const std::vector< int >& slicesPerThread ) const
{
	ConstPointer( T2 ) _in = in;
	int threads = slicesPerThread.size();
	if( threads<1 ) return Multiply( in , out );
	std::vector< int > startEntry;

	startEntry.resize( entriesPerLine.size() );
	startEntry[0] = 0;
	for( int i=1 ; i<entriesPerLine.size() ; i++ ) startEntry[i] = startEntry[i-1] + entriesPerLine[i-1];
#pragma omp parallel for num_threads( threads )
	for( int s=0 ; s<threads ; s++ )
	{
		int startSlice = 0 , endSlice;
		for( int i=0 ; i<s ; i++ ) startSlice += slicesPerThread[i];
		endSlice = startSlice + slicesPerThread[s];
		for( int i=startSlice*linesPerSlice ; i<endSlice*linesPerSlice ; i++ )
		{
			int start = startEntry[i];
			int stop  = startEntry[i] + entriesPerLine[i];
			for( size_t j=start ; j<stop ; j++ )
			{
				T2 temp;
				temp *= 0;
				ConstPointer( T2 ) in = _in;
				const_iterator e = end( j );
				for( const_iterator iter = begin( j ) ; iter!=e ; iter++ ) temp += in[ iter->N ] * iter->Value;
				out[j] += temp;
			}
		}
	}
}

template< class T , class const_iterator >
template< class T2 >
void SparseMatrixInterface< T , const_iterator >::MultiplyTransposeParallel( 
	ConstPointer( T2 ) in , Pointer( T2 ) out ,
	const std::vector< int >& entriesPerSlice , const std::vector< int >& slicesPerThread , int sliceDependence ) const
{
	T2* _out = &out[0];
	int threads = slicesPerThread.size();
	if( threads<1 ) return MultiplyTranspose( in , out );
	std::vector< int > startEntry;

	startEntry.resize( entriesPerSlice.size() );
	startEntry[0] = 0;
	for( int i=1 ; i<entriesPerSlice.size() ; i++ ) startEntry[i] = startEntry[i-1] + entriesPerSlice[i-1];

#pragma omp parallel for num_threads( threads )
	for( int s=0 ; s<threads ; s++ )
	{
		int startSlice = 0 , endSlice;
		for( int i=0 ; i<s ; i++ ) startSlice += slicesPerThread[i];
		endSlice = startSlice + slicesPerThread[s];
		startSlice += sliceDependence;
		for( int i=startSlice ; i<endSlice ; i++ )
		{
			int start = startEntry[i];
			int stop = startEntry[i] + entriesPerSlice[i];
			for( int j=start ; j<stop ; j++ )
			{
				T2* out = _out;
				const_iterator e = end( j );
				for( const_iterator iter = begin( j ) ; iter!=e ; iter++ ) out[ iter->N ] += in[ j ] * iter->Value;
			}
		}
	}
	for( int s=0 ; s<threads ; s++ )
	{
		int startSlice = 0 , endSlice;
		for( int i=0 ; i<s ; i++ ) startSlice += slicesPerThread[i];
		endSlice = startSlice + sliceDependence;
		for( int i=startSlice ; i<endSlice ; i++ )
		{
			size_t start = startEntry[i];
			size_t stop = startEntry[i] + entriesPerSlice[i];
			for( int j=start ; j<stop ; j++ )
			{
				T2* out = _out;
				const_iterator e = end( j );
				for( const_iterator iter = begin( j ) ; iter!=e ; iter++ ) out[ iter->N ] += in[ j ] * iter->Value;
			}
		}
	}
}
template< class T , class const_iterator >
template< class T2 >
void SparseMatrixInterface< T , const_iterator >::MultiplyTransposeParallel( 
	ConstPointer( T2 ) in , Pointer( T2 ) out ,
	const std::vector< int >& entriesPerLine , int linesPerSlice , const std::vector< int >& slicesPerThread , int sliceDependence ) const
{
	T2* _out = &out[0];
	int threads = slicesPerThread.size();
	if( threads<1 ) return MultiplyTranspose( in , out );
	std::vector< int > startEntry;

	startEntry.resize( entriesPerLine.size() );
	startEntry[0] = 0;
	for( int i=1 ; i<entriesPerLine.size() ; i++ ) startEntry[i] = startEntry[i-1] + entriesPerLine[i-1];

#pragma omp parallel for num_threads( threads )
	for( int s=0 ; s<threads ; s++ )
	{
		int startSlice = 0 , endSlice;
		for( int i=0 ; i<s ; i++ ) startSlice += slicesPerThread[i];
		endSlice = startSlice + slicesPerThread[s];
		startSlice += sliceDependence;
		for( int i=startSlice*linesPerSlice ; i<endSlice*linesPerSlice ; i++ )
		{
			int start = startEntry[i];
			int stop = startEntry[i] + entriesPerLine[i];
			for( int j=start ; j<stop ; j++ )
			{
				T2* out = _out;
				const_iterator e = end( j );
				for( const_iterator iter = begin( j ) ; iter!=e ; iter++ ) out[ iter->N ] += in[ j ] * iter->Value;
			}
		}
	}
	for( int s=0 ; s<threads ; s++ )
	{
		int startSlice = 0 , endSlice;
		for( int i=0 ; i<s ; i++ ) startSlice += slicesPerThread[i];
		endSlice = startSlice + sliceDependence;
		for( int i=startSlice*linesPerSlice ; i<endSlice*linesPerSlice ; i++ )
		{
			size_t start = startEntry[i];
			size_t stop = startEntry[i] + entriesPerLine[i];
			for( int j=start ; j<stop ; j++ )
			{
				T2* out = _out;
				const_iterator e = end( j );
				for( const_iterator iter = begin( j ) ; iter!=e ; iter++ ) out[ iter->N ] += in[ j ] * iter->Value;
			}
		}
	}
}

template< class T , class const_iterator >
template< class T2 >
void SparseMatrixInterface< T , const_iterator >::BMinusMXParallel( 
	ConstPointer( T2 ) B , ConstPointer( T2 ) X , Pointer( T2 ) out ,
	const std::vector< int >& entriesPerSlice , const std::vector< int >& slicesPerThread ) const
{
	ConstPointer( T2 ) x = &x[0];
	int threads = slicesPerThread.size();
	if( threads<1 ) return BMinusMX( B , X , out );
	std::vector< int > startEntry;

	startEntry.resize( entriesPerSlice.size() );
	startEntry[0] = 0;
	for( int i=1 ; i<entriesPerSlice.size() ; i++ ) startEntry[i] = startEntry[i-1] + entriesPerSlice[i-1];

#pragma omp parallel for num_threads( threads )
	for( int s=0 ; s<threads ; s++ )
	{
		int startSlice = 0 , endSlice;
		for( int i=0 ; i<s ; i++ ) startSlice += slicesPerThread[i];
		endSlice = startSlice + slicesPerThread[s];
		for( int i=startSlice ; i<endSlice ; i++ )
		{
			int start = startEntry[i];
			int stop = startEntry[i] + entriesPerSlice[i];
			for( size_t j=start ; j<stop ; j++ )
			{
				T2 temp;
				temp *= 0;
				ConstPointer( T2 ) _x = x;
				const_iterator e = end( j );
				for( const_iterator iter = begin( j ) ; iter!=e ; iter++ ) temp += _x[ iter->N ] * iter->Value;
				out[j] += B[j] - temp;
			}
		}
	}
}

#ifndef USE_EIGEN
#define USE_EIGEN 0
#endif // USE_EIGEN

#ifndef USE_CHOLMOD
#define USE_CHOLMOD 0
#endif // USE_CHOLMOD

#if USE_CHOLMOD
#pragma message( "[WARNING] Need to explicitly exclude VCOMP.lib" )
#include <Cholmod/cholmod.h>
#pragma comment( lib , "CHOLMOD.lib" )
#define DLONG
#ifdef DLONG
typedef long long SOLVER_LONG;
#define CHOLMOD( name ) cholmod_l_ ## name
#else // !DLONG
typedef       int SOLVER_LONG;
#define CHOLMOD( name ) cholmod_ ## name
#endif // DLONG
#endif // USE_CHOLMOD

#include <Util/SparseMatrix.h>
double                        SquareNorm( const double* values , int dim ){ double norm2 = 0 ; for( int i=0 ; i<dim ; i++ ) norm2 += values[i] * values[i] ; return norm2; }
double                        SquareNorm( const  float* values , int dim ){ double norm2 = 0 ; for( int i=0 ; i<dim ; i++ ) norm2 += values[i] * values[i] ; return norm2; }
template< class Type > double SquareNorm( const   Type* values , int dim ){ double norm2 = 0 ; for( int i=0 ; i<dim ; i++ ) norm2 += values[dim].squareNorm()  ; return norm2 ; }

double                        SquareDifference( const double* values1 , const double* values2 , int dim ){ double norm2 = 0 ; for( int i=0 ; i<dim ; i++ ) norm2 += ( values1[i] - values2[i] ) * ( values1[i] - values2[i] ) ; return norm2; }
double                        SquareDifference( const  float* values1 , const  float* values2 , int dim ){ double norm2 = 0 ; for( int i=0 ; i<dim ; i++ ) norm2 += ( values1[i] - values2[i] ) * ( values1[i] - values2[i] ) ; return norm2; }
template< class Type > double SquareDifference( const   Type* values1 , const   Type* values2 , int dim ){ double norm2 = 0 ; for( int i=0 ; i<dim ; i++ ) norm2 += ( values1[dim] - values2[dim] ).squareNorm()  ; return norm2 ; }


// This is the conjugate gradients solver.
// The assumption is that the class SPDOperator defines a method operator()( const Real* , Real* ) which corresponds to applying a symmetric positive-definite operator.
template< class Real >
struct CGScratch
{
	Real *r , *d , *q;
	CGScratch( void ) : r(NULL) , d(NULL) , q(NULL) , _dim(0){ ; }
	CGScratch( int dim ) : r(NULL) , d(NULL) , q(NULL){ resize(dim); }
	~CGScratch( void ){ resize(0); }
	void resize( int dim )
	{
		if( dim!=_dim )
		{
			if( r ) delete[] r ; r = NULL;
			if( d ) delete[] d ; d = NULL;
			if( q ) delete[] q ; q = NULL;
			if( dim ) r = new Real[dim] , d = new Real[dim] , q = new Real[dim];
			_dim = dim;
		}
	}
protected:
	int _dim;
};
template< class Real >
struct PreconditionedCGScratch : public CGScratch< Real >
{
	Real *s;
	PreconditionedCGScratch( void ) : CGScratch< Real >() , s(NULL){ ; }
	PreconditionedCGScratch( int dim ) : CGScratch< Real >() { resize(dim); }
	~PreconditionedCGScratch( void ){ resize(0); }
	void resize( int dim )
	{
		if( dim!=CGScratch< Real >::_dim )
		{
			if( s ) delete[] s; s = NULL;
			if( dim ) s = new Real[dim];
		}
		CGScratch< Real >::resize( dim );
	}
};
template< class Real >
struct DiagonalPreconditioner
{
	Real* iDiagonal;
	DiagonalPreconditioner( void ) : iDiagonal(NULL) , _dim(0){ ; }
	~DiagonalPreconditioner( void ){ if( iDiagonal ) delete[] iDiagonal ; iDiagonal = NULL; }
	void set( const SparseMatrix< Real , int >& M )
	{
		if( _dim!=M.Rows() )
		{
			_dim = (int)M.Rows();
			if( iDiagonal ) delete[] iDiagonal , iDiagonal = NULL;
			if( _dim>0 ) iDiagonal = new Real[_dim];
		}
		memset( iDiagonal , 0 , sizeof(Real)*_dim );
#pragma omp parallel for
		for( int i=0 ; i<M.Rows() ; i++ )
		{
			for( int j=0 ; j<M.rowSizes[i] ; j++ ) if( M[i][j].N==i ) iDiagonal[i] += M[i][j].Value;
			iDiagonal[i] = (Real)1./iDiagonal[i];
		}
	}
	void operator()( const Real* in , Real* out ) const
	{
#pragma omp parallel for
		for( int i=0 ; i<_dim ; i++ ) out[i] = in[i] * iDiagonal[i];
	}
protected:
	int _dim;
};

template< class Real , class SPDOperator >
int SolveCG( SPDOperator& L , int iters , int dim , const Real* b , Real* x , CGScratch< Real >* scratch=NULL , double eps=1e-8 , int threads=1 , bool verbose=false )
{
	eps *= eps;
	Real *r , *d , *q;
	if( scratch ) r = scratch->r    , d = scratch->d    , q = scratch->q;
	else          r = new Real[dim] , d = new Real[dim] , q = new Real[dim];
	memset( r , 0 , sizeof(Real)*dim ) , memset( d , 0 , sizeof(Real)*dim ) , memset( q , 0 , sizeof(Real)*dim );
	double delta_new = 0 , delta_0;

	L( x , r );
#pragma omp parallel for num_threads( threads ) reduction( + : delta_new )
	for( int i=0 ; i<dim ; i++ ) d[i] = r[i] = b[i] - r[i] , delta_new += r[i] * r[i];

	delta_0 = delta_new;
	if( delta_new<eps )
	{
		if( !scratch ) delete[] r , delete[] d , delete[] q;
		return 0;
	}

	int ii;
	for( ii=0 ; ii<iters && delta_new>eps*delta_0 ; ii++ )
	{
		L( d , q );
        double dDotQ = 0;
#pragma omp parallel for num_threads( threads ) reduction( + : dDotQ )
		for( int i=0 ; i<dim ; i++ ) dDotQ += d[i] * q[i];
		Real alpha = Real( delta_new / dDotQ );

		double delta_old = delta_new;
		delta_new = 0;

		const int RESET_COUNT = 50;
		if( (ii%RESET_COUNT)==(RESET_COUNT-1) )
		{
#pragma omp parallel for num_threads( threads )
			for( int i=0 ; i<dim ; i++ ) x[i] += d[i] * alpha;
			L( x , r );
#pragma omp parallel for num_threads( threads ) reduction ( + : delta_new )
			for( int i=0 ; i<dim ; i++ ) r[i] = b[i] - r[i] , delta_new += r[i] * r[i];
		}
		else
#pragma omp parallel for num_threads( threads ) reduction( + : delta_new )
			for( int i=0 ; i<dim ; i++ ) r[i] -= q[i] * alpha , delta_new += r[i] * r[i] , x[i] += d[i] * alpha;

		Real beta = Real( delta_new / delta_old );
#pragma omp parallel for num_threads( threads )
		for( int i=0 ; i<dim ; i++ ) d[i] = r[i] + d[i] * beta;
	}
	if( verbose )
	{
		L( x , r );
#pragma omp parallel for num_threads( threads )
		for( int i=0 ; i<dim ; i++ ) r[i] -= b[i];
		printf( "CG: %d %g -> %g\n" , ii , SquareNorm( b , dim ) , SquareNorm( r , dim ) );
	}
	if( !scratch ) delete[] r , delete[] d , delete[] q;
	return ii;
}
template< class Real , class SPDOperator , class SPDPreconditioner >
int SolvePreconditionedCG( SPDOperator& L , SPDPreconditioner& Pinverse , int iters , int dim , const Real* b , Real* x , PreconditionedCGScratch< Real >* scratch=NULL , double eps=1e-8 , int threads=1 , bool verbose=false )
{
	eps *= eps;
	Real *r , *d , *q , *s;
	if( scratch ) r = scratch->r    , d = scratch->d    , q = scratch->q , s = scratch->s;
	else          r = new Real[dim] , d = new Real[dim] , q = new Real[dim] , s = new Real[dim];
	memset( r , 0 , sizeof(Real)*dim ) , memset( d , 0 , sizeof(Real)*dim ) , memset( q , 0 , sizeof(Real)*dim ) , memset( s , 0 , sizeof(Real)*dim );
	double delta_new = 0 , delta_0;

	L( x , r );
#pragma omp parallel for num_threads( threads )
	for( int i=0 ; i<dim ; i++ ) r[i] = b[i] - r[i];
	Pinverse( r , d );
#pragma omp parallel for num_threads( threads ) reduction( + : delta_new )
	for( int i=0 ; i<dim ; i++ ) delta_new += r[i] * d[i];

	delta_0 = delta_new;
	if( delta_new<eps )
	{
		if( !scratch ) delete[] r , delete[] d , delete[] q;
		return 0;
	}
	int ii;
	for( ii=0 ; ii<iters && delta_new>eps*delta_0 ; ii++ )
	{
		L( d , q );
        double dDotQ = 0;
#pragma omp parallel for num_threads( threads ) reduction( + : dDotQ )
		for( int i=0 ; i<dim ; i++ ) dDotQ += d[i] * q[i];
		Real alpha = Real( delta_new / dDotQ );

		const int RESET_COUNT = 50;
#pragma omp parallel for num_threads( threads )
		for( int i=0 ; i<dim ; i++ ) x[i] += d[i] * alpha;
		if( (ii%RESET_COUNT)==(RESET_COUNT-1) )
		{
			L( x , r );
#pragma omp parallel for num_threads( threads )
			for( int i=0 ; i<dim ; i++ ) r[i] = b[i] - r[i];
		}
		else
#pragma omp parallel for num_threads( threads ) reduction( + : delta_new )
			for( int i=0 ; i<dim ; i++ ) r[i] -= q[i] * alpha;
		Pinverse( r , s );

		double delta_old = delta_new;
		delta_new = 0;
#pragma omp parallel for num_threads( threads ) reduction( + : delta_new )
		for( int i=0 ; i<dim ; i++ ) delta_new += r[i] * s[i];

		Real beta = Real( delta_new / delta_old );
#pragma omp parallel for num_threads( threads )
		for( int i=0 ; i<dim ; i++ ) d[i] = s[i] + d[i] * beta;
	}
	if( verbose )
	{
		L( x , r );
#pragma omp parallel for num_threads( threads )
		for( int i=0 ; i<dim ; i++ ) r[i] -= b[i];
		printf( "PCCG: %d %g -> %g\n" , ii , SquareNorm( b , dim ) , SquareNorm( r , dim ) );
	}
	if( !scratch ) delete[] r , delete[] d , delete[] q , delete[] s;
	return ii;
}


#if USE_EIGEN
#include <Eigen/Sparse>

class EigenCholeskySolver
{
	typedef Eigen::SimplicialLLT< Eigen::SparseMatrix< double > > Eigen_Solver;
	typedef Eigen::VectorXd                                       Eigen_Vector;
	Eigen_Solver _solver;
	Eigen_Vector _eigenB , _eigenX;
public:
	template< class Real >
	EigenCholeskySolver( const SparseMatrix< Real , int >& M )
	{
		Eigen::SparseMatrix< double > eigenM( int( M.Rows() ) , int( M.Rows() ) );
		std::vector< Eigen::Triplet<double> > triplets;
		triplets.reserve( M.Entries() );
		for( int i=0 ; i<M.Rows() ; i++ ) for( int j=0 ; j<M.RowSize(i) ; j++ ) triplets.push_back( Eigen::Triplet< double >( i , M[i][j].N , M[i][j].Value ) );
		eigenM.setFromTriplets( triplets.begin() , triplets.end() );
		_solver.analyzePattern( eigenM );
		_solver.factorize( eigenM );
		if( _solver.info()!=Eigen::Success ) fprintf( stderr , "[ERROR] EigenSovler::EigenCholeskySolver Failed to factorize matrix\n" ) , exit(0);
		_eigenB.resize( M.Rows() ) , _eigenX.resize( M.Rows() );
	}
	template< class Real >
	void solve( const Real* b , Real* x )
	{
#pragma omp parallel for
		for( int i=0 ; i<_eigenB.size() ; i++ ) _eigenB[i] = b[i];
		_eigenX = _solver.solve( _eigenB );
#pragma omp parallel for
		for( int i=0 ; i<_eigenX.size() ; i++ ) x[i] = _eigenX[i];
	}
	size_t dimension( void ) const { return _eigenB.size(); }
	template< class Real >
	static void Solve( const SparseMatrix< Real , int >& M , const Real* b , Real* x ){ EigenCholeskySolver solver( M ) ; solver.solve( b , x ); }
};
class EigenCGSolver
{
	Eigen::ConjugateGradient< Eigen::SparseMatrix< double > > _solver;
	Eigen::VectorXd _eigenB , _eigenX;
	Eigen::SparseMatrix< double > _eigenM;
public:
	template< class Real >
	EigenCGSolver( const SparseMatrix< Real , int >& M )
	{
		_eigenM.resize( (int)M.Rows() , (int)M.Rows() );
		std::vector< Eigen::Triplet<double> > triplets;
		triplets.reserve( M.Entries() );
		for( int i=0 ; i<M.Rows() ; i++ ) for( int j=0 ; j<M.RowSize(i) ; j++ ) triplets.push_back( Eigen::Triplet< double >( i , M[i][j].N , M[i][j].Value ) );
		_eigenM.setFromTriplets( triplets.begin() , triplets.end() );
		_solver.compute( _eigenM );
		_solver.analyzePattern( _eigenM );
		if( _solver.info()!=Eigen::Success ) fprintf( stderr , "[ERROR] EigenSovler::EigenCGSolver Failed to factorize matrix\n" ) , exit(0);
		_eigenB.resize( M.Rows() ) , _eigenX.resize( M.Rows() );
	}
	template< class Real >
	void solve( const Real* b , Real* x , int iters )
	{
		_solver.setMaxIterations( iters );
#pragma omp parallel for
		for( int i=0 ; i<_eigenB.size() ; i++ ) _eigenB[i] = b[i] , _eigenX[i] = x[i];
		_eigenX = _solver.solveWithGuess( _eigenB , _eigenX );
#pragma omp parallel for
		for( int i=0 ; i<_eigenX.size() ; i++ ) x[i] = _eigenX[i];
	}
	size_t dimension( void ) const { return _eigenB.size(); }
	template< class Real >
	static void Solve( const SparseMatrix< Real , int >& M , const Real* b , Real* x , int iters ){ EigenCGSolver solver( M ) ; solver._solver.setMaxIterations( iters ) ; solver.solve( b , x ); }
};
#endif // USE_EIGEN

#if USE_CHOLMOD
class CholmodSolver
{
	const static bool LOWER_TRIANGULAR = true;
	int dim;
	cholmod_factor* cholmod_L;
	cholmod_dense*  cholmod_b;
	cholmod_sparse* cholmod_M;
	std::vector< bool > flaggedValues;
	template< class Real > void   _init( const SparseMatrix< Real , int >& M );
	template< class Real > bool _update( const SparseMatrix< Real , int >& M );
public:
	static cholmod_common cholmod_C;
	static bool cholmod_C_set;

	template< class Real >
	CholmodSolver( const SparseMatrix< Real , int >& M );
	~CholmodSolver( void );

	template< class Real > void solve( ConstPointer( Real ) b , Pointer( Real ) x );
	int nonZeros( void ) const;

};
bool CholmodSolver::cholmod_C_set = false;
cholmod_common CholmodSolver::cholmod_C;

template< class Real > CholmodSolver::CholmodSolver( const SparseMatrix< Real , int >& M ){ _init( M ) , _update( M ); }
template< class Real >
void CholmodSolver::_init( const SparseMatrix< Real , int >& M )
{
	{
		if( !cholmod_C_set ) CHOLMOD(start)( &cholmod_C );
		cholmod_C_set = true;
	}
	dim = M.rows;

	int maxEntries;
	if( LOWER_TRIANGULAR )
	{
		maxEntries = (int)( ( M.Entries()-M.rows ) / 2 + M.rows );
		cholmod_M = CHOLMOD(allocate_sparse)( dim , dim , maxEntries , 0 , 1 , -1 , CHOLMOD_REAL , &cholmod_C );
	}
	else
	{
		maxEntries = (int)M.Entries();
		cholmod_M = CHOLMOD(allocate_sparse)( dim , dim , maxEntries , 0 , 1 ,  0 , CHOLMOD_REAL , &cholmod_C );
	}
	cholmod_M->i = malloc( sizeof( SOLVER_LONG ) * maxEntries );
	cholmod_M->x = malloc( sizeof( double ) * maxEntries );

	SOLVER_LONG *_p = (SOLVER_LONG*)cholmod_M->p;
	SOLVER_LONG *_i = (SOLVER_LONG*)cholmod_M->i;

	int off = 0;
	dim = 0;

	for( int i=0 ; i<M.rows ; i++ )
	{
		_p[dim++] = off;
		for( int j=0 ; j<M.rowSizes[i] ; j++ ) if( !LOWER_TRIANGULAR || M[i][j].N>=i ) _i[off++] = M[i][j].N;
	}
	_p[dim] = off;

	cholmod_L = CHOLMOD(analyze)( cholmod_M , &cholmod_C );
	cholmod_b = CHOLMOD(allocate_dense)( dim , 1 , dim , cholmod_M->xtype , &cholmod_C );
}
template< class Real >
bool CholmodSolver::_update( const SparseMatrix< Real , int >& M )
{
	double *_x = (double*)cholmod_M->x;
	int off = 0;

	SOLVER_LONG *_p = (SOLVER_LONG*)cholmod_M->p;
#pragma omp parallel for
	for( int i=0 ; i<M.rows ; i++ )
	{
		int off = (int)_p[i];
		for( int j=0 ; j<M.rowSizes[i] ; j++ ) if( !LOWER_TRIANGULAR || M[i][j].N>=i ) _x[off++] = double( M[i][j].Value );
	}

	cholmod_C.print = 0;
	CHOLMOD(factorize)( cholmod_M , cholmod_L , &cholmod_C );
	if( cholmod_C.status==CHOLMOD_NOT_POSDEF )
	{
		fprintf( stderr , "[WARNING] Matrix not positive-definite\n" );
		return false;
	}
	else if( cholmod_C.status==CHOLMOD_OUT_OF_MEMORY )
	{
		fprintf( stderr , "[WARNING] CHOLMOD ran out of memory\n" );
		return false;
	}
	else if( cholmod_C.status!=CHOLMOD_OK )
	{
		fprintf( stderr , "[WARNING] CHOLMOD status not OK: %d\n" , cholmod_C.status );
		return false;
	}
	return true;
}
CholmodSolver::~CholmodSolver( void )
{
	if( cholmod_L ) CHOLMOD(free_factor)( &cholmod_L , &cholmod_C ) , cholmod_L = NULL;
	if( cholmod_b ) CHOLMOD(free_dense )( &cholmod_b , &cholmod_C ) , cholmod_b = NULL;
	if( cholmod_M ) CHOLMOD(free_sparse)( &cholmod_M , &cholmod_C ) , cholmod_M = NULL;
}

template< class Real >
void CholmodSolver::solve( ConstPointer( Real ) b , Pointer( Real ) x )
{
	double* _b = (double*)cholmod_b->x;
	for( int i=0 ; i<dim ; i++ ) _b[i] = (double)b[i];

	cholmod_dense* cholmod_x = CHOLMOD(solve)( CHOLMOD_A , cholmod_L , cholmod_b , &cholmod_C );
	double* _x = (double*)cholmod_x->x;
	for( int i=0 ; i<dim ; i++ ) x[i] = (Real)_x[i];

	CHOLMOD(free_dense)( &cholmod_x , &cholmod_C );
}
int CholmodSolver::nonZeros( void ) const
{
	long long nz = 0;
	if( cholmod_L->xtype != CHOLMOD_PATTERN && !(cholmod_L->is_super ) ) for( int i=0 ; i<cholmod_L->n ; i++ ) nz += ((SOLVER_LONG*)cholmod_L->nz)[i];
	bool examine_super = false;
	if( cholmod_L->xtype != CHOLMOD_PATTERN ) examine_super = true ;
	else                                      examine_super = ( ((int*)cholmod_L->s)[0] != (-1));
	if( examine_super )
	{
		/* check and print each supernode */
		for (int s = 0 ; s < cholmod_L->nsuper ; s++)
		{
			int k1 = ((int*)cholmod_L->super) [s] ;
			int k2 = ((int*)cholmod_L->super) [s+1] ;
			int psi = ((int*)cholmod_L->pi)[s] ;
			int psend = ((int*)cholmod_L->pi)[s+1] ;
			int nsrow = psend - psi ;
			int nscol = k2 - k1 ;
			nz += nscol * nsrow - (nscol*nscol - nscol)/2 ;
		}
	}
	return (int)nz;
}
#endif // USE_CHOLMOD
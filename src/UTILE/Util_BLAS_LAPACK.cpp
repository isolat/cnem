///////////////////////////////////////////////////////////////////////////////
//                                                                           // 
// illoul_amran@yahoo.fr                                                     // 
//                                                                           // 
// ILLOUL Lounes, LMSP, ENSAM Paris.                                         //
///////////////////////////////////////////////////////////////////////////////

//---------------------------------------------------------------------------//
// fonction utile BLAS et LAPACK 

#include "Util_BLAS_LAPACK.h"

//---------------------------------------------------------------------------//

int Valeurs_Vecteurs_Propres_MKL_LAPACK_33(double* A,double* V,double* W) 
{
    //-----------------------------------------------------------------------//
    /*

    void DSYEV (char *jobz,char *uplo,int *n,double *a,int *lda,double *w,double *work,int *lwork,int *info)

    Input Parameters
    
    jobz CHARACTER*1. Must be 'N' or 'V'.
    If jobz = 'N', then only eigenvalues are computed.
    If jobz = 'V', then eigenvalues and eigenvectors are
    computed.
    uplo CHARACTER*1. Must be 'U' or 'L'.
    If uplo = 'U', a stores the upper triangular part of A.
    If uplo = 'L', a stores the lower triangular part of A.
    n INTEGER. The order of the matrix A (n = 0).
    a, work REAL for ssyev
    DOUBLE PRECISION for dsyev
    Arrays:
    a(lda,*) is an array containing either upper or lower
    triangular part of the symmetric matrix A, as specified by
    uplo.
    The second dimension of a must be at least max(1, n).
    work is a workspace array, its dimension max(1, lwork).
    lda INTEGER. The first dimension of the array a.
    Must be at least max(1, n).
    lwork INTEGER.
    The dimension of the array work.
    Constraint: lwork = max(1, 3n-1).
    If lwork = -1, then a workspace query is assumed; the
    routine only calculates the optimal size of the work array,
    returns this value as the first entry of the work array, and
    no error message related to lwork is issued by xerbla.
    See Application Notes for the suggested value of lwork.

    Output Parameters
    On exit, if jobz = 'V', then if info = 0, array a contains
    the orthonormal eigenvectors of the matrix A.
    a
    If jobz = 'N', then on exit the lower triangle
    (if uplo = 'L') or the upper triangle (if uplo = 'U') of
    A, including the diagonal, is overwritten.
    w REAL for ssyev
    DOUBLE PRECISION for dsyev
    Array, DIMENSION at least max(1, n).
    If info = 0, contains the eigenvalues of the matrix A in
    ascending order.
    On exit, if lwork > 0, then work(1) returns the required
    minimal size of lwork.
    work(1)
    info INTEGER.
    If info = 0, the execution is successful.
    If info = -i, the i-th parameter had an illegal value.
    If info = i, then the algorithm failed to converge; i
    indicates the number of elements of an intermediate
    tridiagonal form which did not converge to zero.

    */
    //-----------------------------------------------------------------------//
    // Vecteur propre en lignes c
    // ...............en colone f
    //-----------------------------------------------------------------------//

    MKL_INT i;
    for(i=0;i<3;i++)
    {
        MKL_INT j;for(j=0;j<=i;j++)V[3*i+j]=A[3*i+j];
    }

    char Jobz='V';
    char Uplo='U';
    MKL_INT N=3;
    MKL_INT Lda=N;
    
    MKL_INT LWork=102;
    double Work[102];
    MKL_INT Info;

    DSYEV(&Jobz,&Uplo,&N,V,&Lda,W,Work,&LWork,&Info);
    //cout<<"Work(0) : "<<Work[0]<<endl;
    
    return Info;   
}

int Inversion_Matrice_MKL_LAPACK_33(double* A,double* Inv_A) 
{
    //-----------------------------------------------------------------------//
    /*
    Input Parameters
    m INTEGER. The number of rows in the matrix A (m = 0).
    n INTEGER. The number of columns in A; n = 0.
    a REAL for sgetrf
    DOUBLE PRECISION for dgetrf
    COMPLEX for cgetrf
    DOUBLE COMPLEX for zgetrf.
    Array, DIMENSION (lda,*). Contains the matrix A. The second
    dimension of a must be at least max(1, n).
    lda INTEGER. The first dimension of array a.
    
    Output Parameters
    a Overwritten by L and U. The unit diagonal elements of L are not stored.
    ipiv INTEGER.
    Array, DIMENSION at least max(1,min(m, n)). The pivot indices: row
    i was interchanged with row ipiv(i).
    info INTEGER. If info=0, the execution is successful.
    If info = -i, the i-th parameter had an illegal value.
    If info = i, uii is 0. The factorization has been completed, but U is
    exactly singular. Division by 0 will occur if you use the factor U for
    solving a system of linear equations.
    
    */
    //-----------------------------------------------------------------------//
    
    MKL_INT M=3;
    MKL_INT N=M;
    MKL_INT LDA=M;

    double Inv_A_F[9];

    MKL_INT i;
    for(i=0;i<M;i++){MKL_INT j;for(j=0;j<N;j++)Inv_A_F[3*i+j]=A[3*j+i];}

    MKL_INT IPIV[3]={0,0,0};
    MKL_INT INFO;

    DGETRF(&M,&N,Inv_A_F,&LDA,IPIV,&INFO);
    
    MKL_INT LWORK=3;
    double WORK[3];
    
    DGETRI(&N,Inv_A_F,&LDA,IPIV,WORK,&LWORK,&INFO);
    //cout<<"WORK(0)="<<WORK[0]<<endl;

    for(i=0;i<M;i++){MKL_INT j;for(j=0;j<N;j++)Inv_A[3*i+j]=Inv_A_F[3*j+i];}

    return INFO;
}

int Inversion_Matrice_MKL_LAPACK_22(double* A,double* Inv_A) 
{    
    MKL_INT M=2;
    MKL_INT N=M;
    MKL_INT LDA=M;

    double Inv_A_F[4];

    MKL_INT i;
    for(i=0;i<M;i++){MKL_INT j;for(j=0;j<N;j++)Inv_A_F[2*i+j]=A[2*j+i];}

    MKL_INT IPIV[2]={0,0};
    MKL_INT INFO;

    DGETRF(&M,&N,Inv_A_F,&LDA,IPIV,&INFO);
    
    MKL_INT LWORK=2;
    double WORK[2];
    
    DGETRI(&N,Inv_A_F,&LDA,IPIV,WORK,&LWORK,&INFO);
    //cout<<"WORK(0)="<<WORK[0]<<endl;

    for(i=0;i<M;i++){MKL_INT j;for(j=0;j<N;j++)Inv_A[2*i+j]=Inv_A_F[2*j+i];}

    return INFO;
}

void Prod_MM_cblas(double *A, double *B, double* C,bool tA,bool tB,MKL_INT m,MKL_INT n,MKL_INT k)
{
    double alpha=1.;
    double beta=0.;
        
    CBLAS_ORDER order=CblasRowMajor;// c,c++ 
    
    CBLAS_TRANSPOSE transA, transB;
    MKL_INT lda, ldb, ldc;
      
    if(tA)
    {
        lda=m;
        transA=CblasTrans;
    }
    else
    {
        lda=k;
        transA=CblasNoTrans;
    }
    
    if(tB)
    {
        ldb=k;
        transB=CblasTrans;
    }
    else
    {
        ldb=n;
        transB=CblasNoTrans;
    }
        
    ldc=n;

    cblas_dgemm(order,transA,transB,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc);
}

void Prod_MV_cblas(double *A, double *b, double* c,bool tA,MKL_INT m,MKL_INT n)
{
    double alpha=1.;
    double beta=0.;
        
    CBLAS_ORDER order=CblasRowMajor;// c,c++ 
    
    CBLAS_TRANSPOSE transA;
    MKL_INT lda=n;
      
    if(tA)
    {
        transA=CblasTrans;
    }
    else
    {
        transA=CblasNoTrans;
    }

    cblas_dgemv(order,transA,m,n,alpha,A,lda,b,1,beta,c,1);
}

void Prod_M3M3(double* M_0,double* M_1,double* P)
{
    long i;
    for(i=0;i<3;i++)
    {
        long j;
        for(j=0;j<3;j++)
        {
            P[3*i+j]=0.;
            long k;
            for(k=0;k<3;k++)
                P[3*i+j]+=M_0[3*i+k]*M_1[3*k+j];
        }
    }
}

void Prod_tM3M3(double* M_0,double* M_1,double* P)
{
    long i;
    for(i=0;i<3;i++)
    {
        long j;
        for(j=0;j<3;j++)
        {
            P[3*i+j]=0.;
            long k;
            for(k=0;k<3;k++)
                P[3*i+j]+=M_0[3*k+i]*M_1[3*k+j];
        }
    }
}

void Prod_M3tM3(double* M_0,double* M_1,double* P)
{
    long i;
    for(i=0;i<3;i++)
    {
        long j;
        for(j=0;j<3;j++)
        {
            P[3*i+j]=0.;
            long k;
            for(k=0;k<3;k++)
                P[3*i+j]+=M_0[3*i+k]*M_1[3*j+k];
        }
    }
}

double Det_M33(double* M)
{
    double Det=M[0]*(M[4]*M[8]-M[5]*M[7])+
               M[1]*(M[5]*M[6]-M[3]*M[8])+
               M[2]*(M[3]*M[7]-M[4]*M[6]);
    return Det;
}

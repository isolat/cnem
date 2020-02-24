///////////////////////////////////////////////////////////////////////////////
//                                                                           // 
// illoul_amran@yahoo.fr                                                     // 
//                                                                           // 
// ILLOUL Lounes, LMSP, ENSAM Paris.                                         //
///////////////////////////////////////////////////////////////////////////////

#pragma once
//---------------------------------------------------------------------------//
// fonction utile BLAS et LAPACK 

#include <mkl_rci.h>
#include <mkl_blas.h>
#include <mkl_spblas.h>
#include <mkl_dss.h>
#include <mkl_cblas.h>
#include <mkl_lapack.h>
#include <mkl_vsl.h>

void Prod_MV_cblas(double *A, double *b, double* c,bool tA,MKL_INT m,MKL_INT n);
void Prod_MM_cblas(double *A, double *B, double* C,bool tA,bool tB,MKL_INT m,MKL_INT n,MKL_INT k);

int Valeurs_Vecteurs_Propres_MKL_LAPACK_33(double* A,double* V,double* W);
int Inversion_Matrice_MKL_LAPACK_33(double* A,double* Inv_A);
int Inversion_Matrice_MKL_LAPACK_22(double* A,double* Inv_A);
void Prod_M3M3(double* M_0,double* M_1,double* P);
void Prod_tM3M3(double* M_0,double* M_1,double* P);
void Prod_M3tM3(double* M_0,double* M_1,double* P);
double Det_M33(double* M);

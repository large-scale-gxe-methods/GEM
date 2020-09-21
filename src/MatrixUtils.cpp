/****************************************************************************
  MatrixUtils.cpp provides routines for matrix or vector multiplication
****************************************************************************/

#include <cstdio>
#include <stdlib.h>

extern "C"{
  // product C= alphaA.B + betaC                                               
 void dgemm_(char* TRANSA, char* TRANSB, const int* M,
             const int* N, const int* K, double* alpha, double* A,
             const int* LDA, double* B, const int* LDB, double* beta,
             double* C, const int* LDC);
  // product Y= alphaA.X + betaY                                               
 void dgemv_(char* TRANS, const int* M, const int* N,
             double* alpha, double* A, const int* LDA, double* X,
             const int* INCX, double* beta, double* C, const int* INCY);

 void dgetrf_(const int* M, const int* N, double* A, const int* LDA, 
              const int* IPIV, const int* INFO);

 void dgetri_(const int* N, double* A, const int* LDA, const int* IPIV,
              double* WORK, const int* LWORK, const int* INFO);

 void daxpy_(const int* N, double* DA, double* DX, const int* INCX, 
             double* DY, const int* INCY);

 void dlacpy_(char* uplo, const int* m, const int* n, double* a,
              const int* lda, double* b, const int* ldb);

} 

void initvec(double* v, int N){
  for(int i= 0; i<N; ++i){
    v[i]= 0.0;
  }
}

void matTmatprod(double* A, double* v, double* u, int Nrow, int Ncol, int NcolB){
  double alpha= 1.0, beta= 0.0;
  char no= 'N', tr= 'T';
  int m= Ncol, n= NcolB, k=Nrow, lda= Ncol, incx= NcolB, incy= Ncol;
  double* tmp= new double[Ncol*NcolB];
  initvec(tmp, Ncol*NcolB);
  dgemm_(&no,&tr,&m,&n,&k,&alpha,A,&lda,v,&incx,&beta,tmp,&incy);
  for(int i= 0; i<Ncol*NcolB; ++i){
    u[i]= tmp[i];
  }
  delete [] tmp;
}

void matNmatNprod(double* A, double* v, double* u, int Nrow, int Ncol, int NcolB){
  double alpha= 1.0, beta= 0.0;
  char no= 'N';
  int m= Nrow, n= NcolB, k=Ncol, lda= Nrow, incx= Ncol, incy= Nrow;
  double* tmp= new double[Nrow*NcolB];
  initvec(tmp, Nrow*NcolB);
  dgemm_(&no,&no,&m,&n,&k,&alpha,A,&lda,v,&incx,&beta,tmp,&incy);
  for(int i= 0; i<Nrow*NcolB; ++i){
    u[i]= tmp[i];
  }
  delete [] tmp;
}

void matmatprod(double* A, double* v, double* u, int Nrow, int Ncol, int NcolB){
  double alpha= 1.0, beta= 0.0;
  char tr= 'T';
  int m= Nrow, n= NcolB, k=Ncol, lda= Ncol, incx= NcolB, incy= Nrow;
  double* tmp= new double[Nrow*NcolB];
  initvec(tmp, Nrow*NcolB);
  dgemm_(&tr,&tr,&m,&n,&k,&alpha,A,&lda,v,&incx,&beta,tmp,&incy);
  for(int i= 0; i<Nrow*NcolB; ++i){
    u[i]= tmp[i];
  }
  delete [] tmp;
}

void matmatTprod(double* A, double* v, double* u, int Nrow, int Ncol, int NcolB){
  double alpha= 1.0, beta= 0.0;
  char no= 'N', tr= 'T';
  int m= Nrow, n= NcolB, k=Ncol, lda= Ncol, incx= Ncol, incy= Nrow;
  double* tmp= new double[Nrow*NcolB];
  initvec(tmp, Nrow*NcolB);
  dgemm_(&tr,&no,&m,&n,&k,&alpha,A,&lda,v,&incx,&beta,tmp,&incy);
  for(int i= 0; i<Nrow*NcolB; ++i){
    u[i]= tmp[i];
  }
  delete [] tmp;
}

void matTvecprod(double* A, double* v, double* u, int Nrow, int Ncol){
  double alpha= 1.0, beta= 0.0;
  char no= 'N', tr= 'T';
  int m= Ncol, n= 1, k=Nrow, lda= Ncol, incx= 1, incy= Ncol;
  double* tmp= new double[Ncol];
  initvec(tmp, Ncol);
  dgemm_(&no,&tr,&m,&n,&k,&alpha,A,&lda,v,&incx,&beta,tmp,&incy);
  for(int i= 0; i<Ncol; ++i){
    u[i]= tmp[i];
  }
  delete [] tmp;
}

void matvecprod(double* A, double* v, double* u, int Nrow, int Ncol){
  double alpha= 1.0, beta= 0.0;
  char tr= 'T';
  int m= Nrow, n= 1, k=Ncol, lda= Ncol, incx= 1, incy= Nrow;
  double* tmp= new double[Nrow];
  initvec(tmp, Nrow);
  dgemm_(&tr,&tr,&m,&n,&k,&alpha,A,&lda,v,&incx,&beta,tmp,&incy);
  for(int i= 0; i<Nrow; ++i){
    u[i]= tmp[i];
  }
  delete [] tmp;
}

void vecmatprod(double* v, double* A, double* u, int N){
  double alpha= 1.0, beta= 0.0;
  char no= 'N', tr= 'T';
  int m= 1, n= N, k= N, lda= 1, ldb= N, ldc= 1;
  double* tmp= new double[N];
  initvec(tmp, N);
  dgemm_(&no,&tr,&m,&n,&k,&alpha,v,&lda,A,&ldb,&beta,tmp,&ldc);
  for(int i= 0; i<N; ++i){ 
    u[i]= tmp[i];
  }
  delete [] tmp;
}

void matInv(double* A, int N) {
  int info = 0;
  int lworkspace = N;
  int *ipiv = new int[N];
  dgetrf_(&N, &N, A, &N, ipiv, &info);
  double* workspace = new double[lworkspace*sizeof(double)];
  dgetri_(&N, A, &N, ipiv, workspace, &lworkspace, &info);
  delete [] ipiv;
  delete [] workspace;
}

void matAdd(double* A, double* B, int N, double alpha) {
  int length = N;
  int incOne = 1;
  daxpy_(&length, &alpha, B, &incOne, A, &incOne); 
}


void subMatrix(double* A, double* u, int M, int N, int NrowA, int NrowB, int start) {
    char uplo = 'A';
    dlacpy_(&uplo, &M, &N, &A[start], &NrowA, u, &NrowB);
}
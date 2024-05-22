#ifndef MATRIXUTILS_H
#define MATRIXUTILS_H

void initvec(double* v, int N);

void matTmatprod(double* A, double* v, double* u, int Nrow, int Ncol, int NcolB);

void matmatprod(double* A, double* v, double* u, int Nrow, int Ncol, int NcolB);

void matNmatNprod(double* A, double* v, double* u, int Nrow, int Ncol, int NcolB);

void matmatTprod(double* A, double* v, double* u, int Nrow, int Ncol, int NcolB);

void matTvecprod(double* A, double* v, double* u, int Nrow, int Ncol);

void matvecprod(double* A, double* v, double* u, int Nrow, int Ncol);

void vecmatprod(double* v, double* A, double* u, int N);

void matInv(double* A, int N);

void matAdd(double* A, double* B, int N, double alpha);

void subMatrix(double* A, double* u, int M, int N, int NrowA, int NrowB, int start);

void matvecSprod(double* A, double* v, double* u, int Nrow, int Ncol, int start);
#endif

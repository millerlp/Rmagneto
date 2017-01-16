/*
	Rmagneto.c, a program to calculate bias corrections for
	magnetometer and accelerometer data using an algorithm derived from
	Li, Q. and J.G. Griffiths (2004). Least squares ellipsoid specific fitting. Geometric Modeling and Processing, 2004. Proceedings, IEEE
	as implemented by Merlin Oz 2013.
	https://sites.google.com/site/sailboatinstruments1/home
	Copyright (C) 2013 www.sailboatinstruments.blogspot.com

 Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

Modified for use with R 3.1.3 by Luke Miller July 2015

*/

#include <R.h>
#include "Rmagneto.h"


////////////////////////////////////////////////////////////////////////////////////
//	The main magneto ellipsoid fitting (calibration) routine. 
//  Original code by Merlin Oz
//	From https://sites.google.com/site/sailboatinstruments1/ 

// For the call below, from R we will expect to receive pointers to 3 input vectors
// X, Y, Z, which for now will be double precision floating point values.
// norm will be a floating point value for the local magnetic field total norm (or 
// acceleration norm of 1000 milli-g). 
// b will be an empty 1x3 vector to hold the offset corrections output data
// Xcorr will be an empty 1x3 vector to hold the coefficients for the X-axis 
// scaling + soft iron correction
// Ycorr will be an empty 1x3 vector to hold the coefficients for the Y-axis 
// scaling + soft iron correction
// Zcorr will be an empty 1x3 vector to hold the coefficients for the Z-axis 
// scaling + soft iron correction

void calibrate(double *X, double *Y, double *Z, int *nrows, double *localnorm, double *b, double *Xcorr, double *Ycorr, double *Zcorr)
{
 int nlines = *nrows; // get a copy of the length of the input vectors X,Y,Z
 double *D, *S, *C, *S11, *S12, *S12t, *S22, *S22_1, *S22a, *S22b, *SS, *E, *U, *SSS;
 double *eigen_real, *eigen_imag, *v1, *v2, *v, *Q, *Q_1, *B, *QB, J, hmb, *SSSS;
 int i, index;
 double maxval, norm, btqb, *eigen_real3, *eigen_imag3, *Dz, *vdz, *SQ, *A_1, hm, norm1, norm2, norm3;

 
 D = (double*)malloc(10 * nlines * sizeof(double)); // allocate a vector to hold data
 for( i = 0; i < nlines; i++)
 {
										// Original code
	D[i] = X[i] * X[i]; 				// D[i] = x * x;
	D[nlines+i] = Y[i] * Y[i]; 			// D[nlines+i] = y * y;
	D[nlines*2+i] = Z[i] * Z[i]; 		// D[nlines*2+i] = z * z;
	D[nlines*3+i] = 2.0 * Y[i] * Z[i];	// D[nlines*3+i] = 2.0 * y * z;
	D[nlines*4+i] = 2.0 * X[i] * Z[i];	// D[nlines*4+i] = 2.0 * x * z;
	D[nlines*5+i] = 2.0 * X[i] * Y[i];	// D[nlines*5+i] = 2.0 * x * y;
	D[nlines*6+i] = 2.0 * X[i]; 		// D[nlines*6+i] = 2.0 * x;
	D[nlines*7+i] = 2.0 * Y[i];			// D[nlines*7+i] = 2.0 * y;
	D[nlines*8+i] = 2.0 * Z[i];			// D[nlines*8+i] = 2.0 * z;
	D[nlines*9+i] = 1.0;				// D[nlines*9+i] = 1.0;
 }

 
 // allocate memory for matrix S
 S = (double*)malloc(10 * 10 * sizeof(double));
 Matrix_x_Its_Transpose(S, D, 10, nlines);
 
 // Create pre-inverted constraint matrix C
 C = (double*)malloc(6 * 6 * sizeof(double));
 C[0] = 0.0; C[1] = 0.5; C[2] = 0.5; C[3] = 0.0;  C[4] = 0.0;  C[5] = 0.0;
 C[6] = 0.5;  C[7] = 0.0; C[8] = 0.5; C[9] = 0.0;  C[10] = 0.0;  C[11] = 0.0;
 C[12] = 0.5;  C[13] = 0.5; C[14] = 0.0; C[15] = 0.0;  C[16] = 0.0;  C[17] = 0.0;
 C[18] = 0.0;  C[19] = 0.0;  C[20] = 0.0;  C[21] = -0.25; C[22] = 0.0;  C[23] = 0.0;
 C[24] = 0.0;  C[25] = 0.0; C[26] = 0.0;  C[27] = 0.0;  C[28] = -0.25; C[29] = 0.0;
 C[30] = 0.0;  C[31] = 0.0; C[32] = 0.0;  C[33] = 0.0;  C[34] = 0.0;  C[35] = -0.25;
  
 S11 = (double*)malloc(6 * 6 * sizeof(double));
 Get_Submatrix(S11, 6, 6, S, 10, 0, 0);
 S12 = (double*)malloc(6 * 4 * sizeof(double));
 Get_Submatrix(S12, 6, 4, S, 10, 0, 6);
 S12t = (double*)malloc(4 * 6 * sizeof(double));
 Get_Submatrix(S12t, 4, 6, S, 10, 6, 0);
 S22 = (double*)malloc(4 * 4 * sizeof(double));
 Get_Submatrix(S22, 4, 4, S, 10, 6, 6);
 
 S22_1 = (double*)malloc(4 * 4 * sizeof(double));
 for(i = 0; i < 16; i++) {
   S22_1[i] = S22[i];
 }
 Choleski_LU_Decomposition(S22_1, 4);
 Choleski_LU_Inverse(S22_1, 4);
 
 // Calculate S22a = S22_1 * S12t   4*6 = 4x4 * 4x6   C = AB
 S22a = (double*)malloc(4 * 6 * sizeof(double));
 Multiply_Matrices(S22a, S22_1, 4, 4, S12t, 6);
 
 // Then calculate S22b = S12 * S22a      ( 6x6 = 6x4 * 4x6)
 S22b = (double*)malloc(6 * 6 * sizeof(double));
 Multiply_Matrices(S22b, S12, 6, 4, S22a, 6);
 
 // Calculate SS = S11 - S22b
 SS = (double*)malloc(6 * 6 * sizeof(double));
 for(i = 0; i < 36; i++) 
 {
   SS[i] = S11[i] - S22b[i];
 }
 E = (double*)malloc(6 * 6 * sizeof(double));
 Multiply_Matrices(E, C, 6, 6, SS, 6);

 SSS = (double*)malloc(6 * 6 * sizeof(double));
 Hessenberg_Form_Elementary(E, SSS, 6);

 eigen_real = (double*)malloc(6 * sizeof(double));
 eigen_imag = (double*)malloc(6 * sizeof(double));

 QR_Hessenberg_Matrix(E, SSS, eigen_real, eigen_imag, 6, 100);
 
 index = 0;
 maxval = eigen_real[0];
 for(i = 1; i < 6; i++)
 {
   if(eigen_real[i] > maxval)
   {
      maxval = eigen_real[i];
      index = i;
   }
 }

 v1 = (double*)malloc(6 * sizeof(double));
 
 v1[0] = SSS[index]; 
 v1[1] = SSS[index+6];
 v1[2] = SSS[index+12];
 v1[3] = SSS[index+18];
 v1[4] = SSS[index+24];
 v1[5] = SSS[index+30];

 // normalize v1
 norm = sqrt(v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2] + v1[3] * v1[3] + v1[4] * v1[4] + v1[5] * v1[5]);
 v1[0] /= norm;
 v1[1] /= norm;
 v1[2] /= norm;
 v1[3] /= norm;
 v1[4] /= norm;
 v1[5] /= norm;

 if(v1[0] < 0.0)
 {
   v1[0] = -v1[0];
   v1[1] = -v1[1];
   v1[2] = -v1[2];
   v1[3] = -v1[3];
   v1[4] = -v1[4];
   v1[5] = -v1[5];
 }

 // Calculate v2 = S22a * v1      ( 4x1 = 4x6 * 6x1)
 v2 = (double*)malloc(4 * sizeof(double));
 Multiply_Matrices(v2, S22a, 4, 6, v1, 1);

 v = (double*)malloc(10 * sizeof(double));

 v[0] = v1[0];
 v[1] = v1[1];
 v[2] = v1[2];
 v[3] = v1[3];
 v[4] = v1[4];
 v[5] = v1[5];
 v[6] = -v2[0];
 v[7] = -v2[1];
 v[8] = -v2[2];
 v[9] = -v2[3];

 Q = (double*)malloc(3 * 3 * sizeof(double));
 
 Q[0] = v[0];
 Q[1] = v[5];
 Q[2] = v[4];
 Q[3] = v[5];
 Q[4] = v[1];
 Q[5] = v[3];
 Q[6] = v[4];
 Q[7] = v[3];
 Q[8] = v[2];

 U = (double*)malloc(3 * sizeof(double));

 U[0] = v[6];
 U[1] = v[7];
 U[2] = v[8];
 Q_1 = (double*)malloc(3 * 3 * sizeof(double));
 for(i = 0; i < 9; i++)
 {
	Q_1[i] = Q[i];
 }
 Choleski_LU_Decomposition(Q_1, 3);
 Choleski_LU_Inverse(Q_1, 3);

 // Calculate B = Q-1 * U   ( 3x1 = 3x3 * 3x1)
 B = (double*)malloc(3 * sizeof(double));
 Multiply_Matrices(B, Q_1, 3, 3, U, 1);
 B[0] = -B[0];     // x-axis combined bias
 B[1] = -B[1];     // y-axis combined bias
 B[2] = -B[2];     // z-axis combined bias
 // The vector B gives the combined bias values to be used to 
 // correct the X, Y, Z axis values. Put them in the vector b
 // so they get returned to R. 
 b[0] = B[0];
 b[1] = B[1];
 b[2] = B[2];
 
// for(i = 0; i < 3; i++)
// {
	// Rprintf("%lf\r\n", B[i]);
// }

//--------------------------------------------------------------
// Now calculate the corrections for the combined scale factors,
// misalignments, and soft iron.
	
 // First calculate QB = Q * B   ( 3x1 = 3x3 * 3x1)
 QB = (double*)malloc(3 * sizeof(double));
 Multiply_Matrices(QB, Q, 3, 3, B, 1);

 // Then calculate btqb = BT * QB    ( 1x1 = 1x3 * 3x1)
 Multiply_Matrices(&btqb, B, 1, 3, QB, 1);

 // Calculate hmb = sqrt(btqb - J).
 J = v[9];
 hmb = sqrt(btqb - J);

 // Calculate SQ, the square root of matrix Q
 SSSS = (double*)malloc(3 * 3 * sizeof(double));
 Hessenberg_Form_Elementary(Q, SSSS, 3);

 eigen_real3 = (double*)malloc(3 * sizeof(double));
 eigen_imag3 = (double*)malloc(3 * sizeof(double));
 QR_Hessenberg_Matrix(Q, SSSS, eigen_real3, eigen_imag3, 3, 100);

 // normalize eigenvectors
 norm1 = sqrt(SSSS[0] * SSSS[0] + SSSS[3] * SSSS[3] + SSSS[6] * SSSS[6]);
 SSSS[0] /= norm1;
 SSSS[3] /= norm1;
 SSSS[6] /= norm1;
 norm2 = sqrt(SSSS[1] * SSSS[1] + SSSS[4] * SSSS[4] + SSSS[7] * SSSS[7]);
 SSSS[1] /= norm2;
 SSSS[4] /= norm2;
 SSSS[7] /= norm2;
 norm3 = sqrt(SSSS[2] * SSSS[2] + SSSS[5] * SSSS[5] + SSSS[8] * SSSS[8]);
 SSSS[2] /= norm3;
 SSSS[5] /= norm3;
 SSSS[8] /= norm3;
 
 Dz = (double*)malloc(3 * 3 * sizeof(double));
 for(i = 0; i < 9; i++)
   Dz[i] = 0.0;
 Dz[0] = sqrt(eigen_real3[0]);
 Dz[4] = sqrt(eigen_real3[1]);
 Dz[8] = sqrt(eigen_real3[2]);

 vdz = (double*)malloc(3 * 3 * sizeof(double));
 Multiply_Matrices(vdz, SSSS, 3, 3, Dz, 3);

 Transpose_Square_Matrix(SSSS, 3);

 SQ = (double*)malloc(3 * 3 * sizeof(double));
 Multiply_Matrices(SQ, vdz, 3, 3, SSSS, 3);

 hm = *localnorm;	// norm of the field (magnetic field total norm, or accel 1000 milli-g
 // Rprintf("Local norm: %d\n",hm);
 // A_1 will hold the correction for the combined scale factors,
 // misalignments, and soft iron. (A^{-1})
 A_1 = (double*)malloc(3 * 3 * sizeof(double));
 for(i = 0; i < 9; i++) 
 {
   A_1[i] = SQ[i] * hm / hmb;
 }
   
 // for(i = 0; i < 3; i++) 
 // {
   // Rprintf("%lf %lf %lf\r\n", A_1[i*3], A_1[i*3+1], A_1[i*3+2]);
 // } 
 
 // Copy the scale factors into their appropriate output vectors. 
 Xcorr[0] = A_1[0];
 Xcorr[1] = A_1[1];
 Xcorr[2] = A_1[2];
 Ycorr[0] = A_1[3];
 Ycorr[1] = A_1[4];
 Ycorr[2] = A_1[5];
 Zcorr[0] = A_1[6];
 Zcorr[1] = A_1[7];
 Zcorr[2] = A_1[8];
   
 
 free(D);
 free(S);
 free(C);
 free(S11);
 free(S12);
 free(S12t);
 free(S22);
 free(S22_1);
 free(S22a);
 free(S22b);
 free(SS);
 free(E);
 free(U);
 free(SSS);
 free(eigen_real);
 free(eigen_imag);
 free(v1);
 free(v2);
 free(v);
 free(Q);
 free(Q_1);
 free(B);
 free(QB);
 free(SSSS);
 free(eigen_real3);
 free(eigen_imag3);
 free(Dz);
 free(vdz);
 free(SQ);
 free(A_1);

}
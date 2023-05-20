/* 2020-02-03:
   e02baf & e02bbf from NAG math library.
   Extracted from A.Lyapustin's gauss.cpp
 
   Prototypes:
       void e02baf(int, double *, double *, double *k, double *c)
	   double e02bbf(int, double *, double *, double)

*/


#include <math.h>
#define max(a,b)    (((a) > (b)) ? (a) : (b))
#define min(a,b)    (((a) < (b)) ? (a) : (b))
//
//************************************************************************************
/*    MARK 11.5(F77) REVISED. (SEPT 1985.)
      ******************************************************
      E01BAF - currently absent.  AN ALGORITHM, WITH CHECKS, TO DETERMINE THE 
      COEFFICIENTS IN THE B-SPLINE REPRESENTATION OF A CUBIC SPLINE WHICH 
      INTERPOLATES (PASSES EXACTLY THROUGH) A GIVEN SET OF POINTS.
 
      INPUT PARAMETERS
         M        THE NUMBER OF DISTINCT POINTS WHICH THE SPLINE IS TO 
                     INTERPOLATE. (M MUST BE AT LEAST 4.)
         X        ARRAY CONTAINING THE DISTINCT VALUES OF THE
                     INDEPENDENT VARIABLE. NB X(I) MUST BE STRICTLY GREATER 
                     THAN X(J) WHENEVER I IS STRICTLY GREATER THAN J.
         Y        ARRAY CONTAINING THE VALUES OF THE DEPENDENT
                     VARIABLE.
         LCK = M + 4.
 
      OUTPUT PARAMETERS
         K        ON SUCCESSFUL EXIT, K CONTAINS THE KNOTS
                     SET UP BY THE ROUTINE. IF THE SPLINE IS
                     TO BE EVALUATED (BY NPL ROUTINE E02BEF,
                     FOR EXAMPLE) THE ARRAY K MUST NOT BE
                     ALTERED BEFORE CALLING THAT ROUTINE.
         C        ON SUCCESSFUL EXIT, C CONTAINS THE B-SPLINE
                     COEFFICIENTS OF THE INTERPOLATING SPLINE.
                     THESE ARE ALSO REQUIRED BY THE EVALUATING
                     ROUTINE E02BEF.
          
      WORKSPACE (AND ASSOCIATED DIMENSION) PARAMETERS
         WRK     WORKSPACE ARRAY, OF LENGTH LWRK /Memory is allocated in E02BAF/
         LWRK    ACTUAL DECLARED DIMENSION OF WRK.
                     MUST BE AT LEAST 6 * M + 16.
****************************************************************************************/ 
/*       NAG LIBRARY SUBROUTINE  E02BAF
 
      E02BAF  COMPUTES A WEIGHTED LEAST-SQUARES APPROXIMATION
      TO AN ARBITRARY SET OF DATA POINTS BY A CUBIC SPLINE
      WITH KNOTS PRESCRIBED BY THE USER.  CUBIC SPLINE
      INTERPOLATION CAN ALSO BE CARRIED OUT.
 
      COX-DE BOOR METHOD FOR EVALUATING B-SPLINES WITH
      ADAPTATION OF GENTLEMAN*S PLANE ROTATION SCHEME FOR
      SOLVING OVER-DETERMINED LINEAR SYSTEMS.
 
      REDESIGNED TO USE CLASSICAL GIVENS ROTATIONS IN
      ORDER TO AVOID THE OCCASIONAL UNDERFLOW (AND HENCE
      OVERFLOW) PROBLEMS EXPERIENCED BY GENTLEMAN*S 3-
      MULTIPLICATION PLANE ROTATION SCHEME
 
      WORK1  AND  WORK2  ARE WORKSPACE AREAS.
      WORK1(R)  CONTAINS THE VALUE OF THE  R TH  DISTINCT DATA
      ABSCISSA AND, SUBSEQUENTLY, FOR  R = 1, 2, 3, 4,  THE
      VALUES OF THE NON-ZERO B-SPLINES FOR EACH SUCCESSIVE
      ABSCISSA VALUE.
      WORK2(L, J)  CONTAINS, FOR  L = 1, 2, 3, 4,  THE VALUE OF
      THE  J TH  ELEMENT IN THE  L TH  DIAGONAL OF THE
      UPPER TRIANGULAR MATRIX OF BANDWIDTH  4  IN THE
      TRIANGULAR SYSTEM DEFINING THE B-SPLINE COEFFICIENTS.
*/ 


#include <stdio.h>
#include <stdlib.h>

void e02baf(int M, double *x, double *y, double *k, double *c)
{
//     .. Parameters ..
      int NCAP7 = M + 4; 
      double *w, *WORK1, **WORK2, SS;   
//    .. Local Scalars ..
      double  ACOL, AROW, CCOL, COSINE, CROW, D, D4, D5, D6, D7, D8, 
              D9, DPRIME, E2, E3, E4, E5, K1, K2, K3, K4, K5, K6, 
              N1, N2, N3, RELEMT, S, SIGMA, SINE, WI, XI;
      int     i, IPLUSJ, IU, j, JOLD, JPLUSL, JREV, l,
              L4, LPLUS1, LPLUSU, NCAP, NCAP3, NCAPM1;
              
      w = new double[M];
      WORK1 = new double[M];
      WORK2 = new double*[4];
      for(i=0; i<4; i++)
          WORK2[i] = new double[NCAP7];
//
//    INITIALISE THE ARRAY OF KNOTS AND THE ARRAY OF WEIGHTS
//
      for(i=0; i<M; i++) w[i] = 1.;
      if(M > 4) 
      {  for(i=4; i<M; i++) k[i] = x[i-2];  }
//
//     TESTS FOR ADEQUACY OF ARRAY LENGTHS AND THAT M IS GREATER
//     THAN 4.
//
      if(M < 4) 
        { printf("\n E02BAF Failure: M < 4");
	      exit(0);  }
//
//     TESTS FOR THE CORRECT ORDERING OF THE X(I)
//
      for(i=1; i<M; i++)
         if( x[i] <= x[i-1] ) 
        {  printf("\n Failure: x[i] <= x[i-1]");
           exit(0);  }

//     .. Executable Statements ..
      NCAP = NCAP7 - 7;
      NCAPM1 = NCAP - 1;
      NCAP3 = NCAP + 3;
//
//     IN ORDER TO DEFINE THE FULL B-SPLINE BASIS, AUGMENT THE
//     PRESCRIBED INTERIOR KNOTS BY KNOTS OF MULTIPLICITY FOUR
//     AT EACH END OF THE DATA RANGE.
//
      for(j=0; j<4; j++)
      {   i = NCAP3 + j;
          k[j] = x[0];
          k[i] = x[M-1];
       } 
//
//     CHECK THAT THE DATA ABSCISSAE ARE ORDERED, THEN FORM THE
//     ARRAY  WORK1  FROM THE ARRAY  X.  THE ARRAY  WORK1  CONTAINS
//     THE SET OF DISTINCT DATA ABSCISSAE.
//
      WORK1[0] = x[0];
      j = 1;
      for(i=1; i<M; i++)
      {  if( x[i] < WORK1[j-1]) 
          { printf("\n Error in E02BAF"); 
	        exit(0); }
         if( x[i] != WORK1[j-1] ) 
          {  WORK1[j] = x[i];
              j++;
           }
       }   
//
//     INITIALISE A BAND TRIANGULAR SYSTEM (I.E. A MATRIX AND A RIGHT HAND SIDE) TO 
//     ZERO. THE PROCESSING OF EACH DATA POINT IN TURN RESULTS IN AN UPDATING OF THIS 
//     SYSTEM. THE SUBSEQUENT SOLUTION OF THE RESULTING BAND TRIANGULAR SYSTEM
//     YIELDS THE COEFFICIENTS OF THE B-SPLINES.
//
      for(i=0; i<NCAP3; i++)
      {   for(l=0; l<4; l++) WORK2[l][i] = 0.0;
          c[i] = 0.0;
      }           
      SIGMA = 0.0;
      j = 0;
      JOLD = 0;
      for(i=1; i<=M; i++) 
      { 
//
//        FOR THE DATA POINT  (X(I), Y(I))  DETERMINE AN INTERVAL 
//        K(J + 3) .LE. X .LT. K(J + 4)  CONTAINING  X(I).  (IN THE CASE 
//        J + 4 .EQ. NCAP  THE SECOND EQUALITY IS RELAXED TO INCLUDE EQUALITY).
//
         WI = w[i-1];
         XI = x[i-1];
         while( (XI >= k[j+4-1]) && (j <= NCAPM1) ) j++;
         if(j != JOLD) 
         {
//        SET CERTAIN CONSTANTS RELATING TO THE INTERVAL
//        K(J + 3) .LE. X .LE. K(J + 4).
//
  	         K1 = k[j];
    	     K2 = k[j+1];
        	 K3 = k[j+2];
         	 K4 = k[j+3];
             K5 = k[j+4];
             K6 = k[j+5];
             D4 = 1.0/(K4-K1);
             D5 = 1.0/(K5-K2);
             D6 = 1.0/(K6-K3);
             D7 = 1.0/(K4-K2);
             D8 = 1.0/(K5-K3);
             D9 = 1.0/(K4-K3);
             JOLD = j;
           } 
//   COMPUTE AND STORE IN  WORK1(L) (L = 1, 2, 3, 4)  THE VALUES OF 
//   THE FOUR NORMALIZED CUBIC B-SPLINES WHICH ARE NON-ZERO AT X=X(I).
//
          E5 = K5 - XI;
          E4 = K4 - XI;
          E3 = XI - K3;
          E2 = XI - K2;
          N1 = WI*D9;
          N2 = E3*N1*D8;
          N1 = E4*N1*D7;
          N3 = E3*N2*D6;
          N2 = (E2*N1+E5*N2)*D5;
          N1 = E4*N1*D4;
          WORK1[3] = E3*N3;
          WORK1[2] = E2*N2 + (K6-XI)*N3;
          WORK1[1] = (XI-K1)*N1 + E5*N2;
          WORK1[0] = E4*N1;
          CROW = y[i-1]*WI;
//
//  ROTATE THIS ROW INTO THE BAND TRIANGULAR SYSTEM USING PLANE ROTATIONS.
//
         for(LPLUS1=1; LPLUS1<=4; LPLUS1++)
         {  l = LPLUS1 - 1;
            RELEMT = WORK1[LPLUS1-1];
            if(RELEMT == 0.0) goto m320;
            JPLUSL = j + l;
            L4 = 4 - l;
            D = WORK2[0][JPLUSL-1];
            if(fabs(RELEMT) >= D) 
              DPRIME = fabs(RELEMT) * sqrt( 1.0 + (D/RELEMT)*(D/RELEMT) );
            else DPRIME = D * sqrt( 1.0 + (RELEMT/D)*(RELEMT/D) );
            WORK2[0][JPLUSL-1] = DPRIME;
            COSINE = D/DPRIME;
            SINE = RELEMT/DPRIME;
            if( L4 >= 2) 
            {   for(IU=2; IU<=L4; IU++)
                {  LPLUSU = l + IU;;
                   ACOL = WORK2[IU-1][JPLUSL-1];
                   AROW = WORK1[LPLUSU-1];
                   WORK2[IU-1][JPLUSL-1] = COSINE*ACOL + SINE*AROW;
                   WORK1[LPLUSU-1] = COSINE*AROW - SINE*ACOL;
                 }
             }    
            CCOL = c[JPLUSL-1];
            c[JPLUSL-1] = COSINE*CCOL + SINE*CROW;
            CROW = COSINE*CROW - SINE*CCOL;
m320:  ; }
         SIGMA += CROW*CROW;
      }
      SS = SIGMA;    
//
//    SOLVE THE BAND TRIANGULAR SYSTEM FOR THE B-SPLINE COEFFICIENTS. 
//    IF A DIAGONAL ELEMENT IS ZERO, AND HENCE THE TRIANGULAR SYSTEM 
//    IS SINGULAR, THE IMPLICATION IS THAT THE SCHOENBERG-WHITNEY 
//    CONDITIONS ARE ONLY JUST SATISFIED. THUS IT IS APPROPRIATE TO EXIT 
//    IN THIS CASE WITH THE SAME VALUE  (IFAIL=5)  OF THE ERROR INDICATOR.
//
      l = -1;
      for(JREV=1; JREV<=NCAP3; JREV++)
      {  j = NCAP3 - JREV + 1;
         D = WORK2[0][j-1];
         if( D == 0.0 ) exit(0);
         if(l < 3) l++;
         S = c[j-1];
         if( l != 0 ) 
         {  for(i=1; i<=l; i++)
            {  IPLUSJ = i + j;
               S -= WORK2[i][j-1]*c[IPLUSJ-1];
            }
         } 
         c[j-1] = S/D;
      }                
      delete(w);  delete(WORK1);  for(i=0; i<4; i++) delete(WORK2[i]);   delete(WORK2);
      return;
}


/*     NAG LIBRARY SUBROUTINE  E02BBF
       E02BBF  EVALUATES A CUBIC SPLINE FROM ITS B-SPLINE REPRESENTATION.
       DE BOOR*S METHOD OF CONVEX COMBINATIONS.
       MARK 11.5(F77) REVISED. (SEPT 1985.)
*/

double e02bbf(int NCAP7, double *K, double *C, double X)
{
      double  C1, C2, C3, E2, E3, E4, E5, S;
      int     J, J1, L;

      if(NCAP7 < 8)
      { printf("\nFailure in E02BBF: NCAP7<8");
        return 0.;
       } 
      if(X < K[4-1] && X > K[NCAP7-3-1]) 
      { printf("\nFailure in E02BBF: x-value is outside the bounds of spline knots");
        return 0.;
       } 
//
//     DETERMINE  J  SUCH THAT  K(J + 3) .LE. X .LE. K(J + 4).
//
      J1 = 0;
      J = NCAP7 - 7;
      L = (J1+J)/2;
      while( (J-J1) > 1) 
      {
         if(X >= K[L+4-1]) J1 = L;
         else            J = L; 
         L = (J1+J)/2;
       }          
//
//     USE THE METHOD OF CONVEX COMBINATIONS TO COMPUTE  S(X).
//
      E2 = X - K[J+1];
      E3 = X - K[J+2];
      E4 = K[J+3] - X;
      E5 = K[J+4] - X;
      C2 = C[J];
      C3 = C[J+1];
      C1 = ((X - K[J])*C2 + E4*C[J-1])/(K[J+3] - K[J]);
      C2 = (E2*C3+E5*C2)/(K[J+4] - K[J+1]);
      C3 = (E3*C[J+2]+(K[J+5] - X)*C3)/(K[J+5] - K[J+2]);
      C1 = (E2*C2+E4*C1)/(K[J+3] - K[J+1]);
      C2 = (E3*C3+E5*C2)/(K[J+4] - K[J+2]);
      S = (E3*C2+E4*C1)/(K[J+3] - K[J+2]);
  
      return S;
}
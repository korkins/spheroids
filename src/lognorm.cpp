#include <math.h> /* log */
//
double lognorm(double const s, double const r0, double const r)
/*--------------------------------------------------------------------------------------------------
TASK:
	To simulate lognormal distribution.
IN:
	s    d   standard deviation, s > 0
	r0   d   mean radius, r0 > 0
	r    d   radius, r > 0
OUT:
	lognorm   d   distribution
NOTE:
    See [1] for definitions.
REFS:
	1. https://en.wikipedia.org/wiki/Log-normal_distribution
--------------------------------------------------------------------------------------------------*/
{
	double const
	    sqrt_2pi = 2.506628274631000502415765284811;
	double
		lnr;
//--------------------------------------------------------------------------------------------------
//
	lnr = log( r/r0 );
	return exp( -0.5*( lnr*lnr/s/s ) )/sqrt_2pi/s;
} // double lognorm(...)
/*--------------------------------------------------------------------------------------------------
23/05/01: cosmetic changes.
23/04/26: first created based on SLOG (DLS) and SLog (SHARM) and tested in
          bimodal mode.
    in:
	    rf = 0.12;
	    sgmf = 0.5;
	    rc = 1.9;
	    sgmc = 0.6;
	    cvf = 0.03;
	    cvc = 0.45;
    out:
	sized[ir] = cvf * lognorm(sgmf, rf, r[ir]) + cvc * lognorm(sgmc, rc, r[ir]);
        r[ 0] =  5.0000000e-02  sd(r) =  5.1681542e-03
        r[ 1] =  6.5603700e-02  sd(r) =  1.1543435e-02
        r[ 2] =  8.6076800e-02  sd(r) =  1.9194996e-02
        r[ 3] =  1.1293900e-01  sd(r) =  2.3765796e-02
        r[ 4] =  1.4818400e-01  sd(r) =  2.1933496e-02
        r[ 5] =  1.9442900e-01  sd(r) =  1.5243585e-02
        r[ 6] =  2.5510500e-01  sd(r) =  8.7806014e-03
        r[ 7] =  3.3471600e-01  sd(r) =  7.4623256e-03
        r[ 8] =  4.3917200e-01  sd(r) =  1.6027865e-02
        r[ 9] =  5.7622600e-01  sd(r) =  4.1605449e-02
        r[10] =  7.5605200e-01  sd(r) =  9.2023323e-02
        r[11] =  9.9199500e-01  sd(r) =  1.6642502e-01
        r[12] =  1.3015700e+00  sd(r) =  2.4527718e-01
        r[13] =  1.7077600e+00  sd(r) =  2.9451528e-01
        r[14] =  2.2407000e+00  sd(r) =  2.8811284e-01
        r[15] =  2.9399600e+00  sd(r) =  2.2962683e-01
        r[16] =  3.8574500e+00  sd(r) =  1.4910262e-01
        r[17] =  5.0612600e+00  sd(r) =  7.8877400e-02
        r[18] =  6.6407400e+00  sd(r) =  3.3995831e-02
        r[19] =  8.7131400e+00  sd(r) =  1.1937184e-02
        r[20] =  1.1432300e+01  sd(r) =  3.4149062e-03
        r[21] =  1.5000000e+01  sd(r) =  7.9591393e-04
//
// Original codes
//
      REAL FUNCTION SLOG(C,S,RM,RR)
C******************************************************************
C**     Lognormal function d(RR)/dlnR                            **
C******************************************************************     
C  C   R  - concentration
C  S   R  - halfwidth
C  RM  R  - mean radius
C  RR  R  - value
C******************************************************************
      PI = ACOS( -1.0 )
      IF(S.LE.1) THEN
      WRITE(*,*) S,'TROUBLES in SLOG:S.LE.1'
	STOP 'STOP: in SLOG:S.LE.1'
	ENDIF
	A1=LOG(RR/RM)
	A2=LOG(S)
      SLOG=C/SQRT(2.0*PI)/A2*EXP(-0.5*((A1*A1)/(A2*A2)))
      RETURN
      END
//--------------------------------------------------------------------------------------------------
//
double SLog(double c, double s, double r, double rr)
{
   double slg, a1, a2;

//		if(s<1)	printf("trouble in SLOG: s<1.\n");
	
	a1 = log(rr/r);
	a2 = log(s);
	slg = c/sqrt(2*_PI)/a2*exp(-0.5*((a1*a1)/(a2*a2)) );
	return slg;
}
--------------------------------------------------------------------------------------------------*/
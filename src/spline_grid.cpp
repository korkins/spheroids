# include <math.h>	/* log, exp */
// NAG:
void e02baf(int M, double *x, double *y, double *k, double *c);
double e02bbf(int NCAP7, double *K, double *C, double X);
//
int spline_grid(int const nx, double const *x, double const *fx, int const nu, double const *u,
                 double *fu)
/*--------------------------------------------------------------------------------------------------
TASK:
	To interpolate function f() from grid x onto grid u: f(x) -> f(u)
IN:
    nx   i       len(x)
	x    d[nx]   original grid, x[0] < x[1] < ... < x[nx-1]
	fx   d[nx]   f(x)
	nu   i       len(u)
	u    d[nx]   user grid within x: x[0] < u[0] < u[1] < ... < u[nu-1] < x[nx-1]
OUT:
    code   i       error code (see table)
	fu     d[nu]   interpolated values
NOTES:
    code values:
		-1: u[0] <= x[0]
		 0: all ok
		+1: u[nu-1] > x[nx-1]
REFS:
    1. NAG e02baf manual
	2. NAG e02bbf manual
--------------------------------------------------------------------------------------------------*/
{
	double const
		fill_value = -999.0;
	int
		code, iu, ix, lck, nxx;
	double
		*krab, *crab, *xx, *fxx;
//--------------------------------------------------------------------------------------------------
//
	code = 0;
	if (u[0] < x[0])
	{
		code = -1;
		for (iu = 0; iu < nu; iu++)
		{
			fu[iu] = fill_value;
		}
	}
	else if (u[nu-1] > x[nx-1])
	{
		code = 1;
		for (iu = 0; iu < nu; iu++)
		{
			fu[iu] = fill_value;
		}
	}
	else
	{
//		NAG: see [1] for 4
		lck = nx + 4;
		krab = new double [lck];
		crab = new double [lck];
//
//		NAG functions do not accept 'const': copy arrays
		nxx = nx;
		xx = new double [nxx];
		fxx = new double [nxx];
		for (ix = 0; ix < nxx; ix++)
		{
			xx[ix] = x[ix];
			fxx[ix] = fx[ix];
		}
//
		e02baf(nxx, xx, fxx, krab, crab);
		for (iu = 0; iu < nu; iu++)
			fu[iu] = e02bbf(lck, krab, crab, u[iu]);
//
		delete[] krab;
		delete[] crab;
		delete[] xx;
		delete[] fxx;
	}
	return code;
} // void spline_grid(...)
/*--------------------------------------------------------------------------------------------------
23/05/01: cosmetic changes.
23/04/28: changed:
		   u[0] <= x[0] -> u[0] < x[0]  &  u[nu-1] >= x[nx-1] -> u[nu-1] > x[nx-1]
23/03/13: first created and tested
          a) x = [0:15:180], f(x) = cos(rad(x)), u = [20:20:160]
			                         spline_grid    python
			 u[iu] =  20.0, fu(u) =   0.93964700    0.9396470027561222
             u[iu] =  40.0, fu(u) =   0.76604302    0.7660430248781633
             u[iu] =  60.0, fu(u) =   0.50000000    0.5000000000000001
             u[iu] =  80.0, fu(u) =   0.17364709    0.1736470861051101
             u[iu] = 100.0, fu(u) =  -0.17364709   -0.1736470861051101
             u[iu] = 120.0, fu(u) =  -0.50000000   -0.4999999999999998
             u[iu] = 140.0, fu(u) =  -0.76604302   -0.7660430248781633
             u[iu] = 160.0, fu(u) =  -0.93964700   -0.9396470027561221
		  
		  b) same as (a) excepot for "random" u:
		     u[iu] =   1.0, fu(u) =   0.99989516    0.9998951565582959
             u[iu] =   5.0, fu(u) =   0.99632186    0.9963218576673675
             u[iu] =  13.0, fu(u) =   0.97440172    0.9744017217963545
             u[iu] =  55.0, fu(u) =   0.57356832    0.5735683162887505
             u[iu] =  77.0, fu(u) =   0.22495109    0.2249510941246620
             u[iu] =  90.0, fu(u) =   0.00000000    0.0
             u[iu] = 130.0, fu(u) =  -0.64277931   -0.6427793098229079
             u[iu] = 179.0, fu(u) =  -0.99989516   -0.9998951565582959

		 c) same as (b) except for u[0] = -1.0 < x[0]:
            code = -1
            u[iu] =  -1.0, fu(u) = -999.00000000
            u[iu] =   5.0, fu(u) = -999.00000000
            u[iu] =  13.0, fu(u) = -999.00000000
            u[iu] =  55.0, fu(u) = -999.00000000
            u[iu] =  77.0, fu(u) = -999.00000000
            u[iu] =  90.0, fu(u) = -999.00000000
            u[iu] = 130.0, fu(u) = -999.00000000
            u[iu] = 179.0, fu(u) = -999.00000000

		c) same as (b) except for u[nu-1] = 181.0 > x[nx-1]:
		    code = 1
            u[iu] =   1.0, fu(u) = -999.00000000
            u[iu] =   5.0, fu(u) = -999.00000000
            u[iu] =  13.0, fu(u) = -999.00000000
            u[iu] =  55.0, fu(u) = -999.00000000
            u[iu] =  77.0, fu(u) = -999.00000000
            u[iu] =  90.0, fu(u) = -999.00000000
            u[iu] = 130.0, fu(u) = -999.00000000
            u[iu] = 181.0, fu(u) = -999.00000000
--------------------------------------------------------------------------------------------------*/
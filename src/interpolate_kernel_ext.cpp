#include "const_param_spheroids.h"
//
void read_fixkernel_ext_bin(char const *fname, double *ufext, double *ufabs, int const nx);
double bilinear(double const x,   double const y,
                double const x1,  double const x2,
                double const y1,  double const y2,
                double const f11, double const f21, double const f12, double const f22);
//
int interpolate_kernel_ext(char const *fname,
	                       int const irefre1, int const irefre2, double const refre1, double const refre2, double const refre,
	                       int const irefim1, int const irefim2, double const refim1, double const refim2, double const refim,
	                       double *extx, double *absx)
/*--------------------------------------------------------------------------------------------------
TASK:
	To interpolate the extinction and absorption kernels over real and imaginary parts of refractive
	index using bilinear interpoaltion.
IN:
    fname     a[]   binary file name with kernel values
	irefre1   i
	irefre2   i
	refre1    d     refre[irefre1] = refre1 <= refre <= refre2 <= refre[irefre2]
	refre2    d
	refre     d
	irefim1   i
	irefim2   i
	refim1    d     refim[irefim1] = refim1 <= refim <= refim2 <= refim[irefim2]
	refim2    d
	refim     d
OUT:
	extx      d[nrgrid_fix]   extinction on rgrid; NOT scaled by wavel_fix/wavel
	absx      d[nrgrid_fix]   absorption on rgrid; NOT scaled by wavel_fix/wavel
	code      i               error code, see TODO list
NOTE:
	TODO list: implement code values:
	    0 - all good (default value);
		1 -
		2 -
		3 -
		4 -
	See [1] for details and "const_param_spheroids.h" for paramters.
REFS:
	1. Dubovik et al, 2006: JGR 111, D11208; https://doi.org/10.1029/2005JD006619
--------------------------------------------------------------------------------------------------*/
{
	int
		code, nx, ir;
	double
		ext11,  ext12,  ext21,  ext22,
		abs11,  abs12,  abs21,  abs22;
	double
		*kext, *kabs;
//--------------------------------------------------------------------------------------------------
	code = -1;
//
	nx = nrgrid_fix * nrefim * nrefre;
	kext = new double [nx];
	kabs = new double [nx];
	read_fixkernel_ext_bin(fname, kext, kabs, nx);
//
	for (ir = 0; ir < nrgrid_fix; ir++)
	{
		ext11 = kext[irefre1 * nrefim * nrgrid_fix + irefim1 * nrgrid_fix + ir];
		ext21 = kext[irefre2 * nrefim * nrgrid_fix + irefim1 * nrgrid_fix + ir];
		ext12 = kext[irefre1 * nrefim * nrgrid_fix + irefim2 * nrgrid_fix + ir];
		ext22 = kext[irefre2 * nrefim * nrgrid_fix + irefim2 * nrgrid_fix + ir];
//		f(x, y): x -> refre, y -> refim
		extx[ir] = bilinear(refre, refim, refre1, refre2, refim1, refim2, ext11, ext21, ext12, ext22);
//
		abs11 = kabs[irefre1 * nrefim * nrgrid_fix + irefim1 * nrgrid_fix + ir];
		abs21 = kabs[irefre2 * nrefim * nrgrid_fix + irefim1 * nrgrid_fix + ir];
		abs12 = kabs[irefre1 * nrefim * nrgrid_fix + irefim2 * nrgrid_fix + ir];
		abs22 = kabs[irefre2 * nrefim * nrgrid_fix + irefim2 * nrgrid_fix + ir];
		absx[ir] = bilinear(refre, refim, refre1, refre2, refim1, refim2, abs11, abs21, abs12, abs22);
	} // for ir
//
	delete[] kext;
	delete[] kabs;
//
	code = 0;
	return code;
} // int interpolate_kernel_ext(...)
/*--------------------------------------------------------------------------------------------------
23/05/12: first created & tested with optichar(...) for 2 different scenarios.
--------------------------------------------------------------------------------------------------*/
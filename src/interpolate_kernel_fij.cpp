#include "const_param_spheroids.h"
//
void read_fixkernel_fij_bin(char const *fname, double *ufij, int const nxij);
double bilinear(double const x,   double const y,
                double const x1,  double const x2,
                double const y1,  double const y2,
                double const f11, double const f21, double const f12, double const f22);
//
int interpolate_kernel_fij(char const *fname,
	                       int const irefre1, int const irefre2, double const refre1, double const refre2, double const refre,
	                       int const irefim1, int const irefim2, double const refim1, double const refim2, double const refim,
	                       double *fijx)
/*--------------------------------------------------------------------------------------------------
TASK:
	To interpolate the phase function kernel from fixed values to that defined by a user.
	Cubic spline and linear interpolations are used for radii and refractive indices, respectively.
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
    fijx      d[nsca_fix * nrgrid_fix]   the phase function kernel for the given refractive index,
	                                        on the fixed kernels r-grid at wavel_fix = 0.340 (um).
	code      i   error code - see table below
NOTE:
	TODO list: implement code values
	           0 - all good (default value);
		       1 -
		       2 -
		       3 -
		       4 -
	See [1] for ATBD and "const_param_spheroids.h" for paramters.
REFS:
	1. Dubovik et al, 2006: JGR 111, D11208; https://doi.org/10.1029/2005JD006619
--------------------------------------------------------------------------------------------------*/
{
	int
		code, nx, ir, isca;
	double
		fk11, fk12, fk21, fk22;
	double
		*fkern;
//--------------------------------------------------------------------------------------------------
//
	code = -1;
//
	nx = nsca_fix * nrgrid_fix * nrefim * nrefre;
	fkern = new double [nx];
	read_fixkernel_fij_bin(fname, fkern, nx);
//
	for (isca = 0; isca < nsca_fix; isca++)
		for (ir = 0; ir < nrgrid_fix; ir++)
		{
			fk11 = fkern[isca * nrefre * nrefim * nrgrid_fix + irefre1 * nrefim * nrgrid_fix + irefim1 * nrgrid_fix + ir];
			fk21 = fkern[isca * nrefre * nrefim * nrgrid_fix + irefre2 * nrefim * nrgrid_fix + irefim1 * nrgrid_fix + ir];
			fk12 = fkern[isca * nrefre * nrefim * nrgrid_fix + irefre1 * nrefim * nrgrid_fix + irefim2 * nrgrid_fix + ir];
			fk22 = fkern[isca * nrefre * nrefim * nrgrid_fix + irefre2 * nrefim * nrgrid_fix + irefim2 * nrgrid_fix + ir];
			fijx[isca * nrgrid_fix + ir] = bilinear(refre, refim, refre1, refre2, refim1, refim2, fk11, fk21, fk12, fk22);
		} // for ir
//
	delete[] fkern;
//
	code = 0;
	return code;
} // int interpolate_kernel_fij(...)
/*--------------------------------------------------------------------------------------------------
23/05/12: first created and tested with optichar(...) in 2 scenarios.
--------------------------------------------------------------------------------------------------*/
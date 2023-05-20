//#include <stdio.h> /* printf */
#include <math.h>  /* cos */
#include "const_param_spheroids.h"
//
int getix(double const x0, double const *x, int const nx, int &ix1, int &ix2);
int interpolate_kernel_ext(char const *fname,
	                       int const irefre1, int const irefre2, double const refre1, double const refre2, double const refre,
	                       int const irefim1, int const irefim2, double const refim1, double const refim2, double const refim,
	                       double *extx, double *absx);
int interpolate_kernel_fij(char const *fname,
	                       int const irefre1, int const irefre2, double const refre1, double const refre2, double const refre,
	                       int const irefim1, int const irefim2, double const refim1, double const refim2, double const refim,
	                       double *fijx);
double lognorm(double const s, double const r0, double const r);
double ddot(double const *x1, double const *x2, int const nx);
// optional steps:
int spline_grid(int const nx, double const *x, double const *fx, int const nu, double const *u, double *fu);
double simpson(double const *f, int const nx, double const dx);
//
int optichar(double const wavel, double const refre, double const refim, double const mie_fraq,
	        double const cvf, double const radf, double const sgmf,
	        double const cvc, double const radc, double const sgmc,
	        double *f11, double *f22, double *f33, double *f44, double *f12, double *f34,
			double &cext, double &ssa, double &conc_ratio)
/*--------------------------------------------------------------------------------------------------
TASK:
	To calculate optical characteristics - normalized phase matrix, extinction, and single
	scattering albedo - for given paramters of bimodal distribution (volumetric),
	refractive index at a given wavelength and spherical (Mie) and spheroidal
	particles mxture.
	For checking purpose, the bimodal distributon is integrated on the tabilated
	and fine grids to get the total concentration (must be close to cvf + cvc).
IN:
	wavel      d   wavelength, microns
	refre      d   refractive index: real part, ...
	refim      d   ... imaginary part; same for Mie & spheroids, for fine & coarse fractions
	mie_fraq   d   fraqtion of spherical vs spheroidal particles; 1.0 - only Mie
	cvf        d   fine mode volumetric concentration, a.u.
	radf       d   fine mode mean radius
	sgmf       d   fine mode lognormal distribution width parameter
	cvc        d   coarse mode volumetric concentration, a.u
	radc       d   coarse mode mean radius
	sgmc       d   coarse mode lognormal distribution width parameter
OUT:
	fij        d[nsca_fix]   normalized phase matrix elements (see NOTEs below)
	cext       d             extinction in a.u. (because of cvf & cvc)
	ssa        d             single scattering albedo
NOTE:
	Normalized phase matrix elements: f22/f11, ..., f34/f11.
	Normalized phase function:

		norm = 1/2 integrate{f11(x) sin(x) dx, x = 0:pi}
		f11 = f11/norm
	
	so that norm ~= 1.0
	Average scattering cosine for normalized phase function:

		g = 1/2 integrate{f11(x) sin(x) cos(x) dx, x = 0:pi} < 1.0
REFS:
	-
--------------------------------------------------------------------------------------------------*/
{	
	int
		ir, ix, nx, isca, code, irefre1, irefre2, irefim1, irefim2;
	double
		cabs, csca, conc, dlnr, dtheta, f, norm,
		refre1, refre2, refim1, refim2,
		spxf, spxc, srd_fraq, x0, x1;
    double
		*sizeparx_fix,
		*extx_mie, *absx_mie, *extx_srd, *absx_srd, *extx, *absx, *sized,
		*f11x_mie, *f22x_mie, *f33x_mie, *f44x_mie, *f12x_mie, *f34x_mie,
		*f11x_srd, *f22x_srd, *f33x_srd, *f44x_srd, *f12x_srd, *f34x_srd,
		*f11x, *f22x, *f33x, *f44x, *f12x, *f34x, *fsin, *fcos, *theta, *theta_fine, *lnf, *lnf_fine;
//--------------------------------------------------------------------------------------------------
//
	code = -1;
//
	spxf = radf/wavel;
	spxc = radc/wavel;
//
	sizeparx_fix = new double [nrgrid_fix];
	extx_mie = new double [nrgrid_fix];
	absx_mie = new double [nrgrid_fix];
	extx_srd = new double [nrgrid_fix];
	absx_srd = new double [nrgrid_fix];
	extx = new double [nrgrid_fix];
	absx = new double [nrgrid_fix];
	sized = new double [nrgrid_fix];
	for (ir = 0; ir < nrgrid_fix; ir++)
	{
		extx_mie[ir] = 0.0;
		absx_mie[ir] = 0.0;
		extx_srd[ir] = 0.0;
		absx_srd[ir] = 0.0;
		extx[ir] = 0.0;
		absx[ir] = 0.0;
		sized[ir] = 0.0;
	} // for ir
//
   	nx = nsca_fix * nrgrid_fix;
	f11x = new double [nx];
	f22x = new double [nx];
	f33x = new double [nx];
	f44x = new double [nx];
	f12x = new double [nx];
	f34x = new double [nx];
	f11x_mie = new double [nx];
	f22x_mie = new double [nx];
	f33x_mie = new double [nx];
	f44x_mie = new double [nx];
	f12x_mie = new double [nx];
	f34x_mie = new double [nx];
	f11x_srd = new double [nx];
	f22x_srd = new double [nx];
	f33x_srd = new double [nx];
	f44x_srd = new double [nx];
	f12x_srd = new double [nx];
	f34x_srd = new double [nx];
//
	for (ix = 0; ix < nx; ix++)
	{
		f11x[ix] = 0.0;
		f22x[ix] = 0.0;
		f33x[ix] = 0.0;
		f44x[ix] = 0.0;
		f12x[ix] = 0.0;
		f34x[ix] = 0.0;
		f11x_mie[ix] = 0.0;
		f22x_mie[ix] = 0.0;
		f33x_mie[ix] = 0.0;
		f44x_mie[ix] = 0.0;
		f12x_mie[ix] = 0.0;
		f34x_mie[ix] = 0.0;
		f11x_srd[ix] = 0.0;
		f22x_srd[ix] = 0.0;
		f33x_srd[ix] = 0.0;
		f44x_srd[ix] = 0.0;
		f12x_srd[ix] = 0.0;
		f34x_srd[ix] = 0.0;
	} // for ix
//
	code = getix(refre, refre_fix, nrefre, irefre1, irefre2);
	refre1 = refre_fix[irefre1];
	refre2 = refre_fix[irefre2];
//
	code = getix(refim, refim_fix, nrefim, irefim1, irefim2);
	refim1 = refim_fix[irefim1];
	refim2 = refim_fix[irefim2];
/*
	for (ir = 0; ir < nrgrid; ir++)
		sizeparx[ir] = rgrid[ir] / wavel;
*/
	for (ir = 0; ir < nrgrid_fix; ir++)
		sizeparx_fix[ir] = rgrid_fix[ir] / wavel_fix;
//
	srd_fraq = 1.0 - mie_fraq;
	if (srd_fraq < tiny)
	{
		code = interpolate_kernel_ext(fname_fixkern_ext_mie,
			                          irefre1, irefre2, refre1, refre2, refre,
		                              irefim1, irefim2, refim1, refim2, refim,
								      extx, absx);
		code = interpolate_kernel_fij(fname_fixkern_f11_mie,
									  irefre1, irefre2, refre1, refre2, refre,
		                              irefim1, irefim2, refim1, refim2, refim,
								      f11x);
		code = interpolate_kernel_fij(fname_fixkern_f22_mie,
									  irefre1, irefre2, refre1, refre2, refre,
		                              irefim1, irefim2, refim1, refim2, refim,
								      f22x);
		code = interpolate_kernel_fij(fname_fixkern_f33_mie,
									  irefre1, irefre2, refre1, refre2, refre,
		                              irefim1, irefim2, refim1, refim2, refim,
								      f33x);
		code = interpolate_kernel_fij(fname_fixkern_f44_mie,
									  irefre1, irefre2, refre1, refre2, refre,
		                              irefim1, irefim2, refim1, refim2, refim,
								      f44x);
		code = interpolate_kernel_fij(fname_fixkern_f12_mie,
									  irefre1, irefre2, refre1, refre2, refre,
		                              irefim1, irefim2, refim1, refim2, refim,
								      f12x);
		code = interpolate_kernel_fij(fname_fixkern_f34_mie,
									  irefre1, irefre2, refre1, refre2, refre,
		                              irefim1, irefim2, refim1, refim2, refim,
								      f34x);
	} // if srd_fraq < tiny
	else if (mie_fraq < tiny)
	{
		code = interpolate_kernel_ext(fname_fixkern_ext_srd,
			                          irefre1, irefre2, refre1, refre2, refre,
		                              irefim1, irefim2, refim1, refim2, refim,
								      extx, absx);
		code = interpolate_kernel_fij(fname_fixkern_f11_srd,
									  irefre1, irefre2, refre1, refre2, refre,
		                              irefim1, irefim2, refim1, refim2, refim,
								      f11x);
		code = interpolate_kernel_fij(fname_fixkern_f22_srd,
									  irefre1, irefre2, refre1, refre2, refre,
		                              irefim1, irefim2, refim1, refim2, refim,
								      f22x);
		code = interpolate_kernel_fij(fname_fixkern_f33_srd,
									  irefre1, irefre2, refre1, refre2, refre,
		                              irefim1, irefim2, refim1, refim2, refim,
								      f33x);
		code = interpolate_kernel_fij(fname_fixkern_f44_srd,
									  irefre1, irefre2, refre1, refre2, refre,
		                              irefim1, irefim2, refim1, refim2, refim,
								      f44x);
		code = interpolate_kernel_fij(fname_fixkern_f12_srd,
									  irefre1, irefre2, refre1, refre2, refre,
		                              irefim1, irefim2, refim1, refim2, refim,
								      f12x);
		code = interpolate_kernel_fij(fname_fixkern_f34_srd,
									  irefre1, irefre2, refre1, refre2, refre,
		                              irefim1, irefim2, refim1, refim2, refim,
								      f34x);
	} // if mie_fraq < tiny
	else // mix
	{
		code = interpolate_kernel_ext(fname_fixkern_ext_mie,
			                          irefre1, irefre2, refre1, refre2, refre,
		                              irefim1, irefim2, refim1, refim2, refim,
								      extx_mie, absx_mie);
		code = interpolate_kernel_fij(fname_fixkern_f11_mie,
									  irefre1, irefre2, refre1, refre2, refre,
		                              irefim1, irefim2, refim1, refim2, refim,
								      f11x_mie);
		code = interpolate_kernel_fij(fname_fixkern_f22_mie,
									  irefre1, irefre2, refre1, refre2, refre,
		                              irefim1, irefim2, refim1, refim2, refim,
								      f22x_mie);
		code = interpolate_kernel_fij(fname_fixkern_f33_mie,
									  irefre1, irefre2, refre1, refre2, refre,
		                              irefim1, irefim2, refim1, refim2, refim,
								      f33x_mie);
		code = interpolate_kernel_fij(fname_fixkern_f44_mie,
									  irefre1, irefre2, refre1, refre2, refre,
		                              irefim1, irefim2, refim1, refim2, refim,
								      f44x_mie);
		code = interpolate_kernel_fij(fname_fixkern_f12_mie,
									  irefre1, irefre2, refre1, refre2, refre,
		                              irefim1, irefim2, refim1, refim2, refim,
								      f12x_mie);
		code = interpolate_kernel_fij(fname_fixkern_f34_mie,
									  irefre1, irefre2, refre1, refre2, refre,
		                              irefim1, irefim2, refim1, refim2, refim,
								      f34x_mie);
//
		code = interpolate_kernel_ext(fname_fixkern_ext_srd,
			                          irefre1, irefre2, refre1, refre2, refre,
		                              irefim1, irefim2, refim1, refim2, refim,
								      extx_srd, absx_srd);
		code = interpolate_kernel_fij(fname_fixkern_f11_srd,
									  irefre1, irefre2, refre1, refre2, refre,
		                              irefim1, irefim2, refim1, refim2, refim,
								      f11x_srd);
		code = interpolate_kernel_fij(fname_fixkern_f22_srd,
									  irefre1, irefre2, refre1, refre2, refre,
		                              irefim1, irefim2, refim1, refim2, refim,
								      f22x_srd);
		code = interpolate_kernel_fij(fname_fixkern_f33_srd,
									  irefre1, irefre2, refre1, refre2, refre,
		                              irefim1, irefim2, refim1, refim2, refim,
								      f33x_srd);
		code = interpolate_kernel_fij(fname_fixkern_f44_srd,
									  irefre1, irefre2, refre1, refre2, refre,
		                              irefim1, irefim2, refim1, refim2, refim,
								      f44x_srd);
		code = interpolate_kernel_fij(fname_fixkern_f12_srd,
									  irefre1, irefre2, refre1, refre2, refre,
		                              irefim1, irefim2, refim1, refim2, refim,
								      f12x_srd);
		code = interpolate_kernel_fij(fname_fixkern_f34_srd,
									  irefre1, irefre2, refre1, refre2, refre,
		                              irefim1, irefim2, refim1, refim2, refim,
								      f34x_srd);
//
		for (ir = 0; ir < nrgrid_fix; ir++)
		{
			extx[ir] = mie_fraq * extx_mie[ir] + srd_fraq * extx_srd[ir];
			absx[ir] = mie_fraq * absx_mie[ir] + srd_fraq * absx_srd[ir];
		} // for ir
//
		for (isca = 0; isca < nsca_fix; isca++)
			for (ir = 0; ir < nrgrid_fix; ir++)
			{
				f11x[isca * nrgrid_fix + ir] = mie_fraq * f11x_mie[isca * nrgrid_fix + ir] + srd_fraq * f11x_srd[isca * nrgrid_fix + ir];
				f22x[isca * nrgrid_fix + ir] = mie_fraq * f22x_mie[isca * nrgrid_fix + ir] + srd_fraq * f22x_srd[isca * nrgrid_fix + ir];
				f33x[isca * nrgrid_fix + ir] = mie_fraq * f33x_mie[isca * nrgrid_fix + ir] + srd_fraq * f33x_srd[isca * nrgrid_fix + ir];
				f44x[isca * nrgrid_fix + ir] = mie_fraq * f44x_mie[isca * nrgrid_fix + ir] + srd_fraq * f44x_srd[isca * nrgrid_fix + ir];
				f12x[isca * nrgrid_fix + ir] = mie_fraq * f12x_mie[isca * nrgrid_fix + ir] + srd_fraq * f12x_srd[isca * nrgrid_fix + ir];
				f34x[isca * nrgrid_fix + ir] = mie_fraq * f34x_mie[isca * nrgrid_fix + ir] + srd_fraq * f34x_srd[isca * nrgrid_fix + ir];
			} // for ir
	} // mie_frac > tiny & srd_frac > tiny
//
//  Note: integration over size parameter x = r/wavel; dlnr = dlnx
	conc = 0.0;
	for (ir = 0; ir < nrgrid_fix; ir++)
	{
		sized[ir] = cvf * lognorm(sgmf, spxf, sizeparx_fix[ir]) + cvc * lognorm(sgmc, spxc, sizeparx_fix[ir]);
		conc += sized[ir];
	} // for ir
	dlnr =  log( sizeparx_fix[nrgrid_fix-1]/sizeparx_fix[0] ) / (nrgrid_fix-1);
	conc_ratio = conc * dlnr/(cvf + cvc);
/*
//  Note: integration over size r - result must be close to integration over size paramter (ideally - the same) 
    conc = 0.0;
	dlnr = log( rmax_grid_fix/rmin_grid_fix ) / (nrgrid_fix-1);
	for (ir = 0; ir < nrgrid_fix; ir++)
	{
		r = exp( log(rmin_grid_fix) + ir * dlnr );
		conc += cvf * lognorm(sgmf, radf, r) + cvc * lognorm(sgmc, radc, r);
	} // for ir
	conc_ratio = conc * dlnr / (cvf + cvc);
*/
	cext = ddot(extx, sized, nrgrid_fix);
	cabs = ddot(absx, sized, nrgrid_fix);
	csca = cext - cabs;
	ssa = 1.0 - cabs/cext;
//
	cext = cext * scalef_um3_um2 * wavel_fix / wavel;
//
	for (isca = 0; isca < nsca_fix; isca++)
	{
		f = ddot( &f11x[isca*nrgrid_fix], sized, nrgrid_fix ); 
		f11[isca] = f / csca;
		f22[isca] =  ddot( &f22x[isca*nrgrid_fix], sized, nrgrid_fix ) / f;
		f33[isca] =  ddot( &f33x[isca*nrgrid_fix], sized, nrgrid_fix ) / f;
		f44[isca] =  ddot( &f44x[isca*nrgrid_fix], sized, nrgrid_fix ) / f;
		f12[isca] = -ddot( &f12x[isca*nrgrid_fix], sized, nrgrid_fix ) / f; // note '-'
		f34[isca] =  ddot( &f34x[isca*nrgrid_fix], sized, nrgrid_fix ) / f;
	} // for isca
/*
//  Smooth f33 & f44 - as in the DLS code
//	f33:
	isca1 = 38;
	isca2 = 50;
	sca1 = sca_fix[isca1];
	f1 = log(f33[isca1]);
	sca2 = sca_fix[isca2];
	f2 = log(f33[isca2]);
	for (isca = isca1; isca < isca2+1; isca++)
		f33[isca] = exp( linear(sca_fix[isca], sca1, sca2, f1, f2) );
//	f44:
	isca1 = 48;
	isca2 = 60;
	sca1 = sca_fix[isca1];
	f1 = log(f44[isca1]);
	sca2 = sca_fix[isca2];
	f2 = log(f44[isca2]);
	for (isca = isca1; isca < isca2+1; isca++)
		f44[isca] = exp( linear(sca_fix[isca], sca1, sca2, f1, f2) );
*/
//  Normalization of f11
	theta_fine = new double [nsca_fine];
	lnf_fine = new double [nsca_fine];
	fcos = new double [nsca_fine];
	fsin = new double [nsca_fine];
//
	dtheta = pi/(nsca_fine-1);
	for (isca = 0; isca < nsca_fine; isca++)
		theta_fine[isca] = isca*dtheta;
	theta = new double [nsca_fix];
	lnf = new double [nsca_fix];
	for (isca = 0; isca < nsca_fix; isca++)
	{
		theta[isca] = sca_fix[isca]*rad;
		lnf[isca] = log(f11[isca]);
	} // for isca
//
	code = spline_grid(nsca_fix, theta, lnf, nsca_fine, theta_fine, lnf_fine);
	for (isca = 0; isca < nsca_fine; isca++)
		fsin[isca] = exp(lnf_fine[isca])*sin(theta_fine[isca]);
//
	norm = simpson(fsin, nsca_fine, dtheta) / 2.0;
	for (isca = 0; isca < nsca_fix; isca++)
	{
		f11[isca] /= norm;
		lnf[isca] = log(f11[isca]);
	} // for isca
//
//	Moments x0 ~= 1.0 and x1 = g (average scattering cosine)
	code = spline_grid(nsca_fix, theta, lnf, nsca_fine, theta_fine, lnf_fine);
	for (isca = 0; isca < nsca_fine; isca++)
	{
		fsin[isca] = exp(lnf_fine[isca])*sin(theta_fine[isca]);
		fcos[isca] = exp(lnf_fine[isca])*sin(theta_fine[isca])*cos(theta_fine[isca]);
	} // for isca
	x0 = simpson(fsin, nsca_fine, dtheta) / 2.0;
	x1 = simpson(fcos, nsca_fine, dtheta) / 2.0;
//
	delete[] sizeparx_fix;
	delete[] extx_mie;
	delete[] absx_mie;
	delete[] extx_srd;
	delete[] absx_srd;
	delete[] extx;
	delete[] absx;
	delete[] sized;
	delete[] f11x;
	delete[] f22x;
	delete[] f33x;
	delete[] f44x;
	delete[] f12x;
	delete[] f34x;
	delete[] f11x_mie;
	delete[] f22x_mie;
	delete[] f33x_mie;
	delete[] f44x_mie;
	delete[] f12x_mie;
	delete[] f34x_mie;
	delete[] f11x_srd;
	delete[] f22x_srd;
	delete[] f33x_srd;
	delete[] f44x_srd;
	delete[] f12x_srd;
	delete[] f34x_srd;
	delete[] theta_fine;
	delete[] lnf;
	delete[] lnf_fine;
	delete[] fcos;
	delete[] fsin;
	delete[] theta;
//
	code = 0;
	return code;
}
/*--------------------------------------------------------------------------------------------------
23/05/13: multiple changes and improvments; tested for 2 cases.
23/05/08: new interpoaltion of fij-kernels: log(fij/f11)+2 for all elements; tested for 2 cases.
23/05/01: minor changes in output variables and comments.
23/04/27: first created; tested vs DLS for 2 different cases at wav = 0.4 and 0.6 (um):
		  Good agreement with origianl DLS - see ./doc/ folder for details.
--------------------------------------------------------------------------------------------------*/
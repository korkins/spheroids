#include <stdio.h>  /* printf */
#include <time.h>   /* clock  */
#include "const_param_spheroids.h"
//
int optichar(double const wavel, double const refre, double const refim, double const mie_fraq,
	        double const cvf, double const radf, double const sgmf,
	        double const cvc, double const radc, double const sgmc,
	        double *f11, double *f22, double *f33, double *f44, double *f12, double *f34,
			double &cext, double &ssa, double &conc_ratio);
//
int main()
{
	int const
		nruns = 50;
	int
		irun, code, isca;
	double
		wavel, refre, refim, time_start, time_end, mie_fraq,
		cvf, cvc, radf, radc, sgmf, sgmc, conc_ratio, cext, ssa,
		cext_avr, ssa_avr;
    double
		*f11, *f22, *f33, *f44, *f12, *f34;
//
	f11 = new double [nsca_fix];
	f22 = new double [nsca_fix];
	f33 = new double [nsca_fix];
	f44 = new double [nsca_fix];
	f12 = new double [nsca_fix];
	f34 = new double [nsca_fix];
	for (isca = 0; isca < nsca_fix; isca++)
	{
		f11[isca] = 0.0;
		f22[isca] = 0.0;
		f33[isca] = 0.0;
		f44[isca] = 0.0;
		f12[isca] = 0.0;
		f34[isca] = 0.0;
	} // for isca
//
//  **************
//	*** Case 1 ***
//  **************
//
	mie_fraq = 0.7;
	wavel = 0.600;
	refre = 1.6;
	refim = 0.001;
//
	radf = 0.3;
	sgmf = 0.4;
	radc = 3.0;
	sgmc = 0.9;
	cvf = 0.5;
	cvc = 1.0;
//
	printf("\nCase 1 inp: wavel = %.3f", wavel);
	printf("\n            mie_fraq = %.3f  refre = %.2f  refim = %.2f", mie_fraq, refre, refim);
	printf("\n            cvf = %.2f  radf = %.2f  sgmf = %.2f", cvf, radf, sgmf);
	printf("\n            cvc = %.2f  radc = %.2f  sgmc = %.2f", cvc, radc, sgmc);
//
	code = optichar(wavel, refre, refim, mie_fraq, cvf, radf, sgmf, cvc, radc, sgmc,
					f11, f22, f33, f44, f12, f34, cext, ssa, conc_ratio);
//
	printf("\nCase 1 out:");
	for (isca = 0; isca < nsca_fix; isca++)
		printf("\n%5.1f  %14.6e  %14.6e  %14.6e  %14.6e  %14.6e  %14.6e",
			       sca_fix[isca], f11[isca], f22[isca], f33[isca], f44[isca], f12[isca], f34[isca]);
	printf("\nconcentration ratio: %9.6f", conc_ratio);
	printf("\next = %8.5f  ssa = %8.5f\n", cext, ssa);
//
//  **************
//	*** Case 2 ***
//  **************
//
	mie_fraq = 0.3;
	wavel = 0.400;
	refre = 1.3;
	refim = 0.1;
//
	radf = 0.3;
	sgmf = 0.4;
	radc = 3.0;
	sgmc = 0.9;
	cvf = 1.5;
	cvc = 1.0;
//
	printf("\nCase 2 inp: wavel = %.3f", wavel);
	printf("\n            mie_fraq = %.3f  refre = %.2f  refim = %.2f", mie_fraq, refre, refim);
	printf("\n            cvf = %.2f  radf = %.2f  sgmf = %.2f", cvf, radf, sgmf);
	printf("\n            cvc = %.2f  radc = %.2f  sgmc = %.2f", cvc, radc, sgmc);
//
	code = optichar(wavel, refre, refim, mie_fraq, cvf, radf, sgmf, cvc, radc, sgmc,
					f11, f22, f33, f44, f12, f34, cext, ssa, conc_ratio);
//
	printf("\nCase 2 out:");
	for (isca = 0; isca < nsca_fix; isca++)
		printf("\n%5.1f  %14.6e  %14.6e  %14.6e  %14.6e  %14.6e  %14.6e",
			       sca_fix[isca], f11[isca], f22[isca], f33[isca], f44[isca], f12[isca], f34[isca]);
	printf("\nconcentration ratio: %9.6f", conc_ratio);
	printf("\next = %8.5f  ssa = %8.5f\n", cext, ssa);
//
//  **************
//	*** Timing ***
//  **************
//
	if (nruns > 1)
	{
		printf("\nCase 1: timing ... ");
		mie_fraq = 0.7;
		wavel = 0.600;
		refre = 1.6;
		refim = 0.001;
//
		radf = 0.3;
		sgmf = 0.4;
		radc = 3.0;
		sgmc = 0.9;
		cvf = 0.5;
		cvc = 1.0;
//
		cext_avr = 0.0;
		ssa_avr = 0.0;
		time_start = (double)clock() /(double)CLOCKS_PER_SEC;
		for (irun = 0; irun < nruns; irun++)
		{
			code = optichar(wavel, refre, refim, mie_fraq, cvf, radf, sgmf, cvc, radc, sgmc,
							f11, f22, f33, f44, f12, f34, cext, ssa, conc_ratio);
			cext_avr += cext;
			ssa_avr += ssa;
		}
		time_end = (double)clock() /(double)CLOCKS_PER_SEC;
//
		printf("done:");
		for (isca = 0; isca < nsca_fix; isca++)
			printf("\n%5.1f  %14.6e  %14.6e  %14.6e  %14.6e  %14.6e  %14.6e",
			       sca_fix[isca], f11[isca], f22[isca], f33[isca], f44[isca], f12[isca], f34[isca]);
		printf("\nconcentration ratio: %9.6f", conc_ratio);
		printf("\next = %8.5f  ssa = %8.5f", cext_avr/nruns, ssa_avr/nruns);
//
		printf("\ntiming: nruns = %i ...", nruns);
		printf("\nruntime per run: %8.3f (sec.)", (time_end - time_start) / nruns);
	} // if nruns > 1
//
	printf("\n");
	return 0;
}
/*------------------------------------------------------------------------------------------------*/
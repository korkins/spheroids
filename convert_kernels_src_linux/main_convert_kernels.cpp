#include <stdio.h> /* printf */
#include <math.h>  /* cos */
#include <time.h>   /* clock */
#include <iostream>
#include "const_param_spheroids.h"
//
void read_Rkext_fix(char const* fname, float* ufext, float* ufabs);
void read_Rkernel_fix(char const* fname, float* ufij);
void write_Rkext_fix_bin(char const *fname, double *ufext, double *ufabs, int const nx);
void write_Rkernel_fix_bin(char const *fname, double *ufij, int const nxij);
void read_fixkernel_ext_bin(char const *fname, double *ufext, double *ufabs, int const nx);
void read_fixkernel_fij_bin(char const *fname, double *ufij, int const nxij);
//
int main()
{
	FILE
		*pfile;
	int const
		str_len_max = 256;
	int
		ix, nx, nxij, iy, ix0, iy0, isangl0, irgrid0, irefim0, irefre0, isangl, irgrid, irefim, irefre;
	double
		time_start, time_end;
	float
		*ufext_srd, *ufabs_srd,
		*uf11_srd, *uf22_srd, *uf33_srd, *uf44_srd, *uf12_srd, *uf34_srd,
		*ufext_mie, *ufabs_mie,
		*uf11_mie, *uf22_mie, *uf33_mie, *uf44_mie, *uf12_mie, *uf34_mie;
	double
		*vfext_srd, *vfabs_srd,
		*vf11_srd, *vf22_srd, *vf33_srd, *vf44_srd, *vf12_srd, *vf34_srd,
		*vfext_mie, *vfabs_mie,
		*vf11_mie, *vf22_mie, *vf33_mie, *vf44_mie, *vf12_mie, *vf34_mie;
	double
		*tfext_srd, *tfabs_srd,
		*tf11_srd, *tf22_srd, *tf33_srd, *tf44_srd, *tf12_srd, *tf34_srd,
		*tfext_mie, *tfabs_mie,
		*tf11_mie, *tf22_mie, *tf33_mie, *tf44_mie, *tf12_mie, *tf34_mie;
//-------------------------------------------------------------------------
	irefre0 = 9;  // 0:21 = 22 total
	irefim0 = 10; // 0:14 = 15 total
	irgrid0 = 11; // 0:33 = 34 total
	isangl0 = 12; // 0:180 = 181 total
//
	nx = nrefre * nrefim * nrgrid_fix;
	nxij = nsca_fix * nrefre * nrefim * nrgrid_fix;
//
//  1. Read ASCII files
	ufext_mie = new float [nx];
	ufabs_mie = new float [nx];
	uf11_mie = new float [nxij];
	uf22_mie = new float [nxij];
	uf33_mie = new float [nxij];
	uf44_mie = new float [nxij];
	uf12_mie = new float [nxij];
	uf34_mie = new float [nxij];
	ufext_srd = new float [nx];
	ufabs_srd = new float [nx];
	uf11_srd = new float [nxij];
	uf22_srd = new float [nxij];
	uf33_srd = new float [nxij];
	uf44_srd = new float [nxij];
	uf12_srd = new float [nxij];
	uf34_srd = new float [nxij];
//
	time_start = (double)clock() /(double)CLOCKS_PER_SEC;
	read_Rkext_fix("Rkext1_s.fix", ufext_mie, ufabs_mie);
	read_Rkernel_fix("Rkernel1.11_s.fix", uf11_mie);
	read_Rkernel_fix("Rkernel1.22_s.fix", uf22_mie);
	read_Rkernel_fix("Rkernel1.33_s.fix", uf33_mie);
	read_Rkernel_fix("Rkernel1.44_s.fix", uf44_mie);
	read_Rkernel_fix("Rkernel1.12_s.fix", uf12_mie);
	read_Rkernel_fix("Rkernel1.34_s.fix", uf34_mie);
//
	read_Rkext_fix("Rkext1.fix", ufext_srd, ufabs_srd);
	read_Rkernel_fix("Rkernel1.11.fix", uf11_srd);
	read_Rkernel_fix("Rkernel1.22.fix", uf22_srd);
	read_Rkernel_fix("Rkernel1.33.fix", uf33_srd);
	read_Rkernel_fix("Rkernel1.44.fix", uf44_srd);
	read_Rkernel_fix("Rkernel1.12.fix", uf12_srd);
	read_Rkernel_fix("Rkernel1.34.fix", uf34_srd);
	time_end = (double)clock() /(double)CLOCKS_PER_SEC;
	printf("\nread ASCII - runtime: %8.3fs\n", time_end - time_start);
//
	printf("   ufext[0] = %12.6e, ufext[nx-1] = %12.6e\n", ufext_srd[0], ufext_srd[nx-1]);
	printf("   ufabs[0] = %12.6e, ufabs[nx-1] = %12.6e\n", ufabs_srd[0], ufabs_srd[nx-1]);
	printf("    uf11[0] = %12.6e,  uf11[nx-1] = %12.6e\n", uf11_srd[0], uf11_srd[nx-1]);
	printf("    uf22[0] = %12.6e,  uf22[nx-1] = %12.6e\n", uf22_srd[0], uf22_srd[nx-1]);
	printf("    uf33[0] = %12.6e,  uf33[nx-1] = %12.6e\n", uf33_srd[0], uf33_srd[nx-1]);
	printf("    uf44[0] = %12.6e,  uf44[nx-1] = %12.6e\n", uf44_srd[0], uf44_srd[nx-1]);
	printf("    uf12[0] = %12.6e,  uf12[nx-1] = %12.6e\n", uf12_srd[0], uf12_srd[nx-1]);
	printf("    uf34[0] = %12.6e,  uf34[nx-1] = %12.6e\n", uf34_srd[0], uf34_srd[nx-1]);
//
//	2. Convert float to binary. For scattering matrix, vfij, make rgrid lead diemansion (originally, in ufij, scattering angle is the lead)
	vfext_mie = new double [nx];
	vfabs_mie = new double [nx];
	vf11_mie = new double [nxij];
	vf22_mie = new double [nxij];
	vf33_mie = new double [nxij];
	vf44_mie = new double [nxij];
	vf12_mie = new double [nxij];
	vf34_mie = new double [nxij];
	vfext_srd = new double [nx];
	vfabs_srd = new double [nx];
	vf11_srd = new double [nxij];
	vf22_srd = new double [nxij];
	vf33_srd = new double [nxij];
	vf44_srd = new double [nxij];
	vf12_srd = new double [nxij];
	vf34_srd = new double [nxij];
//
	for (ix = 0; ix < nx; ix++)
	{
		vfext_mie[ix] = (double) ufext_mie[ix];
		vfabs_mie[ix] = (double) ufabs_mie[ix];
		vfext_srd[ix] = (double) ufext_srd[ix];
		vfabs_srd[ix] = (double) ufabs_srd[ix];
	}
//
	for (irefre = 0; irefre < nrefre; irefre++)
		for (irefim = 0; irefim < nrefim; irefim++)
			for (irgrid = 0; irgrid < nrgrid_fix; irgrid++)
				for (isangl = 0; isangl < nsca_fix; isangl++)
				{
					ix = isangl * nrefre * nrefim * nrgrid_fix + irefre * nrefim * nrgrid_fix + irefim * nrgrid_fix + irgrid;
					iy = irefre * nrefim * nrgrid_fix * nsca_fix + irefim * nrgrid_fix * nsca_fix + irgrid * nsca_fix + isangl;
//
					vf11_mie[ix] = (double) uf11_mie[iy];
					vf22_mie[ix] = (double) uf22_mie[iy];
					vf33_mie[ix] = (double) uf33_mie[iy];
					vf44_mie[ix] = (double) uf44_mie[iy];
					vf12_mie[ix] = (double) uf12_mie[iy];
					vf34_mie[ix] = (double) uf34_mie[iy];
//
					vf11_srd[ix] = (double) uf11_srd[iy];
					vf22_srd[ix] = (double) uf22_srd[iy];
					vf33_srd[ix] = (double) uf33_srd[iy];
					vf44_srd[ix] = (double) uf44_srd[iy];
					vf12_srd[ix] = (double) uf12_srd[iy];
					vf34_srd[ix] = (double) uf34_srd[iy];
				}
//
//	3. Save to binary
	write_Rkext_fix_bin("Rkext1_s.fix.bin", vfext_mie, vfabs_mie, nx);
	write_Rkext_fix_bin("Rkext1.fix.bin", vfext_srd, vfabs_srd, nx);
	write_Rkernel_fix_bin("Rkernel1.11_s.fix.bin", vf11_mie, nxij);
	write_Rkernel_fix_bin("Rkernel1.11.fix.bin", vf11_srd, nxij);
	write_Rkernel_fix_bin("Rkernel1.22_s.fix.bin", vf22_mie, nxij);
	write_Rkernel_fix_bin("Rkernel1.22.fix.bin", vf22_srd, nxij);
	write_Rkernel_fix_bin("Rkernel1.33_s.fix.bin", vf33_mie, nxij);
	write_Rkernel_fix_bin("Rkernel1.33.fix.bin", vf33_srd, nxij);
	write_Rkernel_fix_bin("Rkernel1.44_s.fix.bin", vf44_mie, nxij);
	write_Rkernel_fix_bin("Rkernel1.44.fix.bin", vf44_srd, nxij);
	write_Rkernel_fix_bin("Rkernel1.12_s.fix.bin", vf12_mie, nxij);
	write_Rkernel_fix_bin("Rkernel1.12.fix.bin", vf12_srd, nxij);
	write_Rkernel_fix_bin("Rkernel1.34_s.fix.bin", vf34_mie, nxij);
	write_Rkernel_fix_bin("Rkernel1.34.fix.bin", vf34_srd, nxij);
//
//  4. Read from binary
	tfext_mie = new double [nx];
	tfabs_mie = new double [nx];
	tf11_mie = new double [nxij];
	tf22_mie = new double [nxij];
	tf33_mie = new double [nxij];
	tf44_mie = new double [nxij];
	tf12_mie = new double [nxij];
	tf34_mie = new double [nxij];
	tfext_srd = new double [nx];
	tfabs_srd = new double [nx];
	tf11_srd = new double [nxij];
	tf22_srd = new double [nxij];
	tf33_srd = new double [nxij];
	tf44_srd = new double [nxij];
	tf12_srd = new double [nxij];
	tf34_srd = new double [nxij];
//
	time_start = (double)clock() /(double)CLOCKS_PER_SEC;
	read_fixkernel_ext_bin("Rkext1_s.fix.bin", tfext_mie, tfabs_mie, nx);
	read_fixkernel_ext_bin("Rkext1.fix.bin", tfext_srd, tfabs_srd, nx);
	read_fixkernel_fij_bin("Rkernel1.11_s.fix.bin", tf11_mie, nxij);
	read_fixkernel_fij_bin("Rkernel1.11.fix.bin", tf11_srd, nxij);
	read_fixkernel_fij_bin("Rkernel1.22_s.fix.bin", tf22_mie, nxij);
	read_fixkernel_fij_bin("Rkernel1.22.fix.bin", tf22_srd, nxij);
	read_fixkernel_fij_bin("Rkernel1.33_s.fix.bin", tf33_mie, nxij);
	read_fixkernel_fij_bin("Rkernel1.33.fix.bin", tf33_srd, nxij);
	read_fixkernel_fij_bin("Rkernel1.44_s.fix.bin", tf44_mie, nxij);
	read_fixkernel_fij_bin("Rkernel1.44.fix.bin", tf44_srd, nxij);
	read_fixkernel_fij_bin("Rkernel1.12_s.fix.bin", tf12_mie, nxij);
	read_fixkernel_fij_bin("Rkernel1.12.fix.bin", tf12_srd, nxij);
	read_fixkernel_fij_bin("Rkernel1.34_s.fix.bin", tf34_mie, nxij);
	read_fixkernel_fij_bin("Rkernel1.34.fix.bin", tf34_srd, nxij);
	time_end = (double)clock() /(double)CLOCKS_PER_SEC;
	printf("\nread BIN - runtime: %8.3fs\n", time_end - time_start);
//
	printf("   tfext[0] = %12.6e, tfext[nx-1] = %12.6e\n", tfext_srd[0], tfext_srd[nx-1]);
	printf("   tfabs[0] = %12.6e, tfabs[nx-1] = %12.6e\n", tfabs_srd[0], tfabs_srd[nx-1]);
	printf("    tf11[0] = %12.6e,  tf11[nx-1] = %12.6e\n", tf11_srd[0], tf11_srd[nx-1]);
	printf("    tf22[0] = %12.6e,  tf22[nx-1] = %12.6e\n", tf22_srd[0], tf22_srd[nx-1]);
	printf("    tf33[0] = %12.6e,  tf33[nx-1] = %12.6e\n", tf33_srd[0], tf33_srd[nx-1]);
	printf("    tf44[0] = %12.6e,  tf44[nx-1] = %12.6e\n", tf44_srd[0], tf44_srd[nx-1]);
	printf("    tf12[0] = %12.6e,  tf12[nx-1] = %12.6e\n", tf12_srd[0], tf12_srd[nx-1]);
	printf("    tf34[0] = %12.6e,  tf34[nx-1] = %12.6e\n", tf34_srd[0], tf34_srd[nx-1]);
//
//  7. Test original vs binary
	printf("\nTEST: nodes and values::");
	printf("\nrefre[%i] = %13.7e", irefre0, refre_fix[irefre0]);
	printf("\nrefim[%i] = %13.7e", irefim0, refim_fix[irefim0]);
	printf("\nrgrid_fix[%i] = %13.7e", irgrid0, rgrid_fix[irgrid0]);
	printf("\nsangl_fix[%i] = %5.1f", isangl0, sca_fix[isangl0]);
//
	ix0 = irefre0 * nrefim * nrgrid_fix + irefim0*nrgrid_fix + irgrid0;
	printf("\n ix0 = %i", ix0);
	printf("\n ufext_mie[ix0] = %12.6e, tfext_mie[ix0] = %12.6e", ufext_mie[ix0], tfext_mie[ix0]);
	printf("\n ufabs_mie[ix0] = %12.6e, tfabs_mie[ix0] = %12.6e", ufabs_mie[ix0], tfabs_mie[ix0]);
	printf("\n ufext_srd[ix0] = %12.6e, tfext_srd[ix0] = %12.6e", ufext_srd[ix0], tfext_srd[ix0]);
	printf("\n ufabs_srd[ix0] = %12.6e, tfabs_srd[ix0] = %12.6e", ufabs_srd[ix0], tfabs_srd[ix0]);
//
	iy0 = irefre0 * nrefim * nrgrid_fix * nsca_fix + irefim0 * nrgrid_fix * nsca_fix + irgrid0 * nsca_fix + isangl0;
	ix0 = isangl0 * nrefre * nrefim * nrgrid_fix + irefre0 * nrefim * nrgrid_fix + irefim0 * nrgrid_fix + irgrid0;
	printf("\n iy0 = %i,  ix0 = %i", iy0, ix0);
	printf("\n uf11_mie[iy0] = %14.6e, tf11_mie[ix0] = %14.6e", uf11_mie[iy0], tf11_mie[ix0]);
	printf("\n uf11_srd[iy0] = %14.6e, tf11_srd[ix0] = %14.6e", uf11_srd[iy0], tf11_srd[ix0]);
	printf("\n uf22_mie[iy0] = %14.6e, tf22_mie[ix0] = %14.6e", uf22_mie[iy0], tf22_mie[ix0]);
	printf("\n uf22_srd[iy0] = %14.6e, tf22_srd[ix0] = %14.6e", uf22_srd[iy0], tf22_srd[ix0]);
	printf("\n uf33_mie[iy0] = %14.6e, tf33_mie[ix0] = %14.6e", uf33_mie[iy0], tf33_mie[ix0]);
	printf("\n uf33_srd[iy0] = %14.6e, tf33_srd[ix0] = %14.6e", uf33_srd[iy0], tf33_srd[ix0]);
	printf("\n uf44_mie[iy0] = %14.6e, tf44_mie[ix0] = %14.6e", uf44_mie[iy0], tf44_mie[ix0]);
	printf("\n uf44_srd[iy0] = %14.6e, tf44_srd[ix0] = %14.6e", uf44_srd[iy0], tf44_srd[ix0]);
	printf("\n uf12_mie[iy0] = %14.6e, tf12_mie[ix0] = %14.6e", uf12_mie[iy0], tf12_mie[ix0]);
	printf("\n uf12_srd[iy0] = %14.6e, tf12_srd[ix0] = %14.6e", uf12_srd[iy0], tf12_srd[ix0]);
	printf("\n uf34_mie[iy0] = %14.6e, tf34_mie[ix0] = %14.6e", uf34_mie[iy0], tf34_mie[ix0]);
	printf("\n uf34_srd[iy0] = %14.6e, tf34_srd[ix0] = %14.6e", uf34_srd[iy0], tf34_srd[ix0]);
//
	printf("\ndone!\n");
	return 0;
//
}
#include <stdio.h>
//#include <iostream>
#include <string.h>     /* strcpy, strcat */
#include "const_param_spheroids.h" // rootdir
//
void write_Rkernel_fix_bin(char const *fname, double *ufij, int const nxij) {
/*------------------------------------------------------------------------------
TASK:
	To read data from 'Rkext.fix' file.
IN:
	fname   'a'     file name
	ufij    d[nx]   extinction
	nxij    i       len(ufij)
OUT:
	-
NOTE:
    FYI:
		nx = nrgrid_fix*nrefim*nrefre;
		see: "const_param_spheroids.h" for 'nrefre', 'nrefim', 'nrgrid_fix';
		lead dimension: nrgrid_fix.
REFS:
	-
------------------------------------------------------------------------------*/
//
	FILE
		*pfile;
	int const
		path_len_max = 256;
	char
		fpath[path_len_max];
//-----------------------------------------------------------------------------
//
	strcpy(fpath, rootdir);
	strcat(fpath, fname);
	pfile = fopen(fpath, "wb");
	fwrite(ufij, nxij*sizeof(double), 1, pfile);
	fclose(pfile);
	printf("\n(in write_Rkernel_fix_bin) bin saved: %s", fpath);
//
} // void write_Rkernel_fix_bin
/*------------------------------------------------------------------------------
23/03/02: first created and tested.
------------------------------------------------------------------------------*/
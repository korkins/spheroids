#include <stdio.h>
//#include <iostream>
#include <cstring>     /* strcpy, strcat */
#include "const_param_spheroids.h" // rootdir
//
void write_Rkext_fix_bin(char const *fname, double *ufext, double *ufabs, int const nx) {
/*------------------------------------------------------------------------------
TASK:
	To read data from 'Rkext.fix' file.
IN:
	fname   'a'     file name
	ufext   f[nx]   extinction
	ufabs   f[nx]   absorption
	nx      i       length of arrays
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
	fwrite(ufext, nx*sizeof(double), 1, pfile);
	fwrite(ufabs, nx*sizeof(double), 1, pfile);
	fclose(pfile);
	printf("\n(in write_Rkext_fix_bin) bin saved: %s", fpath);
//
} // void write_Rkext_fix_bin
/*------------------------------------------------------------------------------
23/03/02: first created and tested.
------------------------------------------------------------------------------*/
#include <stdio.h>
//#include <iostream>
#include <cstring>     /* strcpy, strcat */
#include "const_param_spheroids.h"
//
void read_Rkernel_fix(char const *fname, float *ufij) {
/*------------------------------------------------------------------------------
TASK:
	To read data from 'Rkext.fix' file.
IN:
	fname   'a'   file name
OUT:
	ufij   d[nsca*nrgrid_fix*nrefim*nrefre]   ij-element of the phase matrix
NOTE:
    Dimensions 'nsca', 'nrefre', 'nrefim', 'nrgrid_fix' are defined in
	"const_param_spheroids.h"; nsca is the lead dimension:

	    for irefre:
	        for irefim:
		        for irgrid:
				    for isca:
			        read [irefre, irefim, irgrid, isca]

	Thus function is slow because fscanf is called too many times
REFS:
	-
------------------------------------------------------------------------------*/
//
	FILE
		*pfile;
	int const
		str_len_max = 256;
	int
		irefre, irefim, ix, nx, ixoffset;
	char
		chr;
	char
		str_buffer[str_len_max], fpath[path_len_max];
//-----------------------------------------------------------------------------
//
	strcpy(fpath, rootdir);
	strcat(fpath, fname);
	pfile = fopen(fpath, "r");
//	printf("\n(in read_Rkernel_fix) opened: %s", fpath);
//
	for (ix = 0; ix < 60; ix++) // skip header: 60 lines (update if '25  a number of grid aspect ratios' changes)
		fgets(str_buffer, str_len_max, pfile);
//
	nx = nrgrid_fix*nsca_fix;
	for (irefre = 0; irefre < nrefre; irefre++)
	{
		for (irefim = 0; irefim < nrefim; irefim++)
		{
			fgets(str_buffer, str_len_max, pfile);
			fgets(str_buffer, str_len_max, pfile);
			ixoffset = (irefre*nrefim+irefim)*nx;
			for (ix = 0; ix < nx; ix++)
			{
				fscanf(pfile, "%f", &ufij[ixoffset+ix]); // this is slow
			} // for ix
			fscanf(pfile, "%c", &chr); // read '\n'
		}// for irefim
	} // for irefre
//
	fclose(pfile);
//	printf("\n(in read_Rkernel_fix) closed: %s", fpath);
//
} // void read_grid_fix
/*------------------------------------------------------------------------------
23/02/18: fscanf ignores whitespaces - read all 'ngrid_fix*nsca_fix' elements at once;
          tested -- ok; however, there is no gain in time vs. previous versions!
22/04/04: float arrays instead of double, fgets() skips lines; tested -- ok;
22/04/01: first created and tested in 'matrix_fix()' -- ok.
*/
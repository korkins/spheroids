#include <stdio.h>
//#include <iostream>
#include <cstring>     /* strcpy, strcat */
#include "const_param_spheroids.h"
//
void read_Rkext_fix(char const *fname, float *ufext, float *ufabs) {
/*------------------------------------------------------------------------------
TASK:
	To read data from 'Rkext.fix' file.
IN:
	fname   'a'   file name
OUT:
	ufext   f[nrgrid_fix*nrefim*nrefre]   extinction
	ufabs   f[nrgrid_fix*nrefim*nrefre]   absorption
NOTE:
    Dimensions 'nrefre', 'nrefim', 'nrgrid_fix' are defined in
	"const_param_spheroids.h"; nrgrid_fix is the lead dimension:

	    for irefre:
	        for irefim:
		        for irgrid:
			        read [irefre, irefim, irgrid]

	Lead dimension: nrgrid_fix

	THINKME: which way to skip lines is better
	a) fgets() or
	b) for (ix = 0; ix < 33; ix++) // skip header: 33 lines
	   {
	       chr = 'X';
		   while(chr != '\n')
		       fscanf(pfile, "%c", &chr);	// read the line, proceed to next in the end
	   } // for ix
	Both are slow, (a) is a bit slower (?) than (b) and must fit the entire line
	including EOL - otherwise it does not jump to the next line. But it is shorter.

	Some acceleration can be acived by skippink the loop over 7 columns:
	%f %f %f ... %f - repeat 7 times, increment 'irgrid' by 7.

REFS:
	-
------------------------------------------------------------------------------*/
//
	FILE
		*pfile;
	int const
		str_len_max = 256;
	int
		irefre, irefim, ir, ix, ixoffset;
	char
		chr;
	char
		str_buffer[str_len_max], fpath[path_len_max];
//-----------------------------------------------------------------------------
//
	strcpy(fpath, rootdir);
	strcat(fpath, fname);
	pfile = fopen(fpath, "r");
//	printf("\n(in read_Rkext_fix) opened: %s", fpath);
//
	for (ix = 0; ix < 33; ix++) // skip header: 33 lines (update if '25  a number of grid aspect ratios' changes)
        fgets(str_buffer, str_len_max, pfile);
//
	for (irefre = 0; irefre < nrefre; irefre++)
	{
		for (irefim = 0; irefim < nrefim; irefim++)
		{
            ixoffset = irefre*nrefim*nrgrid_fix+irefim*nrgrid_fix;
//
//          EXTINCTION part: skip 3 lines
			fgets(str_buffer, str_len_max, pfile);
			fgets(str_buffer, str_len_max, pfile);
			fgets(str_buffer, str_len_max, pfile);
			for (ir = 0; ir < nrgrid_fix; ir++)
			{
				fscanf(pfile, "%f", &ufext[ixoffset+ir]);
			} // for ir
			fscanf(pfile, "%c", &chr); // read '\n'
//
//          ABSORPTION part: skip 1 line
			fgets(str_buffer, str_len_max, pfile);
			for (ir = 0; ir < nrgrid_fix; ir++)
			{
				fscanf(pfile, "%f", &ufabs[ixoffset+ir]);
			} // for ir
			fscanf(pfile, "%c", &chr); // read '\n'
		}// for irefim
	} // for irefre
//
	fclose(pfile);
//	printf("\n(in read_Rkext_fix) closed: %s", fpath);
//
} // void read_Rkext_fix
/*------------------------------------------------------------------------------
23/02/18: fscanf ignores whitespaces - read all 'ngrid_fix*nsca_fix' elements at once;
          tested -- ok; however, there is no gain in time vs. previous versions!
22/04/04: float arrays instead of double, fgets() skips lines; tested -- ok;
22/04/01: first created and tested in 'matrix_fix()' -- ok.
------------------------------------------------------------------------------*/
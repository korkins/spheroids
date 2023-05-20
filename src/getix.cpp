//
int getix(double const x0, double const *x, int const nx, int &ix1, int &ix2)
/*--------------------------------------------------------------------------------------------------
TASK:
	To find inidices, ix1 ancd ix2, so that x0 belongs to x[ix1]:x[ix2].
IN:
	x0   d       value
	x    d[nx]   array
	nx   i       len(x)
OUT:
	ix1, ix2   i   indices
	code       i   error code
NOTE:
    1) 'x' is sorted in ascensing order: x[0] < x[1] < ... x[nx-1].
	2) code:
	    0 - all is good (default value)
	   -1 - x0 < x[0],    ix1 = ix2 = 0
	   +1 - x0 > x[nx-1], ix1 = ix2 = nx-1
REFS:
	-
--------------------------------------------------------------------------------------------------*/
{
	int
		code, ix;
	double
		xi, xmin, xmax;
//--------------------------------------------------------------------------------------------------
//
	code = 0;
	xmin = x[0];
	xmax = x[nx-1];
	if (x0 < xmin)
	{
		ix1 = 0;
		ix2 = 0;
		code = -1;
	}
	else if (x0 == xmin)
	{
		ix1 = 0;
		ix2 = 1;
		code = 0;
	}
	else if (x0 == xmax)
	{
		ix1 = nx-2;
		ix2 = nx-1;
		code = 0;
	}
	else if (x0 > xmax)
	{
		ix1 = nx-1;
		ix2 = nx-1;
		code = 1;
	}
	else
	{
		ix = 0;
		xi = x[0];
		while (xi < x0)
		{
			ix += 1;
			xi = x[ix];
		}
		ix1 = ix-1;
		ix2 = ix;
		code = 0;
	}
	return code;
} // int getix(...)
/*--------------------------------------------------------------------------------------------------
23/05/01: cosmetic changes.
23/03/05: first created and tested using 'refre' & 'refim' from spheroids.
--------------------------------------------------------------------------------------------------*/
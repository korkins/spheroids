//
double bilinear(double const x,   double const y,
                double const x1,  double const x2,
                double const y1,  double const y2,
                double const f11, double const f21, double const f12, double const f22)
/*--------------------------------------------------------------------------------------------------
TASK:
	To perform bilinear interpolation.
IN:
    x, y             d   point of interpolation
	x1, x2, y1, y2   d   points where f(x,y) is given: x1 < x2, y1 < y2
	f11 ... f22      d   values of function: fij = f(xi, yj), e.g. f12 = f(x1, y2)
OUT:
	bilinear         d   result
NOTE:
    The function assumes x & y are sufficiently different: dx & dy /= 0.0.
    The function was not tested for x2 < x1 and/or y2 < y1.
REFS:
	1. https://en.wikipedia.org/wiki/Bilinear_interpolation
--------------------------------------------------------------------------------------------------*/
{
    double
        u, f1, f2;
//--------------------------------------------------------------------------------------------------
//
    u = (x - x1) / (x2 - x1);
    f1 = f11 + u * (f21 - f11);
    f2 = f12 + u * (f22 - f12);
    return f1 + (f2 - f1) * (y - y1) / (y2 - y1);
} // double bilinear(...)
/*--------------------------------------------------------------------------------------------------
23/05/01: cosmetic changes
23/03/10: first created from bilinear.f90 and tested:
          a) https://java2blog.com/bilinear-interpolation-python/
             f11: f(0, 1) = 12
             f12: f(0, 3) = -4
             f21: f(4, 1) =  0
             f22: f(4, 3) =  8

             f(1,2) = 4 - ok.

          b) Check (a) using Python:
             x = np.array([[  0.0, 4.0 ], [  0.0, 4.0 ]])
             y = np.array([[  1.0, 1.0 ], [  3.0, 3.0 ]])
             z = np.array([[ 12.0, 0.0 ], [ -4.0, 8.0 ]])
             x0 = 1.0
             y0 = 2.0
             cf = intrp.interp2d(x, y, z, kind='linear')
             z0 = cf(x0, y0)
             print(z0)

             z0 = 4.0 - ok.

          c) f(x, y) = exp(-a*(x2+y2)), a = 0.5
             x1 = 0.1; y1 = 0.3; f11 = 0.95122942;
	         x1 = 0.1; y2 = 1.5; f12 = 0.32303326;
	         x2 = 1.0; y1 = 0.3; f21 = 0.57984178;
	         x2 = 1.0; y2 = 1.5; f22 = 0.19691168;

	         x = 0.7;
	         y = 0.8;

             z0  = 0.51001872 - ok.
             zpy = 0.51001872 - python's 'from scipy import interpolate as intrp'
             z   = 0.56836015 - exact analytical calculation

             Accuracy of the interpolation increases as the function flattens: a = 0.05, 0.005...
--------------------------------------------------------------------------------------------------*/
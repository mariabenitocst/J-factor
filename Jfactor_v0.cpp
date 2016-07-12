/****************************************************************************************

Computes the J-factor for a generalized NFW profile in a given region of interest (ROI).

The J-factor is the integral of the density profile's square over the line of sight.
It is a 3-dimensional integral computed via interpolation tables and integrated over
the longitude, latitude of the ROI and s or line of sight up to a maximum value s_max. 

The ROI is characterized by the longitude (l) and latitude (b) angles.

The generalized NFW dark matter density profile is characterized by 3 free parameters:
(1) the local dark matter density rho_0
(2) the inner slope alpha
(3) the scale radius

INPUTS: (QUIERO QUE MI PROGRAMA SEA ITERACTIVO Y ME LOS PIDA)
(1) R0 or Sun's galactocentric distance (kpc)
(2) s_max (kpc)
(3) alpha (dimentionless)
(4) rho0 (GeV/cm^3)
(5) rs (kpc)
(6) bmin (degrees)
(7) bmax (degrees)
(8) lmin (degrees)
(9) lmax (degrees)

OUTPUT:
(1) J-factor value (GeV/cm^5)


HAVE TO EXPLAIN THE METHOD OF INTERPOLATION TABLES


SELF-CONSISTENT:

(A) Calore + JCAP use R0 = 8.5 kpc. In all my derivations from the rotation curve I use 
    R0 = 8 kpc. The idea is that at the end of the day, we will study the dependence of the
    results with the local velocity and galactocentric distance (e.g. v0 and R0)

FUTURE TESTS:

(A) Precision: total/relative errors
(B) maximum s_max for which we saturate the integral. I am starting with a initial value 
    of s_max = 200 kpc

*****************************************************************************************/
#include <cmath>
#include <iostream>
#include <cstdio>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h> // I DO NOT KNOW IF
#include <gsl/gsl_interp.h> // I NEED BOTH OR JUST ONE OF THEM...lets see during coding
using namespace std;

struct params {double R0; double alpha; double rho0; double rs; double b; double l;};

/**************************************************
      
     		Interpolation function

***************************************************/


double  f_interpolation(double x, void*parameters)
{
gsl_spline *spline=(gsl_spline*) parameters;
gsl_interp_accel *acc=gsl_interp_accel_alloc();
double f=gsl_spline_eval(spline, x, acc);
gsl_interp_accel_free (acc);
return f;
}

/**************************************************
      
	Integration over the line of sight

***************************************************/

double lof(double s, void * parameters)
{
params p=*(params *) parameters;
double r=sqrt(pow(R0,2)+pow(s,2)-2*R0*s*cos(p.b)*cos(p.l));
double rhos=p.rho0*pow(p.R0/p.rs,p.alpha)*pow(1+p.R0/p.rs,3.-p.alpha);
return pow(rhos*(pow(r/p.rs, -p.alpha)*pow(1+r/p.rs, -3+p.alpha), 2);
}

double integration_lof(double R0, double s_max, double alpha, double rho0, double rs, \
double b, double l)
{
params p={R0, alpha, rho0, rs, b, l};
gsl_function F; F.function=&lof; F.params=&p;
gsl_integration_cquad_workspace* w=gsl_integration_cquad_workspace_alloc(10000);
double result, error;
size_t neval;
gsl_integration_cquad(&F, 0, s_max, 1e-9, 1e-3, w, &result, &error, &neval);
// I HAVE TO CHECK THE PRECISION FOR s_max, and the relative/absolute errors
gsl_integration_cquad_workspace_free(w);
return result;
}

/**************************************************
      
	Integration over the latitude b

***************************************************/

void ingerpolation_b(gsl_spline *spline,double R0, double s_max, double alpha, double rho0, \
double rs, double bmin, double bmax, double l)
{
int j=200;
double step=(bmax-bmin)/(j-1);
double b_values[j];
double integration_s[j];
for (int i=0; i<j; i++)
	{
        b_values[i]=bmin+i*step;
        integration_s[i]=integration_lof(R0, s_max, alpha, rho0, rs, bmin+k*step, l);
        }
gsl_spline *spline2=gsl_spline_alloc(gsl_interp_cspline, j);
gsl_spline_init(spline2, integration_theta_values, integration_R_phi0, j);
*spline=*spline2;

}


double integrate_b(double R0, double s_max, double alpha, double rho0, double rs, double bmin \
double bmax, double l)
{
gsl_spline spline;
interpolation_ext_theta_real(&spline, phi, l, m, r, x_0, y_0, z_0, alpha, r_max, r_0, rho_0);
gsl_function F;F.function=&f_interpolation;F.params=&spline;
double result, error;size_t neval;
gsl_integration_cquad_workspace *w=gsl_integration_cquad_workspace_alloc(1000);
gsl_integration_cquad(&F, bmin, bmax, 1e-3, 0, w, &result, &error, &neval);
// I HAVE TO CHECK THE PRECISION FOR s_max, and the relative/absolute errors
gsl_integration_cquad_workspace_free(w);
return result;

}

/**************************************************
      
	Integration over the longitude l

***************************************************/

void interpolation_l(gsl_spline *spline, double R0, double s_max, double alpha, double rho0 \
double rs, double bmin, double bmax, double lmin, double lmax)
{
int j=200;	     // number of integrations over b ans s. The longitude space is divided 
		     // into j points
double step=(lmax-lmin)/(j-1);
double l_values[j]; // array which contains the longitude's points at which we integrate over
                    // b and s 
double integration_s_b[j]; // values of the integration over s and b for each longitude point 
for (int k=0; k<j; k++)
	{
        l_values[k]=min_phi+k*step_phi;
        integration_s_b[k]=integration_b(R0, s_max, alpha, rho0, rs, bmin, bmax, lmin+k*step);
        }
gsl_spline *spline2=gsl_spline_alloc(gsl_interp_cspline, j);
gsl_spline_init(spline2, l_values, integration_s_b, j);
*spline=*spline2;
}

double integrate_l(double R0, double s_max, double alpha, double rho0, double rs, double bmin \
double bmax, double lmin, double lmax)
{
gsl_spline spline;
interpolation_int_phi_real(&spline, l, m, r, x_0, y_0, z_0, alpha, r_max, r_0, rho_0);
gsl_function F;F.function=&f_interpolation;F.params=&spline;
double result, error;size_t neval;
gsl_integration_cquad_workspace *w=gsl_integration_cquad_workspace_alloc(1000);
gsl_integration_cquad(&F, lmin, lmax, 1e-3, 0, w, &result, &error, &neval);
// I HAVE TO CHECK THE PRECISION FOR s_max, and the relative/absolute errors
gsl_integration_cquad_workspace_free(w);
return (l*pow(r, l))*2*result;
}


/**************************************************
      
		MAIN FUNCTION

***************************************************/


int main()
{
// INPUT PARAMETERS:
double R0=8.  // (kpc)In Calore + JCAP 2015 uses R0=8.5 kpc. This is something I need to take 
	      //into account in order to be self-consistent.
double s_max= 200. // kpc			
double alpha = 1.; double rho0 = 0.4;
double rs = 20.;
// Change degrees to radians :


return 0;
}

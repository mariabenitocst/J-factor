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

 TESTS:

(A) Precision: total/relative errors
(B) maximum s_max for which we saturate the integral. I am starting with a initial value 
    of s_max = 200 kpc
(C) I obtain the same results for the average geometrical J-factor in 0904.3830

:: Note:: It is important to notice that tests (A) and (B) should be done each time 
          the J-factor value is computed


NOTE - ERRORS:

In funtion integrate_l(), when integration through cquad method I was obtaining the 
following error:
>> gsl: interp.c:150: ERROR: interpolation error
>> Default GSL error handler invoked.
>> Abort trap: 6
I change the integration method to quad and obtain the error:
>> gsl: qag.c:119: ERROR: iteration limit exceeds available workspace
>> Default GSL error handler invoked.
>> Abort trap: 6
This error disappears when I increase the available workspace.

The same happen, when increasing the number of points in the interpolation (j) in
the function integrate_b()

I do not know what is causing this error. When I plot the interpolation of b (before 
integration in l), the curve is 
relatively smooth and there is no-oscillation.


¡¡¡¡¡¡¡¡¡QUESTIONS!!!!!!!!

(1) Is the factor sin(M_PI/2 -b) ok?

*****************************************************************************************/

#include <cmath>
#include <iostream>
#include <fstream>
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
double r=sqrt(pow(p.R0,2)+pow(s,2)-2*p.R0*s*cos(p.b)*cos(p.l));
double rhos=p.rho0*pow(p.R0/p.rs,p.alpha)*pow(1+p.R0/p.rs,3.-p.alpha);
return pow(rhos*pow(r/p.rs, -p.alpha)*pow(1+r/p.rs, -3+p.alpha), 2)*sin(M_PI/2.-p.b);
// The factor sin(M_PI/2 - b) is due to the integral over the solid angle (I am not sure
// about this factor. I have obtain it by analogy with the solid angle in spherical coordinates
}

double integration_lof(double R0, double s_max, double alpha, double rho0, double rs, \
double b, double l)
{
params p={R0, alpha, rho0, rs, b, l};
gsl_function F; F.function=&lof; F.params=&p;
gsl_integration_cquad_workspace* w=gsl_integration_cquad_workspace_alloc(10000);
double result, error;
size_t neval;
gsl_integration_cquad(&F, 0, s_max, 1e-12, 1e-6, w, &result, &error, &neval);
gsl_integration_cquad_workspace_free(w);
return result;
}

/**************************************************
      
	Integration over the latitude b

***************************************************/


double integrate_b(double R0, double s_max, double alpha, double rho0, double rs, double bmin, \
double bmax, double l)
{
int j=400;
gsl_spline * spline = gsl_spline_alloc(gsl_interp_cspline, j);
double step=(bmax-bmin)/(j-1);
double b_values[j];
double integration_s[j];
for (int i=0; i<j; i++)
	{
        b_values[i]=bmin+i*step;
        integration_s[i]=integration_lof(R0, s_max, alpha, rho0, rs, bmin+i*step, l);
        }
gsl_spline_init(spline, b_values, integration_s, j);

gsl_function F; F.function=&f_interpolation; F.params=spline;

double result, error; size_t limit = sizeof(double);
gsl_integration_workspace * w = gsl_integration_workspace_alloc(10000000);
gsl_integration_qag(&F, bmin, bmax, 1e-6, 1e-6, limit, 0, w, &result, &error);
gsl_integration_workspace_free(w);

/* double result, error;size_t neval;
gsl_integration_cquad_workspace *w=gsl_integration_cquad_workspace_alloc(1000);
gsl_integration_cquad(&F, bmin, bmax, 1e-3, 0, w, &result, &error, &neval);
gsl_integration_cquad_workspace_free(w);*/
gsl_spline_free(spline);
return result;
}

double integrate_l(double R0, double s_max, double alpha, double rho0, double rs, double bmin, \
double bmax, double lmin, double lmax)
{
int j=400;  /* number of integrations over b ans s. The longitude space is divided 
               into j points*/	
gsl_spline * spline = gsl_spline_alloc(gsl_interp_cspline, j);
double step=(lmax-lmin)/(j-1);
double l_values[j]; // array which contains the longitude's points at which we integrate over
                    // b and s 
double integration_s_b[j]; // values of the integration over s and b for each longitude point 
for (int i=0; i<j; i++)
	{
        l_values[i]=lmin+i*step;
        integration_s_b[i]=integrate_b(R0, s_max, alpha, rho0, rs, bmin, bmax, lmin+i*step);
        }
gsl_spline_init(spline, l_values, integration_s_b, j);
/*
ofstream B("interpolationB.dat");
if (B.is_open())
{
	B<<"# L_values     integration_s_b_values"<<endl;
	for (int i=0; i<j; i++)
		B<<l_values[i]<<"   "<<integration_s_b[i]<<endl;
}
B.close();
return 8.;*/
/*
gsl_function F; F.function=&f_interpolation; F.params=spline;
double result, error;size_t neval;
gsl_integration_cquad_workspace *w=gsl_integration_cquad_workspace_alloc(100000);
gsl_integration_cquad(&F, lmin, lmax, 1e-3, 0, w, &result, &error, &neval);
gsl_integration_cquad_workspace_free(w);
*/
gsl_function G; G.function=&f_interpolation; G.params=spline;
double result, error; size_t limit = sizeof(double);
gsl_integration_workspace * w = gsl_integration_workspace_alloc(10000000);
gsl_integration_qag(&G, lmin, lmax, 1e-6, 1e-6, limit, 0, w, &result, &error);
gsl_integration_workspace_free(w);
gsl_spline_free(spline);
// 1 kpc = 3.0857e21 cm 
return result*(3.0857e21);
}


/**************************************************
      
		MAIN FUNCTION

***************************************************/

int main()
{
double R0=8.;  // (kpc)In Calore + JCAP 2015 uses R0=8.5 kpc. This is something I need to take 
	      //into account in order to be self-consistent.
// Region Of Interest (ROI) (see 1409.0042)
double bmin = 2*M_PI/180.; double bmax = 20*M_PI/180.; 
double lmin = 0*M_PI/180.; double lmax = 20.*M_PI/180.; // (I should multiply the result by 4)

double smax = 400.;

// Morphology AiX			
//double alpha = 1.006; double rho0 = 0.397; double rs = 63.540;
//double smax= 400.; // kpc
// AjX
//double alpha =0.835; double rho0 = 0.460; double rs = 10.006;
//double smax=400.;
//AkX
//double alpha = 1.769; double rho0 = 0.212; double rs = 99.926;
//double smax = 400.;
//AlX
//double alpha = 0.951; double rho0 = 0.427; double rs = 10.002;
//double smax = 400.;
//AmX
//double alpha = 0.907; double rs = 10.011; double rho0 = 0.399;
//double smax = 400.;

//BiX
//double alpha = 0.862; double rs = 35.488; double rho0 = 0.399;
//double smax = 400.;
//BjX
//double alpha = 0.808; double rs = 10.001; double rho0 = 0.458;
//double smax = 400.; 
//BkX
//double alpha = 1.756; double rs = 99.584; double rho0 = 0.212;
//double smax = 400.;
/* BlX */ //double alpha = 0.919; double rs = 10.002; double rho0 = 0.427;
/* BmX */ //double alpha = 0.873; double rs = 10.;    double rho0 = 0.399;
/* CiX */ //double alpha = 0.423; double rs = 13.858; double rho0 = 0.398;
/* CjX */ //double alpha = 0.697; double rs = 10.004; double rho0 = 0.464;
/* CkX */ //double alpha = 1.730; double rs = 99.724; double rho0 = 0.216;
/* ClX */ //double alpha = 0.792; double rs = 10.001; double rho0 = 0.435; 
/* CmX */ //double alpha = 0.647; double rs = 10.002; double rho0 = 0.412;
/* DiX */ //double alpha = 0.002; double rs = 10.602; double rho0 = 0.402;
/* DjX */ //double alpha = 0.719; double rs = 10.005; double rho0 = 0.464;
/* DkX */ //double alpha = 1.621; double rs = 99.626; double rho0 = 0.221;
/* DlX */ //double alpha = 0.835; double rs = 10.006; double rho0 = 0.431;
/* DmX */ //double alpha = 0.757; double rs = 10.001; double rho0 = 0.405; 
/* EiX */ //double alpha = 0.996; double rs = 50.827; double rho0 = 0.397;
/* EjX */ //double alpha = 0.864; double rs = 10.005; double rho0 = 0.463;
/* EkX */ //double alpha = 1.774; double rs = 99.700; double rho0 = 0.213;
/* EmX */ //double alpha = 0.985; double rs = 10.004; double rho0 = 0.429;
/* ElX */ //double alpha = 0.942; double rs = 10.006; double rho0 = 0.401;
/* FiX */ //double alpha = 1.130; double rs = 99.497; double rho0 = 0.393;
/* FjX */ 
/* FkX */ //double alpha = 1.817; double rs = 99.611; double rho0 = 0.205;
/* FlX */ double alpha = 1.139; double rs = 10.006; double rho0 = 0.415;
/* FmX */

double Jfactor =  integrate_l(R0, smax, alpha, rho0, rs, bmin, bmax, lmin, lmax);
// The units of the J-factor are GeV/cm⁵

cout<<" The Jfactor is "<<4*Jfactor<<endl;

return 0;
}

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
#include <vector>
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
vector<double> alpha, rs, rho0;
vector<char> model;
ifstream E("/Users/mariajosebenitocastano/Desktop/DATA/BAYESIAN/MAPestimationsR8v214.dat");
string line;
while ( getline(E, line))
{
	if (line.empty() || line[0] == '#')
		continue;
	double a=0., r=0., rho=0.;
	char reference[10];
	sscanf(line.c_str(), "%s %lf %lf %lf", reference, &a, &r, &rho);
	model.push_back(*reference);	     
	alpha.push_back(a);
        rs.push_back(r);
	rho0.push_back(rho);
}
ofstream J("JfactorR8v214.dat");
if (J.is_open())
{
	J<<"# alpha   rs    rho0   J-factor (GeV/cm^3)"<<endl;
	for (int i = 0; i < int(alpha.size()); i++)
	{
	double Jfactor = integrate_l(R0, smax, alpha[i], rho0[i], rs[i], bmin, bmax, lmin, lmax);
	J<<model[i]<<"   "<<alpha[i]<<"   "<<rs[i]<<"   "<<rho0[i]<<"   "<<4*Jfactor<<endl;
	}
}
J.close();

// The units of the J-factor are GeV/cm⁵

return 0;
}

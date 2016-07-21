/*
 * =====================================================================================
 *
 *       Filename:  sphericalBulge.cpp
 *
 *    Description:  Studying the MOG RC for a potential with constant density up to 2.5 kpc;
 *
 *        Version:  1.0
 *        Created:  16-06-2015 12:12:06
 *       Revision:  none
 *       Compiler:  g++
 *
 *         Author:  MJ Benito, 
 *   Organization:  IFT-UNESP
 *
 * =====================================================================================
 */

//#include "mpi.h"
#include "spherical.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstdio>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_sf_bessel.h>
#include <vector>
#include <iomanip>
using namespace std;

//Defining constructors and destructors:
sphericalBulge::sphericalBulge(double R_sun, double alpha, double mu)
{
	set_variables(R_sun, alpha, mu);
}

sphericalBulge::~sphericalBulge()
{
	//do nothing
}

void sphericalBulge::set_variables(double R_s, double a, double m)
{
	R_sun = R_s;
	rmax = 2.5; Mtotal = 4.6e9;
	rho0 = Mtotal / (4./3.*M_PI*pow(rmax, 3)); // [rho0] = M_sun/kpc^3 (ARE THIS UNITS OK??)
	alpha = a; mu = m;
}


/***************************************************************************************
*******     INTERPOLATION FUNCTION, equal for all the cases               ************* 
***************************************************************************************/


double  f_interpolation(double x, void*parameters)
{
gsl_spline *spline=(gsl_spline*) parameters;
gsl_interp_accel *acc=gsl_interp_accel_alloc();
double f=gsl_spline_eval(spline, x, acc);
//gsl_spline_free (spline);
gsl_interp_accel_free (acc);
return f;
}


struct params_R {double theta; double phi; int l; int m; double alpha; double mu;};
/* The radial gradient has two integrals corresponding to the internal and external potential. Each of those is, furthermore, divided into two terms corresponding to the real and imaginary part of the spherical harmonic*/

/***************************************************************************
		   
				INTEGRATES 

***************************************************************************/

//  EXTERNAL POTENTIAL:
double integrate_ext_R_real(double a, void*parameters)
{
	params_R p=*(params_R*) parameters;

	double fr;
	if (a < 2.5) 	fr = 1.;
	else		fr = 0.;
	
	return pow(a, p.l+2)*sin(p.theta)*gsl_sf_legendre_sphPlm(p.l, p.m, cos(p.theta))*cos(-p.m*p.phi)*fr;
}

double integrate_ext_R_realMOG(double a, void*parameters)
{
	params_R p=*(params_R*) parameters;
  
	double fr;
	if (a < 2.5) 	fr = 1.;
	else		fr = 0.;

	return pow(a, 2)*sin(p.theta)*gsl_sf_legendre_sphPlm(p.l, p.m, cos(p.theta))*cos(-p.m*p.phi)* \
		gsl_sf_bessel_il_scaled(p.l, p.mu*a)*fr / exp(-p.mu*a);
}




//  INTERNAL POTENTIAL:

double integrate_int_R_real(double a, void*parameters)
{
	params_R p=*(params_R*) parameters;
	double fr;
	if (a < 2.5) 	fr = 1.;
	else		fr = 0.;

	return pow(a, 1-p.l)*sin(p.theta)*gsl_sf_legendre_sphPlm(p.l, p.m, cos(p.theta))*cos(-p.m*p.phi)*fr;
}

double integrate_int_R_realMOG(double a, void*parameters)
{
	params_R p=*(params_R*) parameters;
	double fr;
	if (a < 2.5) 	fr = 1.;
	else		fr = 0.;

	return pow(a, 2)*sin(p.theta)*gsl_sf_legendre_sphPlm(p.l, p.m, cos(p.theta))*cos(-p.m*p.phi)* \
		gsl_sf_bessel_kl_scaled(p.l, p.mu*a)*fr / exp(p.mu*a);
}


/*********************************************************************
**************   INTEGRATION IN R ***********************************
*********************************************************************/

//EXTERNAL POTENTIAL: 

double integration_ext_R_real(double theta,double phi,int l,int m,double r, double alpha, double mu, double rho0)
{
//cout<<"Entrando en integration_ext_R_real..."<<endl;
	params_R parameters={theta, phi, l, m, alpha, mu};
	gsl_function F; F.function=&integrate_ext_R_real; F.params=&parameters;
	double result, error; size_t neval;
	gsl_integration_cquad_workspace * w= gsl_integration_cquad_workspace_alloc (1000);
	gsl_integration_cquad(&F, 0., r, 1e-7, 0, w, &result, &error, &neval);
	gsl_integration_cquad_workspace_free (w);

	return rho0*(1+alpha)*result;
}

double integration_ext_R_realMOG(double theta,double phi,int l,int m,double r, double alpha, double mu, double rho0)
{	
	params_R parameters={theta, phi, l, m, alpha, mu};
	gsl_function F; F.function=&integrate_ext_R_realMOG; F.params=&parameters;	

	double result, error; size_t neval;
	gsl_integration_cquad_workspace * w= gsl_integration_cquad_workspace_alloc (1000);
	gsl_integration_cquad(&F, 0., r, 1e-7, 0, w, &result, &error, &neval);
	gsl_integration_cquad_workspace_free (w);

//cout<< " I am inside the integration over r and result = "<<result<<endl;
	return rho0*alpha*pow(mu, 2)*result;
}



//INTERNAL POTENTIAL:

double integration_int_R_real(double theta,double phi,int l,int m,double r, double alpha, double mu, double rho0)
{
	//cout<<"Entrando en integration_int_R_real..."<<endl;
	params_R parameters={theta, phi, l, m, alpha, mu};
	gsl_function F; F.function=&integrate_int_R_real; F.params=&parameters;
	double result, error;
//	if(r<4.)
//		{
		size_t neval;
		gsl_integration_cquad_workspace * w= gsl_integration_cquad_workspace_alloc (1000);
		gsl_integration_cquad(&F, r, 2.5, 1e-10, 0, w, &result, &error, &neval);
		gsl_integration_cquad_workspace_free (w);
//		}
/*	else
		{
		size_t limit=sizeof(double);
		gsl_integration_workspace * w= gsl_integration_workspace_alloc (1000);
		gsl_integration_qagiu (&F, r, 1e-3, 0, limit, w, &result, &error);
		gsl_integration_workspace_free (w);*/
//		}

	return rho0*(1+alpha)*result;
}

double integration_int_R_realMOG(double theta,double phi,int l,int m,double r, double alpha, double mu, double rho0)
{
	params_R parameters={theta, phi, l, m, alpha, mu};

	gsl_function F; F.function=&integrate_int_R_realMOG; F.params=&parameters;
	double result, error;
//	if(r<4.)
//		{
		size_t neval;
		gsl_integration_cquad_workspace * w= gsl_integration_cquad_workspace_alloc (1000);
		gsl_integration_cquad(&F, r, 2.5, 1e-10, 0, w, &result, &error, &neval);
		gsl_integration_cquad_workspace_free (w);
//		}
/*	else
		{
		size_t limitG=sizeof(double);
		gsl_integration_workspace * w= gsl_integration_workspace_alloc (1000);
		gsl_integration_qagiu (&F, r, 1e-3, 0, limit, w, &result, &error);
		gsl_integration_workspace_free (w);*/
//		}


	return rho0*alpha*pow(mu,2)*result;
}


/***************************************************************************************

  			 INTEGRATING ON THETA the interpolation       

***************************************************************************************/

//EXTERNAL POTENTIAL:

double integration_ext_theta_real(double phi,int l,int m,double r,double alpha, double mu, double rho0)
{
	int j=30;
	gsl_spline *spline=gsl_spline_alloc(gsl_interp_cspline, j);

	double min_theta=0;double max_theta=M_PI;double step_theta=(max_theta-min_theta)/(j-1);
	double integration_R_phi0[j];double integration_theta_values[j];
	for (int i=0; i<j; i++){
	        integration_theta_values[i]=min_theta+i*step_theta;
	        integration_R_phi0[i]=integration_ext_R_real(min_theta+i*step_theta, phi, l, m, r, alpha, mu, rho0);
	        }
	gsl_spline_init(spline, integration_theta_values, integration_R_phi0, j);

	gsl_function F; F.function=&f_interpolation; F.params=spline;
	double result, error; size_t neval;
	gsl_integration_cquad_workspace *w=gsl_integration_cquad_workspace_alloc(1000);
	gsl_integration_cquad(&F, 0., M_PI, 1e-3, 0, w, &result, &error, &neval);
	gsl_integration_cquad_workspace_free(w);
	gsl_spline_free(spline);
//cout<<" integration in theta = "<<result<<endl;
	return result;
}

double integration_ext_theta_realMOG(double phi,int l,int m,double r,double alpha, double mu, double rho0)
{
	int j=30;
	gsl_spline *spline=gsl_spline_alloc(gsl_interp_cspline, j);

	double min_theta=0;double max_theta=M_PI;double step_theta=(max_theta-min_theta)/(j-1);
	double integration_R_phi0[j];double integration_theta_values[j];
	for (int i=0; i<j; i++){
	        integration_theta_values[i]=min_theta+i*step_theta;
	        integration_R_phi0[i]=integration_ext_R_realMOG(min_theta+i*step_theta, phi, l, m, r, alpha, mu, rho0);
	        }
	gsl_spline_init(spline, integration_theta_values, integration_R_phi0, j);

	gsl_function F; F.function=&f_interpolation; F.params=spline;
	double result, error; size_t neval;
	gsl_integration_cquad_workspace *w=gsl_integration_cquad_workspace_alloc(1000);
	gsl_integration_cquad(&F, 0., M_PI, 1e-3, 0, w, &result, &error, &neval);
	gsl_integration_cquad_workspace_free(w);
	gsl_spline_free(spline);
	return result;
}


//INTERNAL POTENTIAL:

double integration_int_theta_real(double phi,int l,int m, double r, double alpha, double mu, double rho0)
{
	int j=30;
	gsl_spline * spline=gsl_spline_alloc(gsl_interp_cspline, j);
	double min_theta=0;double max_theta=M_PI;double step_theta=(max_theta-min_theta)/(j-1);
	double integration_R_phi0[j];double integration_theta_values[j];
	for (int i=0; i<j; i++){
	        integration_theta_values[i]=min_theta+i*step_theta;
	        integration_R_phi0[i]=integration_int_R_real(min_theta+i*step_theta, phi, l, m, r, alpha, mu, rho0);
	        }
	gsl_spline_init(spline, integration_theta_values, integration_R_phi0, j);
	
	gsl_function F; F.function=&f_interpolation; F.params=spline;
	double result, error;size_t neval;
	gsl_integration_cquad_workspace *w=gsl_integration_cquad_workspace_alloc(1000);
	gsl_integration_cquad(&F, 0., M_PI, 1e-3, 0, w, &result, &error, &neval);
	gsl_integration_cquad_workspace_free(w);
	gsl_spline_free(spline);
	return result;
}

double integration_int_theta_realMOG(double phi,int l,int m, double r, double alpha, double mu, double rho0)
{
	int j=30;
	gsl_spline * spline=gsl_spline_alloc(gsl_interp_cspline, j);
	double min_theta=0;double max_theta=M_PI;double step_theta=(max_theta-min_theta)/(j-1);
	double integration_R_phi0[j];double integration_theta_values[j];
	for (int i=0; i<j; i++){
	        integration_theta_values[i]=min_theta+i*step_theta;
	        integration_R_phi0[i]=integration_int_R_realMOG(min_theta+i*step_theta, phi, l, m, r, alpha, mu, rho0);
	        }
	gsl_spline_init(spline, integration_theta_values, integration_R_phi0, j);
	
	gsl_function F; F.function=&f_interpolation; F.params=spline;
	double result, error;size_t neval;
	gsl_integration_cquad_workspace *w=gsl_integration_cquad_workspace_alloc(1000);
	gsl_integration_cquad(&F, 0., M_PI, 1e-3, 0, w, &result, &error, &neval);
	gsl_integration_cquad_workspace_free(w);
	gsl_spline_free(spline);
	return result;
}




/******************************************************************************************
******************************************************************************************/

// EXTERNAL POTENTIAL:

double ext_integral_real(int l,int m, double r, double alpha, double mu, double rho0)
{
	int j=200;
	gsl_spline * spline=gsl_spline_alloc(gsl_interp_cspline, j);
	double min_phi=0; double max_phi=M_PI;
	double step_phi=(max_phi-min_phi)/(j-1);
	double integration_R_theta[j];double integration_phi_values[j];
	for (int k=0; k<j; k++)
		{
        	integration_phi_values[k]=min_phi+k*step_phi;
        	integration_R_theta[k]=integration_ext_theta_real(min_phi+k*step_phi, l, m, r, alpha, mu, rho0);
       		}
	gsl_spline_init(spline, integration_phi_values, integration_R_theta, j);

	gsl_function F; F.function= &f_interpolation; F.params=spline;
	double result, error;size_t neval;
	gsl_integration_cquad_workspace *w=gsl_integration_cquad_workspace_alloc(1000);
	gsl_integration_cquad(&F, 0., M_PI, 1e-10, 0, w, &result, &error, &neval);
	gsl_integration_cquad_workspace_free(w);
	gsl_spline_free(spline);
//cout<<" l = "<<l<<"  m = "<<m<<"  r = "<<"   alpha "<<alpha<<"  mu "<<mu<<" rho0 "<<rho0<<endl;
//cout<<" integration in phi = "<<result<<endl;
	return (-(l+1)*pow(r,-l-1))*2*result;
}

double ext_integral_realMOG(int l,int m, double r, double alpha, double mu, double rho0)
{
	int j=200;
	gsl_spline * spline=gsl_spline_alloc(gsl_interp_cspline, j);
	double min_phi=0; double max_phi=M_PI;
	double step_phi=(max_phi-min_phi)/(j-1);
	double integration_R_theta[j];double integration_phi_values[j];
	for (int k=0; k<j; k++)
		{
        	integration_phi_values[k]=min_phi+k*step_phi;
        	integration_R_theta[k]=integration_ext_theta_realMOG(min_phi+k*step_phi, l, m, r, alpha, mu, rho0);
       		}
	gsl_spline_init(spline, integration_phi_values, integration_R_theta, j);

	gsl_function F; F.function= &f_interpolation; F.params=spline;
	double result, error;size_t neval;
	gsl_integration_cquad_workspace *w=gsl_integration_cquad_workspace_alloc(1000);
	gsl_integration_cquad(&F, 0., M_PI, 1e-10, 0, w, &result, &error, &neval);
	gsl_integration_cquad_workspace_free(w);
	gsl_spline_free(spline);
//cout<<" l = "<<l<<"  m = "<<m<<"  r = "<<"   alpha "<<alpha<<"  mu "<<mu<<" rho0 "<<rho0<<endl;
//cout<<" integration in phi = "<<result<<endl;
	return -r*2*result;
}



// INTERNAL POTENTIAL:

double int_integral_real(int l,int m, double r, double alpha, double mu, double rho0)
{
	int j=200;
	gsl_spline * spline=gsl_spline_alloc(gsl_interp_cspline, j);
	double min_phi=0; double max_phi=M_PI;
	double step_phi=(max_phi-min_phi)/(j-1);
	double integration_R_theta[j];double integration_phi_values[j];
	for (int k=0; k<j; k++)
		{
	        integration_phi_values[k]=min_phi+k*step_phi;
	        integration_R_theta[k]=integration_int_theta_real(min_phi+k*step_phi, l, m, r, alpha, mu, rho0);
	        }
	gsl_spline_init(spline, integration_phi_values, integration_R_theta, j);
	
	gsl_function F; F.function= &f_interpolation; F.params=spline;
	double result, error;size_t neval;
	gsl_integration_cquad_workspace *w=gsl_integration_cquad_workspace_alloc(1000);
	gsl_integration_cquad(&F, 0., M_PI, 1e-3, 0, w, &result, &error, &neval);
	gsl_integration_cquad_workspace_free(w);
	gsl_spline_free(spline);

//cout<<" l = "<<l<<"  m = "<<m<<"  r = "<<"   alpha "<<alpha<<"  mu "<<mu<<" rho0 "<<rho0<<endl;
//cout<<" integration in phi = "<<result<<endl;

	return (l*pow(r, l))*2*result;
}

double int_integral_realMOG(int l,int m, double r, double alpha, double mu, double rho0)
{
	int j=200;
	gsl_spline * spline=gsl_spline_alloc(gsl_interp_cspline, j);
	double min_phi=0; double max_phi=M_PI;
	double step_phi=(max_phi-min_phi)/(j-1);
	double integration_R_theta[j];double integration_phi_values[j];
	for (int k=0; k<j; k++)
		{
	        integration_phi_values[k]=min_phi+k*step_phi;
	        integration_R_theta[k]=integration_int_theta_realMOG(min_phi+k*step_phi, l, m, r, alpha, mu, rho0);
	        }
	gsl_spline_init(spline, integration_phi_values, integration_R_theta, j);
	
	gsl_function F; F.function= &f_interpolation; F.params=spline;
	double result, error;size_t neval;
	gsl_integration_cquad_workspace *w=gsl_integration_cquad_workspace_alloc(1000);
	gsl_integration_cquad(&F, 0., M_PI, 1e-3, 0, w, &result, &error, &neval);
	gsl_integration_cquad_workspace_free(w);
	gsl_spline_free(spline);
	return -r*2*result;
}



/******************************************************************

	     PUTTING EVERYTHING TOGETHER  

******************************************************************/

double sphericalBulge::velocity(double r, double theta, double phi)
{
int l_max=2;
double result_integral=0;
for(int l=0; l<(l_max+1); l+=2)
	{
        for(int m=0; m<(l+1); m+=2)
		{
                double result=0;
                gsl_sf_result result_A_factor, irr_l, irr_l1, reg_l, reg_l1;
                gsl_sf_legendre_sphPlm_e (l, m, cos(theta), &result_A_factor);
		gsl_sf_bessel_kl_scaled_e(l+1, mu*r, &irr_l1);
		gsl_sf_bessel_kl_scaled_e(l, mu*r, &irr_l);
		gsl_sf_bessel_il_scaled_e(l+1, mu*r, &reg_l1);
		gsl_sf_bessel_il_scaled_e(l, mu*r, &reg_l);
                double factor_A = -((4*M_PI)/(2*l+1))*result_A_factor.val;
		double factor_ext = (-irr_l1.val + l/(mu*r)*irr_l.val)*(2/M_PI)/exp(mu*r);
//		double factor_ext = 1.;
		double factor_int= (reg_l1.val + l/(mu*r)*reg_l.val)*(2/M_PI)/exp(-mu*r);
//cout<<" factor_ext = "<<factor_ext<<endl;
//cout<<" factor_int = "<<factor_int<<endl;
		if (r < 2.5)
		{
                if(m==0)
			{

                        result=(factor_ext * ext_integral_realMOG(l, m, r, alpha, mu, rho0) + factor_int * int_integral_realMOG(l, m, r, alpha, mu, rho0)) * factor_A;
//			result += (ext_integral_real(l, m, r, alpha, mu, rho0) + int_integral_real(l, m, r, alpha, mu, rho0)) * factor_A;
/*
                        result=(factor_ext * ext_integral_realMOG(l, m, r, alpha, mu, rho0) + \
				ext_integral_real(l, m, r, alpha, mu, rho0) + \
				factor_int * int_integral_realMOG(l, m, r, alpha, mu, rho0) + \
				int_integral_real(l, m, r, alpha, mu, rho0)) * factor_A;
*/
                        }
                else
			{

                        result=(cos(m*phi)+cos(-m*phi))* factor_ext * ext_integral_realMOG(l, m, r, alpha, mu, rho0)* factor_A;
                        result+=(cos(m*phi)+cos(-m*phi))* factor_int * int_integral_realMOG(l, m, r, alpha, mu, rho0) * factor_A;	
//                        result+=(cos(m*phi)+cos(-m*phi))* ext_integral_real(l, m, r, alpha, mu, rho0)* factor_A;
//                        result+=(cos(m*phi)+cos(-m*phi))* int_integral_real(l, m, r, alpha, mu, rho0) * factor_A;		
/*
                        result=(cos(m*phi)+cos(-m*phi))* (ext_integral_real(l, m, r, alpha, mu, rho0) + \
			       factor_ext * ext_integral_realMOG(l, m, r, alpha, mu, rho0))* factor_A;

                        result+=(cos(m*phi)+cos(-m*phi))* (int_integral_real(l, m, r, alpha, mu, rho0) + \
			       factor_int * int_integral_realMOG(l, m, r, alpha, mu, rho0)) * factor_A;
*/

                        }
                result_integral+=result;
                }
		if (r >= 2.5)
			{
                	if(m==0)
				{
                        	result=factor_ext * ext_integral_realMOG(l, m, r, alpha, mu, rho0)*factor_A;

//                        	result=(factor_ext * ext_integral_realMOG(l, m, r, alpha, mu, rho0) + ext_integral_real(l, m, r, alpha, mu, rho0))*factor_A;
                        	}
                	else
				{

                        	result=(cos(m*phi)+cos(-m*phi))* factor_ext * ext_integral_realMOG(l, m, r, alpha, mu, rho0) * factor_A;
//                        	result=(cos(m*phi)+cos(-m*phi))* (ext_integral_real(l, m, r, alpha, mu, rho0) + factor_ext * ext_integral_realMOG(l, m, r, alpha, mu, rho0)) * factor_A;
                        	}
                	result_integral+=result;
                	}

		}
        }
double G=4.51836e-39;
// G = 4.51..e-39 Kpc^3 / M_sun * s^2
	if (result_integral < 0)
		return -sqrt(abs(G*result_integral))*3.08567758e16;
	else
		return sqrt(G*result_integral)*3.08567758e16;
}


void sphericalBulge::rc(double rmin, double rmax, std::ofstream &outputfile)
{
	double theta=M_PI/2;int N_phi=4;//N_phi defines how many phi directions to average the rotation curve
     if(outputfile.is_open())
	{
           outputfile<<"#v (km/s)       r (kpc)"<<endl;
//	   if (rmin <= 20. && rmax <= 20.)
//		{
        	for(double r=rmin ; r < rmax ;r+=1.)
			{
			outputfile.precision(11); 
			double velocity_phi=0;
			for(double phi=0;phi<M_PI;phi+=2*M_PI/N_phi)
			velocity_phi+=2*pow(velocity(r,theta,phi),2);
			double velocity_=sqrt(velocity_phi/N_phi);
//			double error=velocity_*uncertainty_rho_0/(2*rho_0);
               		
			outputfile<<velocity_<<"   "<<r<<"\n"; 
			}
	}
outputfile.close();
}


int main()
{
	double R0=8.;
	sphericalBulge bulge(R0, 8.89, 0.042);
	// The values of alpha = 8.89 and mu = 0.042 are taken from 1306.6383
	double theta = M_PI/2.; double phi = 0.;

	cout<<bulge.velocity(2., theta, phi)<<endl;
	cout<<bulge.velocity(10., theta, phi)<<endl;

	std::ofstream RC("exponential.dat");
	bulge.rc(1., 14., RC);


/*	int myid, numprocess, source, dest, tag;
	double step, rmin, rmax;
	std::ostringstream filename;
	std::ifstream in;
	string line;
	vector <double> x, y, delta_y;
	int master = 0;
	double buffer[2];			// this is the message send by master containing the values rmin and rmax for 
						// which the corresponding processor is going to compute the rotation curve 
						// buffer[0] = rmin; buffer[1] = rmax (the rmax is the rmin in the next call)
	MPI::Status status;
	
	MPI::Init();
	myid = MPI::COMM_WORLD.Get_rank();
	numprocess = MPI::COMM_WORLD.Get_size();

	if (myid == master)
		{
		step = 2.;	
		rmin = 1.; rmax = 1.;
		for (dest = 1; dest < numprocess; dest ++)
			{
			rmax += step;
			buffer[0] = rmin; buffer[1] = rmax;
			rmin = rmax; 
			MPI::COMM_WORLD.Send(&buffer, 2, MPI::DOUBLE, dest, dest);
			}
		for (source = 1; source < numprocess; source ++)
			{
			MPI::COMM_WORLD.Recv(MPI_BOTTOM, 0, MPI::INT, source, 0);
			filename << "sphericalrc_" << source <<".dat";
			in.open(filename.str().c_str()); 
			while (getline(in , line))
				{
				if (line.empty() || line[0] == '#')
					continue;
				double r=0., v=0.;
				sscanf(line.c_str(), "%lf  %lf", &v, &r);
				y.push_back(v);
				x.push_back(r);
				}
			in.close();
			in.clear();
			filename.str(""); // reset the string to be empty
			filename.clear(); // clear any error flags that may be set
			}
		ofstream RC("spherical_rc.dat");
		if (RC.is_open())
			{
			RC.precision(11);
			RC << "# v (km/s)     error (km/s)     R (kpc)"<<endl;
			for (int i = 0; i < int(x.size()); i++)
				RC<< y[i] <<"   "<< x[i]<<endl;
			}
		RC.close();
		}
	else
		{
		MPI::COMM_WORLD.Recv(&buffer, 2, MPI::DOUBLE, master, MPI_ANY_TAG, status);

		tag = status.Get_tag();
		filename << "sphericalrc_" << tag <<".dat";
		std::ofstream out(filename.str().c_str());
		bulge.rc(buffer[0], buffer[1], out);
		MPI::COMM_WORLD.Send(MPI_BOTTOM, 0, MPI::INT, master, 0);					     		
		}	
	MPI::Finalize();*/
	
	return 0; 
}

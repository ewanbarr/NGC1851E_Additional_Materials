#include "cpgplot.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define AMAX 600000 /* dimensions for sub-arrays being displayed */ 
#define MAX 1000
#define MAX2 1800
#define MAX2H 498
#define NOMMAX 1800 
#define MAPMAX 2
#define MPBIN 400 
#define D1NCONT 7
#define D2NCONT 3
#define NCOSI 40
#define HNCOSI 20
#define NM 40
#define HNM 20
#define NH 40
#define HNH 20
#define SIGMAS 3
#define DAY 86400
#define YEAR 365.25
#define PI 3.141592653589793
#define TSUN 4.925490947
#define one3rd 0.33333333333
#define two3rd 0.66666666666
#define four3rd 1.33333333333
#define five3rd 1.66666666666
#define LIM 1E-10
#define DE 1E-5
#define MARGIN 1.0
#define SIZE0 6.5
#define SIZE1 3.5
#define SIZE2 2.5
#define GAP 0.3
#define PMARGIN 1.5

/* *************************************************************************************** */
/* *************************************************************************************** */
/* *************************************************************************************** */


    /* Welcome to the Bayesian Analysis Program. P. Freire, 2009. */

/* issues to address (order of decreasing priority):
     

   e) When smoothing chi2 maps, use real 2-D interpolation.

\
*/

/* first: define auxiliary functions */

float fm(float x, float s, float f, float m1)
{
  float r;
  r = x*x*x*s*s*s - f*x*x - 2*m1*f*x - f*m1*m1;
  return r;
}
 
/* define function that returns the array levels */

float levels(int kx, int ky, double total, float target, float vmax, float v[AMAX])
{
    int l;
    double x_left, x_right, x_midpoint, f_left, f_right, f_midpoint;
    double level, count;
    
    x_left = 0; f_left = 100 - target;

    x_right = vmax; f_right = 0 - target;

    while ( (fabs(x_left - x_right) > 0.00001 ))
      {
	/* calculate where midpoint is */
	x_midpoint = (x_left + x_right) / 2.0;
	
	/* OK, we have a function here for which we want to find a zero */
	
	/* start contour _counter_ from zero */
	count = 0;
	
	/* add all values above the treshold */
	for (l = 0; l < kx*ky; l++)
	    if (v[l] > x_midpoint)
		count = count + v[l];
	
	/* transform this into a percentage above D2target */
	f_midpoint = 100 * count/total - target;
	
	
	/* use bisection method to find where contour level should be */
	if ((f_left * f_midpoint) > 0 )
	  {
	    x_left = x_midpoint;
	    f_left = f_midpoint;
	  }
	else
	  {
	    x_right = x_midpoint;
	    f_right = f_midpoint;
	  }
      };
    
    level = x_midpoint;
    
    printf("\n Contour with height %f contains %f percent of total probability.\n", level, 100 * count/total);
    
    return level;    
}

/* *************************************************************************************** */
/* *************************************************************************************** */
/* *************************************************************************************** */

int main( )
{
    
  /* SCIENTIFIC FLAGS */
  
  /* do we know these values? Is it a simulation (ksim)?*/
  int komdt, kgamma, kpbdot, ksini, kmc, kh3, kh4, kxi;
  int kinetix, kr, kmp, ksim, ksimcount, kinom;
  
  /* SCIENTIFIC PARAMETERS */
  
  /* Keplerian  parameters */
  double f, ecc, pb, nb, asini;
  
  /* post-keplerian parameters (and uncertainties below) */
  double omdotobs, gammaobs, pbdotobs, siniobs, mcobs, h3obs, h4obs, xiobs;
  double domdotobs, dgammaobs, dpbdotobs, dsiniobs, dmcobs, dh3obs, dh4obs, dxiobs;
  
  /* positional and kinematic parameters */
  double raj, decj, pmobs, papmobs, pxobs, xdotobs, dxdotobs;
  
  /* mass ratio (if determined) */
  float mrobs, dmrobs;
  
  /* pulsar mass */
  float mpobs, dmpobs;
  
  /* input special points */
  int ssymbol[20], scolor[20];
  float sm1[20], sm2[20], ssini[20], sf[20];
  
  /* auxiliary variables */
  float omgrid, cosigrid, mtgrid, h3grid;
  float ra_min, ra_max, dec_min, dec_max, raaux, decaux, raaux_before, decaux_before, ra_diff_min, dec_diff_min;
  float sini, siniaux, sinimin, omaux, omaux_before, cosiaux, cosiaux_before, mtaux, mtaux_before, h3aux, h3aux_before;
  float om_min, om_max, cosi_min, cosi_max, mt_min, mt_max, h3_min, h3_max, om_diff_min, cosi_diff_min, mt_diff_min, h3_diff_min;
  float m2aux, m1aux, h4aux, xiaux;
  float omdotk, omdotaux, gammak, gammaux, pbdotaux, pbdotk, xdotaux, mt_aux;
  float mtobs, mtobsl, mtobsu, mtkinl, mtkinu; 
  float gamma_max, q_max, mp_max, mc_max, mtot_max, mc_min;
  float mp_gamma, mc_gamma, mc_gammau, mc_gammal, sini_gamma, sini_gammau, sini_gammal, i_gamma, i_gammau, i_gammal;
  float dmc, mr, M, incl1, incl2, cosimin;
  
  /* for displaying science lines */
  /* These should have more explicit names - code should migrate to these */
  float mp_sci[MAX+1], mc_sci[MAX+1], cosi_sci[MAX+1];
  /* These old names should preferentially be abandoned */
  float x_sci[MAX], y_sci[MAX], xy_sci[MAX], xy_sci_n[MAX], yx_sci[MAX], yx_sci_n[MAX], dx_sci, dy_sci, y_sci2[MAX2], xy_sci2[MAX2];
  float sigma[SIGMAS];
  sigma[0] = 0;
  sigma[1] = 1;
  sigma[2] = -1;

  /* END OF SCIENTIFIC PARAMETERS */

  /* number of arrays and types */
  int map_number, map_counter;
  int map_type[MAPMAX];

  /* files to be read and written */
  FILE *S;
  FILE *Fh;
  FILE *F;
  FILE *COS;
  FILE *CM;
  FILE *OM;
  FILE *PM;
  FILE *H3;
  FILE *H4;
  FILE *COSI_M2_1;
  FILE *COSI_M2_2;
  FILE *M1_M2_1;
  FILE *M1_M2_2;
  
 
  /* chi2 inputs */
  float v, dv, vmin, vaux;
  float vh, vhmin, vhaux;

  /* array calculation variables */
  int calculate_cosim2[MAPMAX], calculate_m1m2[MAPMAX], calculate_h3h4[MAPMAX], calculate_cosih3[MAPMAX], calculate_cosiom[MAPMAX], calculate_radec[MAPMAX];

  /* array counters */
  int i, j, k, l, t, lmax[MAPMAX], mcos, m, n, imp, imt, jmedian;
  int a;

  /*2-D array dimensions */
  int n_cosi[MAPMAX], n_m2[MAPMAX], n_mt[MAPMAX], n_h[MAPMAX], n_h3[MAPMAX], n_om[MAPMAX], n_ra[MAPMAX], n_dec[MAPMAX];
  float i_cosi[MAPMAX], i_m2[MAPMAX], i_mt[MAPMAX], i_h[MAPMAX], i_h3[MAPMAX], i_om[MAPMAX], i_ra[MAPMAX], i_dec[MAPMAX];
  float s_cosi[MAPMAX], s_m2[MAPMAX], s_mt[MAPMAX], s_h[MAPMAX], s_h3[MAPMAX], s_om[MAPMAX], s_ra[MAPMAX], s_dec[MAPMAX], dcosi[MAPMAX], dm2[MAPMAX], dmt[MAPMAX], dh3[MAPMAX];

  /* 2-D arrays */
  float cosim2[MAPMAX][AMAX], m1m2[MAPMAX][AMAX];
  float cosiom[MAPMAX][AMAX], h3h4[MAPMAX][AMAX], radec[MAPMAX][AMAX];
  /* matrix variables for displaying each array */
  float tr_cosim2[MAPMAX][6], tr_cosih3[MAPMAX][6], tr_m1m2[MAPMAX][6];
  float tr_cosiom[MAPMAX][6], tr_h3h4[MAPMAX][6];
  /* totals for each array */
  double cosim2_total[MAPMAX], m1m2_total[MAPMAX], cosiom_total[MAPMAX], h3h4_total[MAPMAX], m1om_total[MAPMAX], radec_total[MAPMAX], total;

  /* 2-D array maxima and minima */
  float cosim2max[MAPMAX], m1m2max[MAPMAX], cosiommax[MAPMAX], h3h4max[MAPMAX], radecmax[MAPMAX];
  /* 2-D array tallies */
  
  /* 1-D arrays */
  float cosi[MAPMAX][MAX], m2[MAPMAX][MAX], m1[MAPMAX][MPBIN], om[MAPMAX][NOMMAX], h3[MAPMAX][MAX], h4[MAPMAX][MAX], mt[MAPMAX][MAX];
  /* 1-D array dimensions, maxima and minima */

  /* for displaying 1-D probability functions */
  float h3x[MAPMAX][MAX], h4x[MAPMAX][MAX];
  float m2x[MAPMAX][MAX], m1x[MAPMAX][MPBIN], mtx[MAPMAX][MPBIN];
  float mcx_m[MAX+2], mpx_m[MAX+2], mtx_m[MAX+2]; 
  float cosix[MAPMAX][MAX], omx[MAPMAX][NOMMAX];
  float rax[MAPMAX][MAX], decx[MAPMAX][NOMMAX];

  float cosimax[MAPMAX], m2max[MAPMAX], ommax[MAPMAX], m1max[MAPMAX], h3max[MAPMAX], h4max[MAPMAX], mtmax[MAPMAX];
  /* 1-D tallies */
  float cumulative;

  /* percentiles of 1-d distributions */
  float D1target[D1NCONT];
  float D1_om[MAPMAX][D1NCONT], D1_omY[MAPMAX][D1NCONT];
  float D1_cosi[MAPMAX][D1NCONT], D1_cosiY[MAPMAX][D1NCONT];
  float D1_m2[MAPMAX][D1NCONT], D1_m2Y[MAPMAX][D1NCONT];
  float D1_m1[MAPMAX][D1NCONT], D1_m1Y[MAPMAX][D1NCONT];
  float D1_mt[MAPMAX][D1NCONT], D1_mtY[MAPMAX][D1NCONT];
 
  /* binning parameters */
  float cosi_l[MAPMAX], cosi_u[MAPMAX], m2_l[MAPMAX], m2_u[MAPMAX], m1_l[MAPMAX], m1_u[MAPMAX], mt_l[MAPMAX], mt_u[MAPMAX], om_l[MAPMAX], om_u[MAPMAX];
 
  /* END OF ARRAY PARAMETERS */

  /* display parameters */

  /* image dimensions */
  float bxl, bxu, byl, byu, xl, xu, yl, yu, xhl, xhu, yhl, yhu, xpl, xpu;
  /* coordinate dimensions */

  /* display type */
  int display_type;
  int display_marginal;
  int display_marginal_percentiles;

  /* graphical contour levels */
  
  float D2target[D2NCONT];
  float level;
   
  /* graphical auxiliary */
  float fp;
  float xg[2], yg[2];
  float mx[4], my[4];

/* define here colors and line widths */
  int background_color = 0;
  int background_fill = 1;
  int excluded_color = 15;
  int excluded_fill = 1;
  int distribution_color[MAPMAX];
  distribution_color[0] = 1;
  distribution_color[1] = 2;
  int distribution_width[MAPMAX];
  distribution_width[0] = 2;
  distribution_width[1] = 2;
  int sini_color, sini_style, sini_width, mc_color, mc_style, mc_width;
  int h3_color, h3_style, h3_width, h4_color, h4_style, h4_width, xi_color, xi_style, xi_width;
  int omdot_color, omdot_style, omdot_width, gamma_color, gamma_style, gamma_width, pbdot_color, pbdot_style, pbdot_width;
  int xdot_color, xdot_style, xdot_width, r_color, r_style, r_width, mp_color, mp_style, mp_width;
  int nlevels = 3;
  

  /* *************************************************************************************** */
  /* *************************************************************************************** */

  /* Let us make some initializations */

  /* first, the binary parameters */

  printf("\n Reading the scientific parameters.\n");
  
  S = fopen("Science.dat","r");
  
  /* READ THE KEPLERIAN PARAMETERS */
  
  fscanf(S, "%lf %lf %lf", &pb, &ecc, &asini);
  printf("\n For Pb = %lf days, e = %lf, x = %lf lt-s ", pb, ecc, asini);
  /* convert pb to seconds */
  pb = DAY * pb;
  /* orbital frequency */
  nb = 2 * PI / pb;
  /* calculate mass function */
  f = pow(asini, 3) * pow(nb, 2) * 1e+6 / TSUN; 
  printf(" f =  %lf solar masses.\n ", f);
  
  /* READ THE POST-KEPLERIAN PARAMETERS */
  
  fscanf(S, "%d", &komdt);
  if (komdt == 1)
    {
      fscanf(S, "%lf %lf %d %d %d", &omdotobs, &domdotobs, &omdot_color, &omdot_style, &omdot_width);
      /* convert degrees per year to radians per second */
      omdotobs = omdotobs * PI / (180.0 * DAY * YEAR);
      /* do the same for the uncertainty */
      domdotobs = domdotobs * PI / (180.0 * DAY * YEAR);
      /* calculate corresponding masses */
      
      omdotk = (1 - ecc * ecc) * pow ( nb , -five3rd) * pow( ( TSUN*1E-6) , -two3rd) / 3;
     
      mtobs = 0;
      mtobsl = 0;
      mtobsu = 0;

      omdotaux = omdotobs;
      if (omdotaux > 0)
	{
	  mtaux = omdotaux * omdotk;
	  mtobs = pow (mtaux,1.5);
	}

      omdotaux = omdotobs - domdotobs;
      if (omdotaux > 0)
	{
	  mtaux = omdotaux * omdotk;
	  mtobsl = pow (mtaux,1.5);
	}

      omdotaux = omdotobs + domdotobs;
      if (omdotaux > 0)
	{
	  mtaux = omdotaux * omdotk;
	  mtobsu = pow (mtaux,1.5);
	}
    
      printf("\n From that and the rate of advance of periastron, the total mass is %f +%f/-%f solar masses.\n", mtobs, mtobsu-mtobs,mtobs-mtobsl);

      /* calculate minimum companion mass and maximum pulsar mass */
      
      mc_min = pow( (mtobs*mtobs * f), one3rd);
      mp_max = mtobs - mc_min;
      printf("\n Given the mass function and the total mass, the maximum pulsar mass is %f and the minimum companion mass is %f solar masses.\n", mp_max, mc_min);

      /* calculate minimum inclination */

      sinimin = pow ( (f / mtobs), one3rd);
      
      printf("The minimum sin i is %f and corresponding inclination is %f degrees.\n", sinimin, asin(sinimin) * 180 / PI);

    };
  
  fscanf(S, "%d", &kgamma);
  if (kgamma == 1)
    {
      fscanf(S, "%lf %lf %d %d %d", &gammaobs, &dgammaobs, &gamma_color, &gamma_style, &gamma_width);
      /* now calculate part of gamma that depends on keplerian parameters */
      gammak = ecc * pow(nb,-one3rd)*pow(( TSUN*1E-6),two3rd);
      
      /* calculate masses compatible with this gamma and total mass */
      
      mc_gamma = (sqrt( mtobs * mtobs + 4 * pow(mtobs, four3rd) * gammaobs / gammak ) - mtobs) / 2.0;
      
      mp_gamma = mtobs - mc_gamma;

      mc_gammau = (sqrt( mtobs * mtobs + 4 * pow(mtobs, four3rd) * (gammaobs + dgammaobs) / gammak ) - mtobs) / 2.0;
      
      mc_gammal = (sqrt( mtobs * mtobs + 4 * pow(mtobs, four3rd) * (gammaobs - dgammaobs) / gammak ) - mtobs) / 2.0;

      printf("\n The companion mass consistent with this gamma and the total mass is %f +%f/-%f solar masses.\n", mc_gamma, mc_gammau-mc_gamma, mc_gamma-mc_gammal);
      printf("\n The pulsar mass consistent with this gamma and the total mass is %f -%f/+%f solar masses.\n", mp_gamma, mc_gammau-mc_gamma, mc_gamma-mc_gammal);
      
      sini_gamma = pow ((f * mtobs * mtobs),one3rd) / mc_gamma;
      sini_gammau = pow ((f * mtobs * mtobs),one3rd) / mc_gammal;
      sini_gammal = pow ((f * mtobs * mtobs),one3rd) / mc_gammau;
      
      i_gamma = asin(sini_gamma) * 180 / PI;
      i_gammau = asin(sini_gammau) * 180 / PI;
      i_gammal = asin(sini_gammal) * 180 / PI;
      
      printf("\n The sine of the orbital inclination is %f +%f/-%f.\n", sini_gamma, sini_gammau - sini_gamma, sini_gamma-sini_gammal);
      
      printf("\n The orbital inclination is %f +%f/-%f degrees.\n", i_gamma, i_gammau - i_gamma, i_gamma - i_gammal);
      
    }
  fscanf(S, "%d", &kpbdot);
  if (kpbdot == 1)
    {
      fscanf(S, "%lf %lf %d %d %d", &pbdotobs, &dpbdotobs, &pbdot_color, &pbdot_style, &pbdot_width);
      pbdotobs = pbdotobs * 1E-12;
      dpbdotobs = dpbdotobs * 1E-12;
      pbdotk = -192 * PI * pow ((nb * 1e-6 * TSUN), five3rd) / 5.0;
      pbdotk = pbdotk * (1 + ecc * ecc * 73.0 / 24.0 +  pow (ecc , 4) * 37.0 / 96.0 );
      pbdotk = pbdotk * pow((1 - ecc * ecc),-3.5);
    }
  fscanf(S, "%d", &ksini);
  if (ksini == 1) fscanf(S, "%lf %lf %d %d %d", &siniobs, &dsiniobs, &sini_color, &sini_style, &sini_width);
  fscanf(S, "%d", &kmc);
  if (kmc == 1) fscanf(S, "%lf %lf %d %d %d", &mcobs, &dmcobs, &mc_color, &mc_style, &mc_width);
  
  /* OR, USING THE NEW, BETTER FORMULATION FOR THE SHAPIRO DELAY */
  
  fscanf(S, "%d", &kh3);
  if (kh3 == 1) fscanf(S, "%lf %lf %d %d %d", &h3obs, &dh3obs, &h3_color, &h3_style, &h3_width);
  fscanf(S, "%d", &kh4);
  if (kh4 == 1) fscanf(S, "%lf %lf %d %d %d", &h4obs, &dh4obs, &h4_color, &h4_style, &h4_width);
  fscanf(S, "%d", &kxi);
  if (kxi == 1) fscanf(S, "%lf %lf %d %d %d", &xiobs, &dxiobs, &xi_color, &xi_style, &xi_width);
  
  /* READ THE KINETIC PARAMETERS */
  
  fscanf(S, "%d", &kinetix);
  if (kinetix == 1)
    {
      /* read proper motion, its position angle and parallax, no errors for these for the time being */
      fscanf(S, "%lf %lf %lf", &pmobs, &papmobs, &pxobs);
      
      /* convert proper motion from milliarcseconds per year to radians per second */
      pmobs = pmobs * PI / (180.0 * 3600000 * DAY * YEAR);
      
      /* read x-dot and error */
      fscanf(S, "%lf %lf %d %d %d, ", &xdotobs, &dxdotobs, &xdot_color, &xdot_style, &xdot_width);
      xdotobs = xdotobs * 1e-12;
      dxdotobs = dxdotobs * 1e-12;
      
      /* determine maximum cos i */
      
      xdotaux = xdotobs + 2 * dxdotobs;
      incl1 = atan( asini * pmobs /  xdotaux);
      
      xdotaux = xdotobs - 2 * dxdotobs;
      incl2 = atan( asini * pmobs /  xdotaux);
      
      if ((incl1 * incl2) < 0)
	{
	  cosimin = 0;
	}
      else
	if ( cos(incl1) < cos (incl2))
	  {
	    cosimin = cos(incl1);
	    printf("\n From the x-dot and the proper motion we deduce a 99.86 percent");
	    printf("\n probability of i < %f degrees, cos (imax) = %f.\n", incl1 * 180/PI, cosimin);
	  }
	else
	  {
	    cosimin = cos(incl2);
	    printf("\n From the x-dot and the proper motion we deduce a 99.86 percent");
	    printf("\n probability of i < %f degrees, cos (imax) = %f.\n", incl2 * 180/PI, cosimin);
	  }
      
      kinom = 0;
      if (komdt == 1)
	{

	  /* given that there is kinetic information, this might be important to correct for the omega-dot */
	  
	  if (ksini == 1)
	    {
	      omdotaux = omdotobs - pmobs / siniobs;
	      if (omdotaux > 0)
		{
		  mtaux = (omdotaux) * omdotk;
		  mtkinl = pow (mtaux,1.5);
		}
	      omdotaux = omdotobs + pmobs / siniobs;
	      if (omdotaux > 0)
		{
		  mtaux = (omdotaux) * omdotk;
		  mtkinu = pow (mtaux,1.5);
		}
	      kinom = 1;

	      printf("\n Taking the uncertainty due to the proper motion, the total mass is %f -> %f -> %f solar masses.\n", mtkinl, mtobs, mtkinu);
	    }
	  
	  if (kxi == 1)
	    {
	      siniaux = 2 * xiobs / (1 + xiobs * xiobs);
	      omdotaux = omdotobs - pmobs / siniaux;
	      if (omdotaux > 0)
		{
		  mtaux = (omdotaux) * omdotk;
		  mtkinl = pow (mtaux,1.5);
		}
	      omdotaux = omdotobs + pmobs / siniaux;
	      if (omdotaux > 0)
		{
		  mtaux = (omdotaux) * omdotk;
		  mtkinu = pow (mtaux,1.5);
		}

	      kinom = 1;

	      printf("\n Taking the uncertainty due to the proper motion, the total mass is %f -> %f -> %f solar masses.\n", mtkinl, mtobs, mtkinu);

	    }
	 
	}
     
    }
  
  /* READ SPECIAL PARAMETER: MASS RATIO (IF THRERE IS A DETECTION OF COMPANION'S ORBITAL VELOCITY) */
  
  fscanf(S, "%d", &kr);
  if (kr == 1) fscanf(S, "%f %f %d %d %d", &mrobs, &dmrobs, &r_color, &r_style, &r_width);
  
  /* READ DERIVED PARAMETER: PULSAR MASS */
  
  fscanf(S, "%d", &kmp);
  if (kmp == 1) fscanf(S, "%f %f %d %d %d", &mpobs, &dmpobs, &mp_color, &mp_style, &mp_width);
  
  /* Now, read special points, with the first integer indicating the number of points to be represented */
  
  /* This could be anything - values used in simulations, other pulsars, etc.
     Note that now you need to specify the PGPLOT symbol to be used */
  
  fscanf(S,"%d", &ksim);
  
  printf("\n ksim = %d", ksim);
  
  ksimcount = 0;
  
  while (ksimcount < ksim)
    {
      fscanf(S,"%f %f %f", &sm1[ksimcount], &sm2[ksimcount], &ssini[ksimcount]);
      fscanf(S,"%d %d", &ssymbol[ksimcount], &scolor[ksimcount]);
      printf("\n Read special point %d of %d: with m1 = %f, m2 = %f, sin i = %f, symbol type = %d and symbol color = %d", ksimcount, ksim, sm1[ksimcount], sm2[ksimcount], ssini[ksimcount], ssymbol[ksimcount], scolor[ksimcount]);
      sf[ksimcount] = pow((sm2[ksimcount]*ssini[ksimcount]), 3)/pow((sm1[ksimcount]+sm2[ksimcount]), 2);
      ksimcount = ksimcount + 1;
    }
  
  fclose(S);
  
  printf("\n Finished reading the measured parameters.\n");

  /* *************************************************************************************** */
  /* *************************************************************************************** */
  /* *************************************************************************************** */

  map_number = 0;
 
  /* Now going to ask whether we load any chi2 maps */

  printf("\n Do you want to make diagrams with just these parameters (0)");
  printf("\n or do you want to add chi2 maps (1 or 2)?\n");
  
  scanf("%d", &map_number);

  if (map_number > MAPMAX)
    {
      printf("\n You introduced a number larger than %d.", MAPMAX);
      printf("\n The program will read %d maps only.", MAPMAX);
      map_number = MAPMAX;
    }

  /* *************************************************************************************** */
  /* *************************************************************************************** */

  /* Now, start reading the array headers - IF NEEDED */

  /* *************************************************************************************** */

  /* Reading first header */

  if (map_number > 0)
    {
      map_counter = 0;

      printf("\n We are now going to read the headers of the arrays you specified.");
      printf("\n The program will automatically marginalize the pdfs along the relevant axes, and output them.");
      printf("\n Later, you will be given the option of making a display of them.");

      /* Define contour level targets for 2-D distributions */
      D2target[0] = 68.2689; /* one sigma */
      D2target[1] = 95.4500; /* two sigma */
      D2target[2] = 99.7300; /* three sigma */
 
      /* first thing: define 1-D probability targets */
      /* D1target[0] = 0.0; lower edge of 3-sigma */
      D1target[0] = 0.001349898; /* lower edge of 3-sigma */
      D1target[1] = 0.022750132; /* lower edge of 2-sigma */
      D1target[2] = 0.158655254; /* lower edge of 1-sigma */
      D1target[3] = 0.5;         /* median */
      jmedian = 3;
      D1target[4] = 0.841344746; /* upper edge of 1-sigma */
      D1target[5] = 0.977249868; /* upper edge of 2-sigma */
      D1target[6] = 0.998650101; /* upper edge of 3-sigma */
      /* D1target[6] = 1.0; /* upper edge of 3-sigma */

      while (map_counter < map_number)
	{

	  /* initialize the logic */
	  map_type[map_counter] = 0;
	  calculate_cosim2[map_counter] = 0;
	  calculate_m1m2[map_counter] = 0;
	  calculate_h3h4[map_counter] = 0;
	  calculate_cosiom[map_counter] = 0;

	  /* dumb thing - open the correct file */
	  if (map_counter == 0) Fh = fopen("F1.hdr","r");
	  if (map_counter == 1) Fh = fopen("F2.hdr","r");
      
	  printf("\n Now loading parameters for array n. %d.\n", map_counter + 1);
  
	  /* First, most important thing: knowing what sort of file this is */
      
	  fscanf(Fh, "%d", &map_type[map_counter]);

	  /*  ********************************* */
	  
	  /* if this is a Monte Carlo cosi - m2 file, read the array dimensions for the masses and inclination */
	  if (map_type[map_counter] == 0)
	    {
	      
	      /* just scan the total number of elements */
	      fscanf(Fh, "%d", &lmax[map_counter]);
	      	      
	    }
	  
	  /*  ********************************* */
	  
	  /* normal cosi - m2 plane map */
	  if (map_type[map_counter] == 1)
	    {
	      
	      /* if this is the case, we want at least calculate the cosi - m2 plane */
	      calculate_cosim2[map_counter] = 1;
	      
	      /* Read the number of elements, start value and element size for the x dimension (cos i) */
	      fscanf(Fh, "%d", &n_cosi[map_counter]);
	      fscanf(Fh, "%f", &i_cosi[map_counter]); tr_cosim2[map_counter][0] = i_cosi[map_counter];
	      tr_cosim2[map_counter][1] = 0;
	      fscanf(Fh, "%f", &s_cosi[map_counter]); tr_cosim2[map_counter][2] = s_cosi[map_counter];
	      /* variable useful for pulsar mass mapping */
	      dcosi[map_counter] = s_cosi[map_counter]/(NCOSI * 1.0);
	      
	      /* Read the number of elements and element size for the Y dimension (m2) */
	      fscanf(Fh, "%d", &n_m2[map_counter]);
	      fscanf(Fh, "%f", &i_m2[map_counter]); tr_cosim2[map_counter][3] = i_m2[map_counter];
	      fscanf(Fh, "%f", &s_m2[map_counter]); tr_cosim2[map_counter][4] = s_m2[map_counter];
	      tr_cosim2[map_counter][5] = 0;
	      /* variable useful for pulsar mass mapping */
	      dm2[map_counter] = s_m2[map_counter]/(NM * 1.0);
	    }
	  
	  /*  ********************************* */
	  
	  /* if this is an orthogonal file, read the array dimensions for h3 */
	  if (map_type[map_counter] == 2)
	    {
	      
	      /* in this case, we want at least display the h3 - h4 plane */
	      calculate_h3h4[map_counter] = 1;
	      
	      /* specify the dimensions for the h3 dimension */
	      fscanf(Fh, "%d", &n_h[map_counter]);
	      fscanf(Fh, "%f", &i_h[map_counter]); tr_h3h4[map_counter][0] = i_h[map_counter];
	      tr_h3h4[map_counter][1] = 0;
	      fscanf(Fh, "%f", &s_h[map_counter]); tr_h3h4[map_counter][2] = s_h[map_counter];
	      
	      /* specify the dimensions for the h4 array - the same as for the h3 array */
	      tr_h3h4[map_counter][3] = tr_h3h4[map_counter][0]; /* starting value is the same as for h3 */
	      tr_h3h4[map_counter][4] = tr_h3h4[map_counter][2];
	      tr_h3h4[map_counter][5] = 0;
	    }
	    
	   
 	  /*  ********************************* */

	  /* normal cosi - h3 plane map */
	  if (map_type[map_counter] == 3)
	    {
	      
	      /* if this is the case, we want at least calculate the cosi - h3 plane */
	      calculate_cosih3[map_counter] = 1;
	      /* Except that at the moment I don't feel like displaying this plane, and go directly to the cosi - m2 plane */
	      
	      /* Read the number of elements, start value and element size for the x dimension (cos i) */
	      fscanf(Fh, "%d", &n_cosi[map_counter]);
	      fscanf(Fh, "%f", &i_cosi[map_counter]); tr_cosih3[map_counter][0] = i_cosi[map_counter];
	      tr_cosih3[map_counter][1] = 0;
	      fscanf(Fh, "%f", &s_cosi[map_counter]); tr_cosih3[map_counter][2] = s_cosi[map_counter];
	      /* variable useful for pulsar mass mapping */
	      dcosi[map_counter] = s_cosi[map_counter]/(NCOSI * 1.0);
	      
	      /* Read the number of elements and element size for the Y dimension (h3) */
	      fscanf(Fh, "%d", &n_h3[map_counter]);
	      fscanf(Fh, "%f", &i_h3[map_counter]); tr_cosih3[map_counter][3] = i_h3[map_counter];
	      fscanf(Fh, "%f", &s_h3[map_counter]); tr_cosih3[map_counter][4] = s_h3[map_counter];
	      tr_cosih3[map_counter][5] = 0;
	    }
	  
	  
	  /*  ********************************* */
	  
	  /* cosi - omega plane map (FOR SYSTEMS WITH VERY GOOD OMEGA-DOT BUT NO OTHER STRONG CONSTRAINTS) */
	  if (map_type[map_counter] == 5)
	    {
	      if (komdt == 1)
		{
		  
		  /* if this is the case, we want at least calculate the cosi - omega plane */
		  calculate_cosiom[map_counter] = 1;
		  calculate_m1m2[map_counter] = 1;
		  
		  /* Read the number of elements and element size for the x dimension (cos i) */
		  fscanf(Fh, "%d", &n_cosi[map_counter]);
		  fscanf(Fh, "%f", &i_cosi[map_counter]); tr_cosiom[map_counter][0] = i_cosi[map_counter];
		  tr_cosiom[map_counter][1] = 0;
		  fscanf(Fh, "%f", &s_cosi[map_counter]); tr_cosiom[map_counter][2] = s_cosi[map_counter];
		  /* variable useful for pulsar mass mapping */
		  dcosi[map_counter] = s_cosi[map_counter]/(NCOSI * 1.0);
		  
		  /* Read the number of elements and element size for the Y dimension (omega) */
		  fscanf(Fh, "%d", &n_om[map_counter]);
		  fscanf(Fh, "%f", &i_om[map_counter]); tr_cosiom[map_counter][3] = i_om[map_counter];
		  fscanf(Fh, "%f", &s_om[map_counter]); tr_cosiom[map_counter][4] = s_om[map_counter];
		  tr_cosiom[map_counter][5] = 0;
		  
		  /* in this case I've decided to make an automatic calculation of the pulsar (and companion) mass vector */
		  
		  /* these are independent of the distribution number */
		  mc_min = pow( (mtobs*mtobs * f), one3rd);
		  mp_max = mtobs - mc_min;
		  /*display variables also set automatically */
		  xpu = mp_max;
		  xpl = 0;
		  /* mass binning */
		  m1_l[map_counter] = xpl;
		  m1_u[map_counter] = xpu;
		  
		  i_m2[map_counter] = mtobs;
		  s_m2[map_counter] = - mp_max / (MPBIN * 1.0);
		  n_m2[map_counter] = MPBIN;
		}
	      else
		{
		  printf("\n There is no omega-dot in your input file (Science.dat). This sort of map does not make sense in that case.\n");
		}
	    }
	  
	  /*  ********************************* */		  
	  
	  /* Now, we have a special case, a map with the position of the pulsar */
	  
	  if (map_type[map_counter] == 10)
	    {
	      
	      /* if this is the case, we want at least calculate the cosi - omega plane */
	      calculate_radec[map_counter] = 1;
	      	      	  
	      /* In this case nothing needs to be done because the array dimensions are set from the data, however, the header should have
		 a number of elements*/
	      
	      /* just scan the total number of elements, the  */
	      fscanf(Fh, "%d %f %f %d %f %f %d", &lmax[map_counter], &ra_min, &ra_diff_min, &n_ra[map_counter], &dec_min, &dec_diff_min, &n_dec[map_counter]);
	  
	      /* let's print what we have found */
	      
	      printf("\n For RA, the minimum is %f, the step is %f, and the number of elements is %d.\n", ra_min, ra_diff_min, n_ra[map_counter]);
	      printf("\n For Dec, the minimum is %f, the step is %f, and the number of elements is %d.\n", dec_min, dec_diff_min, n_dec[map_counter]);
	      	      
	      /* This imaging stuff is not really needed now */
	      
	      /* establish approximate dimensions of array */
	      i_ra[map_counter] = ra_min;
	      s_ra[map_counter] = ra_diff_min;
		  		  
	      i_dec[map_counter] = dec_min;
	      s_dec[map_counter] = dec_diff_min;
	  
	    }
	  
	  
	  /* Now, we have the case of 3-D maps. These are based on cos i and Omega. The third quantity could be either the Total mass or h3 */
	  if ((map_type[map_counter] == 11)||(map_type[map_counter] == 12))
	    {
	      
	      /* if this is the case, we want at least calculate the cosi - omega plane */
	      calculate_cosiom[map_counter] = 1;
	      
	      
	      /* In this case nothing needs to be done because the array dimensions are set from the data, however, the header should have
		 a number of elements*/
	      
	      
	      /* just scan the total number of elements */
	      fscanf(Fh, "%d", &lmax[map_counter]);
	      
	      
	      /* set the minimum at a very high level */
	      vmin = 1e9;
	      
	      om_min = 360;
	      cosi_min = 2;
	      mt_min = 10;
	      h3_min = 1e9;
	      
	      om_diff_min = 1e9;
	      cosi_diff_min = 1e9;
	      mt_diff_min = 1e9;
	      h3_diff_min = 1e9;
	      
	      om_max = -360;
	      cosi_max = -2;
	      mt_max = 0;
	      h3_max = 0;
	      
	      omaux_before = -1e-9;
	      cosiaux_before = -1e9;
	      mtaux_before = -1e-9;
	      h3aux_before = -1e9;
	      
	      /* we can use this opportunity to do an early load of the super array. Let's open that file */
	      
	      if (map_counter == 0) F = fopen("F1.dat","r");
	      if (map_counter == 1) F = fopen("F2.dat","r");
	      
	      
	      /* cycling through the cos i coordinate */
	      for (i = 0 ; i < lmax[map_counter]; i++)
		{
		  
		  /* read all important variables */
		  if (map_type[map_counter] == 11)
		    {
		      fscanf(F, "%e %e %e %e", &omaux, &cosiaux, &mtaux, &v);		      
		      
		      if (mtaux < mt_min) mt_min = mtaux;
		      if (mtaux > mt_max) mt_max = mtaux;
		      
		      if ( ( fabs(mtaux - mtaux_before) < mt_diff_min) && ( fabs(mtaux - mtaux_before) > 1e-6)  )
			mt_diff_min = fabs(mtaux - mtaux_before);
		      
		      mtaux_before = mtaux;
		    }
		  
		  /* read all important variables */
		  if (map_type[map_counter] == 12)
		    {
		       
		      fscanf(F, "%e %e %e %e", &omaux, &cosiaux, &h3aux, &v);
		      
		      if (h3aux < h3_min) h3_min = h3aux;
		      if (h3aux > h3_max) h3_max = h3aux;
		      
		      if ( (fabs((h3aux - h3aux_before)) < h3_diff_min) && (fabs(h3aux - h3aux_before) > 1e-6)  )
			h3_diff_min = fabs(h3aux - h3aux_before);
		      
		      h3aux_before = h3aux;
		    }
		  	  
		  if (omaux < om_min) om_min = omaux;
		  if (omaux > om_max) om_max = omaux;
		  
		  /* printf("\n %f", om_max); */
		  
		   
		  if ( (fabs(omaux - omaux_before) < om_diff_min) && (fabs(omaux - omaux_before) > 1e-6)  )
		    om_diff_min = fabs(omaux - omaux_before);
		  
		  omaux_before = omaux;
		  
		  
		  if (cosiaux < cosi_min) cosi_min = cosiaux;
		  if (cosiaux > cosi_max) cosi_max = cosiaux;
		  
		  if ( (fabs(cosiaux - cosiaux_before) < cosi_diff_min) && (fabs(cosiaux - cosiaux_before) > 1e-6)  )
		    cosi_diff_min = fabs(cosiaux - cosiaux_before);
		  
		  cosiaux_before = cosiaux;
		  
		  /* update minima and maxima for both sorts of  */
		  
		  if (v < vmin) vmin = v; /* update minimum */
		  
		}
	      
	      fclose(F);
	      
	      /* let's print what we have found */
	      
	      printf("\n For Omega, the minimum is %f, the maximum is %f and the step is %f.\n", om_min, om_max, om_diff_min);
	      printf("\n For Cos i, the minimum is %f, the maximum is %f and the step is %f.\n", cosi_min, cosi_max, cosi_diff_min);
	      
	      if (map_type[map_counter] == 11)
		printf("\n For Mt, the minimum is %f, the maximum is %f and the step is %f.\n", mt_min, mt_max, mt_diff_min);
	      
	      if (map_type[map_counter] == 12)
		printf("\n For h3, the minimum is %f, the maximum is %f and the step is %f.\n", h3_min, h3_max, h3_diff_min);
	      
	      
	      /* establish approximate dimensions of array */
	      n_cosi[map_counter] = floor( (cosi_max - cosi_min) / cosi_diff_min) + 1;
	      i_cosi[map_counter] = cosi_min;
	      s_cosi[map_counter] = cosi_diff_min;
		  
	      /* establish display matrix  */
	      
	      tr_cosiom[map_counter][0] = i_cosi[map_counter];
	      tr_cosiom[map_counter][1] = 0;
	      tr_cosiom[map_counter][2] = s_cosi[map_counter];
	      
	      /* variables useful for smoothing */
	      dcosi[map_counter] = s_cosi[map_counter]/(NCOSI * 1.0);

              if (map_type[map_counter] == 11)
		dmt[map_counter] = mt_diff_min/(NM * 1.0);
	      if (map_type[map_counter] == 12)
		dh3[map_counter] = h3_diff_min/(NM * 1.0);
	      
	      n_om[map_counter] = floor( (om_max - om_min) / om_diff_min) + 1;
	      i_om[map_counter] = om_min;
	      s_om[map_counter] = om_diff_min;
	      
	      
	      /* define the matrix, Y dimension (omega) */
	      tr_cosiom[map_counter][3] = i_om[map_counter];
	      tr_cosiom[map_counter][4] = s_om[map_counter];
	      tr_cosiom[map_counter][5] = 0;

	      if (map_type[map_counter] == 11)
		{
		  n_mt[map_counter] = floor( (mt_max - mt_min) / mt_diff_min) + 1;
		  i_mt[map_counter] = mt_min;
		  /* s_mt[map_counter] = (mt_max - mt_min) / n_mt[map_counter]; */
		  s_mt[map_counter] = mt_diff_min;

		}
	    }
	  	  
	  /* ******* CALCULATE MC - COSI diagram? ******* */

	  if ((map_type[map_counter] == 0)||(map_type[map_counter] == 2))
	    {
	      printf("\n Do you want to make a cos i - m2 map? (1 = Yes, 0 = No)\n");
	      scanf("%d", &calculate_cosim2[map_counter]);

	      /* in this case we want at least a mass-cos i plot */
	      if (calculate_cosim2[map_counter] == 1)
		{
		  /* get dimensions for m2 - cosi array */
		  printf("\n How many elements in the cosi array? (here you should know roughly the range of inclinations involved)\n");
		  scanf("%d", &n_cosi[map_counter]);
		  printf("\n What is the initial cosi?\n");
		  scanf("%f", &i_cosi[map_counter]); tr_cosim2[map_counter][0] = i_cosi[map_counter];
		  tr_cosim2[map_counter][1] = 0;
		  printf("\n What is the pixel width?\n");
		  scanf("%f", &s_cosi[map_counter]); tr_cosim2[map_counter][2] = s_cosi[map_counter];
		  
		  /* Read the number of elements and size for the Y dimension (m2) */
		  printf("\n How many elements in the m2 array? (here you should know roughly the range of companion masses involved)\n");
		  scanf("%d", &n_m2[map_counter]);
		  printf("\n What is the starting m2?\n");
		  scanf("%f", &i_m2[map_counter]); tr_cosim2[map_counter][3] = i_m2[map_counter];
		  printf("\n What is the pixel width?\n");
		  scanf("%f", &s_m2[map_counter]); tr_cosim2[map_counter][4] = s_m2[map_counter];
		  tr_cosim2[map_counter][5] = 0;

		  /* values for mass binning */
		  m2_l[map_counter] = i_m2[map_counter];
		  m2_u[map_counter] = m2_l[map_counter] + n_m2[map_counter] * s_m2[map_counter];

		  cosi_l[map_counter] = i_cosi[map_counter];
		  cosi_u[map_counter] = cosi_l[map_counter] + n_cosi[map_counter] * s_cosi[map_counter];
		  
		}
	    }


	  if ((map_type[map_counter] == 3)||(map_type[map_counter] == 11)||(map_type[map_counter] == 12))
	    {
	      printf("\n Do you want to make a cos i - m2 map? (1 = Yes, 0 = No)\n");
	      scanf("%d", &calculate_cosim2[map_counter]);

	      /* in this case we want at least a mass-cos i plot */
	      if (calculate_cosim2[map_counter] == 1)
		{

		  /* get dimensions for m2 - cosi array - these are already specified by known dimensions of this array */
		  tr_cosim2[map_counter][0] = i_cosi[map_counter];
		  tr_cosim2[map_counter][1] = 0;
		  tr_cosim2[map_counter][2] = s_cosi[map_counter];
		  
		  /* Read the number of elements and size for the Y dimension (m2) */
		  printf("\n How many elements in the m2 array? (here you should know roughly the range of companion masses involved)\n");
		  scanf("%d", &n_m2[map_counter]);
		  printf("\n What is the starting m2?\n");
		  scanf("%f", &i_m2[map_counter]); tr_cosim2[map_counter][3] = i_m2[map_counter];
		  printf("\n What is the pixel width?\n");
		  scanf("%f", &s_m2[map_counter]); tr_cosim2[map_counter][4] = s_m2[map_counter];
		  tr_cosim2[map_counter][5] = 0;

		  /* values for mass binning */
		  m2_l[map_counter] = i_m2[map_counter];
		  m2_u[map_counter] = m2_l[map_counter] + n_m2[map_counter] * s_m2[map_counter];

		}
	    }

	  /*  ******* CALCULATE PSR MASSES? ******* */
	  
	  /* for any of the previously loaded types of maps, we could decide to calculate distributions for the pulsar mass */

	  if ((map_type[map_counter] == 0)||(map_type[map_counter] == 1)||(map_type[map_counter] == 2)||(map_type[map_counter] == 3)||(map_type[map_counter] == 11)||(map_type[map_counter] == 12))
	    { 
	      printf("\n Do you want to make a mass-mass map?\n");
	      
	      scanf("%d", &calculate_m1m2[map_counter]);
	      
	      if (calculate_m1m2[map_counter] == 1)
		{
		  printf("\n What are the lower and upper limits for the pulsar mass vector (calculations only)? \n");
		  scanf("%f %f", &m1_l[map_counter], &m1_u[map_counter]);
		  
		  /* define print vector for mass-mass diagram */
		  tr_m1m2[map_counter][0] = m1_l[map_counter];
		  tr_m1m2[map_counter][1] = 0;
		  tr_m1m2[map_counter][2] = (m1_u[map_counter] - m1_l[map_counter]) / (MPBIN * 1.0);
		  
		  /* if m2 array has already been defined, then vertical part is almost the same, since m-m has mc in the vertical axis.
		     The only difference is that we run through the x-axis first */
		  if (calculate_cosim2[map_counter] == 1)
		    {
		      tr_m1m2[map_counter][3] = tr_cosim2[map_counter][3];
		      tr_m1m2[map_counter][4] = tr_cosim2[map_counter][4];
		      tr_m1m2[map_counter][5] = tr_cosim2[map_counter][5];
		    }
		  else
		    {
		      /* Read the number of elements and size for the Y dimension (m2) */
		      printf("\n How many elements in the m2 array? (here you should know roughly the range of companion masses involved)\n");

		      scanf("%d", &n_m2[map_counter]);
		      printf("\n What is the starting m2?\n");
		      scanf("%f", &i_m2[map_counter]); tr_m1m2[map_counter][3] = i_m2[map_counter];
		      printf("\n What is the pixel width?\n");
		      scanf("%f", &s_m2[map_counter]); tr_m1m2[map_counter][4] = s_m2[map_counter];
		      tr_m1m2[map_counter][5] = 0;
		      
		      /* values for mass binning */
		      /* m1 limits already defined above */
		      m2_l[map_counter] = i_m2[map_counter];
		      m2_u[map_counter] = m2_l[map_counter] + i_m2[map_counter] * s_m2[map_counter];
		    }
		}
	    }

	  /* if distribution is of type 11, we'll calculate the distribution of total masses */

	  if (map_type[map_counter] == 11)
	    {
	      mt_l[map_counter] = i_mt[map_counter];
	      mt_u[map_counter] = mt_l[map_counter] + n_mt[map_counter] * s_mt[map_counter];
	    }


	  /*  ********************************* */
	  
	  /* we're done with the header. */

	  fclose(Fh);

	  /* now, let's initialize the array variables */

	  if ((calculate_cosim2[map_counter] == 1)||(calculate_cosiom[map_counter] == 1))
	    {
	      /* in these cases we need cosine arrays */
	      for (i = 0 ; i < n_cosi[map_counter]; i++)
		{
		  cosi[map_counter][i] = 0; /* load value */
		  cosix[map_counter][i] = i_cosi[map_counter] + i * s_cosi[map_counter]; /* load x-vector */
		}
	    }
	    
	  if ((calculate_cosih3[map_counter] == 1))
	    {
	      /* in these cases we need h3ne arrays */
	      for (i = 0 ; i < n_h3[map_counter]; i++)
		{
		  h3[map_counter][i] = 0; /* load value */
		  h3x[map_counter][i] = i_h3[map_counter] + i * s_h3[map_counter]; /* load x-vector */
		}
	    }  


	  if ((calculate_m1m2[map_counter] == 1)||(calculate_cosim2[map_counter] == 1))
	    {
	      /* in these cases we need m2 arrays */
	      for (i = 0 ; i < n_m2[map_counter]; i++)
		{
		  m2[map_counter][i] = 0; /* load value */
		  m2x[map_counter][i] = i_m2[map_counter] + i * s_m2[map_counter]; /* load x-vector */
		}
	    }

	  if ((calculate_m1m2[map_counter] == 1)||(calculate_cosiom[map_counter] == 1))
	    {
	      /* in these cases we need m1 arrays */
	      for (i = 0 ; i < MPBIN; i++)
		{
		  m1[map_counter][i] = 0; /* load value */
		  m1x[map_counter][i] = m1_l[map_counter] + i *  (m1_u[map_counter] - m1_l[map_counter]) / (MPBIN * 1.0); /* load x-vector */
		}
	    }

	  if (calculate_cosiom[map_counter] == 1)
	    {
	      /* in these cases we need omega arrays */
	      for (i = 0 ; i < n_om[map_counter]; i++)
		{
		  om[map_counter][i] = 0; /* load value */
		  omx[map_counter][i] = i_om[map_counter] + i * s_om[map_counter]; /* load x-vector */
		}
	    }

	  /* exception here, tied to the case of having an input array of type 11 */
	  if (map_type[map_counter] == 11)
	    {
	      /* in these cases we need mt arrays */
	      for (i = 0 ; i < n_mt[map_counter]; i++)
		{
		  mt[map_counter][i] = 0; /* load value */
		  mtx[map_counter][i] = i_mt[map_counter] + i * s_mt[map_counter]; /* load x-vector */
		}

	    }

	  if (calculate_radec[map_counter] == 1)
	    {
	      /* in these cases we need ra and dec arrays */
	      for (i = 0 ; i < n_ra[map_counter]; i++)
		{
		   rax[map_counter][i] = i_ra[map_counter] + i * s_ra[map_counter];
		};
		
              for (i = 0 ; i < n_dec[map_counter]; i++)
	        {		  
		   decx[map_counter][i] = i_dec[map_counter] + i * s_dec[map_counter];   
		};
		
		
	    }
	    

	  /* initialize 2-D arrays */

	  /* (it does not really matter whether we will use these or not */
      
	  for (l = 0 ; l < AMAX ; l++)
	    {
	      cosim2[map_counter][l] = 0;    
	      m1m2[map_counter][l] = 0;
	      h3h4[map_counter][l] = 0;
	      cosiom[map_counter][l] = 0;
	      radec[map_counter][l] = 0;
	    }

	  /* initialize total probability counters and maxima */
	  cosim2_total[map_counter] = 0.0;
	  m1m2_total[map_counter] = 0.0;
          m1om_total[map_counter] = 0.0;
	  cosiom_total[map_counter] = 0.0;
	  h3h4_total[map_counter] = 0.0;
	  radec_total[map_counter] = 0.0;

	  cosim2max[map_counter] = 0;
	  m1m2max[map_counter] = 0;
	  cosiommax[map_counter] = 0;
	  h3h4max[map_counter] = 0;
	  radecmax[map_counter] = 0;
	  
	  cosimax[map_counter] = 0;
	  ommax[map_counter] = 0;
	  m2max[map_counter] = 0;
	  m1max[map_counter] = 0;
	  mtmax[map_counter] = 0;
	  h3max[map_counter] = 0;
	  h4max[map_counter] = 0;
	  
 
	  /* *************************************************************************************** */
	  
	  /* LOAD ARRAY LOOP  */

	  /* cosi - m2 plane map */

	  if (map_type[map_counter] == 1)
	    {
	      /* FIRST: read the map */

	      vmin = 1E+9;
	      
	      /* dumb thing: open the correct file */
	      if (map_counter == 0) F = fopen("F1.dat","r");
	      if (map_counter == 1) F = fopen("F2.dat","r");
 	      
	      /* cycling through the cos i coordinate */
	      for (i = 0 ; i < n_cosi[map_counter]; i++)
		{
		  /* cycling through the omega coordinate */
		  for (j = 0 ; j < n_m2[map_counter]; j++)
		    {
		      /* read chi-2 */
		      fscanf(F, "%e", &v);
		      if (v < vmin) vmin = v; /* update minimum */
		    }
		}
	      fclose(F);

	      /* Second: Calculate the probabilities */

	      printf("\n Read the chi2 map. Now will calculate the probabilities. \n");
	  
	      /* open file again, since we're going to need the values again */
	      if (map_counter == 0) F = fopen("F1.dat","r");
	      if (map_counter == 1) F = fopen("F2.dat","r");

	      /* cycling through the cos i coordinate */
	      for (i = 0 ; i < n_cosi[map_counter]; i++)
		{

		  /* cycling through the y coordinate */
		  for (j = 0 ; j < n_m2[map_counter]; j++)
		    {

		      /* have to read array element here, not in if statement that follows,
			 otherwise there are gaps in the reading */
		      fscanf(F, "%e", &v);
		    
		      /* calculate contribution to 2-D pdf (for cosi - m2 plane) */
		      vaux = exp ((vmin - v)/2);
		      /* printf("\n i = %d, j = %d, v = %f, vmin = %f, vaux = %f.\n", i, j, v, vmin, vaux);
			 printf("\n cosi = %f, sini = %f, m2 = %f, m1 = %f.\n", cosiaux, siniaux, m2aux, m1aux); */
		      
		      /* update counters */
		      
		      /* make a condition for this: whether the pulsar mass is within a NS range */
		      
		     
		      cosiaux = cosix[map_counter][i];
		      siniaux = sqrt(1 -  cosiaux * cosiaux);
		      m2aux = m2x[map_counter][j];
		      m1aux = sqrt( pow( (m2aux*siniaux), 3 ) / (1.0 * f ) ) - m2aux;

		      if( (m1aux > 1.17) && (m1aux < 3.2))
		      	{
			  cosim2[map_counter][i*n_m2[map_counter] + j] = vaux; /* calculating it this way, not as above, because each point should load _once_ */
			  if (cosim2[map_counter][i*n_m2[map_counter] + j] > cosim2max[map_counter])
			    cosim2max[map_counter] = cosim2[map_counter][i*n_m2[map_counter] + j]; /* update maximum of 2-D pdf */
		      
			  cosi[map_counter][i] = cosi[map_counter][i] + vaux; /*update counters for cosi side distribution */
			  if (cosi[map_counter][i] > cosimax[map_counter])
			    cosimax[map_counter] = cosi[map_counter][i]; /* update maximum of 1-D pdf*/
		      
			  m2[map_counter][j] = m2[map_counter][j] + vaux; /*update counters for m2 side distribution */
			  if (m2[map_counter][j] > m2max[map_counter])
			    m2max[map_counter] = m2[map_counter][j]; /* update maximum of 1-D pdf */
			  
			  cosim2_total[map_counter] = cosim2_total[map_counter] + vaux; /* total distribution */
			}
			  
		      /* update pulsar mass distribution - we will need a finer cosi grid, particularly for the low inclinations where steps in pulsar mass are huge*/


		      for (imp = -HNCOSI; imp < HNCOSI; imp++)
			{
			  /* vaux = cosim2[map_counter][i*n_m2[map_counter] + j] + (1.0 * imp)/(1.0 * NCOSI) * dv; */
			  cosiaux = cosix[map_counter][i] + dcosi[map_counter] * imp;
			  siniaux = sqrt(1 -  cosiaux * cosiaux);
			  
			  for (imt = -HNM; imt < HNM; imt++)
			    {
			      /* now: estimate masses */
			      m2aux = m2x[map_counter][j] + dm2[map_counter] * imt;
			      
			      m1aux = sqrt( pow( (m2aux*siniaux), 3 ) / (1.0 * f ) ) - m2aux;

			      if ((m1aux > 1.17)&&(m1aux < 3.2))
				{
				  m = floor( MPBIN * (m1aux - m1_l[map_counter])/(m1_u[map_counter] - m1_l[map_counter]) );
				  /* if (( MPBIN * (m1aux - m1_l[map_counter])/(m1_u[map_counter] - m1_l[map_counter]) - m) > 0.5) m = m+1; */
				  
				  if ((m > -1) && (m < MPBIN))
				    {
				      /* printf("\n i = %d, j = %d, m = %d, v = %f, vmin = %f, vaux = %f.\n", i, j, m, v, vmin, vaux);
					 printf("\n cosi = %f, sini = %f, m2 = %f, m1 = %f.\n", cosiaux, siniaux, m2aux, m1aux); */
				      
				      /* update pulsar mass array */
				      m1m2[map_counter][m * n_m2[map_counter] + j] = m1m2[map_counter][m * n_m2[map_counter] + j] + vaux;
				      /* update its maximum */
				      if (m1m2[map_counter][m * n_m2[map_counter] + j] >  m1m2max[map_counter])
					m1m2max[map_counter] = m1m2[map_counter][m * n_m2[map_counter] + j];
				      
				      /* update lateral arrays */
				      m1[map_counter][m] = m1[map_counter][m] + vaux;
				      if ( m1[map_counter][m] > m1max[map_counter] ) m1max[map_counter] = m1[map_counter][m];
				      
				      /* update counter of total pulsar mass probability */
				      m1m2_total[map_counter] = m1m2_total[map_counter] + vaux;
				    }
				  
				} /* close negatime mass condition */
			      
			    } /* close imt loop */

			} /* close imp loop */
		      
		    } /* close j loop */
		  
		} /* close i loop */
	      	      
	      fclose(F);
	    }



	  /* *************************************************************************************** */
	  
	  /* LOAD ARRAY LOOP  */

	  /* cosi - h3 plane map */

	  if (map_type[map_counter] == 3)
	    {
	      /* FIRST: read the map */

	      vmin = 1E+9;
	      
	      /* dumb thing: open the correct file */
	      if (map_counter == 0) F = fopen("F1.dat","r");
	      if (map_counter == 1) F = fopen("F2.dat","r");
 	      
	      /* cycling through the cos i coordinate */
	      for (i = 0 ; i < n_cosi[map_counter]; i++)
		{
		  /* cycling through the omega coordinate */
		  for (j = 0 ; j < n_h3[map_counter]; j++)
		    {
		      /* read chi-2 */
		      fscanf(F, "%e", &v);
		      if (v < vmin) vmin = v; /* update minimum */
		    }
		}
	      fclose(F);

	      /* Second: Calculate the probabilities */
	      
	      printf(" n_m2 = %d, m2_low = %f, m2_high = %f. \n", n_m2[map_counter], m2_l[map_counter], m2_u[map_counter]);

	      printf("\n Read the chi2 map. Now will calculate the probabilities. \n");
	  
	      /* open file again, since we're going to need the values again */
	      if (map_counter == 0) F = fopen("F1.dat","r");
	      if (map_counter == 1) F = fopen("F2.dat","r");

	      /* cycling through the cos i coordinate */
	      for (i = 0 ; i < n_cosi[map_counter]; i++)
		{

		  cosiaux = cosix[map_counter][i];
		  siniaux = sqrt(1 - cosiaux * cosiaux);
		  xiaux = siniaux / (1 + cosiaux);
		  
		  /* cycling through the y coordinate */
		  for (j = 0 ; j < n_h3[map_counter]; j++)
		    {
		      /* now: estimate masses */
		      
		      h3aux = h3x[map_counter][j];
		      m2aux = h3aux / (TSUN * pow(xiaux,3)); 
		      
		      m1aux = sqrt( pow( (m2aux*siniaux), 3 ) / (1.0 * f ) ) - m2aux;

		      /* have to read array element here, not in if statement that follows,
			 otherwise there are gaps in the reading */
		      fscanf(F, "%e", &v);
		      
		      /* If pulsar mass is larger than zero, update counters */
		      
		      if (m1aux > 0)
			{

			  /* calculate contribution to 2-D pdf (for cosi - m2 plane) */
			  vaux = exp ((vmin - v)/2);
			  /* printf("\n i = %d, j = %d, v = %f, vmin = %f, vaux = %f.\n", i, j, v, vmin, vaux);
			     printf("\n cosi = %f, sini = %f, m2 = %f, m1 = %f.\n", cosiaux, siniaux, m2aux, m1aux); */
			  

			  /* update pulsar mass distribution - we will need a finer cosi grid, particularly for the low inclinations where steps in pulsar mass are huge*/

			  /*start preparing next loop */
			  if (i > 1)
			    {
			      dv = cosim2[map_counter][i*n_m2[map_counter] + j] - cosim2[map_counter][(i-1)*n_m2[map_counter] + j];
			      
			    }			  
			  else
			    {
			      dv = 0;
			    }

			  for (imp = -HNCOSI; imp < HNCOSI; imp++)
			    {
			      /*vaux = cosim2[map_counter][i*n_m2[map_counter] + j] + (1.0 * imp)/(1.0 * NCOSI) * dv; */
			      cosiaux = cosix[map_counter][i] + dcosi[map_counter] * imp;
			      siniaux = sqrt(1 -  cosiaux * cosiaux);
			      xiaux = siniaux / (1 + cosiaux);
			      /* now: estimate masses */

         		      m2aux = h3aux / (TSUN * pow(xiaux,3)); 
			      
			      /* let's calculate where this bin should fall in the m2 range */			  
			       n = floor( n_m2[map_counter] * (m2aux - m2_l[map_counter])/(m2_u[map_counter] - m2_l[map_counter]) );
			  
			  
			      if ((n > -1) && (n < n_m2[map_counter]))
			        {			  
			          /* let's update that bin */

			          /* update the cosine vectors */
			          cosim2[map_counter][i*n_m2[map_counter] + n] = cosim2[map_counter][i*n_m2[map_counter] + n] + vaux; 
			  
 			          if (cosim2[map_counter][i*n_m2[map_counter] + n] > cosim2max[map_counter])
			              cosim2max[map_counter] = cosim2[map_counter][i*n_m2[map_counter] + n]; /* update maximum of 2-D pdf */
			  
			          cosi[map_counter][i] = cosi[map_counter][i] + vaux; /*update counters for cosi side distribution */
			          if (cosi[map_counter][i] > cosimax[map_counter])
			              cosimax[map_counter] = cosi[map_counter][i]; /* update maximum of 1-D pdf*/
				  
			  
			          m2[map_counter][n] = m2[map_counter][n] + vaux; /*update counters for m2 side distribution */
			          if (m2[map_counter][n] > m2max[map_counter])
			          m2max[map_counter] = m2[map_counter][n]; /* update maximum of 1-D pdf */
			  
			          cosim2_total[map_counter] = cosim2_total[map_counter] + vaux; /* total distribution */
			        }

		              m1aux = sqrt( pow( (m2aux*siniaux), 3 ) / (1.0 * f ) ) - m2aux;

			      m = floor( MPBIN * (m1aux - m1_l[map_counter])/(m1_u[map_counter] - m1_l[map_counter]) );
			      if (( MPBIN * (m1aux - m1_l[map_counter])/(m1_u[map_counter] - m1_l[map_counter]) - m) > 0.5) m = m+1;
			      
			      if ((m > -1) && (m < MPBIN))
				{
				  /* printf("\n i = %d, j = %d, m = %d, v = %f, vmin = %f, vaux = %f.\n", i, j, m, v, vmin, vaux);
				     printf("\n cosi = %f, sini = %f, m2 = %f, m1 = %f.\n", cosiaux, siniaux, m2aux, m1aux); */
				  
				  /* update pulsar mass array */
				  m1m2[map_counter][m * n_m2[map_counter] + n] = m1m2[map_counter][m * n_m2[map_counter] + n] + vaux;
				  /* update its maximum */
				  if (m1m2[map_counter][m * n_m2[map_counter] + n] >  m1m2max[map_counter])
				    m1m2max[map_counter] = m1m2[map_counter][m * n_m2[map_counter] + n];
				  
				  /* update lateral arrays */
				  m1[map_counter][m] = m1[map_counter][m] + vaux;
				  if ( m1[map_counter][m] > m1max[map_counter] ) m1max[map_counter] = m1[map_counter][m];
				  
				  /* update counter of total pulsar mass probability */
				  m1m2_total[map_counter] = m1m2_total[map_counter] + vaux;
				}
			      
			    }		     			    	
		       
			}
		      
		    } /* close j loop */
		  
		} /* close i loop */
	      	      
	      fclose(F);
	    }

	  /* *************************************************************************************** */

	  /* cosi - omega plane map (FOR SYSTEMS WITH VERY GOOD OMEGA-DOT BUT NO OTHER STRONG CONSTRAINTS) */
	  if (map_type[map_counter] == 5)
	    {
	      /* FIRST: read the map */

	      vmin = 1E+9;
	      
	      /* dumb thing: open the correct file */
	      if (map_counter == 0) F = fopen("F1.dat","r");
	      if (map_counter == 1) F = fopen("F2.dat","r");
 	      
	      /* cycling through the cos i coordinate */
	      for (i = 0 ; i < n_cosi[map_counter]; i++)
		{
		  /* cycling through the omega coordinate */
		  for (j = 0 ; j < n_om[map_counter]; j++)
		    {
		      /* read chi-2 */
		      fscanf(F, "%e", &v);
		      if (v < vmin) vmin = v; /* update minimum */
		    }
		}
	      fclose(F);

	      /* Second: Calculate the probabilities */

	      printf("\n Read the chi2 map. Now will calculate the probabilities. \n");
	  
	      /* open file again, since we're going to need the values again */
	      if (map_counter == 0) F = fopen("F1.dat","r");
	      if (map_counter == 1) F = fopen("F2.dat","r");

              /* cycling through the y coordinate */
	      for (j = 0 ; j < n_om[map_counter]; j++)
		{

	          /* cycling through the cos i coordinate */
	          for (i = 0 ; i < n_cosi[map_counter]; i++)
		    {

		      cosiaux = cosix[map_counter][i];
		      siniaux = sqrt(1 - cosiaux * cosiaux);
		      /* now: estimate masses */
		      m2aux = pow( (mtobs*mtobs * f), one3rd) / siniaux;
		      m1aux = mtobs - m2aux;

		      /* have to read array element here, not in if statement that follows,
			 otherwise there are gaps in the reading */
		      fscanf(F, "%e", &v);
		      
		      /* If pulsar mass is larger than zero, update counters */
		      
		      if (m1aux > 0)
			{
			  /* calculate contribution to 2-D pdf (for cosi - m2 plane) */
			  vaux = exp ((vmin - v)/2);
			  /* printf("\n i = %d, j = %d, v = %f, vmin = %f, vaux = %f.\n", i, j, v, vmin, vaux); */
			  			  
			  /* update counters */
			  /* cosiom[map_counter][i*n_om[map_counter] + j] = cosiom[map_counter][i*n_om[map_counter] + j] + vaux; update 2-D counter */
			  cosiom[map_counter][i*n_om[map_counter] + j] = vaux; /* calculating it this way, not as above, because each point should load _once_ */
			  if (cosiom[map_counter][i*n_om[map_counter] + j] > cosiommax[map_counter])
			    cosiommax[map_counter] = cosiom[map_counter][i*n_om[map_counter] + j]; /* update maximum of 2-D pdf */
			  
			  cosi[map_counter][i] = cosi[map_counter][i] + vaux; /*update counters for cosi side distribution */
			  if (cosi[map_counter][i] > cosimax[map_counter])
			    cosimax[map_counter] = cosi[map_counter][i]; /* update maximum of 1-D pdf*/
			  
			  om[map_counter][j] = om[map_counter][j] + vaux; /*update counters for om side distribution */
			  if (om[map_counter][j] > ommax[map_counter])
			    ommax[map_counter] = om[map_counter][j]; /* update maximum of 1-D pdf */
			  
			  cosiom_total[map_counter] = cosiom_total[map_counter] + vaux; /* total distribution */
			  
			  /* update pulsar mass distribution - we will need a finer cosi grid, particularly for the low inclinations where steps in pulsar mass are huge*/

			  /* Value to add is an interpolation between current cosine and previous */
			  /* first calculate the difference between prob. density at this cosine and previous cosine at same omega */
			  if (i > 0)
			    {
			      dv = cosiom[map_counter][i*n_om[map_counter] + j] - cosiom[map_counter][(i-1)*n_om[map_counter] + j];
			      dv = dv / 2.0;
			    }
			  else dv = 0;

			  for (imp = -HNCOSI; imp < HNCOSI; imp++)
			    {
			      /* vaux = cosiom[map_counter][i*n_om[map_counter] + j] + (1.0 * imp)/(1.0 * NCOSI) * dv; */
			      
			      /* Let's vary the cosine by 1/NCOSI of its com1uted interval, and load it to the aux. variable cosi */
			      cosiaux = cosix[map_counter][i-1] + dcosi[map_counter] * imp;
			      siniaux = sqrt(1 -  cosiaux * cosiaux);
			      /* now: estimate masses */
			      m2aux = pow( (mtobs*mtobs * f), one3rd) / siniaux;
			      m1aux = mtobs - m2aux;
			      
			      m = floor( MPBIN * (m1aux - m1_l[map_counter])/(m1_u[map_counter] - m1_l[map_counter]) );
			      
			      if ((m > -1) && (m < MPBIN))
				{
				  
				  /* update pulsar mass array */
				  m1[map_counter][m] = m1[map_counter][m] + vaux;
				  /* update companion mass array */
				  m2[map_counter][m] = m2[map_counter][m] + vaux;
				  
				  if ( m1[map_counter][m] > m1max[map_counter] ) m1max[map_counter] = m1[map_counter][m];
				  
				  /* update counter of total pulsar mass probability for this map */
				  m1om_total[map_counter] = m1om_total[map_counter] + vaux;
				}
			    }
			}	  
		      else
			{
			  cosiom[map_counter][i*n_om[map_counter] + j] = 0;
			  /* no counter then needs updating */
			}
		      
		    } /* close j loop */
		  
		} /* close i loop */
	      m1m2_total[map_counter] = m1om_total[map_counter];
	      
	      fclose(F);
	    }
	  

	  /* *************************************************************************************** */

	  /* RA - Dec plane map  */
	  
	  if (map_type[map_counter] == 10)
	    {
	      /* FIRST: read the map */

	      vmin = 1E+9;
	      
	      /* dumb thing: open the correct file */
	      if (map_counter == 0) F = fopen("F1.dat","r");
	      if (map_counter == 1) F = fopen("F2.dat","r");
 	      
	      /* cycling through the cos i coordinate */
	      for (i = 0 ; i < n_ra[map_counter]; i++)
		{
		  /* cycling through the omega coordinate */
		  for (j = 0 ; j < n_dec[map_counter]; j++)
		    {
		      /* read chi-2 */
		      fscanf(F, "%e %e %e", &raaux, &decaux, &v);
		      if (v < vmin) vmin = v; /* update minimum */
		    }
		}
	      fclose(F);

	      /* Second: Calculate the probabilities */

	      printf("\n Read the chi2 map. Now will calculate the probabilities. \n");
	  
	      /* open file again, since we're going to need the values again */
	      if (map_counter == 0) F = fopen("F1.dat","r");
	      if (map_counter == 1) F = fopen("F2.dat","r");

	      /* cycling through the cos i coordinate */
	      for (i = 0 ; i < n_ra[map_counter]; i++)
		{

		  /* cycling through the y coordinate */
		  for (j = 0 ; j < n_dec[map_counter]; j++)
		    {

		      /* have to read array element here, not in if statement that follows,
			 otherwise there are gaps in the reading */
		      fscanf(F, "%e  %e %e",  &raaux, &decaux, &v);
		    
		      /* calculate contribution to 2-D pdf (for cosi - m2 plane) */
		      vaux = exp ((vmin - v)/2);
		      /* printf("\n i = %d, j = %d, v = %f, vmin = %f, vaux = %f.\n", i, j, v, vmin, vaux);
			 printf("\n cosi = %f, sini = %f, m2 = %f, m1 = %f.\n", cosiaux, siniaux, m2aux, m1aux); */
		      
		      /* update counters - no need for side distributions yet*/

		      radec[map_counter][i*n_dec[map_counter] + j] = vaux; /* calculating it this way, not as above, because each point should load _once_ */
		      if ( radec[map_counter][i*n_dec[map_counter] + j] > radecmax[map_counter])
			   radecmax[map_counter] = radec[map_counter][i*n_dec[map_counter] + j]; /* update maximum of 2-D pdf */
		      
	              radec_total[map_counter] = radec_total[map_counter] + vaux; /* total distribution */
			
		      
		    } /* close j loop */
		  
		} /* close i loop */
	      	      
	      fclose(F);
	    }




	  
	  /* *************************************************************************************** */
	  
	  /* Let's calculate all the planes (FOR SYSTEMS WITH 3-D mapping) */
	  
	  /* first case (11): systems where the 3rd dimension is MTOT */
	  /* second case (12): systems where the 3rd dimension is h3 */
	  if ((map_type[map_counter] == 11)||(map_type[map_counter] == 12))
	    {
	      
	      /* first step was already done above - we found, among other things, the minimum v */
	      
	      /* Second: Calculate the probabilities */
	      
	      printf("\n Will now re-read the chi2 maps and calculate the probabilities. \n");
	      
	      /* open file again, since we're going to need the values again */
	      
	      if (map_counter == 0) F = fopen("F1.dat","r");
	      if (map_counter == 1) F = fopen("F2.dat","r");
	      
              /* now, just one single loop through all elements in the array */
              for (l = 0 ; l < lmax[map_counter]; l++)
		{
		  if (map_type[map_counter] == 11)
		    fscanf(F, "%e %e %e %e", &omgrid, &cosigrid, &mtgrid, &v);		  
		  if (map_type[map_counter] == 12)
		    fscanf(F, "%e %e %e %e", &omgrid, &cosigrid, &h3grid, &v);		  
		  
		  
		  /* calculate contribution to 2-D pdf (for cosi - m2 plane) */
		  vaux = exp ((vmin - v)/2);
		  
		  /* Do the calculations for Omega bin here - this is not being smoothed */
		  omaux = omgrid;
		  j = floor( (omaux - i_om[map_counter]) / s_om[map_counter] );
		  
		  /* let's do a fine grid of cos i and mtot, and store everything after that*/
		      
		  for (imp = -HNCOSI; imp < HNCOSI; imp++)
		    {
		      
		      cosiaux = cosigrid + dcosi[map_counter] * imp;
		      siniaux = sqrt(1 - cosiaux * cosiaux);
		      
		      /* calculate the bin for cos i here */
		      i = floor( (cosiaux - i_cosi[map_counter]) / s_cosi[map_counter] );
		      
		      /* Now, we have to check whether user wants any mass maps. If so, then go a bit firther down to do the smoothing in mass too */
		      
		      for (imt = -HNM; imt < HNM; imt++)
			{
			  
			  /* now, do some relevant calculations for the masses */
			  if (map_type[map_counter] == 11)
			    {
			      mtaux = mtgrid + dmt[map_counter] * imt;
			      m2aux = pow( (mtaux*mtaux * f), one3rd) / siniaux;
			      m1aux = mtaux - m2aux;
			    }
			  
			  if (map_type[map_counter] == 12)
			    {
			      h3aux = h3grid + dh3[map_counter] * imt;
			       /* No need to convert h3 to seconds, as TSUN is already in microseconds */
			      m2aux = h3aux * pow( ((1 + sqrt(1 - siniaux*siniaux)) / siniaux), 3) / TSUN;
			      m1aux = sqrt( pow( (m2aux*siniaux), 3 ) / (1.0 * f ) ) - m2aux;
			      mtaux = m1aux + m2aux;    
			    }
			  			  
			  /* now, if pulsar mass is positive, update counters */
			  
			  if (m1aux > 0)
			    {
			      
			      /* calculate the bins for total mass, companion mass and pulsar mass where this should be updated */
			       
			       n = floor( (m2aux - i_m2[map_counter]) / s_m2[map_counter] );
			       m = floor( MPBIN * (m1aux - m1_l[map_counter])/(m1_u[map_counter] - m1_l[map_counter]) );
			       
			       
			       if (map_type[map_counter] == 11)
				 t = floor( (mtaux - i_mt[map_counter]) / s_mt[map_counter] );
			       
			       /* Now, we must specify all the conditions necessary for updating the grid */
			       if ((j > -1) && (j < (n_om[map_counter])) && (i > -1) && (i < (n_cosi[map_counter])))
				 {
				   
				   /* We could in principle have done this in the loop above. */
				   /* HOWEVER, we want to make this conditional on the pulsar mass being positive */
				   
				   /* update cosi - omega maps and side arrays */
				   cosiom[map_counter][i*n_om[map_counter] + j] =  cosiom[map_counter][i*n_om[map_counter] + j] + vaux; 
				  
				   if (cosiom[map_counter][i*n_om[map_counter] + j] > cosiommax[map_counter])
				     cosiommax[map_counter] = cosiom[map_counter][i*n_om[map_counter] + j]; 
				   
				   cosi[map_counter][i] = cosi[map_counter][i] + vaux; /*update counters for cosi side distribution */
				   
				   if (cosi[map_counter][i] > cosimax[map_counter])
				     cosimax[map_counter] = cosi[map_counter][i]; /* update maximum of 1-D pdf*/
				   
				   om[map_counter][j] = om[map_counter][j] + vaux; /*update counters for om side distribution */
				   
				   if (om[map_counter][j] > ommax[map_counter])
				     ommax[map_counter] = om[map_counter][j]; /* update maximum of 1-D pdf */
				   
				   /* update total count for this map*/
				   cosiom_total[map_counter] = cosiom_total[map_counter] + vaux;

				   if ( (calculate_cosim2[map_counter] == 1) && (n > -1) && (n < (n_m2[map_counter])))
				     {
				       /* update cosi - m2 maps and side arrays. No need to update the cosi map array here, as it was already loaded above */
				       
				       cosim2[map_counter][i*n_m2[map_counter] + n] = cosim2[map_counter][i*n_m2[map_counter] + n] + vaux; 
				       
				       if (cosim2[map_counter][i*n_m2[map_counter] + n] > cosim2max[map_counter])
					 cosim2max[map_counter] = cosim2[map_counter][i*n_m2[map_counter] + n]; /* update maximum of 2-D pdf */
				       
				       m2[map_counter][n] = m2[map_counter][n] + vaux; /*update counters for m2 side distribution */
				       if (m2[map_counter][n] > m2max[map_counter])
					 m2max[map_counter] = m2[map_counter][n]; /* update maximum of 1-D pdf */
				       
				       /* update total distributioncount for this map*/
				       cosim2_total[map_counter] = cosim2_total[map_counter] + vaux; /* total distribution */

				       if ( (calculate_m1m2[map_counter] == 1) && (m > -1) && (m < MPBIN) )
					 {
					   /* update mass-mass array. No need to update the m2 array here, as it was already loaded above */
					   m1m2[map_counter][m * n_m2[map_counter] + n] = m1m2[map_counter][m * n_m2[map_counter] + n] + vaux;
					   
					   if (m1m2[map_counter][m * n_m2[map_counter] + n] >  m1m2max[map_counter])
					     m1m2max[map_counter] = m1m2[map_counter][m * n_m2[map_counter] + n];
					   
					   m1[map_counter][m] = m1[map_counter][m] + vaux; /*update counters for m1 side distribution */
					   if ( m1[map_counter][m] > m1max[map_counter] ) m1max[map_counter] = m1[map_counter][m];
					   
					   /* update total count for this map */
					   m1m2_total[map_counter] = m1m2_total[map_counter] + vaux;
					 }
				     }
				   
				   /* update side distribution for mt */
				   if ( (map_type[map_counter] == 11) && (t > -1) && (t < n_mt[map_counter]) )
				     {
				       mt[map_counter][t] = mt[map_counter][t] + vaux;
				   
				       if (mt[map_counter][t] > mtmax[map_counter])
					 mtmax[map_counter] = mt[map_counter][t]; /* update maximum of 1-D pdf*/
				     }
				 }
			       
			     }
			   
			 }
			 
		    }
		  
		}
	      
	      fclose(F);
	    }
	  
          /* Finished reading the maps */
	  
	  /* *************************************************************************************** */

	  /* now it's time to normalize the maps and write them */

	  /* in which cases do we need to write the distribution of cosines? */
	  if((calculate_cosim2[map_counter] == 1)||(calculate_cosiom[map_counter] == 1))
	    {
	      if (calculate_cosiom[map_counter] == 1) total = cosiom_total[map_counter];
	      if (calculate_cosim2[map_counter] == 1) total = cosim2_total[map_counter];
	      

	      /* normalize the distribution */
	      cumulative = 0;
	      /* start target counter */
	      j = 0;
	      
	      /* open files */
	      if (map_counter == 0) PM = fopen("pdf_cosi_1.dat","w");
	      if (map_counter == 1) PM = fopen("pdf_cosi_2.dat","w");
	      
	      for (i = 0 ; i < n_cosi[map_counter] ; i++)
		{
		  cumulative = cumulative + cosi[map_counter][i];
		  
		  if ( ( (cumulative/total) > D1target[j]) && (j < D1NCONT) )
		    {
		      /* since last iteration, where mc_cumulative/total was still below target,
			 mc_cumulative/total increased by m2[i]/total
			 mc has increased by one pixel width, tr[4].
			 So, where was the line crossed?
			 Lets do linear interpolation */
		      
		      fp = (cumulative - D1target[j] * total) / cosi[map_counter][i];
		      D1_cosi[map_counter][j] = cosix[map_counter][i] - s_cosi[map_counter] * fp;
		      
		      /* determine the height in the probability _density_ function where this occurs */
		      D1_cosiY[map_counter][j] = cosi[map_counter][i] * (1 - fp ) / cosimax[map_counter] + fp * cosi[map_counter][i-1];
		      /* this works because previous entry has been divided by mpmax */
		      		      
		      /* since we crossed this target, we move to the next target */
		      j = j + 1;
		    }
		  
		  cosi[map_counter][i] = cosi[map_counter][i]/cosimax[map_counter];
		  fprintf(PM," %f  %f  %f \n", cosix[map_counter][i], cosi[map_counter][i], cumulative/cosiom_total[map_counter]);
		}

	      printf("\n For distribution %d,\n \n cos i = %f / +%f / -%f to %f percent C. L.;\n", map_counter +1, D1_cosi[map_counter][jmedian], D1_cosi[map_counter][jmedian+1] - D1_cosi[map_counter][jmedian],  D1_cosi[map_counter][jmedian] - D1_cosi[map_counter][jmedian -1], 100*(D1target[jmedian +1] - D1target[jmedian -1]));
	      
              printf(" cos i = %f / +%f / -%f to %f percent C. L.; \n", D1_cosi[map_counter][jmedian], D1_cosi[map_counter][jmedian+2] - D1_cosi[map_counter][jmedian],  D1_cosi[map_counter][jmedian] - D1_cosi[map_counter][jmedian -2], 100*(D1target[jmedian +2] - D1target[jmedian -2]));

              printf("\n i = %f / +%f / -%f degrees to %f percent C. L. \n", 180 * acos(D1_cosi[map_counter][jmedian]) / PI, (acos(D1_cosi[map_counter][jmedian-1]) - acos(D1_cosi[map_counter][jmedian]) )* 180 / PI, ( acos(D1_cosi[map_counter][jmedian]) - acos(D1_cosi[map_counter][jmedian + 1])) * 180 / PI, 100*(D1target[jmedian +1] - D1target[jmedian -1]));
	      
              printf(" i = %f / +%f / -%f degrees to %f percent C. L. \n", 180 * acos(D1_cosi[map_counter][jmedian]) / PI, (acos(D1_cosi[map_counter][jmedian-2]) - acos(D1_cosi[map_counter][jmedian]) )* 180 / PI, ( acos(D1_cosi[map_counter][jmedian]) - acos(D1_cosi[map_counter][jmedian +2])) * 180 / PI, 100*(D1target[jmedian +2] - D1target[jmedian -2]));

	      fclose(PM);
	    }

	  if(calculate_cosim2[map_counter] == 1)
	    {
	      total = cosim2_total[map_counter];

	      /* normalize the distribution */
	      cumulative = 0;
	      /* start target counter */
	      j = 0;
	      
	      /* open files */
	      if (map_counter == 0) PM = fopen("pdf_mc_1.dat","w");
	      if (map_counter == 1) PM = fopen("pdf_mc_2.dat","w");
	      
	      for (i = 0 ; i < n_m2[map_counter]; i++)
		{
		  cumulative = cumulative + m2[map_counter][i];
		  
		  if ( ( (cumulative/total) > D1target[j]) && (j < D1NCONT) )
		    {
		      /* since last iteration, where mc_cumulative/total was still below target,
			 mc_cumulative/total increased by m2[i]/total
			 mc has increased by one pixel width, tr[4].
			 So, where was the line crossed?
			 Lets do linear interpolation */
		      
		      fp = (cumulative - D1target[j] * total) / m2[map_counter][i];
		      D1_m2[map_counter][j] = m2x[map_counter][i] - s_m2[map_counter] * fp;
		      
		      /* determine the height in the probability _density_ function where this occurs */
		      D1_m2Y[map_counter][j] = m2[map_counter][i] * (1 - fp ) / m2max[map_counter] + fp * m2[map_counter][i-1];
		      /* this works because previous entry has been divided by mpmax */
		     
		      
		      /* since we crossed this target, we move to the next target */
		      j = j + 1;
		    }
		  
		  m2[map_counter][i] = m2[map_counter][i]/m2max[map_counter];
		  fprintf(PM," %f  %f  %f \n", m2x[map_counter][i], m2[map_counter][i], cumulative/cosiom_total[map_counter]);
		}

	      printf("\n m2 = %f / +%f / -%f solar masses to %f percent C. L.;", D1_m2[map_counter][jmedian], D1_m2[map_counter][jmedian+1] - D1_m2[map_counter][jmedian],  D1_m2[map_counter][jmedian] - D1_m2[map_counter][jmedian -1], 100*(D1target[jmedian +1] - D1target[jmedian -1]));
	      
	      printf("\n m2 = %f / +%f / -%f solar masses to %f percent C. L.;\n", D1_m2[map_counter][jmedian], D1_m2[map_counter][jmedian+2] - D1_m2[map_counter][jmedian],  D1_m2[map_counter][jmedian] - D1_m2[map_counter][jmedian -2], 100*(D1target[jmedian +2] - D1target[jmedian -2]));
	      
	      
	      
	      fclose(PM);
	    }

	  if((calculate_m1m2[map_counter] == 1)||(calculate_cosiom[map_counter] == 1))
	    {
	      if (calculate_cosiom[map_counter] == 1) total = m1om_total[map_counter];
	      if (calculate_m1m2[map_counter] == 1) total = m1m2_total[map_counter];
		    
	      /* normalize distribution and write it to file */
	      
	      /* start normalization counter */
	      cumulative = 0;
	      /* start target counter */
	      j = 0;
	      
	      /* open files */
	      if (map_counter == 0) PM = fopen("pdf_mp_1.dat","w");
	      if (map_counter == 1) PM = fopen("pdf_mp_2.dat","w");
	      
	      for (i = 0 ; i < MPBIN ; i++)
		{
		  cumulative = cumulative + m1[map_counter][i];
		  
		  if ( ( (cumulative/total) > D1target[j]) && (j < D1NCONT) )
		    {
		      /* since last iteration, where mc_cumulative/total was still below target,
			 mc_cumulative/total increased by m2[i]/total
			 mc has increased by one pixel width, tr[4].
			 So, where was the line crossed?
			 Lets do linear interpolation */
		      
		      fp = (cumulative - D1target[j] * total) / m1[map_counter][i];
		      D1_m1[map_counter][j] = m1x[map_counter][i] - (m1_u[map_counter] - m1_l[map_counter]) / (MPBIN * 1.0) * fp;
		      
		      /* determine the height in the probability _densitypulsar mass_ function where this occurs */
		      D1_m1Y[map_counter][j] = m1[map_counter][i] * (1 - fp ) / m1max[map_counter] + fp * m1[map_counter][i-1];
		      /* this works because previous entry has been divided by mpmax */
		      
		      /* since we crossed this target, we move to the next target */
		      j = j + 1;
		    }
		  
		  m1[map_counter][i] = m1[map_counter][i]/m1max[map_counter];
		  fprintf(PM," %f  %f  %f \n", m1x[map_counter][i], m1[map_counter][i], cumulative/m1m2_total[map_counter]);
		}
	      fclose(PM);

	      if (calculate_m1m2[map_counter] == 1)
		{
		  printf("\n m1 = %f / +%f / -%f solar masses to %f percent C. L.;", D1_m1[map_counter][jmedian], D1_m1[map_counter][jmedian+1] - D1_m1[map_counter][jmedian],  D1_m1[map_counter][jmedian] - D1_m1[map_counter][jmedian -1], 100*(D1target[jmedian +1] - D1target[jmedian -1]));

		  printf("\n m1 = %f / +%f / -%f solar masses to %f percent C. L.; \n", D1_m1[map_counter][jmedian], D1_m1[map_counter][jmedian+2] - D1_m1[map_counter][jmedian],  D1_m1[map_counter][jmedian] - D1_m1[map_counter][jmedian -2], 100*(D1target[jmedian+2] - D1target[jmedian-2]));	      

		}
	    }

	  if(calculate_cosiom[map_counter] == 1)
	    {
	      total = cosiom_total[map_counter];

	      /* normalize the distribution */
	      cumulative = 0;
	      /* start target counter */
	      j = 0;
	      
	      /* open files */
	      if (map_counter == 0) PM = fopen("pdf_om_1.dat","w");
	      if (map_counter == 1) PM = fopen("pdf_om_2.dat","w");
	      
	      for (i = 0 ; i < n_om[map_counter]; i++)
		{
		  cumulative = cumulative + om[map_counter][i];
		  
		  if ( ( (cumulative/total) > D1target[j]) && (j < D1NCONT) )
		    {
		      /* since last iteration, where mc_cumulative/total was still below target,
			 mc_cumulative/total increased by om[i]/total
			 mc has increased by one pixel width, tr[4].
			 So, where was the line crossed?
			 Lets do linear interpolation */
		      
		      fp = (cumulative - D1target[j] * total) / om[map_counter][i];
		      D1_om[map_counter][j] = omx[map_counter][i] - s_om[map_counter] * fp;
		      
		      /* determine the height in the probability _density_ function where this occurs */
		      D1_omY[map_counter][j] = om[map_counter][i] * (1 - fp ) / ommax[map_counter] + fp * om[map_counter][i-1];
		      /* this works because previous entry has been divided by mpmax */
		      
		      /* since we crossed this target, we move to the next target */
		      j = j + 1;
		    }
		  
		  om[map_counter][i] = om[map_counter][i]/ommax[map_counter];
		  fprintf(PM," %f  %f  %f \n", omx[map_counter][i], om[map_counter][i], cumulative/cosiom_total[map_counter]);
		}
	      
	      fclose(PM);

	      printf("\n Omega = %f / +%f / -%f degrees to %f percent C. L.;\n", D1_om[map_counter][jmedian], D1_om[map_counter][jmedian+1] - D1_om[map_counter][jmedian],  D1_om[map_counter][jmedian] - D1_om[map_counter][jmedian -1], 100*(D1target[jmedian +1] - D1target[jmedian -1]));
	      
	      printf(" Omega = %f / +%f / -%f degrees to %f percent C. L.;\n", D1_om[map_counter][jmedian], D1_om[map_counter][jmedian+2] - D1_om[map_counter][jmedian],  D1_om[map_counter][jmedian] - D1_om[map_counter][jmedian -2], 100*(D1target[jmedian +2] - D1target[jmedian -2]));

	    }

	  /* now, the exception here: if distribution is of type 11, process total mass distribution */

	  if(map_type[map_counter] == 11)
	    {
	      total = cosiom_total[map_counter];

	      /* normalize the distribution */
	      cumulative = 0;
	      /* start target counter */
	      j = 0;
	      
	      /* open files */
	      if (map_counter == 0) PM = fopen("pdf_mt_1.dat","w");
	      if (map_counter == 1) PM = fopen("pdf_mt_2.dat","w");
	      
	      for (i = 0 ; i < n_mt[map_counter]; i++)
		{
		  cumulative = cumulative + mt[map_counter][i];
		  
		  if ( ( (cumulative/total) > D1target[j]) && (j < D1NCONT) )
		    {
		      /* since last iteration, where mc_cumulative/total was still below target,
			 mc_cumulative/total increased by om[i]/total
			 mc has increased by one pixel width, tr[4].
			 So, where was the line crossed?
			 Lets do linear interpolation */
		      
		      fp = (cumulative - D1target[j] * total) / mt[map_counter][i];
		      D1_mt[map_counter][j] = mtx[map_counter][i] - s_mt[map_counter] * fp;
		      
		      /* determine the height in the probability _density_ function where this occurs */
		      D1_mtY[map_counter][j] = mt[map_counter][i] * (1 - fp ) / mtmax[map_counter] + fp * mt[map_counter][i-1];
		      /* this works because previous entry has been divided by mpmax */
		      
		      /* since we crossed this target, we move to the next target */
		      j = j + 1;
		    }
		  
		  mt[map_counter][i] = mt[map_counter][i]/mtmax[map_counter];
		  fprintf(PM," %f  %f  %f \n", mtx[map_counter][i], mt[map_counter][i], cumulative/m1m2_total[map_counter]);
		}
	      
	      fclose(PM);

	      printf("\n Mt = %f / +%f / -%f solar masses to %f percent C. L.;", D1_mt[map_counter][jmedian], D1_mt[map_counter][jmedian+1] - D1_mt[map_counter][jmedian],  D1_mt[map_counter][jmedian] - D1_mt[map_counter][jmedian -1], 100*(D1target[jmedian +1] - D1target[jmedian -1]));
	      
	      printf("\n Mt = %f / +%f / -%f solar masses to %f percent C. L.;\n", D1_mt[map_counter][jmedian], D1_mt[map_counter][jmedian+2] - D1_mt[map_counter][jmedian],  D1_mt[map_counter][jmedian] - D1_mt[map_counter][jmedian -2], 100*(D1target[jmedian +2] - D1target[jmedian -2]));
	      

	    } 
	    
	  /* in which cases do we need to write the distribution of cosines? */
	  
	  /* let's now write the ra_dec distribution, not display it yet at this stage */
	  
	  if(calculate_radec[map_counter] == 1)
	    {
	      radec_total[map_counter];

	      /* start target counter */
	      j = 0;
	      
	      /* open files */
	      if (map_counter == 0) PM = fopen("pdf_ra_dec_1.dat","w");
	      if (map_counter == 1) PM = fopen("pdf_ra_dec_2.dat","w");
	      
	      /* find levels */
              for ( k = 0 ; k < nlevels ; k++)
	        { 
	           /* find countour level */
	           level = levels(n_ra[map_counter], n_dec[map_counter], radec_total[map_counter], D2target[k], radecmax[map_counter], radec[map_counter]);

		 }
		 
	      /* let us now write the whole pdf in the external file */     
	      /* cycling through the cos i coordinate */
	      for (i = 0 ; i < n_ra[map_counter]; i++)
		{
		  /* cycling through the y coordinate */
		  for (j = 0 ; j < n_dec[map_counter]; j++)
		    {
                         fprintf(PM," %f  %f  %f \n", rax[map_counter][i], decx[map_counter][j], radec[map_counter][i*n_dec[map_counter] + j]);
		      
		    } /* close j loop */
		  
		} /* close i loop */

	      fclose(PM);
	    }

	    
	      
	  /* Increase the count to denote the next */
	  map_counter = map_counter + 1;
	}
    }

  /* *************************************************************************************** */
  /* *************************************************************************************** */

  /* MAKE PLOTS  */

  /*  ************************************************************************************************** */

  /* Now: what is the display type? This question arises even if we don't load any maps */

  display_type = 10;
  display_marginal = 0;
  display_marginal_percentiles = 0;

  while(display_type != 0)
    {
      printf("\n What type is this figure?");
      printf("\n 0: no (more) figures");
      printf("\n 1: cosi - m2 display");
      printf("\n 2: m1 - m2 display");
      printf("\n 12: cosi - m2 AND m1_m2 display");
      printf("\n 3: orthogonal display");
      printf("\n 5: cosi - om display\n");
      
      scanf("%d", &display_type);

      /* *************************************************************************************** */
      /* *************************************************************************************** */
      
	  
      /* if there is a display, let's ask about its dimensions */
      if (display_type == 1)
	{
	  printf("\n What are the lower and upper limits for display of cos i?\n");
	  scanf(" %f %f", &xl, &xu);
	  printf("\n What are the lower and upper limits for display of mc (solar masses)?\n");
	  scanf(" %f %f", &yl, &yu);
	  
	  /* if maps have been loaded, then we have the option of displaying the 1-d pdfs */
	  if ((calculate_cosim2[0] == 1)||(calculate_cosim2[1] == 1))
	    {
		  printf("\n Display the 1-D distributions?\n");
		  scanf("%d", &display_marginal);
	    }
	  
	  bxl = MARGIN;
	  byl = MARGIN;
	  if (display_marginal == 1)
	    {
	      bxu = MARGIN + SIZE1; byu = MARGIN + SIZE1;
	    }
	  else
	    {
	      bxu = MARGIN + SIZE0; byu = MARGIN + SIZE0;
	    }	  
	}
      
      /* *************************************************************************************** */
      
      if (display_type == 2)
	{
	  printf("\n What are the lower and upper limits for display of mp (solar masses)?\n");
	  scanf(" %f %f", &xpl, &xpu);
	  printf("\n What are the lower and upper limits for display of mc (solar masses)?\n");
	  scanf(" %f %f", &yl, &yu);
	  
	  /* if maps have been loaded, then we have the option of displaying the 1-d pdfs */
	  if ((calculate_m1m2[0] == 1)||(calculate_m1m2[1] == 1))
	    {
	      printf("\n Display the 1-D distributions?\n");
	      scanf("%d", &display_marginal);
	    }
	  
	  bxl = MARGIN;
	  byl = MARGIN;
	  if (display_marginal == 1)
	    {
	      bxu = MARGIN + SIZE1; byu = MARGIN + SIZE1;
	    }
	  else
	    {
	      bxu = MARGIN + SIZE0; byu = MARGIN + SIZE0;
	    }
	  
	  /* This is only in the case of an inset */  
	  /* bxu = MARGIN + 3.0; byu = MARGIN + 3.0; */
	  
	}

      /* *************************************************************************************** */
      
      if (display_type == 12)
	{
	  /* Now: display the cosi - m2 plane, with addition of m1 - m2 plane */

	  printf("\n What are the lower and upper limits for display of cos i?\n");
	  scanf(" %f %f", &xl, &xu);
	  printf("\n What are the lower and upper limits for display of mp (solar masses)?\n");
	  scanf(" %f %f", &xpl, &xpu);
	  printf("\n What are the lower and upper limits for display of mc (solar masses)?\n");
	  scanf(" %f %f", &yl, &yu);
	  /* if maps have been loaded, then we have the option of displaying the 1-d pdfs */

	  /* any arrays calculated for this plane? If so, ask about displaying sidebars */

	  if ((calculate_cosim2[0] == 1)||(calculate_cosim2[1] == 1)||(calculate_m1m2[0] == 1)||(calculate_m1m2[1] == 1))
	    {
	      printf("\n Display the 1-D distributions?\n");
	      scanf("%d", &display_marginal);
	    }

	  /* here always define a small map for each */
	  bxl = MARGIN; bxu = MARGIN + SIZE1; byl = MARGIN; byu = MARGIN + SIZE1;

	}


      /* *************************************************************************************** */
      
      if (display_type == 3)
	{
	  /* Now: display the orthometric array */

	  printf("\n What are the lower and upper limits for display of h3 and h4 (microseconds)?\n");
	  scanf(" %f %f", &xhl, &xhu);
	  yhl = xhl;
	  yhu = xhu;
	  /* if maps have been loaded, then we have the option of displaying the 1-d pdfs */
	  if ((calculate_h3h4[0] == 1)||(calculate_h3h4[1] == 1))
	    {
	      printf("\n Display the 1-D distributions?\n");
	      scanf("%d", &display_marginal);
	    }

	  bxl = MARGIN; 
	  byl = MARGIN;

	  if (display_marginal == 1)
	    {
	      bxu = MARGIN + SIZE1; byu = MARGIN + SIZE1;
	    }
	  else
	    {
	      bxu = MARGIN + SIZE0; byu = MARGIN + SIZE0;
	    }
	}

      /* *************************************************************************************** */
      
      if (display_type == 5)
	{
	  /* Now: display the cosi - om plane (cosiom array) */

	  /* printf("\n What are the lower and upper limits for display of cos i?\n");
	  scanf(" %f %f", &xl, &xu); */
	  
	  xl = -1;
	  xu = 1;
	  
	  /* printf("\n What are the lower and upper limits for display of omega (degrees)?\n");
	  scanf(" %f %f", &yl, &yu); */
	  yl = 0;
	  yu = 360;

	  if ((calculate_cosiom[0] == 1)||(calculate_cosiom[1] == 1))
	    {
	      printf("\n Display the 1-D distributions?\n");
	      scanf("%d", &display_marginal);
	    }

	  bxl = MARGIN;
	  byl = MARGIN;
	  if (display_marginal == 1)
	    {
	      bxu = MARGIN + SIZE1; byu = MARGIN + SIZE1;
	    }
	  else
	    {
	      bxu = MARGIN + SIZE0; byu = MARGIN + SIZE0;
	    }
	}

      /* *************************************************************************************** */
 	      
      if (display_type > 0)
	{
	      
	  if (display_marginal == 1)
	    {
	      printf("\n Display the percentiles for the 1-D distributions?\n");
	      scanf("%d", &display_marginal_percentiles);
	    }
	  
	  printf("\n Please introduce the file_name.ps/graphics device below.\n");
	  if(cpgbeg(0, "?", 1, 1) != 1)
	    return EXIT_FAILURE;
	  
	  /* we have all that is needed to start opening graphics */
	  cpgsch(1.2);
	  cpgsci(1);
	  cpgvsiz(bxl, bxu, byl, byu);
	}

      /* *************************************************************************************** */

      if ((display_type == 1)||(display_type == 12))
	{
      
	  /* open the window */
	  cpgswin(xl, xu, yl, yu);

	  /* specify arrays for science lines IN THE MC - COS I PLANE*/
	  dx_sci = (xu - xl)/(1.0 * MAX);
	  dy_sci = (yu - yl)/(1.0 * MAX);
	  for (i = 0; i < MAX; i++)
	    {
	      x_sci[i] = xl + i * dx_sci;
	      y_sci[i] = yl + i * dy_sci;
	    }
	  
	  /* display science lines here */
	  cpgsci(1);

	  /* first, we're going to highlight the areas excluded by the x-dot measurement */
	  if (kinetix == 1)
	    {
	      
	      /* cpgsci(xdot_color);
	      cpgsls(xdot_style);
	      cpgslw(xdot_width);
	      
	      yg[0] = yl;
	      yg[1] = yu;

	      xg[0] = cosimin;
	      xg[1] = cosimin;
	      cpgline(2, xg, yg);

	      xg[0] = -cosimin;
	      xg[1] = -cosimin;
	      cpgline(2, xg, yg);
	    
	      cpgsci(1);
	      cpgsls(1);
	      cpgslw(2); */
	      
	      
	      /* Now: grey-out excluded areas */
	      /* set fill color */
	      cpgsci(xdot_color);
	      /* set fill style */
	      cpgsfs(excluded_fill);
	  
	      /* first corner of polygon */
	      mx[0] = -cosimin; my[0] = 0;
	      /* second corner of polygon */
	      mx[1] = -cosimin; my[1] = yu;
	      /* third corner of polygon */
	      mx[2] = cosimin; my[2] = yu;
	      /* fourth corner of polygon */
	      mx[3] = cosimin ; my[3] = 0;
	  
	      cpgpoly(4, mx, my);
	  

              cpgsci(1);
	      cpgsls(1);
              cpgslw(2);	    
	      
	    }
	  


	  /* sini lines */
	  if (ksini == 1)
	    {
	      cpgsci(sini_color);
	      cpgsls(sini_style);
	      cpgslw(sini_width);
	      
	      yg[0] = yl;
	      yg[1] = yu;
	      for (j = 0; j < SIGMAS; j++)
		{
		  siniaux = siniobs+sigma[j] * dsiniobs;
		  if ((siniaux > 0)&&(siniaux < 1))
		    {
		      xg[0] = sqrt(1 - siniaux * siniaux);
		      xg[1] = xg[0];
		      cpgline(2, xg, yg);

		      /* do the negative cos i lines */
		      xg[0] = -sqrt(1 - siniaux * siniaux);
		      xg[1] = xg[0];
		      cpgline(2, xg, yg);

		    }
		}
	      cpgsci(1);
	      cpgsls(1);
	      cpgslw(2);	     
	    }

	  /* companion mass lines */
	  if (kmc == 1)
	    {
	      cpgsci(mc_color);
	      cpgsls(mc_style);
	      cpgslw(mc_width);
	      
	      xg[0] = xl;
	      xg[1] = xu;
	      for (j = 0; j < SIGMAS; j++)
		{
		  m2aux = mcobs + sigma[j] * dmcobs;
		  yg[0] = m2aux;
		  yg[1] = m2aux;
		  cpgline(2, xg, yg);
		}
	      cpgsci(1);
	      cpgsls(1);
	      cpgslw(2);
	    }
	  
	  /* tempo2 prediction, with error bars */
	  if (kmp == 1)
	    {
	      cpgsci(mp_color);
	      cpgsls(mp_style);
	      cpgslw(mp_width);
	      
	      for (j = 0; j < SIGMAS; j++)
		{
		  m1aux = mpobs + sigma[j] * dmpobs;
		  for (i = 0; i < MAX; i++)
		    {
		      sini = f * pow((m1aux + y_sci[i]),2) / pow(y_sci[i],3);
		      sini = pow(sini, 0.333333333 );
		      if (sini < 1) xy_sci[i] = sqrt(1 - sini * sini);
		      else xy_sci[i] = 0;    
		    }
		  cpgline(MAX, xy_sci, y_sci);
		}
	      cpgsci(1);
	      cpgsls(1);
	      cpgslw(2);
	    }
	  
	  /* total mass lines */
	  if (komdt == 1)
	    {
	      cpgsci(omdot_color);
	      cpgsls(omdot_style);
	      cpgslw(omdot_width);
	      
	      mtaux = mtobs;
	      for (i = 0; i < MAX; i++)
		{
		  sini = pow( (f * mtaux * mtaux), 0.333333333 )/ y_sci[i];
		  if (sini < 1) xy_sci[i] = sqrt(1 - sini * sini);
		  else xy_sci[i] = 0;
		  xy_sci_n[i] = -xy_sci[i];
		}
	      cpgline(MAX, xy_sci, y_sci);
              cpgline(MAX, xy_sci_n, y_sci);
	      mtaux = mtobsl;
	      for (i = 0; i < MAX; i++)
		{
		  sini = pow( (f * mtaux * mtaux), 0.333333333 )/ y_sci[i];
		  if (sini < 1) xy_sci[i] = sqrt(1 - sini * sini);
		  else xy_sci[i] = 0;    
		  xy_sci_n[i] = -xy_sci[i];
		}
	      cpgline(MAX, xy_sci, y_sci);
	      cpgline(MAX, xy_sci_n, y_sci);
	      mtaux = mtobsu;
	      for (i = 0; i < MAX; i++)
		{
		  sini = pow( (f * mtaux * mtaux), 0.333333333 )/ y_sci[i];
		  if (sini < 1) xy_sci[i] = sqrt(1 - sini * sini);
		  else xy_sci[i] = 0;    
		  xy_sci_n[i] = -xy_sci[i];
		}
	      cpgline(MAX, xy_sci, y_sci);
	      cpgline(MAX, xy_sci_n, y_sci);
	      if (kinom == 1)
		{
		  cpgsci(omdot_color);
                  cpgsls(omdot_style + 1);
		  mtaux = mtkinl;
		  for (i = 0; i < MAX; i++)
		    {
		      sini = pow( (f * mtaux * mtaux), 0.333333333 )/ y_sci[i];
		      if (sini < 1) xy_sci[i] = sqrt(1 - sini * sini);
		      else xy_sci[i] = 0;    
		      xy_sci_n[i] = -xy_sci[i];
		    }
		  cpgline(MAX, xy_sci, y_sci);
		  cpgline(MAX, xy_sci_n, y_sci);
		  mtaux = mtkinu;
		  for (i = 0; i < MAX; i++)
		    {
		      sini = pow( (f * mtaux * mtaux), 0.333333333 )/ y_sci[i];
		      if (sini < 1) xy_sci[i] = sqrt(1 - sini * sini);
		      else xy_sci[i] = 0;    
		      xy_sci_n[i] = -xy_sci[i];
		    }
		  cpgline(MAX, xy_sci, y_sci);
		  cpgline(MAX, xy_sci_n, y_sci);
		}
	      cpgsci(1);
	      cpgsls(1);
	      cpgslw(2);
	    }
	  
	  if (kgamma == 1)
	    {
	      cpgsci(gamma_color);
	      cpgsls(gamma_style);
	      cpgslw(gamma_width);
	      
	      for (j = 0; j < SIGMAS; j++)
		{
		  gammaux = gammaobs + sigma[j] * dgammaobs;
		  
		  /* calculate the smallest total mass the system can have */
		  /* calculate the value of m2 for m1 = 0 */
		  m2aux = gammaux / (2 * gammak);
		  m2aux = pow(m2aux,1.5);
		  		  
		  /* let's now calculate the largest masses compatible with this gamma */
		  /* maximum mass ratio compatible with a gamma */
      		  q_max =  (gammaux) * pb / (ecc * asini * asini * 2 * PI) - 2;
      
                  /* maximum companion mass compatible with this gamma */
                  mc_max = pow((q_max + 1),2.0) * f;
      
                  /* maximum pulsar mass compatible with gamma */
                  mp_max = q_max * mc_max;
      
                  /* maximum total mass */
                  mtot_max = mp_max + mc_max;

		  for (i = 0; i < MAX; i++)
		    {
		       /* cycle through the total masses */
		       mt_aux = m2aux + i * (mtot_max - m2aux) / (MAX * 1.0);
		       
		       /* for each total mass calculate companion mass */
		       mc_gamma = (sqrt( mt_aux * mt_aux + 4 * pow(mt_aux, four3rd) * gammaux / gammak ) - mt_aux) / 2.0;
      
      	               /* calculate pulsar mass */
                       mp_gamma = mt_aux - mc_gamma;
		       
		       /* calculate cos i */
		       siniaux = pow( (f * mt_aux * mt_aux), one3rd) / mc_gamma;
		       
		       cosiaux = sqrt(1 - siniaux*siniaux); 
		       
		       /* load to vectors */
		       mc_sci[i] = mc_gamma;
		       cosi_sci[i] = cosiaux;
		    }
		  /* now add las points in the chain */
		  mc_sci[MAX] = mc_max;
		  cosi_sci[MAX] = 0;
		    
		  cpgline(MAX + 1, cosi_sci, mc_sci);
		  
		  for (i = 0; i < (MAX + 1); i++)
		    {
		       cosi_sci[i] = -cosi_sci[i];
		    }
		  cpgline(MAX + 1, cosi_sci, mc_sci);
		    
		}
	      cpgsci(1);
	      cpgsls(1);
	      cpgslw(2);
	    }
	  
	  if (kpbdot == 1)
	    {
	      
	      cpgsci(pbdot_color);
	      cpgsls(pbdot_style);
	      cpgslw(pbdot_width);
	      
	      pbdotaux = pbdotobs;
	      yx_sci[0] = 100;
	      xy_sci[0] = 1;
	      for (j = 0; j < SIGMAS; j++)
		{
		  pbdotaux = pbdotobs + sigma[j] * dpbdotobs;
		  for (i = 1; i < MAX; i++)
		    {
		      mr = 0.025 * i;
		      m2aux = pbdotaux * pow((1 + mr), one3rd)/(mr * pbdotk);
		      m2aux = pow(m2aux, 0.6);
		      yx_sci[i] = m2aux;
		      m1aux = m2aux * mr;
		      siniaux = f * pow ((m1aux + m2aux), 2);
		      siniaux = pow(siniaux, one3rd)/ m2aux;
		      xy_sci[i] = sqrt(1 - siniaux * siniaux);
		    }
		  cpgline(MAX, xy_sci, yx_sci);
		}
	      cpgsci(1);
	      cpgsls(1);
	      cpgslw(2);
	    }
	  
	  /* h3 lines */
	  if (kh3 == 1)
	    {
	      cpgsci(h3_color);
	      cpgsls(h3_style);
	      cpgslw(h3_width);
	      
	      for (j = 0; j < SIGMAS; j++)
		{
		  h3aux = h3obs + sigma[j] * dh3obs;
		  if (h3aux > 0)
		    {
		      for (i = 0; i < MAX; i++)
			{
			  sini = sqrt(1 - x_sci[i]*x_sci[i]);
			  xy_sci[i] = h3aux * pow( ((1 + fabs(x_sci[i]))/ sini), 3)/ TSUN;
			  xy_sci_n[i] = -x_sci[i];
			}
		      cpgline(MAX, x_sci, xy_sci);
		      cpgline(MAX, xy_sci_n, xy_sci);
		    }
		}
	      cpgsci(1);
	      cpgsls(1);
	      cpgslw(2);
	    }
	  
	  /* h4 lines */
	  if (kh4 == 1)
	    {
	      cpgsci(h4_color);
	      cpgsls(h4_style);
	      cpgslw(h4_width);
	      
	      for (j = 0; j < SIGMAS; j++)
		{
		  h4aux = h4obs + sigma[j] * dh4obs;
		  if (h4aux > 0)
		    {
		      for (i = 0; i < MAX; i++)
			{
			  sini = sqrt(1 - x_sci[i]*x_sci[i]);
			  xy_sci[i] = h4aux * pow( ((1 + x_sci[i])/ sini), 4)/ TSUN;
			}
		      cpgline(MAX, x_sci, xy_sci);
		    }
		}
	      cpgsci(1);
	      cpgsls(1);
	      cpgslw(2);
	    }
	  
	  /* xi lines */
	  if (kxi == 1)
	    {
	      cpgsci(xi_color);
	      cpgsls(xi_style);
	      cpgslw(xi_width);
	      
	      yg[0] = yl;
	      yg[1] = yu;
	      for (j = 0; j < SIGMAS; j++)
		{
		  xiaux = xiobs + sigma[j] * dxiobs;
		  if ((xiaux > 0)&&(xiaux < 1))
		    {
		      cosiaux = (1 - xiaux*xiaux)/(1 + xiaux*xiaux);
		      xg[0] = cosiaux;
		      xg[1] = xg[0];
		      cpgline(2, xg, yg);
                      /* do the negative */
		      xg[0] = -cosiaux;
		      xg[1] = xg[0];
		      cpgline(2, xg, yg);		      
		    }
		}
	      cpgsci(1);
	      cpgsls(1);
	      cpgslw(2);
	    }
	  
	  /* MAKE AN X-DOT LINE HERE AS WELL */
	  
	 
	  /* print this line in gray */

	  if (xl < 0)
	    {
	      /* display i = 90 degrees */
	      cpgsci(excluded_color);
	      cpgsls(1);
	      cpgslw(2);
	      
	      xg[0] = 0;
	      yg[0] = yl;
	      xg[1] = 0;
	      yg[1] = yu;
	      cpgline(2, xg, yg); 
	    }	      
	   

	  if (kr == 1)
	    {
	      cpgsci(r_color);
	      cpgsls(r_style);
	      cpgslw(r_width);
	      
	      for (j = 0 ; j < SIGMAS; j++)
		{
		  mr = mrobs + sigma[j] * dmrobs;
		  for (i = 0; i < MAX; i++)
		    {
		      sini = pow( (  f*(1 + mr)*(1 + mr)/y_sci[i] ), 0.333333333 );
		      if (sini < 1) xy_sci[i] = sqrt(1 - sini * sini);
		      else xy_sci[i] = 0;    
		    }
		  cpgline(MAX, xy_sci, y_sci);
		}
	      cpgsci(1);
	      cpgsls(1);
	      cpgslw(2);
	    }
	  
	  /* put in the contours (over the science lines) */
	  
	  /* replace i < 1 by i < 1 */

	  if (map_number > 0)
	    {
	      map_counter = 0;
	      while (map_counter < map_number)
		{
		  
		  /* change display origin slightly */
		  if (calculate_cosim2[map_counter] == 1)
		    {
		      tr_cosim2[map_counter][0] = tr_cosim2[map_counter][0] - tr_cosim2[map_counter][2];
		      tr_cosim2[map_counter][3] = tr_cosim2[map_counter][3] - tr_cosim2[map_counter][4];	      
		      cpgsci(distribution_color[map_counter]);
		      cpgslw(distribution_width[map_counter]);
	     
		      for ( k = 0 ; k < nlevels ; k++)
			{ 
			  /* find countour level */
			  level = levels(n_m2[map_counter], n_cosi[map_counter], cosim2_total[map_counter], D2target[k], cosim2max[map_counter], cosim2[map_counter]);
		  
			  /* display it */
			  cpgcont(cosim2[map_counter], n_m2[map_counter], n_cosi[map_counter], 1, n_m2[map_counter], 1, n_cosi[map_counter], &level, -1, tr_cosim2[map_counter]);

			}
		    }

		  /* print map (for Vivek) */
		  if (map_counter == 0)
		    {
		      COSI_M2_1 = fopen("cosi_m2_1.dat","w");		  
		       
		      for (i = 0 ; i < n_m2[map_counter]*n_cosi[map_counter]; i++)
			{
			  fprintf(COSI_M2_1," %f \n", cosim2[map_counter][i]);
			}
		      printf("Array dimensions: %d x %d \n ", n_m2[map_counter], n_cosi[map_counter]);
		      fclose(COSI_M2_1);
		    }

		  if (map_counter == 1)
		    {
		      COSI_M2_2 = fopen("cosi_m2_2.dat","w");		  
		       
		      for (i = 0 ; i < n_m2[map_counter]*n_cosi[map_counter]; i++)
			{
			  fprintf(COSI_M2_2," %f \n", cosim2[map_counter][i]);
			}

		      printf("Array dimensions: %d x %d \n ", n_m2[map_counter], n_cosi[map_counter]);
		      fclose(COSI_M2_2);
		    }
 		      
		  map_counter = map_counter + 1;
		  
		}
	    }

	  cpgslw(2);
	  cpgsci(1);
	
	  /* display i = 90 degrees */
	  cpgsci(excluded_color);
	  cpgsls(1);
	  cpgslw(2);

	  xg[0] = 0;
	  yg[0] = 0;
	  xg[1] = 0;
	  yg[1] = yu;
	  cpgline(2, xg, yg);

	
	  /* finally: if this is a simulation, mark the values used for the simulation */
      
          ksimcount = 0;
	  
	  while (ksimcount < ksim)
	    {
	      cpgsci(scolor[ksimcount]);
	      xg[0] = sqrt(1 - ssini[ksimcount]*ssini[ksimcount]);
	      yg[0] = sm2[ksimcount];
	      cpgpt(1, xg, yg, ssymbol[ksimcount]);
	      cpgsci(1);
	      ksimcount = ksimcount + 1;
	    }
	  

	  /* exclude areas given bu mp < 0 */
	  /* set fill color */
	  cpgsci(excluded_color);
	  
	  /* set fill style */
	  cpgsfs(excluded_fill);
	  
	  /* start at the smallest companion mass */
	  
	  /* first corner of polygon - start at the lowest possible companion mass*/
	  for (i = 0 ; i < MPBIN; i++)
	    {
	      mcx_m[i] = f + i * (yu - f)/ (1.0 * MPBIN);
	      siniaux = pow( (f/mcx_m[i]), one3rd);
	      mpx_m[i] = sqrt( 1 - siniaux * siniaux);
	    }
	  /* last corners */
	  
	  mcx_m[MPBIN] = yu;
	  mpx_m[MPBIN] = mpx_m[i-1];
	  
	  mcx_m[MPBIN+1] = yu;
	  mpx_m[MPBIN+1] = xu;
	  
	  mcx_m[MPBIN+2] = yl;
	  mpx_m[MPBIN+2] = xu;
	  
	  mcx_m[MPBIN+3] = yl;
	  mpx_m[MPBIN+3] = xl;
	  
	  /* Draw the polygon  */
	  cpgpoly(MPBIN+4, mpx_m, mcx_m);

	  for (i = 0 ; i < MPBIN + 4; i++)
	    {
	      mpx_m[i] = -mpx_m[i];
	    }
	  cpgpoly(MPBIN+4, mpx_m, mcx_m);
	  
	  /* back to black */
	  cpgsci(1);
	  cpgsls(1);
	  cpgslw(2);
	  
	  /* last: labels and box marks */
	  
	  cpgbox("BCTNS1", 0.0, 0, "BCTNS1", 0.0, 0);
	  cpglab("cos i", "Companion Mass (M\\d\\(2281)\\u)","");
	  
	}

      /* *************************************************************************************** */

      if ((display_type == 2) || (display_type == 12))
	{

	  /* If display type is 12, open an extra window */
	    
	  if (display_type == 12) cpgvsiz(bxu + GAP, bxu + GAP + SIZE1, byl, byu);
	     
	  cpgswin(xpl, xpu, yl, yu);
	  
	  
	  /* make a white window in the background (important for insets) */
	  /* first corner of polygon */
	  mx[0] = xpl; my[0] = yl;
	  /* second corner of polygon */
	  mx[1] = xpl; my[1] = yu;
	  /* third corner of polygon */
	  mx[2] = xpu; my[2] = yu;
	  /* fourth corner of polygon */
	  mx[3] = xpu; my[3] = yl;
	 
	  /* set fill color grey */
	  cpgsci(background_color);
	  /* set fill style */
	  cpgsfs(background_fill);
	  
	  cpgsci(0);
	  cpgpoly(4, mx, my);


	  /* now: put in the limit derived from the x-dot */


	  if (kinetix == 1)
	    {

              /* calculate the limiting sin i */
              siniaux = sqrt( 1 - cosimin * cosimin);

	      /* these could be a line
	      cpgsci(xdot_color);
	      cpgsls(xdot_style);
	      cpgslw(xdot_width); 
	      
	      for (i = 0; i < MAX; i++)
		{
		  xy_sci[i] = sqrt (pow( (y_sci[i] * siniaux), 3)/ f) - y_sci[i];
		  if (xy_sci[i] < 0) xy_sci[i] = 0;
		}
	      cpgline(MAX, xy_sci, y_sci);

	      cpgsci(1);
	      cpgsls(1);
	      cpgslw(2); 	*/
	      
	      /* or they could be an exclusion zone, like the region with sini > 1 */
	  
	      /* set fill color */
	      cpgsci(xdot_color);
	      /* set fill style */
	      cpgsfs(excluded_fill);
	  
	      /* start at the smallest companion mass */
	  
	      for (i = 0 ; i < MPBIN; i++)
	        {
	          mcx_m[i] = yl + i * (yu - yl)/ (1.0 * MPBIN);
	          mpx_m[i] = sqrt( pow( (mcx_m[i] * siniaux),3)/ f ) - mcx_m[i];
	        }
	      /* last corners */
	  
	      mcx_m[MPBIN] = yu;
	      mpx_m[MPBIN] = mpx_m[i-1];
	  
  	      mcx_m[MPBIN+1] = yu;
	      mpx_m[MPBIN+1] = xpu;
	  
	      mcx_m[MPBIN+2] = yl;
	      mpx_m[MPBIN+2] = xpu;
	  
	      /* Draw the polygon  */
	      cpgpoly(MPBIN+3, mpx_m, mcx_m);
	      
	      	
	    }


	  /* specify arrays for science lines IN THE MC - MP PLANE*/
	  dx_sci = (xpu - xpl)/(1.0 * MAX);
	  dy_sci = (yu - yl)/(1.0 * MAX);

	  for (i = 0; i < MAX; i++)
	    {
	      x_sci[i] = xpl + i * dx_sci;
	      y_sci[i] = yl + i * dy_sci;	      
	    }


	  /* second: display science lines here */
	  
	  /* sini lines */
	  
	  if (ksini == 1)
	    {
	      
	      cpgsci(sini_color);
	      cpgsls(sini_style);
	      cpgslw(sini_width);
	      
	      for (j = 0; j < SIGMAS; j++)
		{
		  sini = siniobs + sigma[j] * dsiniobs;
		  if ((sini > 0)&&(sini < 1))
		    {
		      for (i = 0; i < MAX; i++)
			{
			  xy_sci[i] = sqrt (pow( (y_sci[i] * sini), 3)/ f) - y_sci[i];
			  if (xy_sci[i] < 0) xy_sci[i] = 0;
			}
		      cpgline(MAX, xy_sci, y_sci);
		    }
		}
	      
	      cpgsci(1);
	      cpgsls(1);
	      cpgslw(2);
	      
	    }
	  
	  /* companion mass lines */
	  if (kmc == 1)
	    {
	      cpgsci(mc_color);
	      cpgsls(mc_style);
	      cpgslw(mc_width);
	      
	      xg[0] = xpl;
	      xg[1] = xpu;
	      for (j = 0; j < SIGMAS; j++)
		{
		  m2aux = mcobs + sigma[j] * dmcobs;;
		  yg[0] = m2aux;
		  yg[1] = m2aux;
		  cpgline(2, xg, yg);
		}
	      
	      cpgsci(1);
	      cpgsls(1);
	      cpgslw(2);
	      
	    };
	  
	  /* pulsar lines */
	  if (kmp == 1)
	    {
	      cpgsci(mp_color);
	      cpgsls(mp_style);
	      cpgslw(mp_width);
	      
	      yg[0] = yl;
	      yg[1] = yu;
	      
	      for (j = 0; j < SIGMAS; j++)
		{
		  m1aux = mpobs + sigma[j] * dmpobs;
		  yg[0] = m1aux;
		  yg[1] = m1aux;
		  cpgline(2, xg, yg);
		}
	      
	      cpgsci(1);
	      cpgsls(1);
	      cpgslw(2);
	    }
	  
	  /* total mass lines */
	  if (komdt == 1)
	    {
	      
	      cpgsci(omdot_color);
	      cpgsls(omdot_style);
	      cpgslw(omdot_width);
	      
	      yg[0] = 0;
	      xg[1] = 0;
	      
	      xg[0] = mtobs;
	      yg[1] = mtobs;
	      cpgline(2, xg, yg);
	      
	      xg[0] = mtobsl;
	      yg[1] = mtobsl;
	      cpgline(2, xg, yg);
	      
	      xg[0] = mtobsu;
	      yg[1] = mtobsu;
	      cpgline(2, xg, yg);
	      
	      if (kinom == 1)
		{
		  cpgsci(omdot_color);
                  cpgsls(omdot_style+1);
		  
		  xg[0] = mtkinu;
		  yg[1] = mtkinu;
		  cpgline(2, xg, yg);
		  
		  xg[0] = mtkinl;
		  yg[1] = mtkinl;
		  cpgline(2, xg, yg);	  
		}
              cpgsci(1);
	      cpgsls(1);
	      cpgslw(2); 
	    }
	  
	  /* gamma lines */
	  if (kgamma == 1)
	    {	      	
	      cpgsci(gamma_color);
	      cpgsls(gamma_style);
	      cpgslw(gamma_width);
	  
	      for (j = 0; j < SIGMAS; j++)
		{
		  gammaux = gammaobs + sigma[j] * dgammaobs;
		  
		  /* calculate the smallest total mass the system can have */
		  /* calculate the value of m2 for m1 = 0 */
		  m2aux = gammaux / (2 * gammak);
		  m2aux = pow(m2aux,1.5);
		  		  
		  /* let's now calculate the largest masses compatible with this gamma */
		  /* maximum mass ratio compatible with a gamma */
      		  q_max =  (gammaux) * pb / (ecc * asini * asini * 2 * PI) - 2;
      
                  /* maximum companion mass compatible with this gamma */
                  mc_max = pow((q_max + 1),2.0) * f;
      
                  /* maximum pulsar mass compatible with gamma */
                  mp_max = q_max * mc_max;
      
                  /* maximum total mass */
                  mtot_max = mp_max + mc_max;

		  for (i = 0; i < MAX; i++)
		    {
		       /* cycle through the total masses */
		       mt_aux = m2aux + i * (mtot_max - m2aux) / (MAX * 1.0);
		       
		       /* for each total mass calculate companion mass */
		       mc_gamma = (sqrt( mt_aux * mt_aux + 4 * pow(mt_aux, four3rd) * gammaux / gammak ) - mt_aux) / 2.0;
      
      	               /* calculate pulsar mass */
                       mp_gamma = mt_aux - mc_gamma;
		       
		       /* calculate cos i */
		       siniaux = pow( (f * mt_aux * mt_aux), one3rd) / mc_gamma;
		       
		       cosiaux = sqrt(1 - siniaux*siniaux); 
		       
		       /* load to vectors */
		       mc_sci[i] = mc_gamma;
		       mp_sci[i] = mp_gamma;
		    }
		    
		  /* now add last points in the chain */
		  mp_sci[MAX] = mp_max;
		  mc_sci[MAX] = mc_max;

		  cpgline(MAX+1, mp_sci, mc_sci);
		
		}
	      cpgsci(1);
	      cpgsls(1);
	      cpgslw(2);
	    }
	  
	  if (kpbdot == 1)
	    {
              cpgsci(pbdot_color);
	      cpgsls(pbdot_style);
	      cpgslw(pbdot_width); 
	      
	      /* ludicrous values for first entry in vector 
		 - not likely to be displayed */
	      yx_sci[0] = 100;
	      xy_sci[0] = 0;
	      
	      for (j = 0; j < SIGMAS; j++)
		{
		  pbdotaux = pbdotobs + sigma[j]*dpbdotobs;
		  
		  for (i = 1; i < MAX; i++)
		    {
		      /* First: cycle through mass ratios */
		      mr = 0.025 * i;
		      m2aux = pbdotaux * pow((1 + mr), one3rd)/(mr * pbdotk);
		      m2aux = pow(m2aux,0.6);
		      yx_sci[i] = m2aux;
		      xy_sci[i] = m2aux * mr;
		    }
		  cpgline(MAX, xy_sci, yx_sci);
		}              
	      cpgsci(1);
	      cpgsls(1);
	      cpgslw(2); 
	    }
	  
	  /* h3 lines */
	  if (kh3 == 1)
	    {
	      
              cpgsci(h3_color);
	      cpgsls(h3_style);
	      cpgslw(h3_width); 
	      
	      for (j = 0; j < SIGMAS; j++)
		{
		  h3aux = h3obs + sigma[j]*dh3obs;
		  if (h3aux > 0)
		    {
		      for (i = 0; i < MAX; i++)
			{
			  
			  xiaux = pow ( (h3aux/(TSUN* y_sci[i])), 0.3333333333333 );
			  
			  if (xiaux > 1) siniaux = 1.05;
			  else  siniaux = 2 * xiaux / (1 + xiaux*xiaux);
			  
			  /* calculate pulsar mass from companion mass, sini and f - see eqs. above */
			  xy_sci[i] = sqrt (pow( (y_sci[i] * siniaux), 3)/ f) - y_sci[i];
			}
		      cpgline(MAX, xy_sci, y_sci);
		    }
		}
	      cpgsci(1);
	      cpgsls(1);
	      cpgslw(2); 
	    }
      
	  /* h4 lines */
	  if (kh4 == 1)
	    {
              cpgsci(h4_color);
	      cpgsls(h4_style);
	      cpgslw(h4_width); 

	      for (j = 0; j < SIGMAS; j++)
		{
		  h4aux = h4obs + sigma[j]*dh4obs;
		  if (h4aux > 0)
		    {
		      for (i = 0; i < MAX; i++)
			{
			  
			  xiaux = pow ( (h4aux/(TSUN* y_sci[i])), 0.25 );
			  
			  if (xiaux > 1) siniaux = 1.05;
			  else  siniaux = 2 * xiaux / (1 + xiaux*xiaux);
			  
			  /* calculate pulsar mass from companion mass, sini and f - see eqs. above */
			  xy_sci[i] = sqrt (pow( (y_sci[i] * siniaux), 3)/ f) - y_sci[i];
			}
		      cpgline(MAX, xy_sci, y_sci);
		    }
		}
	      cpgsci(1);
	      cpgsls(1);
	      cpgslw(2); 
	    }
	  
	  if (kxi == 1)
	    {
	      cpgsci(xi_color);
	      cpgsls(xi_style);
	      cpgslw(xi_width); 
	      
	      for (j = 0; j < SIGMAS; j++)
		{
		  xiaux = xiobs + sigma[j]*dxiobs;
		  if ((xiaux > 0) && (xiaux < 1.0))
		    {
		      siniaux = 2 * xiaux / (1 + xiaux * xiaux);
		      for (i = 0; i < MAX; i++)
			{
			  xy_sci[i] = sqrt (pow( (y_sci[i] * siniaux), 3)/ f) - y_sci[i];
			  if (xy_sci[i] < 0) xy_sci[i] = 0;
			}
		      cpgline(MAX, xy_sci, y_sci);
		    }
		}
	      cpgsci(1);
	      cpgsls(1);
	      cpgslw(2); 
	    }
	    
	  /* display the mass ratio */
	  if (kr == 1)
	    {
	      cpgsci(r_color);
	      cpgsls(r_style);
	      cpgslw(r_width); 

	      xg[0] = 0;
	      yg[0] = 0;

	      xg[1] = xpu;

	      for (j = 0; j < SIGMAS; j++)
		{
		  mr = mrobs + sigma[j] * dmrobs;
		  yg[1] = xpu / mr;
		  cpgline(2, xg, yg);
		}
	      cpgsci(1);
	      cpgsls(1);
	      cpgslw(2); 		
	    }

	  if (map_number > 0)
	    {
	      map_counter = 0;
	      while (map_counter < map_number)
		{
		  
		  /* change display origin slightly */
		  if (calculate_m1m2[map_counter] == 1)
		    { 
		      cpgsci(distribution_color[map_counter]);
		      cpgslw(distribution_width[map_counter]);
		      for ( k = 0 ; k < nlevels ; k++)
			{ 
			  /* find countour level */
			  level = levels(n_m2[map_counter], MPBIN, m1m2_total[map_counter], D2target[k], m1m2max[map_counter], m1m2[map_counter]);
		  
			  /* display it */
			  cpgcont(m1m2[map_counter], n_m2[map_counter], MPBIN, 1, n_m2[map_counter], 1, MPBIN, &level, -1, tr_m1m2[map_counter]);

			}
		    }

		  /* print map (for Vivek) */
		  if (map_counter == 0)
		    {
		      M1_M2_1 = fopen("m1_m2_1.dat","w");		  
		       
		      for (i = 0 ; i < n_m2[map_counter]*MPBIN; i++)
			{
			  fprintf(M1_M2_1," %f \n", m1m2[map_counter][i]);
			}
		      printf("Array dimensions: %d x %d \n ", n_m2[map_counter], MPBIN);
		      fclose(M1_M2_1);
		    }

		  if (map_counter == 1)
		    {
		      M1_M2_2 = fopen("m1_m2_2.dat","w");		  
		       
		      for (i = 0 ; i < n_m2[map_counter]*MPBIN; i++)
			{
			  fprintf(M1_M2_2," %f \n", m1m2[map_counter][i]);
			}
		      printf("Array dimensions: %d x %d \n ", n_m2[map_counter], MPBIN);
		      fclose(M1_M2_2);
		    }
		  
		  map_counter = map_counter + 1;
		  
		}
	    }
	      
	  cpgslw(2);
	  cpgsci(1);
	  
	  
	  /* fourth: highlight the region excluded by knowledge of the mass function */
	  
	  /* set fill color grey */
	  cpgsci(excluded_color);
	  /* set fill style */
	  cpgsfs(excluded_fill);
	  
	  /* start at the smallest companion mass */
	  
	  for (i = 0 ; i < MPBIN; i++)
	    {
	      mcx_m[i] = yl + i * (yu - yl)/ (1.0 * MPBIN);
	      mpx_m[i] = sqrt( pow(mcx_m[i],3)/ f ) - mcx_m[i];
	    }
	  /* last corners */
	  
	  mcx_m[MPBIN] = yu;
	  mpx_m[MPBIN] = mpx_m[i-1];
	  
	  mcx_m[MPBIN+1] = yu;
	  mpx_m[MPBIN+1] = xpu;
	  
	  mcx_m[MPBIN+2] = yl;
	  mpx_m[MPBIN+2] = xpu;
	  
	  /* Draw the polygon  */
	  cpgpoly(MPBIN+3, mpx_m, mcx_m);


	  
	  
	  /* ******************** */
	  
	  
	  /* display the special symbols */

	  ksimcount = 0;
	  while (ksimcount < ksim)
	    {
	      cpgsci(scolor[ksimcount]);
	      
	      xg[0] = sm1[ksimcount];
	      yg[0] = sm2[ksimcount];
	      cpgpt(1, xg, yg, ssymbol[ksimcount]);
	      ksimcount = ksimcount + 1;
	    }
	  
	  /* Last: labels and box marks */
	  
	  cpgsci(1);
	  cpgsls(1);
	  cpgslw(2);

	  if (display_type == 2){
	         /* comment out in case of inset */
		 cpglab("Pulsar Mass (M\\d\\(2281)\\u)", "Companion Mass (M\\d\\(2281)\\u)","");
		 cpgbox("BCTNS1", 0.0, 0, "BCTNS1", 0.0, 0);
             }
	  if (display_type == 12){
                 cpglab("Pulsar Mass (M\\d\\(2281)\\u)", "", "");
	         cpgbox("BCTNS1", 0.0, 0, "BCTS", 0.0, 0);
	     } 
	}

      /* *************************************************************************************** */

      if (display_type == 5)
	{
	

	  /* specify arrays for science lines */
	  dx_sci = (xu - xl)/(1.0 * MAX);
	  dy_sci = (yu - yl)/(1.0 * MAX);
	  for (i = 0; i < MAX; i++)
	    {
	      x_sci[i] = xl + i * dx_sci;
	      y_sci[i] = yl + i * dy_sci;
	    }
	  	  
	  cpgswin(xl, xu, yl, yu);
	  
	  
	  /* display the gray levels of the 2-D pdf */
	  if (map_number > 0)
	    {
	      map_counter = 0;
	      while (map_counter < map_number)
		{
		  if (calculate_cosiom[map_counter] == 1)
		    {
		      /* third: put in the contour levels (over the science lines)  */
		      cpgsci(distribution_color[map_counter]);
		      cpgslw(distribution_width[map_counter]);
		      /* move one bin down */
		      tr_cosiom[map_counter][0] = tr_cosiom[map_counter][0] - tr_cosiom[map_counter][2];
		      tr_cosiom[map_counter][3] = tr_cosiom[map_counter][3] - tr_cosiom[map_counter][4];
		      
		      /* add an option here to display a greyscale */
		      cpggray(cosiom[map_counter], n_om[map_counter], n_cosi[map_counter], 1, n_om[map_counter], 1, n_cosi[map_counter], cosiommax[map_counter], 0, tr_cosiom[map_counter]);
		      
		    }
		  map_counter = map_counter + 1;
		}
	    }
	  
	  /* display science lines for this plot */

	  if (ksini == 1)
	    {
	      cpgsci(sini_color);
	      cpgsls(sini_style);
	      cpgslw(sini_width);
	      
	      yg[0] = 0;
	      yg[1] = yu;
	      
	      for (j = 0 ; j < SIGMAS; j++)
		{
		  siniaux = siniobs + sigma[j] * dsiniobs;
		  if ((siniaux > 0)&&(siniaux < 1.0))
		    {
		      /*first solution */
		      xg[0] = sqrt( 1 - siniaux * siniaux);
		      xg[1] = xg[0];
		      cpgline(2, xg, yg);
		      /* second solution */
		      xg[0] = -xg[0];
		      xg[1] = xg[0];
		      cpgline(2, xg, yg);
		    }
		}
	      
	      cpgsci(1);
	      cpgsls(1);
	      cpgslw(2);
	    }
	  
	  /* new variables over the old */
	  	  
	  /* xi lines */
	  
	  if (kxi == 1)
	    {
	      cpgsci(xi_color);
	      cpgsls(xi_style);
	      cpgslw(xi_width);

	      yg[0] = 0;
	      yg[1] = yu;
	      
	      for (j = 0 ; j < SIGMAS; j++)
		{
		  xiaux = xiobs + sigma[j] * dxiobs;
		  /*first solution */
		  if ((xiaux > 0)&&(xiaux < 1.0))
		    {
		      xg[0] = (1 - xiaux*xiaux) / (1 + xiaux*xiaux);
		      xg[1] = xg[0];
		      cpgline(2, xg, yg);
		      /* second solution */
		      xg[0] = -xg[0];
		      xg[1] = xg[0];
		      cpgline(2, xg, yg);
		    }
		}

	      cpgsci(1);
	      cpgsls(1);
	      cpgslw(2); 
	    }

	  /* now, the biggy: the x-dot lines */

	  if (kinetix == 1)
	    {
	      cpgsci(xdot_color);
	      cpgsls(2);
	      cpgslw(xdot_width);

	      /* first: direction of proper motion */
	      xg[0] = -1;
	      yg[0] = papmobs;
	      xg[1] = 1;
	      yg[1] = papmobs;
	      cpgline(2, xg, yg);

	      /* now: style that has been saved for the rest of the lines*/
	      cpgsls(xdot_style);

	      /* first half */
	      for (i = 0; i < (MAX2 -2); i++) y_sci2[i] = (i+1) * 180 / (1.0 * MAX2) + papmobs;
	      for (j = 0; j < SIGMAS; j++)
		{
		  xdotaux = (xdotobs + sigma[j] * dxdotobs);
		  
		  for (i = 0; i < (MAX2 -2); i++)
		    {
		      cosiaux = atan2( asini * pmobs * sin ((y_sci2[i] - papmobs) * PI / 180.0), xdotaux);
		      if ( cosiaux < 0) xy_sci2[i] = cos(cosiaux);
		      else xy_sci2[i] = -cos(cosiaux);
		    }
  		  cpgline(MAX2 - 2, xy_sci2, y_sci2);

		}
	      /* second half */
	      for (i = 0; i < (MAX2 -2); i++) y_sci2[i] = (i+1) * 180 / (1.0 * MAX2) + 180 + papmobs;
	      for (j = 0; j < SIGMAS; j++)
		{
		  xdotaux = (xdotobs + sigma[j] * dxdotobs);
		  for (i = 0; i < (MAX2 -2); i++)
		    {
		      cosiaux = atan2( asini * pmobs * sin ((y_sci2[i] - papmobs) * PI / 180.0), xdotaux);
		      if ( cosiaux < 0) xy_sci2[i] = cos(cosiaux);
		      else xy_sci2[i] = -cos(cosiaux);
		    }
                  /* do the end points here */
		  cpgline(MAX2 - 2, xy_sci2, y_sci2);
		}
	      
	      /* third half */
	      for (i = 0; i < (MAX2 -2); i++) y_sci2[i] = (i+1) * 180 / (1.0 * MAX2) -180 + papmobs;
	      for (j = 0; j < SIGMAS; j++)
		{
		  xdotaux = (xdotobs + sigma[j] * dxdotobs);
		  for (i = 0; i < (MAX2 -2); i++)
		    {
		      cosiaux = atan2( asini * pmobs * sin ((y_sci2[i] - papmobs) * PI / 180.0), xdotaux);
		      if ( cosiaux < 0) xy_sci2[i] = cos(cosiaux);
		      else xy_sci2[i] = - cos(cosiaux);
		    }
                  /* do the end points here */
		  cpgline(MAX2 - 2, xy_sci2, y_sci2);
		}

	      /* fourht half */
	      for (i = 0; i < (MAX2 -2); i++) y_sci2[i] = (i+1) * 180 / (1.0 * MAX2) - 360 + papmobs;
	      for (j = 0; j < SIGMAS; j++)
		{
		  xdotaux = (xdotobs + sigma[j] * dxdotobs);
		  for (i = 0; i < (MAX2 -2); i++)
		    {
		      cosiaux = atan2( asini * pmobs * sin ((y_sci2[i] - papmobs) * PI / 180.0), xdotaux);
		      if ( cosiaux < 0) xy_sci2[i] = cos(cosiaux);
		      else xy_sci2[i] = - cos(cosiaux);
		    }
                  /* do the end points here */
		  cpgline(MAX2 - 2, xy_sci2, y_sci2);
		}
	      cpgsci(1);
	      cpgsls(1);
	      cpgslw(2);
	    }

	  /* now display the contour plots */
	  /* if (map_number > 0)
	    {
	      map_counter = 0;
	      while (map_counter < map_number)
		{
		  if (calculate_cosiom[map_counter] == 1)
		    {
		      /* third: put in the contour levels (over the science lines)   
		      cpgsci(distribution_color[map_counter]);
		      cpgslw(distribution_width[map_counter]);

		       move one bin down 

		      tr_cosiom[map_counter][0] = tr_cosiom[map_counter][0] - tr_cosiom[map_counter][2];
		       tr_cosiom[map_counter][3] = tr_cosiom[map_counter][3] - tr_cosiom[map_counter][4];

		      for ( k = 0 ; k < nlevels ; k++)
			{	  
			  find levels
			  level = levels(n_cosi[map_counter], n_om[map_counter], cosiom_total[map_counter], D2target[k], cosiommax[map_counter], cosiom[map_counter]);
			  /* display levels
			  cpgcont(cosiom[map_counter], n_om[map_counter], n_cosi[map_counter], 1, n_om[map_counter], 1, n_cosi[map_counter], &level, -1, tr_cosiom[map_counter]);
			}
		    }
		  map_counter = map_counter + 1;
		}
	    } */

	  cpgsci(1);
	  cpgslw(2);

	  /* Now: grey-out excluded areas */
	  /* set fill color */
	  cpgsci(excluded_color);
	  /* set fill style */
	  cpgsfs(excluded_fill);
	  
	  /* first corner of polygon */
	  mx[0] = -1; my[0] = 0;
	  /* second corner of polygon */
	  mx[1] = -1; my[1] = 360;
	  /* third corner of polygon */
	  mx[2] = - sqrt(1 - sinimin*sinimin); my[2] = 360;
	  /* fourth corner of polygon */
	  mx[3] = mx[2] ; my[3] = 0;
	  
	  cpgpoly(4, mx, my);
	  
	  /* second excluded area */

	  /* first corner of polygon */
	  mx[0] = 1; my[0] = 0;
	  /* second corner of polygon */
	  mx[1] = 1; my[1] = 360;
	  /* third corner of polygon */
	  mx[2] = sqrt(1 - sinimin*sinimin); my[2] = 360;
	  /* fourth corner of polygon */
	  mx[3] = mx[2] ; my[3] = 0;
	  
	  cpgpoly(4, mx, my);
	  
	  
	  /* display i = 90 degrees */
	  cpgsci(excluded_color);
	  cpgsls(1);
	  cpgslw(2);

	  xg[0] = 0;
	  yg[0] = 0;
	  xg[1] = 0;
	  yg[1] = 360;
	  cpgline(2, xg, yg);

	  
	  /* change color back to black */
	  cpgsci(1);
	  cpgsls(1);
	  cpgslw(2);
  	  
	  /* finally: if this is a simulation, mark the values used for the simulation */
      
          ksimcount = 0;
	  while (ksimcount < ksim)
	    {
	      cpgsci(scolor[ksimcount]);
	      
	      /* need a bit of theory here */
	      /* first calculate simulated xi from sini used in simulation */
	      xiaux = ssini[ksimcount] / (1 + sqrt( 1 - ssini[ksimcount]*ssini[ksimcount]));
	      
	      /*second: calculate h3 */
	      xg[0] = sm2[ksimcount] * TSUN * pow(xiaux, 3);
	      
	      /* third: calculate h4 */
	      yg[0] = xg[0] * xiaux;
	      
	      /*plot point */
	      cpgpt(1, xg, yg, ssymbol[ksimcount]);
	      
	      ksimcount = ksimcount + 1;
	      
	    }

	  /* *************************************************************************** */

	  /* almost last thing: box and labels */
	  cpgslw(2);
	  cpgsci(1);
	  cpgsls(1);

	  cpgbox("BCTNS1", 0.0, 0, "BCTNS1", 90.0, 9);
	  cpglab("Cos i", "\\(0550) (\\u\\(0200)\\d)", "");
	  
	  /* *************************************************************************** */

	  /* last thing: box with alternative world coordinates, if this is the only plot */

	  if (display_marginal == 0)
	    {
	      cpgmtxt("T", 2.0, 0.5, 0.5, "Orbital Inclination (\\(2218))");
	      
	      cpgswin(180, 0, yl, yu);
	      cpgbox("CTMS1", 90, 1, "", 0.0, 0);
	      
	      /* mark 60 degrees */
	      
	      cpgswin(150, 30, yl, yu);
	      cpgbox("CTMS1", 60, 1, "", 0.0, 0);
	    }
	}

      /* *********************************************************************************************** */

      /* let's now do the lateral windows for these cases */

      if (display_marginal == 1)
	{
	  /* Now: display cosi distribution */
	  
	  cpgsci(1);
	  cpgsls(1);
	  cpgsch(0.8);
	  
	  /* Open first window and print first distribution */
	  
	  if ((display_type == 1) || (display_type == 12) || (display_type == 5))
	    {

	      /* in any case, we're writing in the first box */
	      cpgvsiz(bxl, bxu, byu + GAP, byu + GAP + PMARGIN);
	      cpgswin(xl, xu, 0, 1.1);

	      /* we're going to display cosi in this window */
	      
	      /* display 1-sigma predictions of tempo2 */
	      

	      /* if there are kinetic constraints on cos i, display them here as well, but only in a conventional plot
                 (with minimum cosine larger than zero) */
	      

	      /* if we're not in a cosi - om plot, display the limits derived from x-dot */
	      if ((display_type == 1) || (display_type == 12))
		if (kinetix == 1)
		  {	    
		  
	            /* display just lines 	  
		  
		  
		    cpgsci(xdot_color);
		    cpgsls(xdot_style);
		    cpgslw(xdot_width); 
		    
		    yg[0] = 0;
		    yg[1] = 1.1;
		    
		    xg[0] = cosimin;
		    xg[1] = cosimin;
		    cpgline(2, xg, yg);
		    
		    xg[0] = -cosimin;
		    xg[1] = -cosimin;
		    cpgline(2, xg, yg); */
		    
		    

	            /* Now: grey-out excluded areas */
	            /* set fill color */
	            cpgsci(xdot_color);
	            /* set fill style */
	            cpgsfs(excluded_fill);
	  
	            /* first corner of polygon */
	            mx[0] = -cosimin; my[0] = 0;
	            /* second corner of polygon */
	            mx[1] = -cosimin; my[1] = 1.1;
	            /* third corner of polygon */
	            mx[2] = cosimin; my[2] = 1.1;
	            /* fourth corner of polygon */
	            mx[3] = cosimin ; my[3] = 0;
	  
	            cpgpoly(4, mx, my);
	  

		    cpgsci(1);
		    cpgsls(1);
		    cpgslw(2);	    
		};
		
	      if (ksini == 1)
		{	      
		  cpgsci(sini_color);
		  cpgsls(sini_style);
		  cpgslw(sini_width); 
		  
		  yg[0] = 0;
		  yg[1] = 1.1;
		  
		  for (j = 0; j < SIGMAS; j++)
		    {
		      siniaux = siniobs + sigma[j]* dsiniobs;
		      if ((siniaux > 0 ) && (siniaux < 1))
			{
			  cosiaux =  sqrt(1 - siniaux * siniaux);
			  xg[0] = cosiaux;
			  xg[1] = cosiaux;
			  cpgline(2, xg, yg);

			  /* let's do the negative of this */
			  xg[0] = -cosiaux;
			  xg[1] = -cosiaux; 
			  cpgline(2, xg, yg);
			}
		    }
		  cpgsci(1);
		  cpgsls(1);
		  cpgslw(2);
		};
	      
	      
	      if (kxi == 1)
		{
		  cpgsci(xi_color);
		  cpgsls(xi_style);
		  cpgslw(xi_width); 
		  
		  yg[0] = 0;
		  yg[1] = 1.1;
		  
		  for (j = 0 ; j < SIGMAS; j++)
		    {

		      xiaux = xiobs + sigma[j] * dxiobs;
		      if ((xiaux > 0) && (xiaux < 1))
			{
			  cosiaux = (1 - xiaux*xiaux)/(1 + xiaux*xiaux);
			  xg[0] = cosiaux;
			  xg[1] = cosiaux; 
			  cpgline(2, xg, yg);
			  /* let's do the negative of this */
			  xg[0] = -cosiaux;
			  xg[1] = -cosiaux; 
			  cpgline(2, xg, yg);
			}
		    }
		  cpgsci(1);
		  cpgsls(1);
		  cpgslw(2);
		}
	      
              /* the pdfs... */
	      
	      if (map_number > 0) 
		{
		  map_counter = 0;
		  
		  /* this side window only opens if we have distributions to project */
		  while (map_counter < map_number)
		    {
		      
		      if ((calculate_cosim2[map_counter] == 1)||(calculate_cosiom[map_counter] == 1))
			{
			  
			  /* display it */
			  cpgsci(distribution_color[map_counter]);
			  cpgslw(distribution_width[map_counter]);
			  
			  cpgline(n_cosi[map_counter], cosix[map_counter], cosi[map_counter]);
			  
			  /* display the error lines */
			  if (display_marginal_percentiles == 1)
			  for (j = 0; j < D1NCONT; j++)
			    {
			      xg[0] = D1_cosi[map_counter][j];
			      xg[1] = D1_cosi[map_counter][j];
			      yg[0] = 0;
			      yg[1] = D1_cosiY[map_counter][j];
			      cpgline(2, xg, yg);
			    }
			}
		      map_counter = map_counter + 1;
		    }
                }

	      /* display i = 90 degrees */
	      cpgsci(excluded_color);
	      cpgsls(1);
	      cpgslw(2);
	      
	      xg[0] = 0;
	      yg[0] = 0;
	      xg[1] = 0;
	      yg[1] = 1.1;
	      cpgline(2, xg, yg);

	      cpgsci(1);
	      cpgslw(2);

	      /* Now: grey-out excluded areas */
	      /* set fill color */
	      cpgsci(excluded_color);
	      /* set fill style */
	      cpgsfs(excluded_fill);
	  
	      /* first corner of polygon */
	      mx[0] = -1; my[0] = 0;
	      /* second corner of polygon */
	      mx[1] = -1; my[1] = 1.1;
	      /* third corner of polygon */
	      mx[2] = - sqrt(1 - sinimin*sinimin); my[2] = 1.1;
	      /* fourth corner of polygon */
	      mx[3] = mx[2] ; my[3] = 0;
	  
	      cpgpoly(4, mx, my);
	  
	      /* second excluded area */

	      /* first corner of polygon */
	      mx[0] = 1; my[0] = 0;
	      /* second corner of polygon */
	      mx[1] = 1; my[1] = 1.1;
	      /* third corner of polygon */
	      mx[2] = sqrt(1 - sinimin*sinimin); my[2] = 1.1;
	      /* fourth corner of polygon */
	      mx[3] = mx[2] ; my[3] = 0;
	  
	      cpgpoly(4, mx, my);

	      cpgslw(2);
	      cpgsci(1);
	      cpgsls(1);

	      /* make the box */

	      cpgbox("BCTS", 0.0, 0, "BCTSN", 1.0, 10);
	      cpglab("", "Probability density", "");
	      /* now, add the outer marks on orbital inclination */

	      if (display_type == 5)
		{
		  /* back to character size for the main plot */
		  cpgsch(1.2);
		  
		  cpgmtxt("T", 2.0, 0.5, 0.5, "Orbital Inclination (\\(2218))");
		  
		  cpgswin(180, 0, 0, 1.1);
		  cpgbox("CTMS1", 90, 1, "", 0.0, 0);
		  
		  /* mark 60 degrees */
		  
		  cpgswin(150, 30, 0, 1.1);
		  cpgbox("CTMS1", 60, 1, "", 0.0, 0);
		  
		  cpgsch(0.8);
		}
	    }

	  /* *************************************************************************** */

	  if ((display_type == 2) || (display_type == 12))
	    {
	      /* now, display pulsar mass distribution */
	      /* in first case, this is the first horizontal box */
	      if (display_type == 2) cpgvsiz(bxl, bxu, byu + GAP, byu + GAP + PMARGIN);
	      /* in second case, this is the second horizontal box */	      
	      if (display_type == 12) cpgvsiz(bxu + GAP, bxu + GAP + SIZE1, byu + GAP, byu + GAP + PMARGIN);

	      cpgswin(xpl, xpu, 0, 1.1);

	      /* display sigma lines predicted by tempo2 */
	      
	      if (kmp == 1)
		{
		  cpgsci(mp_color);
	          cpgsls(mp_style);
	          cpgslw(mp_width); 
		  
		  yg[0] = 0;
		  yg[1] = 1.1;
		  
		  for (j = 0; j < SIGMAS ; j++)
		    {
		      m1aux = mpobs+sigma[j]*dmpobs;
		      xg[0] = m1aux;
		      xg[1] = m1aux;
		      cpgline(2, xg, yg);
		    }
		  
		  cpgsls(1);
		  cpgsci(1);
		  cpgslw(2);
		}

	      map_counter = 0;
	      
	      /* this side window only opens if we have distributions to project */
	      while (map_counter < map_number)
		{
		  
		  if (calculate_m1m2[map_counter] == 1)
		    {
		      		      
		      /* display it */
		      cpgsci(distribution_color[map_counter]);
		      cpgslw(distribution_width[map_counter]);
		      
		      cpgline(MPBIN, m1x[map_counter], m1[map_counter]);
		      
		      /* display the error lines */
		      if (display_marginal_percentiles == 1)
		      for (j = 0; j < D1NCONT; j++)
			{
			  xg[0] = D1_m1[map_counter][j];
			  xg[1] = D1_m1[map_counter][j];
			  yg[0] = 0;
			  yg[1] = D1_m1Y[map_counter][j];
			  cpgline(2, xg, yg);
			}
		    }
		  map_counter = map_counter + 1;

		}
	      cpgsci(1);
	      cpgslw(2);
	      cpgsls(1);

	      cpgbox("BCTS", 0.0, 0, "BCTS", 1.0, 10);
	      if ((display_type == 2)) cpglab("", "Probability density", "");

	    }

	  /* *************************************************************************** */

	  if ((display_type == 1) || (display_type == 2) || (display_type == 12))
	    {
	      /* display mc distribution - on its side, always */
 
	      if (display_type == 12) cpgvsiz(bxu + GAP + SIZE1 + GAP, bxu + GAP + SIZE1 + GAP + PMARGIN, byl, byu);
	      else cpgvsiz(bxu + GAP, bxu + GAP + PMARGIN, byl, byu);
	      
	      cpgswin(1.1, 0, yl, yu);
	      
	      /* display sigma lines predicted by tempo2 */
	      
	      if (kmc == 1)
		{
		  
		  cpgsci(mc_color);
		  cpgsls(mc_style);
		  cpgslw(mc_width); 
		  
		  /* notice - x-y coordinates changed their roles, because this plot is rotated 90 degrees */
		  xg[0] = 0;
		  xg[1] = 1.1;
		  
		  for (j = 0; j < SIGMAS; j++)
		    {
		      m2aux = mcobs + sigma[j]*dmcobs;
		      yg[0] = m2aux;
		      yg[1] = m2aux;
		      cpgline(2, xg, yg);
		    }	      
		  cpgsls(1);
		  cpgslw(2);
		  cpgsci(1);
		  
		}

	      map_counter = 0;
	      
	      /* this side window only opens if we have distributions to project */
	      while (map_counter < map_number)
		{
		  
		  if ((calculate_cosim2[map_counter] == 1)||(calculate_m1m2[map_counter] ==1))
		    {
		      		      
		      /* display it */
		      cpgsci(distribution_color[map_counter]);
		      cpgslw(distribution_width[map_counter]);
		      
		      if (map_type[map_counter] == 5)
		      	 cpgline(n_m2[map_counter], m1[map_counter], m2x[map_counter]);
		      else
		         cpgline(n_m2[map_counter], m2[map_counter], m2x[map_counter]);
		      
		      /* display the error lines */
	              if (display_marginal_percentiles == 1)
		      for (j = 0; j < D1NCONT; j++)
			{
			  xg[0] = D1_m2[map_counter][j];
			  xg[1] = D1_m2[map_counter][j];
			  yg[0] = 0;
			  yg[1] = D1_m2Y[map_counter][j];
			  cpgline(2, yg, xg);
			}
		    }
		  map_counter = map_counter + 1;
		}

	      cpgslw(2);
	      cpgsci(1);
	      cpgsls(1);
	      
	      cpgbox("BCTSN", 1.0, 10, "BCTS", 0.0, 0);
	      cpglab("Probability density", "", "");

	    }

	  /* *************************************************************************** */

	  if (display_type == 5)
	    {
	      /* display omega distribution - on its side */
	      cpgvsiz(bxu + GAP, bxu + GAP + PMARGIN, byl, byu); 
	      cpgswin(1.1, 0, yl, yu);
	      
	      /* display sigma lines predicted by tempo2 */
	      
	      map_counter = 0;
	      
	      /* this side window only opens if we have distributions to project */
	      while (map_counter < map_number)
		{
		  
		  if (calculate_cosiom[map_counter] == 1)
		    {
		      		      
		      /* display it */
		      cpgsci(distribution_color[map_counter]);
		      cpgslw(distribution_width[map_counter]);
		      
		      cpgline(n_om[map_counter], om[map_counter], omx[map_counter]);
		      
		      /* display the error lines */
		      if (display_marginal_percentiles == 1)
		      for (j = 0; j < D1NCONT; j++)
			{
			  xg[0] = D1_om[map_counter][j];
			  xg[1] = D1_om[map_counter][j];
			  yg[0] = 0;
			  yg[1] = D1_omY[map_counter][j];
			  cpgline(2, yg, xg);
			}
		    }
		  map_counter = map_counter + 1;
		}

	      if (kinetix == 1)
		{
		  cpgsci(xdot_color);
		  cpgsls(2);
		  cpgslw(xdot_width);
		  
		  /* first: direction of proper motion */
		  xg[0] = 0;
		  yg[0] = papmobs;
		  xg[1] = 1.1;
		  yg[1] = papmobs;
		  cpgline(2, xg, yg);
		}

	      cpgslw(2);
	      cpgsci(1);
	      cpgsls(1);
	      
	      /* frame the plot */
	      cpgbox("BCTSN", 1.0, 10, "BCTS", 90.0, 9);
	      cpglab("Probability density", "", "");

	    }
	}

      /* close this graph */
      cpgend();
    }


  /* *************************************************************************************** */
    
  /* in the case of display 5, make second plot for component masses */
  
  if ((display_type == 5)&&(map_number > 0)&&(map_type[map_counter] == 5))
    {
      map_counter = 0;
      while (map_counter < map_number)
	{
	  
	  if (calculate_cosiom[map_counter] == 1)
	    {
	      
	      printf("\n\n Please introduce the file_name.ps/graphics device below for the mass distributions.\n");
	      if(cpgbeg(0, "?", 1, 1) != 1)
		return EXIT_FAILURE;
	      
	      /* we have all that is needed to start opening graphics */
	      cpgsch(1.2);
	      cpgsci(1);
	      
	      bxl = MARGIN;
	      byl = MARGIN;
	      bxu = MARGIN + SIZE0;
	      byu = MARGIN + SIZE0;
	      
	      cpgvsiz(bxl, bxu, byl, byu);
	      cpgswin(xpl, xpu, 0, 1.1);
	      
	      /* now display the pdf for the pulsar mass */
	      
	      map_counter = 0;
	      while (map_counter < map_number)
		{
		  /* display it */
		  cpgsci(distribution_color[map_counter]);
		  cpgslw(distribution_width[map_counter]);
		      
		  cpgline(MPBIN, m1x[map_counter], m1[map_counter]);
		  
		  /* display the error lines */
		  if (display_marginal_percentiles == 1)
		  for (j = 0; j < D1NCONT; j++)
		    {
		      xg[0] = D1_m1[map_counter][j];
		      xg[1] = D1_m1[map_counter][j];
		      yg[0] = 0;
		      yg[1] = D1_m1Y[map_counter][j];
		      cpgline(2, xg, yg);
		    }
		  
		  map_counter = map_counter +1 ;
		}
		
              /* frame the box */
	      cpgsci(1);
	      cpgsls(1);
	      cpgslw(2);
	      cpgbox("BTNSN", 1.0, 0, "BTNS", 0.0, 0);
	      cpgmtxt("L", 3.0, 0.5, 0.5, "Probability density");
	      cpgmtxt("B", 3.0, 0.5, 0.5, "Pulsar Mass (M\\d\\(2281)\\u)");

              /* box with alternative world coordinates */

	      cpgswin(mtobs, mc_min, 0, 1.1);
	      cpgbox("CMTS1", 0.0, 0, "CTMS1", 0.0, 0);
	      cpgmtxt("T", 2.0, 0.5, 0.5, "Companion mass (M\\d\\(2281)\\u)");

		
	    }
	}
      cpgend();
    }

  return EXIT_SUCCESS;

  /* *************************************************************************************** */

}

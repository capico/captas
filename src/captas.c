/*

 Computer Aided Pressure Transient Analysis (and Simulation?) � CAPTA(S?)
     Copyright (C) 2013  Carlos E. Pico

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <malloc.h>
#include <float.h>
#include <getopt.h>
#include <gsl/gsl_cdf.h>

#include "cminpack/cminpack.h"
#include "iniparser/iniparser.h"
#include "gnuplot_i/gnuplot_i.h"
#include "stehfest/stehfest.h"
#include "utils/utils.h"

#include "captas.h"

#define VERSION "0.0"


/*****************************************************************************/
int main(int argc, char *argv[])
{
	/* processing times */
	clock_t
		start,
		finish;

	double duration = 0.0;

	int
		i = 0,
		j = 0,
		ans = 0,
		tmp;

	/* parameters for the fitting function */
	int
		info  = 0,  // exit code
		njev  = 0,  // calls to jacobian
		nfev  = 0,  // calls to residual
		mode  = 0,  // scaling
		totalfcneval = 0,
		totaljaceval = 0;

	double
		tol,
		sqrttol,
		fnorm;

	/* unknown parameters, residuals and jacobian */
	double
		*x,		// unknown parameters
		*ci,    // confidence intervals
		*fvec;	// residual

	double
        **fjac, // jacobian
        **corr; // correlation matrix


	/* data and parameters size */
	int
		m  = 0, /* number of observations */
		n  = 0; /* number of parameters */

	int nevents = 0; // number of production events

	/* x-y tables for pressures and derivatives */
	lmdatatype data;
	lmdatatype dataSol;

	/* */
	modelparameters par;

	double **table; /* table for external (file) data */

	/* gnuplot handles */
	gnuplot_ctrl
		*hcartesian,
		*hsemilog,
		*hloglog;
	/*************************************************************************/
	/*************************************************************************/
	/*************************************************************************/

	parse_args(argc, argv, &par);
	init_parameters(&par);
	read_inifile(&par);
	set_parameters(&par);
	/*************************************************************************/

	m       = par.m;
	n       = par.n;
	nevents = par.nevents;

    /* memory allocation */
    /* observations */
    if((data.t    = (double*)calloc(m, sizeof(double))) == NULL || // time
       (data.p    = (double*)calloc(m, sizeof(double))) == NULL || // pressure
       (data.dp   = (double*)calloc(m, sizeof(double))) == NULL || // delta p = (p@dt=0 - p)
       (data.dt   = (double*)calloc(m, sizeof(double))) == NULL || // elapsed time
       (data.dte  = (double*)calloc(m, sizeof(double))) == NULL || // multirate equivalent-time
       (data.derl = (double*)calloc(m, sizeof(double))) == NULL || // logarithmic derivative
       (data.derg = (double*)calloc(m, sizeof(double))) == NULL){  // derivative group dp/(2*derl)

        printf("\n memory allocation failure \n");
        exit(EXIT_FAILURE);
    }

	/* solution */
    if((dataSol.t    = (double*)calloc(m, sizeof(double))) == NULL ||
       (dataSol.p    = (double*)calloc(m, sizeof(double))) == NULL ||
       (dataSol.dp   = (double*)calloc(m, sizeof(double))) == NULL ||
       (dataSol.dt   = (double*)calloc(m, sizeof(double))) == NULL ||
       (dataSol.dte  = (double*)calloc(m, sizeof(double))) == NULL ||
       (dataSol.derl = (double*)calloc(m, sizeof(double))) == NULL ||
       (dataSol.derg = (double*)calloc(m, sizeof(double))) == NULL){

        printf("\n memory allocation failure \n");
        exit(EXIT_FAILURE);
    }

	/* nonlinear regression parameters, residuals, etc */
    if((x    = (double*)calloc(n, sizeof(double))) == NULL ||		// unknown parameters
       (ci   = (double*)calloc(n, sizeof(double))) == NULL ||		// confidence intervals
       (fvec = (double*)calloc(m, sizeof(double))) == NULL ||		// residuals
       (fjac = (double **)calloc(n, sizeof(double *))) == NULL ||	// jacobian
       (corr = (double **)calloc(n, sizeof(double *))) == NULL){	// correlation matrix

	    printf("\n memory allocation failure \n");
        exit(EXIT_FAILURE);
	}

	/* matrix for jacobian and correlation */
	for(j = 0; j < n; j++) {
		if((fjac[j] = (double *)calloc(m, sizeof(double))) == NULL ||
           (corr[j] = (double *)calloc(n, sizeof(double))) == NULL) {
			printf("\n memory allocation failure \n");
			exit(EXIT_FAILURE);
		}
	}

	/* convergence criteria */
	tol     = sqrt(MACHEPS);
	sqrttol = sqrt(tol);

	/****************** pressure data ****************/
	/*************************************************/
	table = (double **)calloc(m, sizeof(double *));
	for(i = 0; i < m; i++) {
		table[i] = (double *)calloc(2, sizeof(double));
	}
	rTableFile(par.pressfile, table, m, 2);
	/**************************************************/
	for(i = 0; i < m; i++) {
		data.t[i]  = table[i][0];           // time
		data.p[i]  = table[i][1];           // pressure
	}
	for(i = 0; i < m; i++) {
		free(table[i]);
	}
	free(table);

	/***************** flow rate history *************/
	/*************************************************/
	table = (double **)calloc(nevents, sizeof(double *));
	for(i = 0; i < nevents; i++) {
		table[i] = (double *)calloc(2, sizeof(double));
	}
	rTableFile(par.ratefile, table, nevents, 2);
	/**************************************************/
	for(i = 0; i < nevents; i++) {
		par.tps[i]  = table[i][0];
		par.qBps[i] = table[i][1]*par.B;
	}
	for(i = 0; i < nevents; i++) {
		free(table[i]);
	}
	free(table);

	/********* copy of the abscissas (time) ***********/
	for(i = 0; i < m; i++) {
		dataSol.t[i] = data.t[i];
	}

	/*************************************************************************/
	/**************************** plot windows *******************************/
	/*************************************************************************/
	if(par.plots){
		hcartesian = gnuplot_init();
		if(par.testtype != MULTIRATE){
		    hsemilog   = gnuplot_init();
		    hloglog    = gnuplot_init();
		}

		if(hcartesian == NULL){
			printf("\n");
			printf("GNUPLOT error\n");
			printf("GNUPLOT is not available in your path.\n");

			par.plots = 0;
		}
	}

	if(par.plots){
		printf("*** plots - data and model with initial guess ***\n");

		gnuplot_set_xlabel(hcartesian, "t [hr]");
		gnuplot_cmd(hcartesian, "set grid");
        gnuplot_cmd(hcartesian, "set terminal wxt enhanced background '#fff6af'");
		if(par.testtype != MULTIRATE){
            gnuplot_set_xlabel(hsemilog, "{/Symbol D}t_e [hr]");
            gnuplot_set_xlabel(hloglog, "{/Symbol D}t [hr]");

            gnuplot_cmd(hsemilog, "set logscale x");
            gnuplot_cmd(hloglog, "set logscale xy");

            gnuplot_cmd(hsemilog, "set grid");
            gnuplot_cmd(hloglog, "set grid");

            gnuplot_cmd(hsemilog,   "set terminal wxt enhanced background '#fff6af'");
            gnuplot_cmd(hloglog,    "set terminal wxt enhanced background '#fff6af'");
		}

		if(par.units == ANP){
			gnuplot_set_ylabel(hcartesian, "p [kgf/cm^2]");
            if(par.testtype != MULTIRATE){
                gnuplot_set_ylabel(hsemilog, "p [kgf/cm^2]");
                gnuplot_set_ylabel(hloglog, "{/Symbol D}p  and  {/Symbol D}p' [kgf/cm^2]");
            }
		}
		else{
			gnuplot_set_ylabel(hcartesian, "p [psi]");
            if(par.testtype != MULTIRATE){
                gnuplot_set_ylabel(hsemilog, "p [psi]");
                gnuplot_set_ylabel(hloglog, "{/Symbol D}p  and  {/Symbol D}p' [psi]");
            }
		}
	}
	/*************************************************************************/


	/*************************************************************************/
	/********************* wellbore pressure (data) **************************/
	/*************************************************************************/
	if(par.plots){

		gnuplot_setstyle(hcartesian, "points");
		gnuplot_plot_xy(hcartesian, data.t, data.p, m, "p (data)");

        if(par.testtype != MULTIRATE){
            calc_plotting_abscissas(&par, &data);
            calc_pressure_drop(&par, &data);
            calc_log_derivative(&par, &data);

            gnuplot_setstyle(hsemilog,   "points");
            gnuplot_setstyle(hloglog,    "points");

            gnuplot_plot_xy(hsemilog, data.dte, data.p, m, "p (data)");
            gnuplot_plot_xy(hloglog, data.dt, data.dp, m, "{/Symbol D}p (data)");
            gnuplot_plot_xy(hloglog, data.dt, data.derl,  m, "{/Symbol D}p' (data)");
            //gnuplot_plot_xy(hloglog, data.dt, data.derg,  m, "{/Symbol D}p / 2*{/Symbol D}p' (data)");
        }
	}
	/*************************************************************************/


	/*************************************************************************/
	/********************** wellbore pressure (model) ************************/
	/*************************************************************************/
	if(par.plots){

		calc_pressure(&par, &dataSol);

		gnuplot_setstyle(hcartesian, "lines");
		gnuplot_plot_xy(hcartesian, dataSol.t, dataSol.p, m, "p (model)");

		if(par.testtype != MULTIRATE){
            calc_plotting_abscissas(&par, &dataSol);
            calc_pressure_drop(&par, &dataSol);
            calc_log_derivative(&par, &dataSol);

            gnuplot_setstyle(hsemilog, "lines");
            gnuplot_setstyle(hloglog, "lines");

            gnuplot_plot_xy(hsemilog, dataSol.dte, dataSol.p, m, "p (model)");
            gnuplot_plot_xy(hloglog, dataSol.dt, dataSol.dp, m, "{/Symbol D}p (model)");
            gnuplot_plot_xy(hloglog, dataSol.dt, dataSol.derl, m, "{/Symbol D}p' (model)");
            //gnuplot_plot_xy(hloglog, dataSol.dt, dataSol.derg,  m, "{/Symbol D}p / 2*{/Symbol D}p' (model)");
		}
	}
	pressEnterToContinue();
	/*************************************************************************/


	/*************************************************************************/
	/*************** fitting by the Levenberg-Marquardt algorithm ************/
	/*************************************************************************/
	setxtomodel(x, n, &par);
	totalfcneval = 0;
	totaljaceval = 0;
	mode = 1;
	do{
		start = clock();

		/* solve the nonlinear least squares fitting problem */
		lm_nlsf(resjacprn, &data, &par, m, n, x, ci, corr, fvec, fjac, tol,
          &info, &nfev, &njev, &fnorm, mode);

		finish = clock();
		duration += (double)(finish-start)/CLOCKS_PER_SEC;

		totalfcneval += nfev;
		totaljaceval += njev;
		printf("\n");
		printf("Time: %f\n", duration);
		printf("Final L2 norm: %.16E \n", fnorm);
		printf("Exit message: %s \n", lm_nlsf_msg[info]);
		printf("Number of calls to function: %d \n", totalfcneval);
		printf("Number of calls to jacobian: %d \n", totaljaceval);
		printf("\n");

		setmodeltox(x, n, &par);

		print_par(stderr, x, ci, corr, n, &par);
		printf("\n\n");

		/*********************** wellbore pressure estimative ****************/
		if(par.plots){
			calc_pressure_drop(&par, &data);

			gnuplot_resetplot(hcartesian);
			gnuplot_setstyle(hcartesian, "points");
			gnuplot_plot_xy(hcartesian, data.t, data.p, m, "p (data)");

			if(par.testtype != MULTIRATE){
                gnuplot_resetplot(hsemilog);
                gnuplot_resetplot(hloglog);

                gnuplot_setstyle(hsemilog, "points");
                gnuplot_setstyle(hloglog, "points");

                gnuplot_plot_xy(hsemilog, data.dte, data.p, m, "p (data)");
                gnuplot_plot_xy(hloglog, data.dt, data.dp, m, "{/Symbol D}p (data)");
                gnuplot_plot_xy(hloglog, data.dt, data.derl,  m, "{/Symbol D}p' (data)");
                //gnuplot_plot_xy(hloglog, data.dt, data.derg,  m, "{/Symbol D}p / 2*{/Symbol D}p' (data)");
			}

			calc_pressure(&par, &dataSol);

			gnuplot_setstyle(hcartesian, "lines");
			gnuplot_plot_xy(hcartesian, dataSol.t, dataSol.p, m, "p (model)");

			if(par.testtype != MULTIRATE){
                calc_pressure_drop(&par, &dataSol);
                calc_log_derivative(&par, &dataSol);

                gnuplot_setstyle(hsemilog, "lines");
                gnuplot_setstyle(hloglog, "lines");

                gnuplot_plot_xy(hsemilog, dataSol.dte, dataSol.p, m, "p (model)");
                gnuplot_plot_xy(hloglog, dataSol.dt, dataSol.dp, m, "{/Symbol D}p (model)");
                gnuplot_plot_xy(hloglog, dataSol.dt, dataSol.derl, m, "{/Symbol D}p' (model)");
                //gnuplot_plot_xy(hloglog, dataSol.dt, dataSol.derg,  m, "{/Symbol D}p / 2*{/Symbol D}p' (model)");
			}
		}
		/*************************************************************************/

		printf("Acceptable Estimative? (1 = yes, 2 = change mode (mode = %d), otherwise = no) ", mode);
		if (scanf("%d", &ans) == 0){
            scanf("%*s");
            ans = 0;
		}
		//info = 0;

		if(ans == 2)
			mode = (mode == 1) ? 2 : 1;
	}
	while(ans != 1);
	/*************************************************************************/
	/*************************************************************************/

	//pressEnterToContinue();

	/*************************************************************************/
	write_outfile(&par, &dataSol, x, ci, corr);
	/*************************************************************************/

	/*************************************************************************/
	if(par.plots){
		gnuplot_close(hcartesian);
		if(par.testtype != MULTIRATE){
			gnuplot_close(hsemilog);
			gnuplot_close(hloglog);
		}
	}
	/*************************************************************************/

	return 0;
}
/*===========================================================================*/
/*****************************************************************************/
/*===========================================================================*/





/**
Process the command line arguments
*/
void parse_args(int argc, char *argv[], modelparameters *par)
{
	int option_index = 0;

	while (( option_index = getopt(argc, argv, "hv")) != -1){
		switch(option_index){
		case 'h':
			help_msg();
			exit(EXIT_SUCCESS);
			break;
		case 'v':
			print_version();
			exit(EXIT_SUCCESS);
			break;
		case '?':         /* Invalid option or missing argument returns '?' */
		default:
			help_msg();
			exit(EXIT_FAILURE); /* helpmsg() exits with error */
		}
	}

	if(argc < 2){
		printf("<inifile> not specified...\n\n");
		help_msg();
		exit(EXIT_FAILURE);
	}
	else{
		strcpy(par->inifile, argv[1]);
		printf("reading inifile: %s\n\n", par->inifile);
	}

	return;
}
/*****************************************************************************/



void help_msg(void)
{
	printf("\nUsage: capta [options] [ini_file]\n");
	printf("\n\
		   Options:\n\
		   -h,      show this help and exit\n\
		   -v,      print version and exit\n");
	return;
}
/*****************************************************************************/



void print_version(void)
{
	printf("captas version %s\n", VERSION);
	return;
}
/*****************************************************************************/



double f_default(void)
{
    printf("\n you are trying to use a model or an analytic derivative that \
           is not available.\n");
    exit(EXIT_FAILURE);
}



/**
this function initializes the model parameters structure with the pointers to
the dwf models and the corresponding analytical derivatives, and some
parameters to the default values.
*/
void init_parameters(modelparameters *par)
{
	int i, j;

	for(i = 0; i < NPARAMETERS; i++){
		par->rpsta[i] = OFF;
		par->rpjac[i] = CENTRAL;
		par->rpdx[i]  = DX_DEFAULT;
	}

	par->parnames[PERMEABILITY]            = "k";
	par->parnames[SKIN_FACTOR]             = "S";
	par->parnames[WELLBORE_STORAGE]        = "C";
	par->parnames[INITIAL_PRESSURE]        = "pi";
	par->parnames[EXTERNAL_RADIUS]         = "re";
	par->parnames[DISTANCE_TO_FAULT]       = "L";
	par->parnames[DISTANCE_TO_FAULT_1]     = "w1";
	par->parnames[DISTANCE_TO_FAULT_2]     = "w2";
	par->parnames[DISTANCE_TO_FAULT_1_X]   = "w1x";
	par->parnames[DISTANCE_TO_FAULT_2_X]   = "w2x";
	par->parnames[DISTANCE_TO_FAULT_1_Y]   = "w1y";
	par->parnames[DISTANCE_TO_FAULT_2_Y]   = "w2y";
	par->parnames[OMEGA]                   = "omega";
	par->parnames[LAMBDA]                  = "lambda";
	par->parnames[FRACTURE_HALF_LENGTH]    = "xf";
	par->parnames[FRACTURE_CONDUCTIVITY]   = "fc";
	par->parnames[PERMEABILITY_RATIO]      = "krat";
	par->parnames[PENETRATION_RATIO]       = "b";
	par->parnames[MIDPOINT_ELEVATION]      = "zw",
	par->parnames[EFFECTIVE_HEAT_CAPACITY] = "cpt",
	par->parnames[INITIAL_TEMPERATURE]     = "Ti",
	par->parnames[COEFFICIENT_JOULE_THOMSON]    = "ejt";

	/******************** pointers to delta_pwf functions ********************/
	for(i = 0; i <  NPMODELS; i++){
		par->dpwffcn[i] = f_default;
	}
	par->dpwffcn[PWF]       = dpwf;
	par->dpwffcn[PWFT]      = dpwf;
	par->dpwffcn[PWFNFOB]   = dpwfnfob;
	par->dpwffcn[PWFCPOB]   = dpwfcpob;
	par->dpwffcn[PWFF]      = dpwff;
	par->dpwffcn[PWFC]      = dpwfc;
	par->dpwffcn[PWFCL]     = dpwfcl;
	par->dpwffcn[PWFRECT]   = dpwfrect;
	par->dpwffcn[PWFDPPSS]  = dpwfdppss;
	par->dpwffcn[PWFDPTSL]  = dpwfdptsl;
	par->dpwffcn[PWFDPTSP]  = dpwfdptsp;
	par->dpwffcn[PWFICF]    = dpwficf;
	par->dpwffcn[PWFFCF]    = dpwffcf;
	par->dpwffcn[PWFLE]     = dpwfle;

	par->dTwffcn[TWF]       = dTwf;

	/*************** pointers to interporosity flow functions ****************/
	for(i = 0; i <  NPMODELS; i++){
		par->f[i] = f_default;
	}
	par->f[PWFDPPSS]  = fpss;
	par->f[PWFDPTSL]  = ftsl;
	par->f[PWFDPTSP]  = ftsp;

	/********* pointers to the derivatives of the delta_pwf functions ********/
	for(i = 0; i <  NMODELS; i++)
	{
		for(j = 0; j < NPARAMETERS; j++) {
			par->dr_dx[i][j] = f_default;
		}
	}
	par->dr_dx[PWF][PERMEABILITY]          = ddpwf_dk;
	par->dr_dx[PWF][WELLBORE_STORAGE]      = ddpwf_dC;
	par->dr_dx[PWF][SKIN_FACTOR]           = ddpwf_dS;

	par->dr_dx[PWFT][PERMEABILITY]         = ddpwft_dk;
	par->dr_dx[PWFT][WELLBORE_STORAGE]     = ddpwft_dC;
	par->dr_dx[PWFT][SKIN_FACTOR]          = ddpwft_dS;

	par->dr_dx[PWFF][PERMEABILITY]         = ddpwff_dk;
	par->dr_dx[PWFF][WELLBORE_STORAGE]     = ddpwff_dC;
	par->dr_dx[PWFF][SKIN_FACTOR]          = ddpwff_dS;
	par->dr_dx[PWFF][DISTANCE_TO_FAULT]    = ddpwff_dL;

	par->dr_dx[PWFC][PERMEABILITY]         = ddpwfc_dk;
	par->dr_dx[PWFC][WELLBORE_STORAGE]     = ddpwfc_dC;
	par->dr_dx[PWFC][SKIN_FACTOR]          = ddpwfc_dS;
	par->dr_dx[PWFC][DISTANCE_TO_FAULT_1]  = ddpwfc_dw1;
	par->dr_dx[PWFC][DISTANCE_TO_FAULT_2]  = ddpwfc_dw2;

	par->dr_dx[PWFCL][PERMEABILITY]        = ddpwfcl_dk;
	par->dr_dx[PWFCL][SKIN_FACTOR]         = ddpwfcl_dS;
	par->dr_dx[PWFCL][DISTANCE_TO_FAULT_1] = ddpwfcl_dw1;
	par->dr_dx[PWFCL][DISTANCE_TO_FAULT_2] = ddpwfcl_dw2;

	par->dr_dx[PWFCPOB][WELLBORE_STORAGE]  = ddpwfcpob_dC;

	par->dr_dx[PWFNFOB][PERMEABILITY]      = ddpwfnfob_dk;
	par->dr_dx[PWFNFOB][WELLBORE_STORAGE]  = ddpwfnfob_dC;
	par->dr_dx[PWFNFOB][SKIN_FACTOR]       = ddpwfnfob_dS;
	par->dr_dx[PWFNFOB][EXTERNAL_RADIUS]   = ddpwfnfob_dre;

	par->dr_dx[PWFDPPSS][WELLBORE_STORAGE] = ddpwfdppss_dC;

	par->dr_dx[PWFDPTSL][WELLBORE_STORAGE] = ddpwfdptsl_dC;

	par->dr_dx[PWFDPTSP][WELLBORE_STORAGE] = ddpwfdptsp_dC;

	par->dr_dx[PWFICF][WELLBORE_STORAGE]   = ddpwficf_dC;
	par->dr_dx[PWFICF][SKIN_FACTOR]        = ddpwficf_dS;

	par->dr_dx[PWFFCF][WELLBORE_STORAGE]   = ddpwffcf_dC;
	par->dr_dx[PWFFCF][SKIN_FACTOR]        = ddpwffcf_dS;

	par->dr_dx[PWFLE][WELLBORE_STORAGE]    = ddpwfle_dC;

	par->dr_dx[PWF][INITIAL_PRESSURE]      =
	par->dr_dx[PWFNFOB][INITIAL_PRESSURE]  =
	par->dr_dx[PWFCPOB][INITIAL_PRESSURE]  =
	par->dr_dx[PWFF][INITIAL_PRESSURE]     =
	par->dr_dx[PWFC][INITIAL_PRESSURE]     =
	par->dr_dx[PWFCL][INITIAL_PRESSURE]    =
	par->dr_dx[PWFRECT][INITIAL_PRESSURE]  =
	par->dr_dx[PWFDPPSS][INITIAL_PRESSURE] =
	par->dr_dx[PWFDPTSL][INITIAL_PRESSURE] =
	par->dr_dx[PWFDPTSP][INITIAL_PRESSURE] =
	par->dr_dx[PWFICF][INITIAL_PRESSURE]   =
	par->dr_dx[PWFFCF][INITIAL_PRESSURE]   =
	par->dr_dx[PWFLE][INITIAL_PRESSURE]    = dr_dpi;

	par->dr_dx[TWF][INITIAL_TEMPERATURE]   = dr_dTi;

	par->dr_dx[PWFT][INITIAL_PRESSURE]     = drt_dpi;
	/*************************************************************************/

	par->C1  = C1_OILFIELD;
	par->C2  = C2_OILFIELD;
	par->C3  = C3_OILFIELD;

	par->nstehfest = NSTEHFEST_DEFAULT;
	par->Lder      = LDER_DEFAULT;

	return;
}



/**
this function reads a list of parameters from a inifile and copies
the parameters to the structure "par".
*/
void read_inifile(modelparameters *par)
{
	char *str;
	char
		*strval = "Test description:",
		*strsta = "Regression parameters:rp_",
		*strjac = "Regression parameters derivatives:jac_",
		*strdx  = "Regression parameters derivatives:dx_";
	char strtmp[LENGTENTRY];
	dictionary *ini;
	int i;

	ini = iniparser_load(par->inifile);
	if(ini == NULL) {
		printf("\n cannot parse ini file: \n");
		exit(EXIT_FAILURE);
	}

	//iniparser_dump_ini(ini, stdout);

	/********************* reading data in ini file ********************/
	par->mode     = iniparser_getint(ini, "Program mode:mode",         REGRESSION);
	par->units    = iniparser_getint(ini, "Units system:units",        OILFIELD);
	par->testtype = iniparser_getint(ini, "Test description:testtype", DRAWDOWN);

	par->phi     = iniparser_getdouble(ini, "Test description:phi",    0.0);
	par->B       = iniparser_getdouble(ini, "Test description:B",      0.0);
	par->mu      = iniparser_getdouble(ini, "Test description:mu",     0.0);
	par->h       = iniparser_getdouble(ini, "Test description:h",      0.0);
	par->rw      = iniparser_getdouble(ini, "Test description:rw",     0.0);
	par->ct      = iniparser_getdouble(ini, "Test description:ct",     0.0);
	par->rhosc   = iniparser_getdouble(ini, "Test description:rhosc",  0.0);
	par->cp      = iniparser_getdouble(ini, "Test description:cp",     0.0);


	/* values */
    for(i = 0; i < NPARAMETERS; i++){
        strcat(strcpy(strtmp, strval), par->parnames[i]);
        par->rpval[i] = iniparser_getdouble(ini, strtmp, 0.0);
	}

	/* status */
    for(i = 0; i < NPARAMETERS; i++){
        strcat(strcpy(strtmp, strsta), par->parnames[i]);
        par->rpsta[i] = iniparser_getboolean(ini, strtmp, OFF);
	}

    /* derivatives */
    for(i = 0; i < NPARAMETERS; i++){
        strcat(strcpy(strtmp, strjac), par->parnames[i]);
        par->rpjac[i] = iniparser_getint(ini, strtmp, CENTRAL);
	}

	/* dx */
    for(i = 0; i < NPARAMETERS; i++){
        strcat(strcpy(strtmp, strdx), par->parnames[i]);
        par->rpdx[i] = iniparser_getdouble(ini, strtmp, DX_DEFAULT);
	}

	par->model     = iniparser_getint(ini, "Regression model:model", PWF);
	par->nstehfest = iniparser_getint(ini, "Stehfest parameters:nstehfest", NSTEHFEST_DEFAULT);
	par->Lder      = iniparser_getdouble(ini, "Smoothing parameters:Lder", LDER_DEFAULT);
	par->plots     = iniparser_getboolean(ini, "Output:plots", OFF);

	str = iniparser_getstring(ini, "Test description:pressfile", NULL);
	strcpy(par->pressfile, str);
	par->presssize = iniparser_getint(ini, "Test description:presssize", 0);

//    str = iniparser_getstring(ini, "Test description:tempfile", NULL);
//	strcpy(par->tempfile, str);
//	par->tempsize = iniparser_getint(ini, "Test description:tempsize", 0);

	str = iniparser_getstring(ini, "Test description:ratefile",  NULL);
	strcpy(par->ratefile, str);
	par->ratesize = iniparser_getint(ini, "Test description:ratesize",   0);

	str = iniparser_getstring(ini, "Output:outfile", NULL);
	strcpy(par->outfile, str);
	/*************************************************************************/

	iniparser_freedict(ini);

	return;
}
/*****************************************************************************/



/**
This function sets some parameters of the structure "par".
*/
void set_parameters(modelparameters *par)
{
	int i = 0, j = 0;

	par->nevents	= par->ratesize;
	par->m			= par->presssize;

    par->n = 0;
    for(i = 0; i < NPARAMETERS; i++){
        par->n += par->rpsta[i];
    }

	if(par->nstehfest < 4 ||
		par->nstehfest > 20 ||
		par->nstehfest%2 != 0){
			printf("Stehfest's N must be even, 4 <= N <= 20\n");
			printf("using default value N = %d\n", NSTEHFEST_DEFAULT);
			par->nstehfest = NSTEHFEST_DEFAULT;
	}
	par->v = NULL;
	par->v = stehfest_init(par->nstehfest, par->v);

    if((par->tps  = (double*)calloc(par->nevents, sizeof(double))) == NULL ||
        (par->qBps = (double*)calloc(par->nevents, sizeof(double))) == NULL ||
        (par->nps  = (int*)calloc(par->nevents, sizeof(int))) == NULL ||
        (par->partype = (int*)calloc(par->n, sizeof(int))) == NULL ||
        (par->jactype = (int*)calloc(par->n, sizeof(int))) == NULL ||
        (par->dx = (double*)calloc(par->n, sizeof(double))) == NULL){

        printf("\n memory allocation failure \n");
        exit(EXIT_FAILURE);
	}

	for(i = 0, j = 0; i < NPARAMETERS; i++){
		if(par->rpsta[i] == ON){
			par->partype[j] = i;
			par->jactype[j] = par->rpjac[i];
			par->dx[j]      = par->rpdx[i];
			j++;
		}
	}

	switch(par->units){
	case ANP:
		par->C1      = C1_ANP;
		par->C2      = C2_ANP;
		par->C3      = C3_ANP;
		break;
	case OILFIELD:
		par->C1      = C1_OILFIELD;
		par->C2      = C2_OILFIELD;
		par->C3      = C3_OILFIELD;
		break;
	default:
		par->C1      = C1_ANP;
		par->C2      = C2_ANP;
		par->C3      = C3_ANP;
	}
	/*************************************************************************/

	switch(par->testtype){
	case MULTIRATE:
	case DRAWDOWN:
	case FALLOFF:
	case MULTIRATE_DRAWDOWN:
		par->sign = +1.0;
		break;
	case BUILDUP:
	case INJECTION:
	case MULTIRATE_BUILDUP:
		par->sign = -1.0;
		break;
	default:
		par->sign = +1.0;
	}

	return;
}
/*****************************************************************************/



/**
this function calculates the elapsed time dt and the multi-rate equivalent time
dte, for the array t of the structure "data", using the production events tps
and qBs in the structure "par".
*/
void calc_plotting_abscissas(modelparameters *par, lmdatatype *data)
{
	int i, j;
	double qBnumer, qBdenom;

	for(i = 0; i < par->m; i++) {

		data->dt[i]  = data->t[i] - par->tps[par->nevents-1];
		data->dte[i] = 0.0;

		if(data->dt[i] > 0.0) {
			for(j = 0; j < par->nevents - 1; j++) {
				if(data->t[i] > par->tps[j]) {
					qBnumer = (j == 0 ? par->qBps[0] : par->qBps[j] - par->qBps[j-1]);
					qBdenom = par->qBps[par->nevents-1] - par->qBps[par->nevents-2];
					data->dte[i] += (qBnumer / qBdenom)
						* log((data->t[i] - par->tps[j])
						/ (par->tps[par->nevents-1] - par->tps[j]));
				}
			}

			if(data->t[i] > par->tps[par->nevents - 1]) {
				data->dte[i] += log(data->dt[i]);
			}

			data->dte[i] = exp(data->dte[i]);
		}
	}

	return;
}
/*****************************************************************************/



/**
this function calculates the pressure drop dp  = p(dt = 0) - p(dt), for the
values of t in the structure "data", using the production events tps and qBs
in the structure "par".
*/
void calc_pressure_drop(modelparameters *par, lmdatatype *data)
{
	int i;

	double pi, p_at_0;  /* pressure at elapsed time zero */

	pi = par->rpval[INITIAL_PRESSURE];

	switch(par->testtype) {
	case MULTIRATE:
	case DRAWDOWN:
	case INJECTION:
		p_at_0 = pi;
		break;
	case BUILDUP:
	case FALLOFF:
	case MULTIRATE_BUILDUP:
	case MULTIRATE_DRAWDOWN:
		p_at_0 = data->p[0];
		break;
	default:
		p_at_0 = pi;
	}

	for(i = 0; i < par->m; i++) {
		data->dp[i] = par->sign*(p_at_0 - data->p[i]);
	}

	return;
}
/*****************************************************************************/



/**
this function calculates the logarithmic derivative of the pressure drop dp
with respect to the multi-rate equivalent time dte in the structure "data",
using the parameters in the structure "par". The logarithmic derivative is
stored in data.derl and the group dp/(2*derl) in data.derg.
*/
void calc_log_derivative(modelparameters *par, lmdatatype *data)
{
	int i = 0, j = 1, k = 1;
	double drigth, dleft;

	while(i < par->m){

		drigth = 0.0;
		dleft  = 0.0;

		/************************ dln(dp)/dln(dt) ****************************/
		if(i == 0){
			drigth = log(data->dte[i+j]/data->dte[i]);
			while(drigth < par->Lder && (i + j) < (par->m - 1)){
				j++;
				drigth = log(data->dte[i+j]/data->dte[i]);
			}
			data->derg[i] = log(data->dp[i+j]/data->dp[i])/drigth;
		}
		else if(i == (par->m - 1)){
			dleft = log(data->dte[i]/data->dte[i-k]);
			while(dleft < par->Lder && (i - k) > 0){
				k++;
				dleft = log(data->dte[i]/data->dte[i-k]);
			}
			data->derg[i] = log(data->dp[i]/data->dp[i-k])/dleft;
		}
		else{

			drigth = log(data->dte[i+j]/data->dte[i]);
			while(drigth < par->Lder && (i + j) < (par->m - 1)){
				j++;
				drigth = log(data->dte[i+j]/data->dte[i]);
			}

			dleft = log(data->dte[i]/data->dte[i-k]);
			while(dleft < par->Lder && (i - k) > 0){
				k++;
				dleft = log(data->dte[i]/data->dte[i-k]);
			}

			data->derg[i] = dleft*log(data->dp[i+j])
				/ ( drigth*log(data->dte[i+j]/data->dte[i-k]) )

				+ log(data->dte[i+j]*data->dte[i-k] /
				(data->dte[i]*data->dte[i]))
				*log(data->dp[i])
				/ ( drigth*dleft )

				- drigth*log(data->dp[i-k])
				/ ( dleft*log(data->dte[i+j]/data->dte[i-k]) );
		}

		/***************** Bourdet logarithmic derivative ********************/
		data->derl[i] = data->derg[i] * data->dp[i];

		/***** group dp/(2*derl) (Chow, Onur, Reynolds, Duong, Ozkan) ********/
		data->derg[i] = 0.5 / data->derg[i];

		i++;
		j = k = 1;
	}

	return;
}
/*****************************************************************************/



/**
This function calculates the pressure, using the selected model and the
parameters in the "par" structure, for the times in the "t" vector of the
"data" structure.
*/
void calc_pressure(modelparameters *par, lmdatatype *data)
{
	int i, j;
	double pi;

	pi = par->rpval[INITIAL_PRESSURE];

	for(i = 0; i < par->m; i++){
		if(data->t[i] >= par->tps[0]){
			par->qB = par->qBps[0];
			data->dp[i] =
				(par->dpwffcn)[par->model](par, data->t[i] - par->tps[0]);
		}
		for(j = 1; j < par->nevents; j++){
			if(data->t[i] > par->tps[j]){
				par->qB = par->qBps[j] - par->qBps[j-1];
				data->dp[i] +=
					(par->dpwffcn)[par->model](par, data->t[i] - par->tps[j]);
			}
		}
		data->p[i]  = pi - data->dp[i];
	}

	return;
}
/*****************************************************************************/



/**
this function called by lmder calculates the residual, the jacobian, or prints
the values of x during an iteration, depending on the value of iflag.
if iflag = RESIDUAL calculate the functions at x and return this vector
in fvec. do not alter fjac.
if iflag = JACOBIAN calculate the jacobian at x and return this matrix
in fjac. do not alter fvec.
the value of iflag should not be changed by fcn unless the user wants to
terminate execution of lmder. in this case set iflag to a negative integer.
*/
void resjacprn(int m, int n, double *x, double *fv, double **fjac,
			   lmdatatype *data, void *par, int *iflag)
{
	switch(*iflag) {
	case RESIDUAL:
		residual(m, n, x, fv,       data, par, iflag);
		break;
	case JACOBIAN:
		jacobian(m, n, x, fv, fjac, data, par, iflag);
		break;
	case PRINT:
		printiter(n, x);
		break;
	default:
		printiter(n, x);
		break;
	}

	return;
}
/*****************************************************************************/



/**
Calculates the differences between pressure observations and model estimative,
with parameters x. r = p - (p_i - dpwf)
*/
void residual(int m, int n, double *x, double *fv, lmdatatype *data,
			  void *par, int *iflag)
{
	int i = 0, j = 0;
	double pi;
	modelparameters *mpar;

	mpar = (modelparameters *)par;

	setmodeltox(x, n, mpar);

	pi = mpar->rpval[INITIAL_PRESSURE];

	for(i = 0; i < m; i++){
		if(data->t[i] >= mpar->tps[0]){
			mpar->qB = mpar->qBps[0];
			fv[i] =
				(mpar->dpwffcn)[mpar->model](mpar, data->t[i] - mpar->tps[0])
				+ data->p[i] - pi;
		}
		for(j = 1; j < mpar->nevents; j++){
			mpar->qB = mpar->qBps[j] - mpar->qBps[j-1];
			if(data->t[i] > mpar->tps[j]){
				fv[i] +=
					(mpar->dpwffcn)[mpar->model](mpar, data->t[i] - mpar->tps[j]);
			}
		}
	}

	return;
}
/*****************************************************************************/



/**
Calculates the Jacobian, with elements dr/dx. dr/dx = d(dpwf)/dx for all para-
meters, except the initial pressure. dr/dp_i = -1.0
*/
void jacobian(int m, int n, double *x, double *fv, double **fjac,
			  lmdatatype *data, void *par, int *iflag)
{
	int i, j, k;
	modelparameters *mpar;

	mpar = (modelparameters *)par;

	setmodeltox(x, n, mpar);

	/* calculates the analytical part of the jacobian */
	for(i = 0; i < m; i++){
		if(data->t[i] >= mpar->tps[0]){
			mpar->qB = mpar->qBps[0];
			for(k = 0; k < n; k++){
				if(mpar->jactype[k] == ANALYTICAL){
					fjac[k][i] =
						(mpar->dr_dx)[mpar->model][mpar->partype[k]]
					(mpar, data->t[i] - mpar->tps[0]);
				}
			}
		}
		for(j = 1; j < mpar->nevents; j++){
			if(data->t[i] > mpar->tps[j]){
				mpar->qB = mpar->qBps[j] - mpar->qBps[j-1];
				for(k = 0; k < n; k++){
					if(mpar->jactype[k] == ANALYTICAL){
						fjac[k][i] +=
							(mpar->dr_dx)[mpar->model][mpar->partype[k]]
						(mpar, data->t[i] - mpar->tps[j]);
					}
				}
			}
		}
	}

	/* calculates the finite difference part of the jacobian */
	for(k = 0; k < n; k++){
		if(mpar->jactype[k] != ANALYTICAL){
			fd_jacobian(m, n, x, k, fv, fjac, data, par, iflag);
		}
	}

	return;
}
/*****************************************************************************/



void printiter(int n, double *x)
{
	int i;

	/*printf("Parameters: ");*/
	for (i = 0; i < n; i++)
		printf("%23.15e ", x[i]);

	printf("\n");

	return;
}
/*****************************************************************************/



/**
This function copies the parameters of the model, in the struct "par", to the
corresponding vector of estimates x.
*/
void setxtomodel(double *x, int n, modelparameters *par)
{
	int i;

	for(i = 0; i < n; i++){
		if(par->model == PWFT) { // for transformed parameters
			switch(par->partype[i]) {
			case SKIN_FACTOR:
				x[i] = log( (par->rpval[par->partype[i]] + 8.0)/(25.0 - par->rpval[par->partype[i]]) );
				break;
			case PERMEABILITY:
			case WELLBORE_STORAGE:
			case INITIAL_PRESSURE:
				x[i] = log( par->rpval[par->partype[i]] );
				break;
			default:
				x[i] = par->rpval[par->partype[i]];
				break;
			}
		}
		else{
			x[i] = par->rpval[par->partype[i]];
		}
	}

	return;
}
/*****************************************************************************/



/**
This function copies the values of the vector of estimates x to the
corresponding parameters in the struct "par".
*/
void setmodeltox(double *x, int n, modelparameters *par)
{
	int i;

	for(i = 0; i < n; i++){
		if(par->model == PWFT){ // for transformed parameters
			switch(par->partype[i]){
			case SKIN_FACTOR:
				par->rpval[par->partype[i]] = (25.0*exp(x[i]) - 8.0)/(1.0 + exp(x[i]));
				break;
			case PERMEABILITY:
			case WELLBORE_STORAGE:
			case INITIAL_PRESSURE:
				par->rpval[par->partype[i]] = exp(x[i]);
				break;
			default:
				par->rpval[par->partype[i]] = x[i];
				break;
			}
		}
		else{
			par->rpval[par->partype[i]] = x[i];
		}
	}

	return;
}
/*****************************************************************************/



/**
jacobian approximation by the finite difference method.
*/
void fd_jacobian(int m, int n, double *x, int k, double *fv, double **fjac,
				 lmdatatype *data, void *par, int *iflag)
{
	int i;
	double eps, h, xi;
	double *fip1, *fim1;
	modelparameters *mpar;

	mpar = (modelparameters *)par;

	if(mpar->jactype[k] != CENTRAL &&
		mpar->jactype[k] != FORWARD &&
		mpar->jactype[k] != BACKWARD){
			printf("\n error on jactype values \n");
			exit(EXIT_FAILURE);
	}

	//eps = sqrt(sqrt(sqrt(MACHEPS)));
	//eps = sqrt(sqrt(MACHEPS));
	//eps = sqrt(MACHEPS);
	//eps = DX_DEFAULT;

	eps = mpar->dx[k];

	xi = x[k];
	if(xi == 0.0)
		h = eps;
	else
		h = eps * fabs(xi);

	if(mpar->jactype[k] == CENTRAL || mpar->jactype[k] == FORWARD){
		if((fip1 = (double*)calloc(m, sizeof(double))) == NULL){
			printf("\n memory allocation failure \n");
			exit(EXIT_FAILURE);
		}
		x[k] = xi + h;
		residual(m, n, x, fip1, data, par, iflag);
	}
	if(mpar->jactype[k] == CENTRAL || mpar->jactype[k] == BACKWARD){
		if((fim1 = (double*)calloc(m, sizeof(double))) == NULL){
			printf("\n memory allocation failure \n");
			exit(EXIT_FAILURE);
		}
		x[k] = xi - h;
		residual(m, n, x, fim1, data, par, iflag);
	}

	for(i = 0; i < m; i++){
		switch(mpar->jactype[k]){
		case CENTRAL :
			fjac[k][i] = (fip1[i] - fim1[i]) / (2.0*h);
			break;
		case FORWARD :
			fjac[k][i] = (fip1[i] -   fv[i]) / h;
			break;
		case BACKWARD :
			fjac[k][i] = (fv[i]   - fim1[i]) / h;
			break;
		}
	}

	x[k] = xi;

	if(mpar->jactype[k] == CENTRAL || mpar->jactype[k] == FORWARD){
		free(fip1);
	}
	if(mpar->jactype[k] == CENTRAL || mpar->jactype[k] == BACKWARD){
		free(fim1);
	}

	return;
}
/*****************************************************************************/



/**
levenberg-marquardt non-linear least squares fitting. driver for the lmder
minpack function, with confidence intervals.
*/
void lm_nlsf(void fcn(), lmdatatype *data, void *modelpar, int m, int n,
			 double *x, double *ci, double **corr, double *fvec, double **fjac,
			 double tol, int *info, int *nfev, int *njev, double *fnorm,
			 int mode)
{
	int i,
		j,
		maxfev;

	int *ipvt;

	double
		ftol,
		xtol,
		gtol,
		factor,
		norm,
		covfac,
		student05;

	double
		*diag,
		*qtf,
		*wa1,
		*wa2,
		*wa3,
		*wa4;

	/* Check input parameters */
	if (n <= 0 || m < n || tol < 0.0) {
		*info = 0;
		return;
	}
	/* Allocate memory for working arrays. */
	if( (ipvt = (int *)calloc(n, sizeof(int))) == NULL ||
		(diag = (double *)calloc(n, sizeof(double))) == NULL ||
		(qtf  = (double *)calloc(n, sizeof(double))) == NULL ||
		(wa1  = (double *)calloc(n, sizeof(double))) == NULL ||
		(wa2  = (double *)calloc(n, sizeof(double))) == NULL ||
		(wa3  = (double *)calloc(n, sizeof(double))) == NULL ||
		(wa4  = (double *)calloc(m, sizeof(double))) == NULL ){
			*info = 9;
			return;
	}

	/*  */
	for(j = 0; j < n; j++){
		diag[j] = 1.0;
	}

	/* Set convergence tolerances */
	ftol = tol;
	xtol = tol;
	gtol = tol;

	maxfev = (n + 2) * 100;
	factor = 100.0;
	*nfev  = 0;
	*njev  = 0;

	/*************************************************************************/
	lmder(fcn, data, modelpar, m, n, x, fvec, ftol, xtol, gtol, maxfev, diag,
		mode, factor, info, nfev, njev, fjac, ipvt, qtf, wa1, wa2, wa3, wa4);
	/*************************************************************************/

	student05 = gsl_cdf_tdist_Pinv(1.0 - 0.5/2.0, m - n);
	covar(n, fjac, ipvt, tol, wa1);
	norm = enorm(m, fvec);
	covfac = norm*norm/(m - n);

	*fnorm = norm;
	for (i = 0; i < n; i++){
		/* 95% confidence interval */
		ci[i] = student05 * sqrt(covfac*fjac[i][i]);

		/* correlation matrix */
		for(j = 0; j < n; j++){
			corr[i][j] = fjac[i][j] / sqrt(fjac[i][i] * fjac[j][j]);
		}
	}

	free(wa4);
	free(wa3);
	free(wa2);
	free(wa1);
	free(qtf);
	free(diag);
	free(ipvt);

	return;
}
/*****************************************************************************/



/**
*/
void print_par(FILE *file, double *x, double *ci, double **corr, int n,
			   modelparameters *par)
{
	int i, j;

	fprintf(file, "Parameters: \n");

	for(i = 0; i < n; i++){
		if(par->model == PWFT){ // for transformed parameters
			switch(par->partype[i]){
			case PERMEABILITY:
			    fprintf(file, "%s \t = % 12g\n", par->parnames[par->partype[i]], par->rpval[par->partype[i]]);
				fprintf(file, "ln(k) \t = % 12g +/- %12g ", x[i], ci[i]);
				fprintf(file, x[i] == 0.0 ? "\n" : "or +/- %4g%%\n" , ci[i]*100.0/fabs(x[i]));
				break;
            case SKIN_FACTOR:
				fprintf(file, "%s \t = % 12g\n", par->parnames[par->partype[i]], par->rpval[par->partype[i]]);
				fprintf(file, "ln(S*) \t = % 12g +/- %12g ", x[i], ci[i]);
				fprintf(file, x[i] == 0.0 ? "\n" : "or +/- %4g%%\n" , ci[i]*100.0/fabs(x[i]));
				break;
			case WELLBORE_STORAGE:
                fprintf(file, "%s \t = % 12g\n", par->parnames[par->partype[i]], par->rpval[par->partype[i]]);
				fprintf(file, "ln(C) \t = % 12g +/- %12g ", x[i], ci[i]);
				fprintf(file, x[i] == 0.0 ? "\n" : "or +/- %4g%%\n" , ci[i]*100.0/fabs(x[i]));
				break;
			case INITIAL_PRESSURE:
				fprintf(file, "%s \t = % 12g\n", par->parnames[par->partype[i]], par->rpval[par->partype[i]]);
				fprintf(file, "ln(pi) \t = % 12g +/- %12g ", x[i], ci[i]);
				fprintf(file, x[i] == 0.0 ? "\n" : "or +/- %4g%%\n" , ci[i]*100.0/fabs(x[i]));
				break;
			default:
			    fprintf(file, "%s \t = % 12g +/- %12g ", par->parnames[par->partype[i]], x[i], ci[i]);
			    fprintf(file, x[i] == 0.0 ? "\n" : "or +/- %4g%%\n" , ci[i]*100.0/fabs(x[i]));
			    break;
			}
		}
		else{
			fprintf(file, "%s \t = % 12g +/- %12g ", par->parnames[par->partype[i]], x[i], ci[i]);
			fprintf(file, x[i] == 0.0 ? "\n" : "or +/- %4g%%\n" , ci[i]*100.0/fabs(x[i]));
		}
	}

	fprintf(file, "Correlation matrix: \n");
	for(i = 0; i < n; i++){
		for(j = 0; j < n; j++){
			fprintf(file, "% .4E    ", corr[i][j]);
		}
		fprintf(file, "\n");
	}

	return;
}
/*****************************************************************************/



/**
*/
void write_outfile(modelparameters *par, lmdatatype *data, double *x,
				   double *ci, double **corr)
{
	int     j, n, m;
	FILE*   outfile;

	n = par->n;
	m = par->m;

	outfile = fopen(par->outfile, "w");
	if (outfile == NULL){
		printf("\n cannot open out file: \n");
		exit(EXIT_FAILURE);
	}

	print_par(outfile, x, ci, corr, n, par);

	fprintf(outfile, "Model (t, p, dt, dte, dp, derlog): \n");
	for (j = 0; j < m; j++){
		fprintf(outfile, "%12g %12g %12g %12g %12g %12g \n",
			data->t[j],
			data->p[j],
			data->dt[j],
			data->dte[j],
			data->dp[j],
			data->derl[j]);
	}

	fclose(outfile);

	return;
}
/*****************************************************************************/

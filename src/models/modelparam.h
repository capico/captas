#ifndef MODELPARAM_H
#define MODELPARAM_H

#define LENGTHFN    64  // file name length

#define NMODELS     14  // number of wellbore/reservoir models
#define NPARAMETERS 19  // total number of parameters

/* reservoir/well-bore parameters for regression */
#define PERMEABILITY                    0
#define SKIN_FACTOR                     1
#define WELLBORE_STORAGE                2
#define INITIAL_PRESSURE                3
#define EXTERNAL_RADIUS                 4
#define DISTANCE_TO_FAULT               5
#define DISTANCE_TO_FAULT_1             6
#define DISTANCE_TO_FAULT_2             7
#define DISTANCE_TO_FAULT_1_X           8
#define DISTANCE_TO_FAULT_2_X           9
#define DISTANCE_TO_FAULT_1_Y           10
#define DISTANCE_TO_FAULT_2_Y           11
#define OMEGA                           12
#define LAMBDA                          13
#define FRACTURE_HALF_LENGTH            14
#define FRACTURE_CONDUCTIVITY           15
#define PERMEABILITY_RATIO              16
#define PENETRATION_RATIO               17
#define MIDPOINT_ELEVATION              18

typedef struct
{
	/* number of events */
	int nevents;

	/* times for production events, flow rates and number of points */
	double *tps;
	double *qBps;
	int    *nps;

	/* test type: 0 full multirate, 1 drawdown, 2 buildup, 3 injection,
	4 falloff, 5 multirate buildup, 6 multirate drawdown, ... */
	int testtype;

	/* program mode: 0 nonlinear regression, 1 analytical simulation */
	int mode;

	/* units system: 0 oilfield, 1 ANP */
	int units;

	double C1,C2, C3;

	/* fluid and reservoir properties */
	double
		B,      // formation volume factor
		qB,     // sand face flow rate
		phi,    // porosity
		mu,     // viscosity
		ct,     // total compressibility
		rw,     // wellbore radius
		h;      // formation thickness

	/* regression parameters status */
	int rpsta[NPARAMETERS];
	/* regression parameters derivatives type */
	int rpjac[NPARAMETERS];
	/* regression parameters values */
	double rpval[NPARAMETERS];
	/* finite difference parameters derivatives dx */
	double rpdx[NPARAMETERS];

	char *parnames[NPARAMETERS];

	/* regression model */
	int model;

	/* for loglog plot functions */
	double sign;

	/* dpwf functions */
	double (*dpwffcn[NMODELS])();

	/* interporosity flow functions for double porosity */
	double (*f[NMODELS])();

	/* dpwf partial derivatives */
	double (*dr_dx[NMODELS][NPARAMETERS])();

	/* io files */
	char inifile[LENGTHFN];
	char outfile[LENGTHFN];
	char pressfile[LENGTHFN];
	char ratefile[LENGTHFN];

	/* lines in table files */
	int
		presssize,
		ratesize;

	int
		m, // data points
		n; // regressed parameters

	int *partype;   /* type of parameter (k, S, re, ...)*/
	int *jactype;   /* analytical or numerical derivatives */
	double *dx;

    /* Stehfest's parameter */
	int nstehfest;
	double *v;

	/* numerical derivative parameter for loglog plot functions */
	double Lder;

	int plots;      /* use gnuplot */
} modelparameters;

#endif

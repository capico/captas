#ifndef MODELPARAM_H
#define MODELPARAM_H

#define LENGTHFN    64  // file name length

#define NMODELS     13  // number of wellbore/reservoir models
#define NPARAMETERS 16  // total number of parameters

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
    qB,
    phi,    // porosity
    mu,     // viscosity
    ct,     // total compressibility
    rw,     // wellbore radius
    h;      // formation thickness

    /* wellbore/reservoir parameters */
    double
    k,      // permeability
    C,      // wellbore storage
    S,      // skin factor
    pi,     // initial pressure
    re,     // external radius
    w1,     // distance to fault 1
    w2,     // distance to fault 2
    L,      // distance to the sealing fault
    w1x,    // distance to fault 1, x direction
    w2x,    // distance to fault 2, x direction
    w1y,    // distance to fault 1, y direction
    w2y,    // distance to fault 2, y direction
    omega,  // storativity ratio
    lambda, // interporosity flow coefficient
    xf,     // fracture half length
    fc;     // fracture conductivity

    /* Stehfest's parameter */
    int nstehfest;
    double *v;

    /* numerical derivative parameter for loglog plot functions */
    double Lder;

    /* regression parameters */
    int
    rp_S,
    rp_k,
    rp_C,
    rp_pi,
    rp_re,
    rp_L,
    rp_w1,
    rp_w2,
    rp_w1x,
    rp_w2x,
    rp_w1y,
    rp_w2y,
    rp_omega,
    rp_lambda,
    rp_xf,
    rp_fc;

    /* type of regression parameters derivatives */
    int
    jac_S,
    jac_k,
    jac_C,
    jac_pi,
    jac_re,
    jac_L,
    jac_w1,
    jac_w2,
    jac_w1x,
    jac_w2x,
    jac_w1y,
    jac_w2y,
    jac_omega,
    jac_lambda,
    jac_xf,
    jac_fc;

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

    int plots;      /* use gnuplot */
} modelparameters;

#endif

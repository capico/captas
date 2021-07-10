#ifndef CAPTA_H
#define CAPTA_H

#define LENGTENTRY	64	// max length of a inifile entry

#define JACOBIAN	2
#define RESIDUAL	1
#define PRINT		0

/* type of derivatives in the jacobian matrix */
#define ANALYTICAL	 0
#define FORWARD		+1
#define BACKWARD	-1
#define CENTERED	 2

/* mode */
#define REGRESSION	 0
#define SIMULATION	 1

/* units system */
#define OILFIELD     0
#define ANP          1

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


/* reservoir/well-bore models */
#define PWF			0
#define PWFT	    1
#define PWFNFOB		2
#define PWFCPOB		3
#define PWFF        4
#define PWFC        5
#define PWFCL       6
#define PWFRECT     7
#define PWFDPPSS    8
#define PWFDPTSL    9
#define PWFDPTSP    10
#define PWFICF      11
#define PWFFCF      12


/* test type */
#define MULTIRATE			0
#define DRAWDOWN			1
#define BUILDUP				2
#define INJECTION			3
#define FALLOFF				4
#define MULTIRATE_BUILDUP	5
#define MULTIRATE_DRAWDOWN	6

/* unit conversion factors */
#define C1_ANP			0.0003484228
#define C2_ANP			19.032784
#define C3_ANP		    0.15915494309189533576888376337251

#define C1_OILFIELD         0.00026367864
#define C2_OILFIELD         141.20546
#define C3_OILFIELD         0.89358869

const char *lm_nlsf_msg[] =
{
    "improper input parameters.",
    "both actual and predicted relative reductions in the sum of squares are at most ftol.",
    "relative error between two consecutive iterates is at most xtol.",
    "actual and predicted relative reductions are at most ftol, and relative error between two consecutive iterates is at most xtol.",
    "the cosine of the angle between fvec and any column of the jacobian is at most gtol in absolute value.",
    "number of calls to fcn with iflag = 1 has reached maxfev.",
    "ftol is too small. no further reduction in the sum of squares is possible.",
    "xtol is too small. no further improvement in the approximate solution x is possible.",
    "gtol is too small. fvec is orthogonal to the columns of the jacobian to machine precision.",
    "memory allocation failure."
};

void fd_jacobian(int m, int n, double *x, int k, double *fv, double **fjac,
                 lmdatatype *data, void *par, int *iflag);

void setxtomodel(double *x, int n, modelparameters *par);

void setmodeltox(double *x, int n, modelparameters *par);

void read_inifile(modelparameters *par);

void set_parameters(modelparameters *par);

void calc_pressure(modelparameters *par, lmdatatype *data);

void calc_plotting_abscissas(modelparameters *par, lmdatatype *data);

void calc_pressure_drop(modelparameters *par, lmdatatype *data);

void calc_log_derivative(modelparameters *par, lmdatatype *data);

void resjacprn(int m, int n, double *x, double *fv, double **fjac,
               lmdatatype *data, void *modelpar, int *iflag);

void residual(int, int, double *, double *, lmdatatype *, void *, int *);

void jacobian(int, int, double *, double *, double **, lmdatatype *, void *,
              int *);

void print_par(FILE *file_name, double *x, double *ci, double **corr,
               int n, modelparameters *par);

void parse_args(int argc, char *argv[], modelparameters *par);

void help_msg(void);

void print_version(void);

void printiter(int, double *);

void lm_nlsf(void fcn(), lmdatatype *data, void *modelpar, int m, int n,
             double *x, double *ci, double **corr, double *fvec, double **fjac,
             double tol, int *info, int *nfev, int *njev, double *fnorm, int mode);

void write_outfile(modelparameters *par, lmdatatype *data, double *x,
                   double *ci, double **corr);

#endif

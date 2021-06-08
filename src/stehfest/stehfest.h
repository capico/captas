#ifndef STEHFEST_H
#define STEHFEST_H

/*---------------------------------------------------------------------------*/
double* stehfest_init(const int nstehfest, double *v);

double stehfest_ilt(double const (*f)(const void *, double), const void *par,
					int nstehfest, const double *const v, double t);

#endif

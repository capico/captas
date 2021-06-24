# Computer Aided Pressure Transient Analysis and Simulation – CAPTAS?

CAPTAS is a program, with a very simple interface, to fit pressure test data to a mathematical model, using the Levenberg-Marquardt method, as part of a well test interpretation.
The program is useful to refine the estimates you have made by, for example, graphical methods, and to determine the confidence intervals on the estimated parameter values.

The input file, with .ini extension, contains sections and, within each section, keys (parameters). Each section begins with a sentence enclosed in square brackets, followed by identifiers for the parameters, the equal sign, and the parameter value, which can be a number or a string.
For example:

[Title]
title  = Simulated Build Up;
caseid = BU test 1;

When a semicolon appears after the value of a parameter, the rest of the line is considered as a comment. Lines starting with the `#` character are also considered comments. There is no case difference in the names of sections or parameters.
The sections and parameters are as follows:

- `Title`: used to identify the test, it accepts two string parameters: `title` and `caseid`.
- `Program mode`: its only parameter, `mode`, determines whether the program is executed to adjust test data (`mode` = 0) or as an analytical test simulator (`mode` = 1, not programmed yet).
- `Units system`: its only parameter, `units`, specifies the system of units used for input and output data. The value 0 corresponds to the oilfield system and 1 to the ANP system.
- `Test description`: it is the section with the well test data and the well/reservoir system data. Its parameters are:
  - `testtype`: identifies the type of test, according to the following values: 0 for a test with multiple flow rates, with pressures in more than one flow period (see example Gringarten_1979); 1 for a production test (initial) with constant rate; 2 for a pressure build up test, after a period of production with constant rate (see example Bourdet_1983_1); 3 for an injection test (initial), with constant flow rate; 4 for a falloff test, after a period of constant injection rate; 5 for a period of zero flow rate, preceded by periods with different flow rates; and 6 for a period with positive flow rate, preceded by periods with different flow rates.
  - `pressfile`: is a string with the name of the text file with the pressure vs time data table, including its relative path. Each line of this file contains the time value, separated from the pressure value by spaces or tabs, in the `units` specified in `Units system`.
  - `presssize`: is the number of lines in the file with the pressure data.
  - `ratefile`: is a string with the name of the file with the flow vs. time data table, including its relative path. Each line of this file contains the time value, separated from the flow rate value by spaces or tabs, in the units specified in `Units system`. For example, the table: 0.0 250.0 12.0 500.0 24.0 0.0 indicates that the well produced with a flow rate of 250.0 between times 0.0 and 12.0, and with a flow rate of 500.0 between times 12.0 and 24.0, after which it was closed. With the exception of the test identified as `testtype` = 0, the pressures in the `pressfile` file must match the last period in the `ratefile` file, which is the period under analysis.
   - `ratesize`: it is the number of lines in the file with the flow rate data.
   - `phi`: porosity
   - `B`: formation volume factor
   - `mu`: viscosity
   - `h`: formation thickness
   - `rw`: well radius
   - `ct`: total compressibility
   - `pi`: initial pressure
   - `S`: skin factor
   - `k`: permeability
   - `C`: wellbore storage coefficient
   - `re`: external radius, in circular reservoir models (optional)
   - `L`: distance to sealing fault, in sealing linear fault models (optional)
   - `w1`: distance to the first boundary, in the channel reservoir model (optional)
   - `w2`: distance to the second boundary, in the channel reservoir model (optional)
   - `w1x`: distance to the first boundary, in the x direction, in the rectangular reservoir model (optional)
   - `w2x`: distance to the second boundary, in the x direction, in the rectangular reservoir model (optional)
   - `w1y`: distance to the first boundary, in the y direction, in the rectangular reservoir model (optional)
   - `w2y`: distance to the second boundary, in the y direction, in the rectangular reservoir model (optional)
 - `Regression model`: its only parameter, `model`, specifies the mathematical model of the well/reservoir system used in the regression. The value 0 corresponds to the infinite homogeneous reservoir; 1 to the infinite homogeneous reservoir, with transformed parameters, as sugested by Dastan and Horne (2011); 2 to the sealed circular homogeneous reservoir; 3 to the homogeneous circular reservoir with constant pressure at the boundary; 4 to the homogeneous reservoir with the well close to a sealing linear fault; 5 to the homogeneous reservoir with the well between two parallel sealing linear faults (channel); 6 to the homogeneous reservoir with the well between two parallel sealing linear faults (channel), without wellbore storage effects; and 7 to the rectangular homogeneous reservoir with impermeable boundaries.
 - `Regression parameters`: in this section the parameters that will be considered unknown in the non-linear regression method are specified. A parameter can be included in the regression method by declaring as 1 a variable formed by the prefix `rp_` and the parameter name, for example, `rp_S` = 1. The parameters that can be included are: `pi`, `S`, `k`, `C`, `re`, `L` , `w1`, `w2`, `w1x`, `w2x`, `w1y` and `w2y`.
 - `Regression parameters derivatives`: in this section you specify the types of derivatives with which the Jacobian is calculated, for the Levenberg-Marquard method. The type of derivative is given by the variable formed by the prefix `jac_` and the name of the parameter, for example, `jac_k` = 0. The options are: 0 for analytic derivatives, 1 for forward finite differences, -1 for backward finite differences and 2 for centered finite differences. The default value is 2. Analytical derivatives are available for all models with the exception of the models for circular homogeneous reservoir with constant pressure at the boundary (3) and rectangular homogeneous reservoir with impermeable boundaries (7).
 - `Derivative parameters`: in this section the `Lder` parameter indicates the differentiation interval in the calculation of the logarithmic derivative (values between 0.0 and 0.3).
 - `Stehfest parameters`: in this section the `nstehfest` parameter is the number of terms in the Stehfest method, used in the numerical inversion of solutions in Laplace space. Allowed values are even numbers between 4 and 20, with 12 being the default.
 - `Output`: this section specifies the output text file, including its relative path, in the `outfile` parameter. When gnuplot is available, the non-zero `plots` parameter causes graphs to be shown with the data and the solution. For all types of tests a Cartesian graph of pressure vs. time is shown. For tests that fit only one flow period (`Test description:testtype` different from zero), a semi-logarithmic plot of pressure vs. the equivalent time and a graph of the pressure drop and its logarithmic derivative with realtion to the equivalent time vs. the elapsed time are shown. The pressure drop is defined as deltap = sign*(p0 – p(t)), with sign = -1 for injection and pressure buildup tests and sign = 1 for the others. p0 is the pressure at the beginning of the period under analysis, that is, p0 = p(dt = 0). When the pressure at the beginning of the test is not found in the first line of the `pressfile` file, it is replaced by whatever value is in place.

## License

This library is distributed under the GNU General Public License, version 3 or later (GPLv3+).

## How to compile and run?

To compile use:

gcc capta.c cminpack/covar.c cminpack/enorm.c cminpack/lmder.c cminpack/lmpar.c cminpack/qrfac.c cminpack/qrsolv.c gnuplot_i/gnuplot_i.c iniparser/dictionary.c iniparser/iniparser.c models/dpwf.c models/dpwfc.c models/dpwfcl.c models/dpwfcpob.c models/dpwff.c models/dpwfnfob.c models/dpwfrect.c models/dpwft.c stehfest/stehfest.c utils/utils.c -lm -lgsl -I/usr/include/gsl -lgslcblas -o captas

### Requirements

The program uses open source libraries: 
 - iniparser (http://ndevilla.free.fr/iniparser/), for processing input files; 
 - gnuplot_i (http://ndevilla.free.fr/gnuplot/), for communicating with the gnuplot graphics program;
 - gls (https://www.gnu.org/software/gsl/), for math functions. 
 - The Levenberg-Marquardt method is based on the implementation of the MINPACK library (https://www.netlib.org/minpack/), whose routines have been translated to the C language.

In order to display the plots, the gnuplot program (http://www.gnuplot.info/) should be installed.

### Plataform

captas is designed to run on Windows (gcc with msys2 and wsl) and Linux platforms.

## References

- D. Bourdet (2002), Well Test Analysis: The Use of Advanced Interpretation Models, Elsevier.
- R. Horne (1995), Modern Well Test Analysis: A Computer-Aided Approach, Petroway.
- A. Dastan and R. N. Horne (2010), A New Look at Nonlinear Regression in Well Test Interpretation, Society of Petroleum Engineers (SPE), SPE-135606-MS
- D. Bourdet et al (1983), A New Set of Type Curves Simplifies Well Test Analysis, World Oil

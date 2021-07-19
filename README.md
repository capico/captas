# Computer Aided Pressure Transient Analysis and Simulation – captas?

**captas** is a program for the interpretation and simulation of pressure tests in wells. The program fits the pressure test data to a mathematical model, using the Levenberg-Marquardt method, as part of a well test interpretation. It is useful to refine the estimates you have made by, for example, graphical methods, and to determine the confidence intervals on the estimated parameter values.
It contains several analytical solutions and their analytical derivatives and should be relatively easy to include others. In case the analytical derivatives for the Jacobian calculation are not available, it is possible to use numerical derivatives by finite differences.

## License

This library is distributed under the GNU General Public License, version 3 or later (GPLv3+).

## How to compile?

### Requirements

The program uses open source libraries:
 - iniparser (http://ndevilla.free.fr/iniparser/), for processing input files;
 - gnuplot_i (http://ndevilla.free.fr/gnuplot/), for communicating with the gnuplot graphics program;
 - gls (https://www.gnu.org/software/gsl/), for math functions.
 - The Levenberg-Marquardt method is based on the implementation of the MINPACK library (https://www.netlib.org/minpack/), whose routines have been translated to the C language.

The source files for iniparser and gnuplot_i are included in the `src` folder. The GSL library (https://www.gnu.org/software/gsl/) is required to compile `captas`. In order to display the plots, the gnuplot program (http://www.gnuplot.info/) must be installed.

### Platform

captas is designed to run on Windows (gcc with msys2 or wsl) and Linux platforms.

To compile use:

    gcc capta.c cminpack/covar.c cminpack/enorm.c cminpack/lmder.c cminpack/lmpar.c cminpack/qrfac.c cminpack/qrsolv.c gnuplot_i/gnuplot_i.c iniparser/dictionary.c  iniparser/iniparser.c models/dpwf.c models/dpwfc.c models/dpwfcl.c models/dpwfcpob.c models/dpwff.c models/dpwfnfob.c models/dpwfrect.c models/dpwft.c models/dpwfdp.c models/dpwfdppss.c models/dpwfdptsl.c models/dpwfdptsp.c models/dpwficf.c models/dpwffcf.c stehfest/stehfest.c utils/utils.c -lm -lgsl -I/usr/include/gsl -lgslcblas -o captas

## How to run?

### Input File

The input file, with `.ini` extension, contains sections and, within each section, keys (parameters). Each section begins with a sentence enclosed in square brackets, followed by identifiers for the parameters, the equal sign, and the parameter value, which can be a number or a string.
For example:

    [Title]
    title  = Simulated Build Up;
    caseid = BU test 1;

    [Units system]
    units  = 0; Oilfield unit system

When a semicolon appears after the value of a parameter, the rest of the line is considered as a comment. Lines starting with the `#` character are also considered comments. There is no case difference in the names of sections or parameters.
The sections and parameters are as follows:

- `Title`: used to identify the test, it accepts two string parameters: `title` and `caseid`.
- `Program mode`: its only parameter, `mode`, determines whether the program is executed to adjust test data (`mode` = 0) or as an analytical test simulator (`mode` = 1, not programmed yet).
- `Units system`: its only parameter, `units`, specifies the system of units used for input and output data. The value 0 corresponds to the oilfield system and 1 to the ANP system.
- `Test description`: it is the section with the well test data and the well/reservoir system data. Its parameters are:
  - `testtype`: identifies the type of test, according to the following values: 0 for a test with multiple flow rates, with pressures in more than one flow period (see example Gringarten_1979); 1 for a drawdown test (initial) with constant flow rate; 2 for a pressure build up test, after a period of production with constant flow rate (see example Bourdet_1983_1); 3 for an injection test (initial), with constant flow rate; 4 for a falloff test, after a period of constant injection flow rate; 5 for a shut-in period, preceded by periods with different flow rates; and 6 for a period with positive flow rate, preceded by periods with different flow rates.
  - `pressfile`: is a string with the name of the text file with the pressure vs time data table, including its relative path. Each line of this file contains the time value, separated from the pressure value by spaces or tabs, in the `units` specified in `Units system`.
  - `presssize`: is the number of lines in the file with the pressure data.
  - `ratefile`: is a string with the name of the file with the flow rate vs. time data table, including its relative path. Each line of this file contains the time value, separated from the flow rate value by spaces or tabs, in the units specified in `Units system`. For example, the table:  
    `0.00 250.0`  
    `12.0 500.0`  
    `24.0 000.0`  
indicates that the well produced with a flow rate of 250.0 between times 0.0 and 12.0, and with a flow rate of 500.0 between times 12.0 and 24.0, after which it was closed. With the exception of the test identified as `testtype` = 0, the pressures in the `pressfile` file must match the last period in the `ratefile` file, which is the period under analysis.
   - `ratesize`: it is the number of lines in the file with the flow rate data.
   - `phi`: porosity (required)
   - `B`: formation volume factor (required)
   - `mu`: viscosity (required)
   - `h`: formation thickness (required)
   - `rw`: well radius (required)
   - `ct`: total compressibility (required)
   - `pi`: initial pressure (required)
   - `k`: permeability (required)
   - `S`: skin factor (optional)
   - `C`: wellbore storage coefficient (optional)
   - `re`: external radius, in circular reservoir models (optional)
   - `L`: distance to sealing fault, in sealing linear fault models (optional)
   - `w1`: distance to the first boundary, in the channel reservoir model (optional)
   - `w2`: distance to the second boundary, in the channel reservoir model (optional)
   - `w1x`: distance to the first boundary, in the x direction, in the rectangular reservoir model (optional)
   - `w2x`: distance to the second boundary, in the x direction, in the rectangular reservoir model (optional)
   - `w1y`: distance to the first boundary, in the y direction, in the rectangular reservoir model (optional)
   - `w2y`: distance to the second boundary, in the y direction, in the rectangular reservoir model (optional)
   - `omega`: storativity ratio, in the double porosity reservoir models (optional)
   - `lambda`: interporosity flow coefficient, in the double porosity reservoir models (optional)
   - `xf`: fracture half length, in the fractured well models (optional)
   - `fc`: fracture conductivity, in the fractured well models (optional)
 - `Regression model`: its only parameter, `model`, specifies the mathematical model of the well/reservoir system used in the regression. The value 0 corresponds to the infinite homogeneous reservoir; 1 to the infinite homogeneous reservoir, with transformed parameters, as suggested by Dastan and Horne (2011); 2 to the sealed circular homogeneous reservoir; 3 to the homogeneous circular reservoir with constant pressure at the boundary; 4 to the homogeneous reservoir with the well close to a sealing linear fault; 5 to the homogeneous reservoir with the well between two parallel sealing linear faults (channel); 6 to the homogeneous reservoir with the well between two parallel sealing linear faults (channel), without wellbore storage effects; 7 to the rectangular homogeneous reservoir with impermeable boundaries; 8 to the infinite reservoir with double porosity with pseudo steady state interporosity flow; 9 to the infinite reservoir with double porosity with transient interporosity flow (slabs); 10 to the infinite reservoir with double porosity with transient interporosity flow (spheres); 11 to the infinite reservoir with a infinite conductivity fractured well; and 12 to the infinite reservoir with a finite conductivity fractured well, with dimensionless fracture conductivity greater than one. Default value is zero.
 - `Regression parameters`: in this section the parameters that will be considered unknown in the non-linear regression method are specified. A parameter can be included in the regression method by declaring as 1 a variable formed by the prefix `rp_` and the parameter name, for example, `rp_S` = 1. The parameters that can be included are: `pi`, `S`, `k`, `C`, `re`, `L` , `w1`, `w2`, `w1x`, `w2x`, `w1y`, `w2y`, `omega`, `lambda`, `xf` and `fc`. Default values are zero.
 - `Regression parameters derivatives`: in this section you specify the types of derivatives with which the Jacobian is calculated, for the Levenberg-Marquard method. The type of derivative is given by the variable formed by the prefix `jac_` and the name of the parameter, for example, `jac_k` = 0. The options are: 0 for analytic derivatives, 1 for forward finite differences, -1 for backward finite differences and 2 for central finite differences. The default value is 2. Analytical derivatives are available for all models with the exception of the models for circular homogeneous reservoir with constant pressure at the boundary (3), rectangular homogeneous reservoir with impermeable boundaries (7), double porosity models (8, 9 and 10) and the fractured well models (11 and 12).
 - `Derivative parameters`: in this section the `Lder` parameter indicates the differentiation interval in the calculation of the logarithmic derivative (values between 0.0 and 0.3). Default value is 0.1.
 - `Stehfest parameters`: in this section the `nstehfest` parameter is the number of terms in the Stehfest method, used in the numerical inversion of the solutions in Laplace space. Allowed values are even numbers between 4 and 20, with 12 being the default.
 - `Output`: this section specifies the output text file, including its relative path, in the `outfile` parameter. When gnuplot is available, the non-zero `plots` parameter causes graphs to be shown with the data and the solution. For all types of tests a Cartesian graph of pressure vs. time is shown. For tests that fit only one flow period (`Test description:testtype` different from zero), a semi-logarithmic plot of pressure vs. the equivalent time and a log-log plot of the pressure drop and its logarithmic derivative, with relation to the equivalent time, vs. the elapsed time are shown. The pressure drop is defined as deltap = sign*(p0 – p(t)), with sign = -1 for injection and pressure buildup tests and sign = 1 for the others. p0 is the pressure at the beginning of the period under analysis, that is, p0 = p(dt = 0). When the pressure at the beginning of the test is not found in the first line of the `pressfile` file, it is replaced by whatever value is in place.

### Example:

To run the example, type in the command line:

`captas ./examples/Bourdet_1983_1/Bourdet_example_1.ini`

The text below corresponds to the `Bourdet_example_1.ini` file for the fitting of the data presented in Bourdet et al (1983). It is a pressure buildup test, after a constant flow rate period, interpreted with the infinite homogeneous reservoir model with a vertical well, with wellbore storage and skin factor effects.

    [Title]
    title  = Bourdet et al (1983)
    caseid = buildup test example 1;

    [Units system]
    units  = 1; ANP unit system

    [Test description]
    testtype  = 2; buildup test
    pressfile = examples/Bourdet_1983_1/Bourdet_example_1_p.dat; [h vs. kgf/cm2]
    presssize = 106;
    ratefile  = examples/Bourdet_1983_1/Bourdet_example_1_q.dat; [h vs. m3/d]
    ratesize  = 2;
    phi       = 0.25;       porosity
    B         = 1.06;       formation volume fator
    mu        = 2.5;        viscosity [cp]
    h         = 32.6136;    formation thickness [m]
    rw        = 0.088392;   wellbore radius [m]
    ct        = 0.00005974;	total compressibility [cm2/kgf]
    pi        = 272.58;     initial pressure [kgf/cm2]
    S         = 0.0;        skin factor
    k         = 100.0;      permeability [md]
    C         = 0.01;       wellbore storage [m3/kgf/cm2]

    [Regression parameters]
    rp_S	= 1;
    rp_k	= 1;
    rp_C	= 1;
    rp_pi	= 1;

    [Regression parameters derivatives]
    jac_S	 = 0;   analytical derivatives
    jac_k	 = 0;
    jac_C	 = 0;
    jac_pi = 0;

    [Regression model]
    model = 0;    infinite homogeneous reservoir, vertical well with wellbore storage and skin factor

    [Derivative parameters]
    Lder	= 0.1;  differentiation interval

    [Stehfest parameters]
    nstehfest = 16;

    [Output]
    plots   = 1;  gnuplot installed
    outfile = examples/Bourdet_1983_1/Bourdet_example_1.out;

If the option plots was activated, after a successful fitting, gnuplot will show the following three windows with the pressure history, semi-logarithmic and diagnostic log-log plots:

Pressure history (pressure vs. time)
![History](/examples/Bourdet_1983_1/history.png)

Semi-logarithmic (pressure vs. equivalent time)
![Semilog](/examples/Bourdet_1983_1/semilog.png)

Loglog (pressure drop and pressure drop logarithmic derivative vs. elapsed time)
![Loglog](/examples/Bourdet_1983_1/loglog.png)

## References

- D. Bourdet (2002), Well Test Analysis: The Use of Advanced Interpretation Models, Elsevier.
- R. Horne (1995), Modern Well Test Analysis: A Computer-Aided Approach, Petroway.
- A. Dastan and R. N. Horne (2010), A New Look at Nonlinear Regression in Well Test Interpretation, Society of Petroleum Engineers (SPE), SPE-135606-MS
- D. Bourdet et al (1983), A New Set of Type Curves Simplifies Well Test Analysis, World Oil

# captas
## Computer Aided Pressure Transient Analysis and Simulation – CAPTAS?

CAPTAS is a program, with a very simple interface, to fit pressure test data to a mathematical model, using the Levenberg-Marquardt method, as part of a well test interpretation.
The program is useful to refine the estimates you have made by, for example, graphical methods, and to determine the confidence intervals on the estimated parameter values.

The input file, with .ini extension, contains sections and, within each section, keys (parameters). Each section begins with a sentence enclosed in square brackets, followed by identifiers for the parameters, the equal sign, and the parameter value, which can be a number or a string.
For example:

[Title]
title  = Simulated Build Up;
caseid = BU test 1;

When a semicolon appears after the value of a parameter, the rest of the line is considered as a comment. Lines starting with the # character are also considered comments. There is no case difference in the names of sections or parameters.
The sections and parameters are as follows:

· Title: used to identify the test, it accepts two string parameters: title and caseid.
· Program mode: its only parameter, mode, determines whether the program is executed to adjust test data (mode = 0) or as an analytical test simulator (mode = 1, not programmed yet).
· Units system: its only parameter, units, specifies the system of units used for input and output data. The value 0 corresponds to the oilfield system and 1 to the ANP system.
· Test description: it is the section with the well test data and the well/reservoir system data. Its parameters are:
· testtype: identifies the type of test, according to the following values: 0 corresponds to a test with multiple flow rates, with pressures in more than one flow period (see example Gringarten_1979); 1 is for a production test (initial) with constant rate; 2 is for a pressure build up test, after a period of production with constant rate (see example Bourdet_1983_1); 3 for an injection test (initial), with constant flow rate; 4 is for a falloff test, after a period of constant injection rate; 5 to a period of zero flow rate, preceded by periods with different flow rates; and 6 is for a period with positive flow rate, preceded by periods with different flow rates.

## Lincense

This library is distributed under the GNU General Public License, version 3 or later (GPLv3+).

## How to install and run it?
 
### Requirements

Gnuplot.

### Plataform

captas is designed to run on Windows and Linux platforms.

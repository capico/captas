# captas
## Computer Aided Pressure Transient Analysis and Simulation – CAPTAS?

CAPTAS is a program, with a very simple interface, to fit pressure test data to a mathematical model, using the Levenberg-Marquardt method, as part of a well test interpretation.
The program is useful to refine the estimates you have made by, for example, graphical methods, and to determine the confidence intervals on the estimated parameter values.

The input file, with .ini extension, contains sections and, within each section, keys (parameters). Each section begins with a sentence enclosed in square brackets, followed by identifiers for the parameters, the equal sign, and the parameter value, which can be a number or a string.
For example:

[Title]
title  = Simulated Build Up;
caseid = BU test 1;

## Lincense

This library is distributed under the GNU General Public License, version 3 or later (GPLv3+).

## How to install and run it?
 
### Requirements

Gnuplot.

### Plataform

captas is designed to run on Windows and Linux platforms.

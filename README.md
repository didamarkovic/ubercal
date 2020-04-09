# Ubercal
Euclid relative calibration from uber-cal simulation

## LICENCE

    Copyright (C) 2015 Will Percival & Dida Markovic

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see http://www.gnu.org/licenses.

## INSTALLATION

The ubercal/README.bsh script compiles and runs all the code. Sometimes this might not work, so try copy-pasting the commands form the README.bsh script into the command line and running them manually.
Before you do that however, you should install Mangle (GitHub version: https://github.com/mollyswanson/mangle) by following the instructions. Place the 'mangle' folder into the main folder of the repo so the rest of the code can find it and so that git can ignore it nicely. Note that Mangle uses Fortran, so... good luck to you!

On my Mac, I use Conda environments to get everything in one place (https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html).

Then configure the Mangle Makefile and uncomment "STATICFLAGS:= -static-libgfortran" (l.15, I think), and comment the following line.

Then compile Mangle with 'make static'.

## FURTHER INFO

The survey and instrument specifications have been taken from the [NISP page](http://www.euclid-ec.org/?page_id=2490) on the Euclid Consortium website, the [Euclid Redbook](https://arxiv.org/abs/1110.3193) and the [Prieto et al., end of Phase B1 study](http://www.congrexprojects.com/custom/icso/2012/papers/FP_ICSO-160.pdf).

If you have access to the Euclid Redmine, you can read a basic description of this on the [Euclid Redmine](http://euclid.roe.ac.uk/projects/gcswg/wiki/Calibration).

The original Ubercal paper is the [Padmanabhan et al, 2007](http://arxiv.org/abs/astro-ph/0703454), however, our algorithm is even simpler, since we are doing only simulations.

Then there is the [Holmes et al, 2012](http://arxiv.org/abs/1203.6255), which has the code on GitHub. Again, this is more complex than our algorithm, because it mainly focuses on constraining the variation in the flat-field across a detector.

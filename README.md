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

The test-calibration/README.bsh script compiles and runs all the code.
Before you do that however, you should install Mangle (GitHub version!) by following the instructions. Place the 'mangle' folder into the main folder of the repo so the rest of the code can find it and so that git can ignore it nicely.

## FURTHER INFO

If you have access to the Euclid Redmine, you can read a basic description of this on the [Euclid Redmine](http://euclid.roe.ac.uk/projects/gcswg/wiki/Calibration).

The original Ubercal paper is the [Padmanabhan et al, 2007](http://arxiv.org/abs/astro-ph/0703454), however, our algorithm is even simpler, since we are doing only simulations.

Then there is the [Holmes et al, 2012](http://arxiv.org/abs/1203.6255), which has the code on GitHub. Again, this is more complex than our algorithm, because it mainly focuses on constraining the variation in the flat-field across a detector.

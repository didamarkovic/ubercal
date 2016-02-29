# Ubercal

## Euclid relative calibration from uber-cal simulation

### Will Percival & Dida Markovic, 2015

### INSTALLATION

The test-calibration/README.bsh script compiles and runs all the code.
Before you do that however, you should install Mangle (GitHub version!) by following the instructions. Place the 'mangle' folder into the main folder of the repo so the rest of the code can find it and so that git can ignore it nicely.

### FURTHER INFO

If you have access to the Euclid Redmine, you can read a basic description of this on the [Euclid Redmine](http://euclid.roe.ac.uk/projects/gcswg/wiki/Calibration).

The original Ubercal paper is the [Padmanabhan et al, 2007](http://arxiv.org/abs/astro-ph/0703454), however, I think our algorithm is even simpler, since we are doing only simulations.

Then there is the [Holmes et al, 2012](http://arxiv.org/abs/1203.6255), which has the code on GitHub. Again, this is more complex than our algorithm, because it mainly focuses on constraining the variation in the flat-field across a detector.

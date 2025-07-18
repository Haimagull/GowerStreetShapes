# GowerStreetShapes

Welcome to the Gower Street Shapes repository. Here you can find all of the codes I made to study Dark Matter halos shape and intrinsic alignment from Gower Street simulations outputs (FOF Outputs).

**Contact :** _Oriane.Laurens[at]etu.univ-grenoble-alpes.fr/oriane.laurens[at]outlook.fr_

## Understanding shapes
- **Halo Shaper** is a little side project that helps visualize the influence of semi-axis ratios parameters on an ellipsoid. It includes a projection code that allows to go from a 3D ellipsoid to a 2D ellipse projected with a "stick" (length = ellipticity module; angle = orientation, half of ellipticity argument).
  
## Halo Shape Distribution 
- The **Shape Distribution** * code allows to plot the distribution of halo shape for a given run and time slice (and corresponding redshift). If multiple time slices are generated for a single run, using ffmpeg helps visualising the evolution with redshift with a video.
- The **Ratios Moments** * code is made to calculate the main first order moments such as mean, variance and skewness with error bars to visualise statistically the evolution across redshift rather than visually with the video that can be hard to use to compare different runs.

## Shape correlations
- **FOF to Covo** is a small script that converts a FOF output file into a csv readable by Covo, a 3D correlation code developped by Kai Hoffman, available on Bitbucket.
- The **Correlation Functions** code is made to plot the correlation functions of chosen parameters from a Covo output.

  _\* = You better use an HPC to run these ones on multiple redshifts. I haven't found a way to optimise/parallelise this efficiently but the processing is very slow as it contains a loop to read each row of a binary file when reading the FOF output._

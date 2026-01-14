# GowerStreetShapes

Welcome to the Gower Street Shapes repository. Here you can find all of the codes I made to study Dark Matter halos shape and intrinsic alignment from Gower Street simulations outputs (FOF Outputs).

**Contact :** _Oriane.Laurens[at]etu.univ-grenoble-alpes.fr/oriane.laurens[at]outlook.fr_

## Understanding shapes
- **Halo Shaper** is a little side project that helps visualize the influence of semi-axis ratios parameters on an ellipsoid. It includes a projection code that allows to go from a 3D ellipsoid to a 2D ellipse projected with a "stick" (length = ellipticity module; angle = orientation, half of ellipticity argument).
  
## Halo Shape Distribution 
- The **Shape** * code generates a csv file with semi axis of halos at each redshift for every chosen runs. It can be used before the two codes below to avoid great computing time when trying to plot. If those are not needed specifically, run a more relevant code.
- The **Shape Distribution** * code allows to plot the distribution of halo shape for a given run and time slice (and corresponding redshift). If multiple time slices are generated for a single run, using ffmpeg helps visualising the evolution with redshift with a video. You can change it to use the previous Shape code's outputs.
- The **Ratios Moments** * code is made to calculate the main first order moments such as mean, variance and skewness with error bars to visualise statistically the evolution across redshift rather than visually with the video that can be hard to use to compare different runs. You can change it to use the previous Shape code's outputs.

## Shape correlations
- **FOF to Covo** is a small script that converts a FOF output file into a csv readable by Covo, a 3D correlation code developped by Kai Hoffman, available on Bitbucket.
- The **Correlation Functions** code is made to plot the correlation functions of chosen parameters from a Covo output.
- **FOF to IACorr** is a script that converts a FOF output file into a csv readable by IACorr, an ellipticity correlation code developped by Elisabeth JD, available on Github (https://github.com/elizabethjg/IACorr).

**Note on covo usage** : To run it, type : ./covo covo.params_parameter-file-ID

  _\* = You better use an HPC to run these ones on multiple redshifts. I haven't found a way to optimise/parallelise this efficiently but the processing is very slow as it contains a loop to read each row of a binary file when reading the FOF output._

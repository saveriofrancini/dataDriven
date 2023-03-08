![banner image](banner/banner.png)

# Data driven

**A tool for mapping forest resources in a completely data-driven, design-based framework exploiting GEE remote sensing data as auxiliary information**

<br><br>

## [Manuscript](https://onlinelibrary.wiley.com/doi/abs/10.1002/env.2750/) 

## Citation

> Di Biase RM, Fattorini L, Franceschi S, Grotti M, Puletti N, Corona P. (2022) "From model selection to maps: A completely design‐based data‐driven inference for mapping forest resources". Environmetrics, John Wiley & Sons, Ltd., vol. 33(7).

<br><br>

## [Google Earth Engine user interface](https://code.earthengine.google.com/78c531d69c1c17ef1c5651bc041ad5cf) 

# Introduction
This R package produces maps of forest attributes within study regions tessellated into spatial units togheter with the estimation of their precision by means of a completely design-based data-driven procedure, exploiting 
GEE remote sensing data as auxiliary information. 
GEE remote sensing data.....scriverei brevemente cosa fa GEE APP.
Based on a sample of spatial units, Datadriven package reconstructs the map of the interest attribute and the map of its precision according to the following steps:
i) Akaike-type criterion to remothe GEE variables poorly correlated with the interest variable
ii) Least-squared method to predict the values of the interest variable within the spatial units
iii) Inverse distance weighting interpolator to interpolate the residuals in non-sampled units, using the leave-one-out cross-validation for the selection of the smoothing parameter
iv) Harmonization of interest attribute to match traditional total estimates and map estimates
v) Bootstrap procedure for map uncertainty
The [manuscript](https://onlinelibrary.wiley.com/doi/abs/10.1002/env.2750/) provides a detailed documentation of the procedure.

The package functions require as input a dbf file containing information about the value of the interest attribute recorded at a sample of spatial units,the first-order inclusion probabilities of spatial units, an ID of the strata if a stratified sampling was adopted to select the sample of units,
and the remote sensing data downloaded with the GEE app. 
The R package includes two principal functions: 
`Main()` which constructs the predicted map of the interest attribute
`bootStrap()` which constructs both the predicted map and the associated error map 



# Applicazione GEE e parametri

# Required Input 
`Main()` function requires the following input parameters:
  - `inFile`: path for the input dbf file
  - `incProb`: name of the first-order inclusion probabiliy column
  - `depVar`: name of the interest variable column
  - `varsToRemove`:  name of the columns of variables not considered in the analysis (if there are any)
  - `coordinates` : the spatial coordinates of the population spatial units 

`bootStrap()` function requires the following input parameters:
  - `inFile`: path for the input dbf file
  - `incProb`: name of the first-order inclusion probabiliy column
  - `depVar`: name of the interest variable column
  - `varsToRemove`:  name of the columns of variables not considered in the analysis (if there are any)
  - `coordinates` : spatial coordinates of the population spatial units 
  - `cluster` : ID strata of the population spatial units. 
                If the population was not partitioned into strata and simple random sampling without replacement (SRSWOR) was adopted to select the sample, set the strata equal to 1
  - `iterations` : number of bootstrap iterations 
  - `param outDir`: directory to save the output shapefile

# Output 
`Main()` function returns a data.frame object of three variables: (1) x coordinate, (2) y coordinate, (3) the predicted value of the interest attribute of the spatial units

`bootStrap()` function returns a shapefile object containing both the predicted map ("Harmon") and the associated error map ("error")


<br><br>

## Install and Test the package

# Install and test the dataDriven package. 

`install.packages("devtools", repos = "http://cran.us.r-project.org")`
`devtools::install_github("saveriofrancini/dataDriven")`
`s <- dataDriven::bootStrap(outDir = getwd())`

## Data driven help

For more information on the package and on functions it includes, use the help.

<br><br>

###### Thanks for stopping by and have a great day!

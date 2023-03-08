![banner image](banner/banner.png)

# dataDriven: design-based data-driven mapping and per-pixel error estimation

**A Google Earth Engine and R tool for mapping forest resources in a completely data-driven, design-based framework exploiting remote sensing data as auxiliary information**

**dataDriven** is a tool to produce maps of environmental attributes within an area of interest (AOI) tessellated in spatial units. In addition, dataDriven produces estimation of map precision. This inference is performed using a completely design-based data-driven procedure, exploiting Sentinel-2 remote sensing data as auxiliary information. Owing to the design-based nature of the procedure, the massive modeling involved in model-based approaches is avoided. A peculiarity of the procedure is to develop not only the map of the attribute of interest but also the map of the estimated precision, which represents an essential requirement to allow the user to make informed considerations on the investigated attribute. The statistical consistency of the design-based data-driven maps has been proven by Di Biase et al. [(2022)](https://onlinelibrary.wiley.com/doi/abs/10.1002/env.2750/).

<br><br>

## [Science](https://onlinelibrary.wiley.com/doi/abs/10.1002/env.2750/) 

> Di Biase RM, Fattorini L, Franceschi S, Grotti M, Puletti N, Corona P. (2022) "From model selection to maps: A completely design‐based data‐driven inference for mapping forest resources". Environmetrics, John Wiley & Sons, Ltd., vol. 33(7).

<br><br>

## Requirements
Our tool should only be used if the sample of spatial units, for which the attribute of interest was recorded, was selected according to one of the following sampling schemes:
- (1) Simple Random Sampling Without Replacement (SRSWoR): the spatial units are selected at random without replacement (any unit can be sampled only once).
- (2) One Per Stratum Stratified Sampling (OPSS): the population of spatial units is partitioned into strata of contiguous units, and a unit is randomly selected within each of them.
- (3) Systematic Sampling (SYS): the population of spatial units is divided into strata of contiguous units, a unit is randomly selected in one stratum and then repeated in remaining ones. SYS can only be adopted when the population consists of a grid of regular polygons partitioned into equal-shaped strata.
Our tool requires as input a shapefile. For each sampled unit, the shapefiles gives the value of the attribute of interest, and an ID of the stratum to which the units belongs.
The shapefile must be uploaded on Google Earth Engine (GEE) following those simple [steps](https://developers.google.com/earth-engine/guides/table_upload#upload-a-shapefile). If you don't have a GEE account, you can register [here](https://code.earthengine.google.com/register).
You can find [here](https://code.earthengine.google.com/?asset=users/saveriofrancini/PRIN/aoi) an example of shapefile we uploaded to GEE. It includes all the variables you'll need for your samples.

Based on the ID variable `cluster` specified in the input shapefile, the dataDriven package automatically identifies which of the three sampling schemes has been adopted.


The dataDriven tool consists of two parts: 
- (1) a [GEE app](https://code.earthengine.google.com/ba0aa6f437687148f01b15af0017ceee)
- (2) an [R package](https://github.com/saveriofrancini/dataDriven)

## [GEE app](https://code.earthengine.google.com/ba0aa6f437687148f01b15af0017ceee)
The GEE app is used to preprocess remote sensing data and to construct cloud-free Sentinel-2 composites of different periods. The main steps of the GEE app are (1) filtering Sentinel-2 imagery, (2) cloud masking, and (3) compositing. The application allows you to filter Sentinel-2 imagery by specifying the period, the years of interest, and the Maximum percentage of clouds in the image. Clouds are masked from all remaining imagery using the [Sentinel-2 cloud probability dataset](https://developers.google.com/earth-engine/datasets/catalog/COPERNICUS_S2_CLOUD_PROBABILITY). The final cloud composites are generated for the selected years as the *medoid* of all remaining valid observations. The objective of the Medoid algorithm is to populate the final image composite with the pixel that has a surface reflectance most similar to the median surface reflectance value for that pixel. For more details on *medoid* compositing algorithm see [Kennedy et al. (2018)](https://www.mdpi.com/2072-4292/10/5/691) and [Flood et al. (2013)](https://doi.org/10.3390/rs5126481).
By clicking the run button, Sentinel-2 medoid predictors are calculated for the selected periods and years and the average is calculated for each sample unit in the shapefile that you uploaded and defined using the "Input shapefile" parameter.
You can finally download the shapefile from the tasks window in the top right.
The figure below provides more information about the GEE app parameters.

![appParameters](banner/appParameters.png)

## [R package](https://github.com/saveriofrancini/dataDriven)
The R package is then used to process the [shapefile](https://github.com/saveriofrancini/dataDriven/tree/master/inst/data) exported from Google Earth Engine.
You can use the [shapefile](https://github.com/saveriofrancini/dataDriven/tree/master/inst/data) we provide inside the package for testing purposes. Note that the R package requires a shapefile similar to that, but it is not mandatory to build it using our GEE app, and you can include your own predictors.

On the basis of a sample of spatial units, the dataDriven package estimates the map of the attribute of interest and the map of its precision according to the following steps: 
i) Akaike-type criterion for the removal of GEE variables poorly correlated with the variable of interest; 
ii) Linear regression to estimate the values of the variable of interest in the sampled spatial units; 
iii) Inverse distance weighting interpolator of residuals to estimate  the values of the variable of interest in non-sampled units, using leave-one-out cross-validation for the selection of the smoothing parameter; 
iv) Harmonization of maps for matching traditional estimates of totals (achieved by the regression estimator) with those achieved from maps; 
v) Bootstrap procedure for map uncertainty.

The [manuscript](https://onlinelibrary.wiley.com/doi/abs/10.1002/env.2750/) provides a detailed documentation of the procedure.

The R package includes two main features: 
- `Main()` which builds a shapefile of the estimated attribute.
- `bootStrap()` which builds both the shapefile of the estimated attribute and the shapefile of the associated error. 

# Output 
`Main()` function returns a data.frame object of three variables: (1) x coordinate, (2) y coordinate, (3) the estimated value of the interest attribute of the spatial units

`bootStrap()` function returns a shapefile object containing both the estimated variable ("Harmon") and the associated error ("error")

Within the package, we provide a [folder](https://github.com/saveriofrancini/dataDriven/tree/master/inst/run) that includes a [code](https://github.com/saveriofrancini/dataDriven/blob/master/inst/run/run.R) for testing the package and the expected outputs.

Install the package:

`install.packages("devtools", repos = "http://cran.us.r-project.org")`

Find the path of the [run folder](https://github.com/saveriofrancini/dataDriven/tree/master/inst/run):

`system.file("data", package = "dataDriven")`

# Acknowledgments
DataDriven was made by Saverio Francini, Agnese Marcelli, Gherardo Chirici, Rosa Maria Di Biase, Lorenzo Fattorini and Piermaria Corona as part of research project MULTIFOR (PRIN 2020 – prot E52THS) funded by the Italian Ministry of University and Research.
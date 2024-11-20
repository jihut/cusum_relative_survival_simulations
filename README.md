# Monitoring time to event in medical registry data using CUSUMs based on excess hazard models

This repository containes the codes for the simulations done in "Monitoring time to event in medical registry data using CUSUMs based on excess hazard models" (Tran, Kvaløy and Körner, 2024). 

The folder `simple_examples_scripts` contains some script files of a few simple examples calculating the proposed CUSUM charts under the proportional alternative. The examples show which objects or quantities are needed depending on the choice of model in order to compute the charts. 

The remaining folders provide the results presented in Chapter 3. 

All the functions relevant for calculating the charts, simulating excess times from a piecewise constant baseline excess hazard etc. are included in the package `CUSUMrelsurv` for simplicity and can be found from https://github.com/jihut/CUSUMrelsurv. To install the package, one can use for instance `remotes` or `devtools` by simply writing the following:

`devtools::install_github("https://github.com/jihut/CUSUMrelsurv")`
`remotes::install_github("jihut/CUSUMrelsurv")`

Note: Remember to set the root directory of this repository as working directory before running any of the files if the repository is not opened as a project in R! 

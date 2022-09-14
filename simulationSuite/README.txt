Simulation suite for saemix, containing datasets under different types of data, models and design (K=100 for each scenario)
- the suite is organised by the type of the response
    - see below for the data folders (cont to RTTE data)
- in each folder, several scenarios are provided
    - the scenarios are labelled as modelDesignFeature, for example emaxRichHigh denotes an Emax model with a rich design and high IIV (see below for a brief description)
    - for each scenario, a subfolder data contains the simulated data (K=100 datasets) 
    - for each scenario, a subfolder results is used store the estimations
        - the format could be date_modeller_version_feature.res
- scripts: an additional folder contains the estimation scripts (on Eco's laptop)
    - organised by saemix version

Available datasets and results

* cont: continuous longitudinal data modelled using structural model f and error model g
    - emaxSparseProp: Emax model, sparse design, proportional error model, corresponding to the model simulated by Plan et al. 2012 [S1P]
        - dataset name: data_pdemaxXXX.tab
        - name of the file containing individual parameters: param_pdemaxXXX.tab
        - results obtained with the CRAN version of saemix 3.1
            - 220913_eco_cran31_defaultTrue.res: true parameters for the CI, default settings
            - 220913_eco_cran31_defaultFalse.res: wrong CI, default settings
            - 220913_eco_cran31_longFalse.res: wrong CI, 5 chains, (800, 300) iterations    - hillSparseProp: Hill model, sparse design, proportional error model, corresponding to the model simulated by Plan et al. 2012 [S3P]
        - dataset name: data_pdemaxXXX.tab
        - name of the file containing individual parameters: param_pdemaxXXX.tab
        - results obtained with the CRAN version of saemix 3.1
            - 220913_eco_cran31_defaultTrue.res: true parameters for the CI, default settings
            - 220913_eco_cran31_defaultFalse.res: wrong CI, default settings
            - 220913_eco_cran31_longFalse.res: wrong CI, 5 chains, (800, 300) iterations    - emaxRichProp: Emax model, rich design, proportional error model, corresponding to the model simulated by Plan et al. 2012 [R1P]
        - dataset name: data_pdhillhighXXX.tab
        - name of the file containing individual parameters: param_pdhillhighXXX.tab
        - results obtained with the CRAN version of saemix 3.1
            - 220913_eco_cran31_defaultTrue.res: true parameters for the CI, default settings
            - 220913_eco_cran31_defaultFalse.res: wrong CI, default settings
            - 220913_eco_cran31_longFalse.res: wrong CI, 5 chains, (800, 300) iterations
    - hillRichProp: Hill model, rich design, proportional error model, corresponding to the model simulated by Plan et al. 2012 [R3P]
        - dataset name: data_pdhillhighXXX.tab
        - name of the file containing individual parameters: param_pdhillhighXXX.tab
        - results obtained with the CRAN version of saemix 3.1
            - 220913_eco_cran31_defaultTrue.res: true parameters for the CI, default settings
            - 220913_eco_cran31_defaultFalse.res: wrong CI, default settings
            - 220913_eco_cran31_longFalse.res: wrong CI, 5 chains, (800, 300) iterations
    - other datasets to consider adding to the suite
        - high residual error (conditional bootstrap paper with Sofia) where saemix doesn't perform well
        - different number of subjects (performance in small samples)
        - a PK example
        - future: joint PK/PD models

* binary: binary longitudinal data (0/1)

* cat: categorical data

* count: count data

* tte: time-to-event data

* rtte: repeated time-to-event data

* joint (starting from 4.0)
    - PK/PD
    - PK/binary
    - PK/count
    - PK/TTE

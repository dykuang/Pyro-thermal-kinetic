# Pyro-thermal-kinetic
This repo contains matlab implementations of some commonly used algorithms in pyrolysis or gasification experiments.

Part of the codes and their modifications were used in papers:

  * [Kinetics and mechanism of CO2 gasification of coal catalyzed by Na2CO3, FeCO3 and Na2CO3–FeCO3](https://www.sciencedirect.com/science/article/abs/pii/S1743967119308207)

  * [Predicting kinetic triplets using a 1d convolutional neural network](https://www.sciencedirect.com/science/article/abs/pii/S0040603118304325)

  * [Characterization of Powder River Basin coal pyrolysis with cost-effective and environmentally-friendly composite Na-Fe catalysts in a thermogravimetric analyzer and a fixed-bed reactor](https://www.sciencedirect.com/science/article/abs/pii/S036031991830569X)

  * [Some Optional Methods of Activation Energy Determination on Pyrolysis](https://link.springer.com/article/10.1134/S0023158419020125)


Currently, functions are not well structured and organised. More descriptions and usage examples will be added .......


# Some Key Functions

* `preprocess.m`: Preprocess data obtained from experiement equipments.  
* `batch_E.m`: Activation energy calculation from different methods.
* `daem_point_fit.m`: DAEM method in a point-wise style.  
* `daem_distribution.m`: Make plots for decomposed composnents of f(E) from DAEM method.  
* `daem_simulate.m`: Make curves simulated with parameters calculated from DAEM method.
* `analysis_plots.m`: Making plots for general thermal dynamic variables.  
* `masterG_modified.m`: Making masterplots for deciding reaction functions.  


# Some sample plots
TG/DTG plots:  
<img src = '/sample plots/cc_tg.png' height="260" width="360"><img src = '/sample plots/cc_dtg.png' height="260" width="360">    

Activation Energy E plots:    
<img src = '/sample plots/FWO_E.png' height="260" width="360">

Thermal dynamical parameter plot:  
<img src = '/sample plots/FWO_dG.png' height="260" width="360"><img src = '/sample plots/cc_ds.png' height="260" width="360"> 

DAEM Decomposition plot:  
<img src = '/sample plots/cc_daem_rec15.png' height="260" width="360">
<img src = '/sample plots/cc_fe.png' height="260" width="360">  

Thermal dynamical parameter plot with estimation uncertainty (DAEM):    
<img src = '/sample plots/cc_E_bar.png' height="260" width="360"><img src = '/sample plots/cc_E.png' height="260" width="360"> 

Master Plot:  
<img src = '/sample plots/master.png' height="260" width="360">

If you feel like the repo is useful, please kindly cite this repo.



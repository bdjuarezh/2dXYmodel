# 2dXYmodel
Analysis of BKT phase transition of the 2d XY model (2d O(2) model) in Julia

This work is about simulate the 2d XY model using Julia language programing as a tool. We use a Wolff multicluster algorithm improve with Hoshen-Kopelman
algorithm for cluster identification. 

In the `/src/` directory the core files of simulations are stored. In the `/script/` directory the conde for diferent kind of simulations are stored. 

In the first step we measure de standar observables, like energy, magnetisation, specific heat and susceptibility. Also we present measurements of 
helicity modulus and vorticity. 

The main secction of this work is the propousal of a new observable that evidence de Berenzisnki-Kosterlitz-Thouless phase transition. This observable 
relies in the clusters and vortex of the system. 

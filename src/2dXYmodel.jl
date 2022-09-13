#module XYmodel
#Base function for model definition
include("base_model.jl")
include("base_topology.jl")
include("base_vortex.jl")
include("base_statistics.jl")
#end #module XYmodel

#module ClustersAlgorithms
#Algorithms for single and multicluster Monte Carlo Updates
include("cluster_multi.jl")
include("cluster_single.jl")
#end # module ClustersAlgorithms

#module ObservablesFunctions
#Function of system observables 
include("observables_cluster.jl")
include("observables_standard.jl")
include("observables_vortex.jl")
#end # module ObservablesFunctions

#module ExpectedValues
#Routines for compute the expected value or system observables
include("expval_observables_cluster.jl")
include("expval_observables_standard.jl")
include("expval_observables_vortex.jl")
#end #module ExpectedValues

#module PlotVortex
#For spin, clusters and vortex visulization
include("plot_vortex.jl")
#end # module PlotVortex

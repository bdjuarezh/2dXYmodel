#=             Modules or measures vorticity in a configuration

This functions work with single and multi cluster algorithm, in both cases we need
the wolff's vector and a list of list of all the clusters in the configuration 

Functions
-svorticity()
    -list_of_labels()
    -list_of_clusters()
    -reference_conf()
    -partial_flip()
    -count_vortex()


CM : CLUSTER MATRIX
LL : list Of Label
LC : list of Cluster
RC : reference Configuration
RCp : reference Configuration partial 
nC : number of Clusters
=#

function reference_configuration(L::Int64,  Λ::Array{Float64, 3}, r::Array{Float64, 1})
    #=
    Flip all spins to only one side of Wolff's line
    =#
    refΛ = deepcopy(Λ)
    for i in 1:L, j in 1:L
        dotWolff = dot_wolff(refΛ, r, i, j)
        if dotWolff < 0 
            refΛ[i, j, 1], refΛ[i, j, 2] = one_flip(refΛ, r, i, j, dotWolff)
        end
    end
    return refΛ
end

function partial_flip(refΛ::Array{Float64, 3}, elementsCluster::Array{Any, 1}, r::Array{Float64, 1})
    #=
    Función que toma la configuración de referencia 
    y voltea (flip) solo los espines de un cluster
    =#
    parRefΛ = deepcopy(refΛ)
    for l in elementsCluster
        i, j = l
        dotWolff = dot_wolff(parRefΛ, r, i, j)
        parRefΛ[i, j, 1], parRefΛ[i, j, 2] = one_flip(parRefΛ, r, i, j, dotWolff)
    end
    return parRefΛ
end

function boundary_sites(adSites::Array{Array{Int64, 1}, 3}, elementsBoundary::Array{Any, 1}, checkDual::Array{Int64, 2}, flag::Int64)
    #=
    Crea una lista de los sitios duales o sitios de plaquetas que pasan por la frontera
    =#
    revSites = []
    for l in elementsBoundary
        ads =  normal_to_dual(adSites, l[1], l[2])
        for s in ads
            if checkDual[s[1], s[2]] != flag
                push!(revSites, s)
                checkDual[s[1], s[2]] = flag
            end
        end
    end
    return revSites
end

function count_vortex(L::Int64, ADS::Array{Array{Int64, 1}, 3}, r::Array{Float64, 1}, listClusters::Array{Array{Any,1},1}, listBoundaries::Array{Array{Any,1},1}, refΛ::Array{Float64, 3})
    #=
    Función que mide la vorticidad de los clusters con 
    ayuda de lo configuración de referencia 
    =#
    
    checkDual = zeros(Int64, L, L) 
    clusterVorticity = 0
    clusterVorticitySquare = 0
    
    
    for i in 1:length(listClusters)
        elementsCluster = listClusters[i]
        elementsBoundary = listBoundaries[i]
        flag = i
        
        if elementsCluster != []
            parRefΛ = partial_flip(refΛ, elementsCluster, r)
            revSites = boundary_sites(ADS, elementsBoundary, checkDual, flag)
            V, A = find_vortex(L, parRefΛ, revSites)
            
            v = length(V); a = length(A)
            clusterVorticity = clusterVorticity + v
            clusterVorticitySquare = clusterVorticitySquare + v*v
            
            if v != a
                println("Error: ", V, " not equal ", A) 
            end
        end
    end
    return clusterVorticity, clusterVorticitySquare
end

#------------------------------------Main Functions-------------------------------#
function cluster_vorticity(L::Int64, ADS::Array{Array{Int64, 1}, 3}, Λ::Array{Float64, 3}, r::Array{Float64, 1}, clusters::Array{Int64, 2}, listClusters::Array{Array{Any,1},1}, listBoundaries::Array{Array{Any,1},1})
    #=
    Calcula la vorticidad cluster de una configuración
    sumando la vorticidad de cada cluster
    =#
    refΛ = reference_configuration(L, Λ, r)
    vorticity, vorticitySquare = count_vortex(L, ADS, r, listClusters, listBoundaries, refΛ)
    return vorticity, vorticitySquare
end

#=
This module implements the Wolff's multi-cluster algorithm.

We divide the function in various sub functions and define varios types of sweeps for diferents 
proposites. 

create() := The first stage is create the link matrix  and then compute the bound between all 
spins with a probability function based in the wolff's algorithm for bound creations. At this time the cluster already exist, but we don't known where there are. 

find() := The next stage is the cluster identification. For this we use the hoshen-kopelman algorithm that identifies the bound with O(L^2) complexity. 

flip() := Then we flip all spins of each cluster with probability p=0.5 respect the Wolff's line. 

Functions
-> step_sw()
     -create()
     -find()
        -hk_algorithm
        -local_bonds()
        -relabel()
     -flip()
=#


#----------------------------MULTI CLUSTER MODULES----------------------------#
#First Big Function 
function create_bonds(L::Int64, T::Array{Array{Int64,1}, 3}, Λ::Array{Float64, 3}, β::Float64, r::Array{Float64, 1}) 
    #=
    This function codifies the vertical and horizontal bonds. First we define a matrix thad represents
    all posible links between spins. Then we compute if there exist a vertical and horizantal 
    bound with a probability function,  if this functions is true then we put a bound that identifies 
    with a 1.
    
    bonds[:,:,1] horizontal, bonds[:,:,2] vertical 
    i : down,  j : right
    we only choice 2 directions down  and right.
    =#
    
    bonds = zeros(Int64, L, L, 2) 
    for j in 1:L, i in 1:L        
        up, left, down, right = neighbors(T, i, j)
            
        siteDot = dot_wolff(Λ, r, i, j) 
        rightDot = dot_wolff(Λ, r, right[1], right[2])
        downDot = dot_wolff(Λ, r, down[1], down[2])
            
        rightΔE = siteDot*rightDot
        downΔE = siteDot*downDot
        
        rightProb, downProb = probability(β, rightΔE), probability(β, downΔE)
            
        if rand() < rightProb    bonds[i, j, 1] = 1    end 
        if rand() < downProb     bonds[i, j, 2] = 1    end
    end
    return bonds
end

#>>>>>>>>>>>>>>>>>>>>Subroutines of function find()<<<<<<<<<<<<<<<<<<<<#
function hoshen_kopelman(labels::Array{Int64, 1}, clusters::Array{Int64, 2}, i::Int64, j::Int64)
    #=  Hoshen-Kopelman algorithm
    This function provides a fast form for find the real index of a cluster. 
    First check if  the value in the array is the same as the index of the array
    (The value is provided with the function call).  If this is true, returns the 
    index, if not, then obtains the value stored in that entry and uptade the value 
    of the  number k, and then continue checking.
    =#
    
    while labels[clusters[i, j]] != clusters[i, j]
        clusters[i, j] = labels[clusters[i, j]]
    end
    return clusters[i, j]
end

# Genera una lista con los sitios vecinos si pertenecen al cluster
function local_bonds(L::Int64, bonds::Array{Int64, 3}, i::Int64, j::Int64) 
    #=
    Make a list with the neighbor sites if belongs to the cluster
    =#
    
    boundedNbrs = Array[]
    #Horizontal bonds
    if (j > 1 && bonds[i, j - 1, 1] == 1) push!(boundedNbrs, [i, j - 1]) end 
    if (j == L && bonds[i, j, 1] == 1)    push!(boundedNbrs, [i, 1])     end
    #vertical bonds
    if (i > 1 && bonds[i - 1, j, 2] == 1) push!(boundedNbrs, [i - 1, j]) end 
    if (i == L &&  bonds[i, j, 2] == 1)   push!(boundedNbrs, [1, j])     end
    return boundedNbrs
end
# Reetiqueta lo sitios que pertenezcan al mismo cluster usando las clases de eq
function relabel(L::Int64, labels::Array{Int64, 1}, clusters::Array{Int64, 2})
    #=
    Relebel the sites that belongs to the same clusters using the equivalence classes
    =#
    for j in 1:L, i in 1:L
        clusters[i, j] = hoshen_kopelman(labels, clusters, i, j)
    end
    return clusters
end
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<#

#Second Big Function
function find_clusters(L::Int64, bonds::Array{Int64, 3})   #find the clusters
    #=
    This functions find all the clusters in the spin configuration. 
    We check the bonds matrix and identifies the bounded spins locally and then global. 
    We construct a cluster matrix where every value denotes the cluster that the spin 
    in that site belongs. We storage the labels in a array of lenght L^2 (that is the maximun
    number of clusters)
    =#
    
    clusters = zeros(Int64, L, L)
    labels   = zeros(Int64, L*L)
    flag     = 1
    
    for j in 1:L, i in 1:L
        boundedNeighbors = local_bonds(L, bonds, i, j)
        boundedNumber = length(boundedNeighbors)
            
        if boundedNumber == 0                   # There are no bonds 
            clusters[i, j] = flag
            labels[flag] = flag
            flag = flag  + 1
        else  
            minFlag = flag                      # Choose the min flag between neighbors
            for nb in boundedNeighbors
                pFlag = hoshen_kopelman(labels, clusters, nb[1], nb[2])
                if minFlag > pFlag    
                    minFlag = pFlag    
                end
            end
                
            clusters[i, j] = minFlag            # Label the current site
                
            for nb in boundedNeighbors          # Relabel the neighbors wiht the min flag and the classes too.
                pFlag_n = clusters[nb[1], nb[2]] 
                clusters[nb[1], nb[2]] = minFlag
                labels[pFlag_n] = minFlag
            end
        end
    end
    clusters = relabel(L, labels, clusters)
    return clusters
end

#Third Big Function 
function flip_spins(L::Int64, T::Array{Array{Int64,1}, 3}, Λ::Array{Float64, 3}, clusters::Array{Int64, 2}, r::Array{Float64, 1}) 
    #=
    Flip the spins of all clusters respect the Wolff's line. In this algorithm we make a 
    sweep flip of al spin respect the cluster that owns. Then we use this sweep for make 
    measures of size, perimeter and number of clusters in the configuration.
    =#
    Vol = L*L
    newSpin           = zeros(Int64, Vol)
    sizeClusters      = zeros(Int64, Vol)
    perimeterClusters = zeros(Int64, Vol)
    numberClusters    = 0  
    listClusters      = [[] for i in 1:Vol]
    listBoundaries      = [[] for i in 1:Vol]
    
    
    for j in 1:L, i in 1:L
        flag = clusters[i, j]
        if newSpin[flag] == 0
            (rand() < 0.5) ? newSpin[flag] = 1 : newSpin[flag] = -1
            numberClusters = numberClusters + 1       #numer of cluster measure  
        end
            
        if newSpin[flag] == -1
            dotWolff = dot_wolff(Λ, r, i, j)
            Λ[i, j, 1], Λ[i, j, 2] = one_flip(Λ, r, i, j, dotWolff)
        end
        
        sizeClusters[flag] = sizeClusters[flag] + 1   #cluster size measure
        push!(listClusters[flag], [i, j])             #list of elements of each cluster
        
       
        nbrs = neighbors(T, i, j)                     
        key = true
        for n in nbrs
            if flag != clusters[n[1], n[2]]
                perimeterClusters[flag] = perimeterClusters[flag] + 1  #cluster perimeter measure
                if key
                    push!(listBoundaries[flag], [i, j])                  #list of elements of boundary of each cluster
                    key = false
                end
            end
        end
        
        if key
            nbrs2 = second_neighbors(T, i, j)
            key2 = true
            for n in nbrs2 
                if flag != clusters[n[1], n[2]] && key2
                    push!(listBoundaries[flag], [i, j])
                    key2 = false
                end
            end
        end
    end
    #ya están listas listClusters  y  listBoundaries 
    return Λ, sizeClusters, perimeterClusters, numberClusters, listClusters, listBoundaries
end 

#Complete Sweep 
function sweep_multicluster(L::Int64, T::Array{Array{Int64,1}, 3}, Λ::Array{Float64, 3}, β::Float64) 
    #=
    Complete Sweep that flip the spins of all clusters. We define randomly a Wolff's vector for this 
    sweep only. Then create the bonds between spins and then we find the clusters of the 
    configuration. Then we flip randomly all the spins of each cluster.
    =#
    
    θ = 2*pi*rand()
    r = Float64[cos(θ), sin(θ)]
    
    bonds = create_bonds(L, T, Λ, β, r)
    clusters = find_clusters(L, bonds)    
    Λ, sizeClusters, perimeterClusters, numberClusters, listClusters, listBoundaries = flip_spins(L, T, Λ, clusters, r) 
    
    return Λ, r, clusters, sizeClusters, perimeterClusters, numberClusters, listClusters, listBoundaries
end
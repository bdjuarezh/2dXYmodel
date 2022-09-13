#=----------------------------SINGLE CLUSTER MODULES----------------------------
This file contains the algorithms for implement the Wolff's single cluster algorithm
We first select a random site and respect the spin of this site we decide the belonging
to the cluster in his neighbors. An then the spins that results bound whit the initial spin
are in the process of found if yours neighbor belong to the cluster too, via the function 
probability

=#

function flip_spins(T::Array{Array{Int64,1}, 3}, Λ::Array{Float64, 3}, cluster::Array{Int64, 2}, sizeCluster::Int64, sitesCluster, newSpins)
    perimeterCluster = 0
    for k in 1:sizeCluster
        i, j = sitesCluster[k]
        Λ[i, j, 1], Λ[i, j, 2] = newSpins[k]
        
        #cluster perimeter measures
        nbrs = neighbors(T, i, j) 
        for n in nbrs
            if cluster[n[1], n[2]] != 1
                perimeterCluster = perimeterCluster + 1
            end
        end
    end
    return Λ, perimeterCluster
end

function sweep_singlecluster(L::Int64, T::Array{Array{Int64,1}, 3}, Λ::Array{Float64, 3}, β::Float64)
    θ = 2*pi*rand()                
    r = Float64[cos(θ), sin(θ)]  
    
    sitesCluster = Array[]
    newSpins     = Array[]
    check        = Array[]
    checkAux     = Array[]
    
    cluster = zeros(Int64, L, L)
    g, h = rand(1:L), rand(1:L)
    
    cluster[g, h] = 1
    push!(sitesCluster, [g, h])
    push!(check, [g, h])
    
    spinDotWolff = dot_wolff(Λ, r, g, h)
    xNew, yNew = one_flip(Λ, r, g, h, spinDotWolff)
    push!(newSpins, [xNew, yNew])
    
    while length(check) > 0
        for site in check
            i, j = site
            neighborsSite = neighbors(T, i, j)
            siteDotWolff = dot_wolff(Λ, r, i, j)
            
            for nSite in neighborsSite
                k, l = nSite
                newDotWolff = dot_wolff(Λ, r, k, l)
                deltaEnergy = siteDotWolff*newDotWolff
                p = probability(β, deltaEnergy)
                
                if cluster[k, l] == 0 && rand() < p
                    cluster[k, l] = 1
                    xNew, yNew = one_flip(Λ, r, k, l, newDotWolff)
                    push!(newSpins, [xNew, yNew])
                    push!(checkAux, [k, l])
                    push!(sitesCluster, [k, l])
                end
            end
        end
        check = copy(checkAux)
        checkAux = Array[]
    end
    sizeCluster = length(sitesCluster)
    Λ, perimeterCluster = flip_spins(T, Λ, cluster, sizeCluster, sitesCluster, newSpins)
    return Λ, r, cluster, sizeCluster, perimeterCluster, sitesCluster
end

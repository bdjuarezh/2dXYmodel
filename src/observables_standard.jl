#=



=#

function thermalisation(L::Int64,  T::Array{Array{Int64,1}, 3}, Λ::Array{Float64, 3}, β::Float64, thermalisationTime::Int64, thermalisationType::String="multicluster")
    #=
    This function only implement the sweep for a number of times with the porpouse of 
    thermalisation of the system. We can choice between multicluster, multicluster_vorticity, and 
    singlecluster sweeps
    =#
    
    if thermalisationType=="multicluster"
        for t in 1:thermalisationTime
            Λ, r, clusters, sizes, perimeters, numbers, listClusters, listBoundaries = sweep_multicluster(L, T, Λ, β)
        end
  
    elseif thermalisationType=="singlecluster"
        for t in 1:thermalisationTime
            Λ, r, cluster, size, perimeter, sites = sweep_singlecluster(L, T, Λ, β)
        end
    end
    return Λ
end

function energy(L::Int64, T::Array{Array{Int64,1}, 3}, Λ::Array{Float64, 3})
    #= Compute the energy of
    a single configuration =#
    E = 0.0
    for j in 1:L, i in 1:L
        up, left, down, right = neighbors(T, i, j)
        sumX = Λ[up[1], up[2], 1] + Λ[right[1], right[2], 1]
        sumY = Λ[up[1], up[2], 2] + Λ[right[1], right[2], 2]
        e = Λ[i, j, 1]*sumX +  Λ[i, j, 2]*sumY
        E = E + e
    end
    return -E
end

function magnetisation(Λ::Array{Float64, 3})
    #= Compute the magnetization 
    of a single configuration =#
    x = sum(Λ[:,:,1])
    y = sum(Λ[:,:,2])
    m = sqrt(x*x + y*y)
    return m
end

function helicity(L::Int64, T::Array{Array{Int64,1}, 3}, Λ::Array{Float64, 3}, β::Float64)
    #= 
    Compute the helicity modulus 
    of a single configurations 
    h = <s_x*s_{x+1}>/V + (β/V) * [<(s_x^1*s_{x+1}^2)> + ...]
    =#
    h1 = 0.0
    h2 = 0.0
    for j in 1:L, i in 1:L
        upNb, leftNb, downNb, rightNb = neighbors(T, i, j)
  
        spinX = Λ[i, j, 1]                     #s_{x}^{(1)}
        spinY = Λ[i, j, 2]                     #s_{x}^{(2)}
        nbX = Λ[rightNb[1], rightNb[2], 1]     #s_{x+1}^{(1)}
        nbY = Λ[rightNb[1], rightNb[2], 2]     #s_{x+1}^{(2)}
            
        dotProduct   = spinX*nbX + spinY*nbY
        crossProduct = spinX*nbY - spinY*nbX
        
        h1 = h1 + dotProduct
        h2 = h2 + crossProduct
    end
    return h1 - β*h2*h2
end

function vorticity(L::Int64, Λ::Array{Float64, 3})
    #=
    Función que mide la vorticidad normal de 
    una configuración contando de manera directa
    la cantidad de vortices
    =#
    vortex, antivortex = find_vortex(L, Λ)
    v = length(vortex)
    return v
end


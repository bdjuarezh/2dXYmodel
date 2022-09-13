#= Mean value of observables
This module contains abreviation of function that compute de expected value of standar observables. 
We must provide the measure time, the time between mesures (inter time) and the confNumber. 


=#

function mean_energy(L::Int64, T::Array{Array{Int64,1}, 3}, Λ::Array{Float64, 3}, β::Float64, measureTime::Int64, interTime::Int64, jkBins::Int64, typeAlgorithm::String="multicluster")
    EE = Float64[]
    EE2 = Float64[]
    if typeAlgorithm=="multicluster"
        for t in 1:measureTime
            Λ, r, clusters, sizes, perimeters, number, listClusters, listBoundaries = sweep_multicluster(L, T, Λ, β)
            if t%interTime == 0
                e = energy(L, T, Λ)
                push!(EE, e)
                push!(EE2, e*e)
            end
        end
        
    elseif typeAlgorithm=="singlecluster"
        for t in 1:measureTime
            Λ, r, cluster, size, perimeter, sites = sweep_singlecluster(L, T, Λ, β)
            if t%interTime == 0
                e = energy(L, T, Λ)
                push!(EE, e)
                push!(EE2, e*e)
            end
        end
    end
    E, errorE = jackknife(EE, jkBins)
    E2, errorE2 = jackknife(EE2, jkBins)
    return E, errorE, E2, errorE2
end

function mean_magnetisation(L::Int64, T::Array{Array{Int64,1}, 3}, Λ::Array{Float64, 3}, β::Float64, measureTime::Int64, interTime::Int64, jkBins::Int64, typeAlgorithm::String="multicluster")
    MM = Float64[]
    MM2 = Float64[]
    if typeAlgorithm=="multicluster"
        
        for t in 1:measureTime
            Λ, r, clusters, sizes, perimeters, number, listClusters, listBoundaries = sweep_multicluster(L, T, Λ, β)
            if t%interTime == 0
                m = magnetisation(Λ)
                push!(MM, m)
                push!(MM2, m*m)
            end
        end
        
    elseif typeAlgorithm=="singlecluster"
        
        for t in 1:measureTime
            Λ, r, cluster, size, perimeter, sites = sweep_singlecluster(L, T, Λ, β)
            if t%interTime == 0
                m = magnetisation(Λ)
                push!(MM, m)
                push!(MM2, m*m)
            end
        end
        
    end
    M, errorM = jackknife(MM, jkBins)
    M2, errorM2 = jackknife(MM2, jkBins)
    return M, errorM, M2, errorM2
end

function suceptibility(L::Int64, T::Array{Array{Int64,1}, 3}, Λ::Array{Float64, 3}, β::Float64, measureTime::Int64, interTime::Int64, jkBins::Int64, typeAlgorithm::String="multicluster")
    χχ = Float64[]
    
    if typeAlgorithm=="multicluster"
        for t in 1:measureTime
            Λ, r, clusters, sizes, perimeters, number, listClusters, listBoundaries = sweep_multicluster(L, T, Λ, β)
            if t%interTime == 0
                m = magnetisation(Λ)
                push!(χχ, m*m)
            end
        end
        
    elseif typeAlgorithm=="singlecluster"
        for t in 1:measureTime
            Λ, r, cluster, size, perimeter, sites = sweep_singlecluster(L, T, Λ, β)
            if t%interTime == 0
                m = magnetisation(Λ)
                push!(χχ, m*m)
            end
        end
    end
    
    χ, errorχ = jackknife(χχ, jkBins)
    return χ, errorχ
end

function helicity_modulus(L::Int64, T::Array{Array{Int64,1}, 3}, Λ::Array{Float64, 3}, β::Float64, measureTime::Int64, interTime::Int64, jkBins::Int64, typeAlgorithm::String="multicluster")
    ΥΥ = Float64[] #Upsilon
    
    if typeAlgorithm=="multicluster"
        
        for t in 1:measureTime
            Λ, r, clusters, sizes, perimeters, number, listClusters, listBoundaries = sweep_multicluster(L, T, Λ, β)
            if t%interTime == 0
                y = helicity(L, T, Λ, β)
                push!(ΥΥ, y)
            end
        end
    
    elseif typeAlgorithm=="singlecluster"
        
        for t in 1:measureTime
            Λ, r, cluster, size, perimeter, sites = sweep_singlecluster(L, T, Λ, β)
            if t%interTime == 0
                y = helicity(L, T, Λ, β)
                push!(ΥΥ, y)
            end
        end
       
    end
    Υ, errorΥ = jackknife(ΥΥ, jkBins)
    return Υ, errorΥ
end

function mean_vorticity(L::Int64, T::Array{Array{Int64,1}, 3}, Λ::Array{Float64, 3}, β::Float64, measureTime::Int64, interTime::Int64, jkBins::Int64, typeAlgorithm::String="multicluster")
    VV = Float64[]
    
    if typeAlgorithm=="multicluster"
        
        for t in 1:measureTime
            Λ, r, clusters, sizes, perimeters, number, listClusters, listBoundaries = sweep_multicluster(L, T, Λ, β)
            if t%interTime == 0
                v = vorticity(L, Λ)
                push!(VV, v) 
            end
        end
        
    elseif typeAlgorithm=="singlecluster"
        
        for t in 1:measureTime
            Λ, r, cluster, size, perimeter, sites = sweep_singlecluster(L, T, Λ, β)
            if t%interTime == 0
                v = vorticity(L, Λ)
                push!(VV, v) 
            end
        end
        
    end
    V, errorV = jackknife(VV, jkBins)
    return V, errorV
end

function vortex_suceptibility(L::Int64, T::Array{Array{Int64,1}, 3}, Λ::Array{Float64, 3}, β::Float64, measureTime::Int64, interTime::Int64, jkBins::Int64, typeAlgorithm::String="multicluster")
    VV = Float64[]
    VV2 = Float64[]
    
    if typeAlgorithm=="multicluster"
        
        for t in 1:measureTime
            Λ, r, clusters, sizes, perimeters, number, listClusters, listBoundaries = sweep_multicluster(L, T, Λ, β)
            if t%interTime == 0
                v = vorticity(L, Λ)
                push!(VV, v) 
                push!(VV2, v*v) 
            end
        end
        
    elseif typeAlgorithm=="singlecluster"
        
        for t in 1:measureTime
            Λ, r, cluster, size, perimeter, sites = sweep_singlecluster(L, T, Λ, β)
            if t%interTime == 0
                v = vorticity(L, Λ)
                push!(VV, v) 
                push!(VV2, v*v) 
            end
        end
        
    end
    
    V, errorV = jackknife(VV, jkBins)
    V2, errorV2 = jackknife(VV2, jkBins)
    return V, errorV, V2, errorV2 
end

function simulate_2dXYmodel_multicluster(L::Int64, T::Array{Array{Int64,1}, 3}, Λ::Array{Float64, 3}, β::Float64, thermTime::Int64,  measureTime::Int64, interTime::Int64, jkBins::Int64)
    Λ = thermalisation(L, T, Λ, β, thermTime)
    #=
    Functions that compute all observables for only one β with the multicluster algorithm. 
    =#
    
    EE  = Float64[]
    EE2 = Float64[]
    MM  = Float64[]
    χχ  = Float64[]
    ΥΥ  = Float64[]
    VV  = Float64[]
    VV2 = Float64[]
    
    for t in 1:measureTime
        Λ, r, clusters, sizes, perimeters, numbers, listClusters, listBoundaries = sweep_multicluster(L, T, Λ, β)
        if t%interTime == 0
            e = energy(L, T, Λ)
            m = magnetisation(Λ)
            y = helicity(L, T, Λ, β)
            v = vorticity(L, Λ)
            
            push!(EE, e)
            push!(EE2, e*e)
            push!(MM, m)
            push!(χχ, m*m)
            push!(ΥΥ, y)
            push!(VV, v)
            push!(VV2, v*v)
        end
    end
    
    E, errorE = jackknife(EE, jkBins)
    E2, errorE2 = jackknife(EE2, jkBins)
    M, errorM = jackknife(MM, jkBins)
    χ, errorχ = jackknife(χχ, jkBins)
    Υ, errorΥ = jackknife(ΥΥ, jkBins)
    V, errorV = jackknife(VV, jkBins)
    V2, errorV2 = jackknife(VV2, jkBins)
    
    return [β, E, errorE, E2, errorE2, M, errorM, χ, errorχ, Υ, errorΥ, V, errorV, V2, errorV2]
end



function simulate_2dXYmodel_singlecluster(L::Int64, T::Array{Array{Int64,1}, 3}, Λ::Array{Float64, 3}, β::Float64, thermTime::Int64,  measureTime::Int64, interTime::Int64, jkBins::Int64)
    Λ = thermalisation(L, T, Λ, β, thermTime, "singlecluster")
    
    EE = Float64[]
    EE2 = Float64[]
    MM = Float64[]
    χχ = Float64[]
    ΥΥ = Float64[]
    VV = Float64[]
    VV2 = Float64[]
    
    for t in 1:measureTime
        Λ, r, cluster, size, perimeter, sites = sweep_singlecluster(L, T, Λ, β)
        if t%interTime == 0
            e = energy(L, T, Λ)
            m = magnetisation(Λ)
            y = helicity(L, T, Λ, β)
            v = vorticity(L, Λ)
            
            push!(EE, e)
            push!(EE2, e*e)
            push!(MM, m)
            push!(χχ, m*m)
            push!(ΥΥ, y)
            push!(VV, v)
            push!(VV2, v*v)
        end
    end
    
    E, errorE = jackknife(EE, jkBins)
    E2, errorE2 = jackknife(EE2, jkBins)
    M, errorM = jackknife(MM, jkBins)
    χ, errorχ = jackknife(χχ, jkBins)
    Υ, errorΥ = jackknife(ΥΥ, jkBins)
    V, errorV = jackknife(VV, jkBins)
    V2, errorV2 = jackknife(VV2, jkBins)
    
    return [β, E, errorE, E2, errorE2, M, errorM, χ, errorχ, Υ, errorΥ, V, errorV, V2, errorV2]
end
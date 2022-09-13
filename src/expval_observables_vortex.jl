#=

=#

function simulate_2dXYmodel_cluster_vorticity(L::Int64, ADS::Array{Array{Int64, 1}, 3}, T::Array{Array{Int64,1}, 3}, Λ::Array{Float64, 3}, β::Float64, thermTime::Int64,  measureTime::Int64, interTime::Int64, jkBins::Int64)
    #=
    Functions that compute all observables for only one β with the multicluster algorithm. 
    =#
    
    percent = Int64(measureTime/10)
    totalPercent = 0
    
    Λ = thermalisation(L, T, Λ, β, thermTime)
    
    VVc  = Float64[]
    VVc2 = Float64[]
    CCv  = Float64[]
    CCv2 = Float64[]
    
    print("        Percentage of measurement: [  ", totalPercent, " ")
    itime = time()
    for t in 1:measureTime
        Λ, r, clusters, sizes, perimeters, numbers, listClusters, listBoundaries = sweep_multicluster(L, T, Λ, β)
        if t%interTime == 0
            vc, vc2 = cluster_vorticity(L, ADS, Λ, r, clusters,  listClusters, listBoundaries)
            cv, cv2 = vc/numbers, vc2/numbers
            
            push!(VVc, vc)
            push!(VVc2, vc2)
            push!(CCv, cv)
            push!(CCv2, cv2)
        end
        
        if t%percent==0  
            totalPercent = totalPercent + 10
            print(" ", totalPercent, " ")
        end
    end
    ftime = time()
    print("]\n        Measure time: ",(ftime - itime)/60 , " minutes \n")
    
    Vc, errorVc = jackknife(VVc, jkBins)
    Vc2, errorVc2 = jackknife(VVc2, jkBins)
    Cv, errorCv = jackknife(CCv, jkBins)
    Cv2, errorCv2 = jackknife(CCv2, jkBins)
    
    
    return [β, Vc, errorVc, Vc2, errorVc2, Cv, errorCv, Cv2, errorCv2]
end

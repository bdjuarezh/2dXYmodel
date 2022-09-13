function histogram_cluster_obs(L::Int64, T::Array{Array{Int64,1}, 3}, Λ::Array{Float64, 3}, β::Float64, thermTime::Int64, measureTime::Int64)
    Vol = L*L
    
    freqCN = zeros(Int64, Vol)
    freqCS = zeros(Int64, Vol)
    freqCP = zeros(Int64, Vol)
    
    sumCPvCS   = zeros(Int64, Vol)

    Λ = thermalisation(L, T, Λ, β, thermTime)
    
    percent = Int64(measureTime/10)
    totalPercent = 0
    
    print("        Percentage of measurement: [  ", totalPercent, " ")
    for t in 1:measureTime
        Λ, r, clusters, sizes, perimeters, numbers, listClusters, listBoundaries = sweep_multicluster(L, T, Λ, β)
        
        freqCN[numbers] = freqCN[numbers] + 1
        
        for i in 1:Vol
            s = sizes[i]
            p = perimeters[i]
            if s != 0          # if s = 0 then  p = 0
                if p!=0
                freqCS[s] = freqCS[s] + 1    
                freqCP[p] = freqCP[p] + 1    
                sumCPvCS[s] = sumCPvCS[s] + p
                else
                    print(i, " ", p, "  ", s, " ", sizes, " " ,perimeters)
                end
            end
        end
        if t%percent==0  
            totalPercent = totalPercent + 10
            print(" ", totalPercent, " ")
        end
    end
    print("]\n")
    
    meanCPvCS = zeros(Float64, Vol) #i:= cluster size, CPvCS[i]:= cluster perimeter
    
    for i in 1:Vol
        if freqCS[i] != 0
            meanCPvCS[i] = sumCPvCS[i]/freqCS[i]
        else 
            meanCPvCS[i] = 0
        end
    end
    return freqCN, freqCS, freqCP, meanCPvCS
end

function expval_cluster_obs(L::Int64, T::Array{Array{Int64,1}, 3}, Λ::Array{Float64, 3}, β::Float64, thermTime::Int64,  measureTime::Int64, interTime::Int64, jkBins::Int64)
    #=
    Functions that compute all cluster observables for only one β with the multicluster algorithm. 
    =#
    Λ = thermalisation(L, T, Λ, β, thermTime)
    
    NN = Float64[]
    NN2 = Float64[]
    SS  = Float64[]
    SS2  = Float64[]
    PP = Float64[]
    PP2 = Float64[]
    
    
    
    for t in 1:measureTime
        Λ, r, clusters, sizes, perimeters, numbers, listClusters, listBoundaries = sweep_multicluster(L, T, Λ, β)
        
        
        
        if t%interTime == 0
            s = sum(sizes)/numbers
            p = sum(perimeters)/numbers
            
            push!(NN, numbers)
            push!(NN2, numbers*numbers)
            push!(SS, s)
            push!(SS2, s*s)
            push!(PP, p)
            push!(PP2, p*p)
        end
    end
    
    N, errorN = jackknife(NN, jkBins)
    N2, errorN2 = jackknife(NN2, jkBins)
    S, errorS = jackknife(SS, jkBins)
    S2, errorS2 = jackknife(SS2, jkBins)
    P, errorP = jackknife(PP, jkBins)
    P2, errorP2 = jackknife(PP2, jkBins)

    return [β, N, errorN, N2, errorN2, S, errorS, S2, errorS2, P, errorP, P2, errorP2]
end 




function simulate_2dXYmodel_cluster_obs(L::Int64, T::Array{Array{Int64,1}, 3}, Λ::Array{Float64, 3}, β::Float64, thermTime::Int64,  measureTime::Int64, interTime::Int64, jkBins::Int64)
    #=
    Functions that compute all cluster observables for only one β with the multicluster algorithm. 
    =#
    Λ = thermalisation(L, T, Λ, β, thermTime)
    
    NN = Float64[]
    NN2 = Float64[]
    SS  = Float64[]
    SS2  = Float64[]
    PP = Float64[]
    PP2 = Float64[]
    
     
    Vol = L*L
    freqCN = zeros(Int64, Vol)
    freqCS = zeros(Int64, Vol)
    freqCP = zeros(Int64, Vol)
    sumCPvCS   = zeros(Int64, Vol)
    
    
    for t in 1:measureTime
        Λ, r, clusters, sizes, perimeters, numbers, listClusters, listBoundaries = sweep_multicluster(L, T, Λ, β)
        
        freqCN[numbers] = freqCN[numbers] + 1
        
        for i in 1:Vol
            s = sizes[i]
            p = perimeters[i]
            if s != 0          # if s = 0 then  p = 0
                if p!=0
                freqCS[s] = freqCS[s] + 1    
                freqCP[p] = freqCP[p] + 1    
                sumCPvCS[s] = sumCPvCS[s] + p
                else
                         print(i, " ", p, "  ", s, " ", sizes, " " ,perimeters)
                end
            end
        end
        
         
        
        if t%interTime == 0
            s = sum(sizes)/numbers
            p = sum(perimeters)/numbers
            
            push!(NN, numbers)
            push!(NN2, numbers*numbers)
            push!(SS, s)
            push!(SS2, s*s)
            push!(PP, p)
            push!(PP2, p*p)
        end
    end
    
    meanCPvCS = zeros(Float64, Vol) #i:= cluster size, CPvCS[i]:= cluster perimeter
    for i in 1:Vol
            if freqCS[i] != 0
            meanCPvCS[i] = sumCPvCS[i]/freqCS[i]
        else 
            meanCPvCS[i] = 0
        end
    end
    
    N, errorN = jackknife(NN, jkBins)
    N2, errorN2 = jackknife(NN2, jkBins)
    S, errorS = jackknife(SS, jkBins)
    S2, errorS2 = jackknife(SS2, jkBins)
    P, errorP = jackknife(PP, jkBins)
    P2, errorP2 = jackknife(PP2, jkBins)

    return [β, N, errorN, N2, errorN2, S, errorS, S2, errorS2, P, errorP, P2, errorP2, freqCN, freqCS, freqCP, meanCPvCS]
end 

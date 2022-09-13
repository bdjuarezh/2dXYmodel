function jackknife(data::Array{Float64,1}, jkBins::Int64)
    #=
    Compute the Jackknife error from a sample of data. 
    
    CONFS: length of array, or number of configurations measured
    BIN_SIZE: number of elements of each bin
    =#
    
    CONF_NUMBER = length(data) 
    BIN_SIZE = Int(CONF_NUMBER/jkBins)
    
    obs = sum(data)/CONF_NUMBER
    obsError = 0.0
    totalSum = 0.0
    
    if CONF_NUMBER%BIN_SIZE == 0
        for i in 1:jkBins
            partialObs = sum(data[(i - 1)*BIN_SIZE + 1:i*BIN_SIZE])/BIN_SIZE
            partialDif = (partialObs - obs)
            totalSum = totalSum + partialDif*partialDif
        end
        obsError = sqrt(totalSum/(jkBins*(jkBins - 1)))
    else
        print("Ajustar el número de bins para que no haya residuo \n
            por ahora el error será cero al estar mal calculado")
    end
    return obs, obsError
end
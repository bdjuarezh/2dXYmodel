#=                    2d XY model definition

This functions define the two dimensional XY model. We include functions that
constructs a square lattice, from a hot start and cold start. We define the probability
function for spin flip and a short function for dot product with the Wolff's vector.

Notation
L := Lengh of a side lattice
V := L*L
Λ := Matrix of spins 
r := Wolff vector
β := 1/T                    
=#



#------------------------------BASE MODULES------------------------------#
function hot_start(L::Int64) 
    #= 
    Make a configuration with side length L 
    and random spin orientation 
    =#
    Λ = zeros(Float64, L, L, 2)  
    for j in 1:L, i in 1:L
        θ = 2*pi*rand() 
        Λ[i, j, 1] = cos(θ)
        Λ[i, j, 2] = sin(θ)
    end
    return Λ
end


function cold_start(L::Int64, θ::Float64) 
    #= 
    Make a configuration with side length L 
    and fixed spin orientation
    =#
    Λ = zeros(Float64, L, L, 2)
    for j in 1:L, i in 1:L
        Λ[i, j, 1] = cos(θ)
        Λ[i, j, 2] = sin(θ)
    end
    return Λ
end


#Compute the probability for a bond
probability(β::Float64, K::Float64) = Float64(1 - exp(2*β*min(0.0, -K)))

# Make the dot product betwen Wolff vector and only one spin 
dot_wolff(Λ::Array{Float64,3}, r::Array{Float64,1}, i::Int64, j::Int64) = Float64(r[1]*Λ[i, j, 1] + r[2]*Λ[i, j, 2])


function one_flip(Λ::Array{Float64,3}, r::Array{Float64,1}, i::Int64, j::Int64, dotWolff::Float64)
    #= 
    Make the flip of only one spin respect the Wolff's line
    =#
    xNewComp = Float64(Λ[i, j, 1] - 2*r[1]*dotWolff)
    yNewComp = Float64(Λ[i, j, 2] - 2*r[2]*dotWolff)
    return xNewComp, yNewComp
end


function normalize_configuration(L::Int64, Λ::Array{Float64,3})
    #= 
    Normalizes all the spins of the configuration Λ 
    =#
    for j in 1:L, i in 1:L
        x = Λ[i, j, 1] 
        y = Λ[i, j, 2] 
    
        norm = sqrt(x*x + y*y)
           
        Λ[i, j, 1] = Float64(Λ[i, j, 1]/norm)
        Λ[i, j, 2] = Float64(Λ[i, j, 2]/norm)
    end
    return Λ
end
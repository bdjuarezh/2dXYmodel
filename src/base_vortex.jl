#=                   Module for find vortex in a configuration
This functions found the vortex and antivortex of a spin configuration and returns
the plaquetes sites where there are.  
find_vortex() by default is electer the choice "dual" that return the sites in the dual lattice that identifies
the plaquetes where there are. If we choice the option "real" then returns the real sites where the V and AV exists. 
We can provide to the function all the conf or only a list of sites in the conf. 


Notation
L := Lengh of a lattice side
V := L*L
Λ := Matrix of spins
Vortex := array with the positions of all vortex
Antivortex := array with the positions of all anti-vortex                   
=#

function find_vortex(L::Int64, Λ::Array{Float64, 3}, typeResult::String="dual")
    #=
    Find vortex and antivortex in a configuration with the standar definition in the with the position
    of the vortex and antivortex in the real lattice
    =#
    δ = 0.000005      
    Vortex = Array[]
    Antivortex = Array[]
    
    if typeResult=="dual"
        for j in 1:L, i in 1:L
            x1, y1 = Λ[i, j, 1], Λ[i, j, 2]                        
            x2, y2 = Λ[i, j%L + 1, 1], Λ[i, j%L + 1, 2]
            x3, y3 = Λ[i%L + 1, j%L + 1, 1], Λ[i%L + 1, j%L + 1, 2]
            x4, y4 = Λ[i%L + 1, j, 1], Λ[i%L + 1, j, 2]
            
            cos1 = (x1*x2 + y1*y2)
            cos2 = (x2*x3 + y2*y3)
            cos3 = (x3*x4 + y3*y4)
            cos4 = (x4*x1 + y4*y1)
            
            if abs(cos1)>1 print("."); cos1=1 end
            if abs(cos2)>1 print("."); cos2=1 end
            if abs(cos3)>1 print("."); cos3=1 end
            if abs(cos4)>1 print("."); cos4=1 end
            
            m1 = acos(cos1) 
            m2 = acos(cos2)
            m3 = acos(cos3)
            m4 = acos(cos4)
            
            s1 = sign(x1*y2 - x2*y1)
            s2 = sign(x2*y3 - x3*y2)
            s3 = sign(x3*y4 - x4*y3)
            s4 = sign(x4*y1 - x1*y4)
            
            θ = s1*m1 + s2*m2 + s3*m3 + s4*m4
            
            if  abs(2*pi + θ) < δ
                push!(Antivortex, [i, j])
            elseif abs(2*pi - θ) < δ
                push!(Vortex, [i, j])
            end
        end
    elseif typeResult=="real"
        for j in 1:L, i in 1:L
            x1, y1 = Λ[i, j, 1], Λ[i, j, 2]                        
            x2, y2 = Λ[i, j%L + 1, 1], Λ[i, j%L + 1, 2]
            x3, y3 = Λ[i%L + 1, j%L + 1, 1], Λ[i%L + 1, j%L + 1, 2]
            x4, y4 = Λ[i%L + 1, j, 1], Λ[i%L + 1, j, 2]
            
            cos1 = (x1*x2 + y1*y2)
            cos2 = (x2*x3 + y2*y3)
            cos3 = (x3*x4 + y3*y4)
            cos4 = (x4*x1 + y4*y1)
            
            if abs(cos1)>1 print("."); cos1=1 end
            if abs(cos2)>1 print("."); cos2=1 end
            if abs(cos3)>1 print("."); cos3=1 end
            if abs(cos4)>1 print("."); cos4=1 end
            
            m1 = acos(cos1) 
            m2 = acos(cos2)
            m3 = acos(cos3)
            m4 = acos(cos4)
            
            s1 = sign(x1*y2 - x2*y1)
            s2 = sign(x2*y3 - x3*y2)
            s3 = sign(x3*y4 - x4*y3)
            s4 = sign(x4*y1 - x1*y4)
            
            θ = s1*m1 + s2*m2 + s3*m3 + s4*m4
            
            if  abs(2*pi + θ) < δ
                push!(Antivortex, [i + 1/2, j + 1/2])
            elseif abs(2*pi - θ) < δ
                push!(Vortex, [i + 1/2, j + 1/2])
       
            end
        end
    end
    return Vortex, Antivortex
end

function find_vortex(L::Int64, Λ::Array{Float64, 3}, boundarySites::Array{Any, 1}, typeResult::String="dual")
    #=
    Find the vortex and antivortex only in the sites (that defines a list of plaquetes)
    that are provide.
    =#
    δ = 0.000005      
    Vortex = Array[] 
    Antivortex = Array[]
    
   if typeResult=="dual"
        
        for site in boundarySites
            i, j = site
            x1, y1 = Λ[i, j, 1], Λ[i, j, 2]                        
            x2, y2 = Λ[i, j%L + 1, 1], Λ[i, j%L + 1, 2]
            x3, y3 = Λ[i%L + 1, j%L + 1, 1], Λ[i%L + 1, j%L + 1, 2]
            x4, y4 = Λ[i%L + 1, j, 1], Λ[i%L + 1, j, 2]
            
            cos1 = (x1*x2 + y1*y2)
            cos2 = (x2*x3 + y2*y3)
            cos3 = (x3*x4 + y3*y4)
            cos4 = (x4*x1 + y4*y1)
            
            if abs(cos1)>1 print("."); cos1=1 end
            if abs(cos2)>1 print("."); cos2=1 end
            if abs(cos3)>1 print("."); cos3=1 end
            if abs(cos4)>1 print("."); cos4=1 end
            
            m1 = acos(cos1) 
            m2 = acos(cos2)
            m3 = acos(cos3)
            m4 = acos(cos4)
            
            s1 = sign(x1*y2 - x2*y1)
            s2 = sign(x2*y3 - x3*y2)
            s3 = sign(x3*y4 - x4*y3)
            s4 = sign(x4*y1 - x1*y4)
            
            θ = s1*m1 + s2*m2 + s3*m3 + s4*m4
        
            if  abs(2*pi + θ) < δ
                push!(Antivortex, [i, j])
                
            elseif abs(2*pi - θ) < δ
                push!(Vortex, [i, j])
            end
        end
        return Vortex, Antivortex
        
    elseif typeResult=="real"
        
        for site in boundarySites
            i, j = site
            x1, y1 = Λ[i, j, 1], Λ[i, j, 2]                        # s1 --- s2
            x2, y2 = Λ[i, j%L + 1, 1], Λ[i, j%L + 1, 2]            #  |     |
            x3, y3 = Λ[i%L + 1, j%L + 1, 1], Λ[i%L + 1, j%L + 1, 2]#  |     |
            x4, y4 = Λ[i%L + 1, j, 1], Λ[i%L + 1, j, 2]            # s4 --- s3
            
            cos1 = (x1*x2 + y1*y2)
            cos2 = (x2*x3 + y2*y3)
            cos3 = (x3*x4 + y3*y4)
            cos4 = (x4*x1 + y4*y1)
            
            if abs(cos1)>1 print("."); cos1=1 end
            if abs(cos2)>1 print("."); cos2=1 end
            if abs(cos3)>1 print("."); cos3=1 end
            if abs(cos4)>1 print("."); cos4=1 end
            
            m1 = acos(cos1) 
            m2 = acos(cos2)
            m3 = acos(cos3)
            m4 = acos(cos4)
            
            s1 = sign(x1*y2 - x2*y1)
            s2 = sign(x2*y3 - x3*y2)
            s3 = sign(x3*y4 - x4*y3)
            s4 = sign(x4*y1 - x1*y4)
            
            θ = s1*m1 + s2*m2 + s3*m3 + s4*m4
        
            if  abs(2*pi + θ) < δ
                push!(Antivortex, [i+1/2, j+1/2])
                
            elseif abs(2*pi - θ) < δ
                push!(Vortex, [i+1/2, j+1/2])
            end
        end
        return Vortex, Antivortex
    end
end
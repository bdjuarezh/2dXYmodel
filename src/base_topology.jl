#=
Pasa de la representación de un índice a la representación de dos índices.
La convención es para un solo índice será n y se recorrerá como:

1  5  9   13              
2  6  10  14              
3  7  11  15    
4  8  12  16.

El caso de dos índices serán (i, j), donde i será el número de renglon y j el número de columna

(1, 1)   (1, 2)   (1, 3)   (1, 4)
(2, 1)   (2, 2)   (2, 3)   (2, 4)
(3, 1)   (3, 2)   (3, 3)   (3, 4)
(4, 1)   (4. 2)   (4, 3)   (4, 4)
  
La regla de equivalencia será que

n = N*(j - 1) + i 

y 

i = n%N 
j = n%N

=#

n_to_ij(n::Int64, L::Int64) = [Int64(mod(n - 1, L) + 1), Int64(floor((n - 1)/L) + 1)]
ij_to_n(i::Int64, j::Int64, L::Int64) = Int64(L*(j - 1) + i)


function topology(L::Int64)#, nIndex::Int64=2) 
    #= 
    Function that constructs the topology of system. 
    If the value of nIndex=2 then make the set of neighbors from the site (i, j),
    where i:=rows/lines,  j:=columns.
    If the value of nIndex=1 then make the set of neighbors from the site n, where 
    n:=unique index from 1 to L^2.
    =#
    
   
        T = Array{Array{Int64,1}, 3}(undef, L, L,8)
        
        for j in 1:L, i in 1:L                                       
            iUp = mod(mod(i, L) - 2, L) + 1   
            jLeft = mod(mod(j, L) - 2, L) + 1
            iDown = i%L + 1      
            jRight = j%L + 1

            
            T[i, j, 1] = [iUp, j]   #              (iUp, j)
            T[i, j, 2] = [i, jLeft] #                  |
            T[i, j, 3] = [iDown, j] # (i, jLeft)----(i, j)----(i, jRight)
            T[i, j, 4] = [i,jRight] #                  |
                                    #              (iDown, j)
            
            
            T[i, j, 5] = [iUp, jLeft]    # (iUp,jLeft)-----(,)-----(iUp,jRight)
            T[i, j, 6] = [iDown, jLeft]  #     |            |            |
            T[i, j, 7] = [iDown, jRight] #    (,)---------(i, j)--------(,)
            T[i, j, 8] = [iUp,jRight]    #     |            |            |
                                         # (iDown,jLeft)---(,)-----(iDown,jRight)
            
        end
        return T
    
     #if nIndex == 2
     #=       
    elseif nIndex == 1
        L2 = L*L
        T = Array{Array{Int64, 1}, 1}(undef, L2)
        
        for n in 1:L2
            i, j = n_to_ij(n, L) 
                                       #              (iUp, j)
            iUp = mod(i - 2, L) + 1    #                  |
            jLeft = mod(j - 2, L) + 1  # (i, jLeft)----(i, j)----(i, jRight)
            iDown = mod(i, L) + 1      #                  |
            jRight = mod(j, L) + 1     #              (iDown, j)
            
            up = ij_to_n(iUp, j, L)
            left = ij_to_n(i, jLeft, L)
            down = ij_to_n(iDown, j, L)
            right = ij_to_n(i, jRight, L)
        
            T[n] = [up, left, down, right]
        end

        return T
    end
    =#

end

#=
Gives us the neighbors of site (i, j); the order is: [[up], [left], [down], [right]] and
each element of the array is another array of lenght 2, witht the index [i,j].
=#
neighbors(T::Array{Array{Int64,1}, 3}, i::Int64, j::Int64) = [T[i, j, 1], T[i, j, 2], T[i, j, 3], T[i, j, 4]] 

second_neighbors(T::Array{Array{Int64,1}, 3}, i::Int64, j::Int64) = [T[i, j, 5], T[i, j, 6], T[i, j, 7], T[i, j, 8]] 

all_neighbors(T::Array{Array{Int64,1}, 3}, i::Int64, j::Int64) = [T[i, j, 1], T[i, j, 2], T[i, j, 3], T[i, j, 4], T[i, j, 5], T[i, j, 6], T[i, j, 7], T[i, j, 8]] 

#=
Another definition of function that gives us the neighbors of a 1index site n.
Returns an arrays with the four neighbors of the site n
=#
neighbors(T::Array{Array{Int64,1}, 1}, n::Int64) = T[n]


#=
funciones que no ayudaran a identificar los sitios de la 
matriz dual que rodeen un elemento de la matriz normal y 
viceversa

(1,1)-------(1,2)-------(1,3)-------
  |           |           |
  |   [1,1]   |   [1,2]   |   [1,3]
  |           |           |
(2,1)-------(2,2)-------(2,3)-------
  |           |           |
  |   [2,1]   |   [2,2]   |   [2,3]
  |           |           |
(3,1)-------(3,2)-------(3,2)-------
  |           |           |
  |   [3,1]   |   [3,2]   |   [3,3]
  |           |           |
=#
function adjacent_dual_sites(L::Int64)
    adSites = Array{Array{Int64,1},3}(undef, (L, L, 4))
    for j in 1:L, i in 1:L                     #j:left to right and i:up to down 
        adSites[i, j, 1] = [mod(mod(i, L) - 2, L) + 1 , mod(mod(j, L) - 2, L) + 1]  #up left
        adSites[i, j, 2] = [mod(mod(i, L) - 2, L) + 1 , j]                          #up right
        adSites[i, j, 3] = [i, j]                                                   #down right
        adSites[i, j, 4] = [i, mod(mod(j, L) - 2, L) + 1]                           #down left
    end
    return adSites
end 

function normal_to_dual(adSites::Array{Array{Int64,1},3}, i::Int64, j::Int64)
    return adSites[i, j, 1], adSites[i, j, 2], adSites[i, j, 3], adSites[i, j, 4]
end


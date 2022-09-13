include("./simplevorticity_MOD.jl")
include("./2indexto1index_MOD.jl")

#=
Este módulo encuentra los pares ligados de vortices. Es decir parejas V-AV 
que comparten los mismos dos clusters. La función principarl es bond_pairs()
y sus subrutinas son find_vortex_simplevorticity() y find_pairs()

CM : CLUSTER MATRIX
LL : list Of Label
LC : list of Cluster
RC : reference Configuration
RCp : reference Configuration partial 
nC : number of Clusters
=#


function control_real_vortex(N, M) #ABSOLUTE NEW PART!!!!!!!!!!!!!!!!!!!!!1
    RMV = [[] for i in 1:N*N]
    l = 1
    
    V, A = find_vortex_dual(N, M)
    
    if length(V)==length(A)
        for i in 1:length(V)
            v = V[i]
            a = A[i]
            nv = ij_to_n(v[1], v[2], N)
            na = ij_to_n(a[1], a[2], N)
                    
            push!(RMV[nv], 1)
            push!(RMV[na], 1)
        end
    else
        print("Algo va mal")
    end
    l = l + 1
        #print(RMV)
    return RMV 
end


function write_vortex(N, r, LL, LC, RC, RMV)
    #=
    Escribe en cada entrada de una matriz corresponde a las plaquetas, las etiquetas 
    de los clusters a los que esta plaqueta contribuye. 
    =#
    MV = [[] for i in 1:N*N]
    l = 1
    
    for elements in LC
        RCp = partial_flip(RC, elements, r)
        V, A = find_vortex_dual(N, RCp)
    
        if length(V)==length(A)
            for i in 1:length(V)
                v = V[i]
                a = A[i]
                nv = ij_to_n(v[1], v[2], N)
                na = ij_to_n(a[1], a[2], N)
                
                #print(RMV[nv])
                if (RMV[nv]==[1]) push!(MV[nv], LL[l]) end  
                if (RMV[na]==[1]) push!(MV[na], LL[l]) end  
                    
            end
        else
            print("Algo va mal")
        end
        l = l + 1
    end
    return MV
end

function count_vortex(N, MV)
    #=
    Check sites of vortex list and find  pairs vortex anti-vortex
    =#
    aux = zeros(Int64, N*N)
    nB = 0
    nTot = 0
    
    for k in 1:N*N
        site = MV[k]
        if length(site)!=0 
            nTot = nTot + 0.5
            if aux[k]==0
                match = 0
                for l in (k+1):N*N
                    if site==MV[l] && match==0 
                        nB = nB + 1.0
                        match = match + 1
                        aux[l] = 1
                    end
                end
            end
        end
    end
    nF = nTot - nB
    if nF<0
        print("Están saliendo números negativos")
        print(nTot, "  ",nB, "\n")
        print(MV, "\n\n\n")
    end
    return nB, nF
end








function bounded_free_vortex(N, M, CM, r, nC)
    #=
    return number of bounded and free vortex
    =#
    RMV = control_real_vortex(N, M) #NEW PART!!!!!!!
    
    LL = list_of_labels(N, CM)
    LC = list_of_clusters(N, CM, LL)
    RC = reference_conf(N, M, r)
    MV = write_vortex(N, r, LL, LC, RC, RMV)
    nB, nF = count_vortex(N, MV)
    
    #=print(CM,"\n")
    print(RMV, "\n")
    print(MV,"\n\n\n")
    print("ligados: ", nB, "\n")
    print("libres: ", nF, "\n")
    =#
    nc = length(LC)
    print(length(LC), "\n")
    return nB, nF, nc #CAMBIO "nc"
end

function ev_bond_pairs(N, G, M, β, tCalc, inter)
    #=
    Returns the expected value of the number of bonded and free vortex
    =#
    Nconf = tCalc/inter

    B = 0
    B2 = 0
    eB = 0
    
    F = 0
    F2 = 0
    eF = 0
    
    percent = tCalc/10
    cont = 0
    print("    Percentage of measurement: ",cont, "%\n")
    
    for t in 1:tCalc
        M, sizes, CM, r, perimeters, nC = step_sw(N, G, M, β)
        if t%inter == 0
            nB, nF = bounded_free_vortex(N, M, CM, r, nC)
            
            B = B + nB
            B2 = B2 + nB*nB
            
            F = F + nF
            F2 = F2 + nF*nF
        end
        if t%percent ==0
            cont = cont + 10
            print("    Percentage of measurement: ",cont, "%\n")
        end
    end
    
    B = B/Nconf # <B>
    B2 = B2/Nconf # <B^2>
    eB = sqrt((B2 - B*B)/Nconf)
    
    F = F/Nconf
    F2 = F2/Nconf
    eF = sqrt((F2 - F*F)/Nconf)
    
    return B, eB, F, eF
end



#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<NEW SECTION 25 feb<<<<<<<<<<<<<<<<<<<<<<<<<<#
function ev_vorticity(N, G, M, β, tCalc, inter)
    #=
    Returns the expected value of the number of bonded and free vortex
    =#
    Nconf = tCalc/inter

    V = 0
    V2 = 0
    V3 = 0
    V4 = 0
    
    #eV = 0
    
    #F = 0
    #F2 = 0
    #eF = 0
    
    percent = tCalc/10
    cont = 0
    print("    Percentage of measurement: ",cont, "%\n")
    
    for t in 1:tCalc
        M, sizes, CM, r, perimeters, nC = step_sw(N, G, M, β)
        if t%inter == 0
            nB, nF, nc = bounded_free_vortex(N, M, CM, r, nC)
            
            v = (nB + nF)/nC
            print("\n",nC, "\n\n")
            
            V = V + v
            V2 = V2 + v*v
            V3 = V3 + v*v*v
            V4 = V4 + v*v*v*v
            
            #F = F + nF
            #F2 = F2 + nF*nF
        end
        if t%percent ==0
            cont = cont + 10
            print("    Percentage of measurement: ",cont, "%\n")
        end
    end
    
    V = V/Nconf # <B>
    V2 = V2/Nconf # <B^2>
    V3 = V3/Nconf # <B^3>
    V4 = V4/Nconf # <B^4>
    
    #eB = sqrt((B2 - B*B)/Nconf)
    
    #F = F/Nconf
    #F2 = F2/Nconf
    #eF = sqrt((F2 - F*F)/Nconf)
    
    return V, V2, V3, V4
end

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>New section 27ene20<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#Función de prueba para ver qué está pasando
function simple_vorticity_tot(N, M, CLUSTERS, r, nClusters)
    listLabel = list_of_labels(N, CLUSTERS)
    listCluster = list_of_clusters(N, CLUSTERS, listLabel)
    refConf = reference_conf(N, M, r)
    totV = 0
    
    for elements in listCluster
        refConfp = partial_flip(refConf, elements, r)
        #plot_vortex(N, refConfp, elements, "refConfp", true, false) #experiment
        V, A = find_vortex(N, refConfp)
        totV = totV + length(V) # + length(A)
    end
    return totV
end


#nueva moficación 4feb20>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
function diff_VcV(N, M, CM, r, nC)
    #=
    return diference between Vc-V for a configuration
    =#
    v, a = find_vortex(N, M) 
    V = length(v) # + length(a) #Estamos contando las plaquetas ocupadas, es lo importante
    
    Vc = simple_vorticity_tot(N, M, CM, r, nC)
    #LL = list_of_labels(N, CM)
    #LC = list_of_clusters(N, CM, LL)
    #RC = reference_conf(N, M, r)
    #MV = write_vortex(N, r, LL, LC, RC)
    #nB, nF = count_vortex(N, MV)
    #Vc = nB + nF
    print("Vc = ", Vc, "\n")
    print("V = ",V, "\n\n")
    return Vc, V
end


function ev_diff_VcV(N, G, M, β, tCalc, inter)
    #=
    Returns the expected value of the differece between Vc-V for a beta
    =#
    Nconf = tCalc/inter
    percent = tCalc/10
    cont = 0
    
    Vc = 0
    Vc2 = 0
    eVc = 0
    
    V = 0
    V2 = 0
    eV = 0
    
    print("    Percentage of measurement: ",cont, "%\n")
    for t in 1:tCalc
        M, sizes, CM, r, perimeters, nC = step_sw(N, G, M, β)
        if t%inter == 0
            vc, v = diff_VcV(N, M, CM, r, nC) #bounded_free_vortex(N, M, CM, r, nC)
            
            Vc = Vc + vc
            Vc2 = Vc2 + vc*vc
            
            V = V + v
            V2 = V2 + v*v
        end
        if t%percent ==0
            cont = cont + 10
            print("    Percentage of measurement: ",cont, "%\n")
        end
    end
    
    Vc = Vc/Nconf # <B>
    Vc2 = Vc2/Nconf # <B^2>
    eVc = sqrt((Vc2 - Vc*Vc)/Nconf)
    
    V = V/Nconf # <B>
    V2 = V2/Nconf # <B^2>
    eV = sqrt((V2 - V*V)/Nconf)
    
    print(Vc, "\n")
    print(V, "\n\n\n\n")
    
    return Vc, eVc, V, eV
end




include("./simplevorticity_MOD.jl")
include("./2indexto1index_MOD.jl")
#=
Este módulo encuentra los pares ligados de vortices. Es decir parejas V-AV 
que comparten los mismos dos clusters. La función principarl es bond_pairs()
y sus subrutinas son find_vortex_simplevorticity() y find_pairs()

CM : CLUSTER
LL : listLabel
LC : listCluster
RC : refConf
RCp : refConfp
nC : nCluster
=#


function find_vortex_simplevorticity(N, M, CM, r, nC)
    #=
    returns VV = a list of list of vortex from each cluster
    AA = a list of list of anti-vortex from each cluster
    LL = a list of labels of all clusters
    =#
    LL = list_of_labels(N, CM)
    LC = list_of_clusters(N, CM, LL)
    RC = reference_conf(N, M, r)
    VV = []
    AA = []
    
    for elements in LC
        RCp = partial_flip(RC, elements, r)
        V, A = find_vortex_dual(N, RCp)
        push!(VV, V)
        push!(AA, A)
    end
    return VV, AA, LL
end


function find_pairs(N, VV, AA, LL)
    #=
    Usa las listas de vórtices y antivórtices junto con la de 
    etiquetas para relacionar los vortices que pertencen a dos 
    cluster distintos y cuando encuentra pares V-AV que pertenecen 
    a los dos mismos clusters suma 1 a la contribución 
    
    MV : Vortex dual matrix
    =#
    MV = [[] for i in 1:N*N];
    MA = [[] for i in 1:N*N];
    
    n = 0
    if length(VV) ==length(AA) == length(LL)
        #print("todo bien\n")
        
        for i in 1:length(VV)
            V = VV[i]
            A = AA[i]
            L = LL[i]
            
            if length(V)==length(A) #prueba
                
                for j in 1:length(V)
                    v = V[j]
                    a = A[j]
                    nv = ij_to_n(v[1], v[2], N)
                    na = ij_to_n(a[1], a[2], N)
                    
                    push!(MV[nv], L)
                    push!(MA[na], L)
                    
                end
            else
                print("la rutina de identificación de vortices no funciona bien")
            end
        end
    else
        print("algo va mal")
    end
    return MV, MA
end

function count_pairs(MV, MA)
    #=
    Check sites of vortex list and find  pairs vortex anti-vortex
    =#
    nB = 0
    nTot = 0
    for k in 1:N*N
        site = MV[k]
        if length(site)!=0
            #print(site)
            nTot = nTot + 0.5
            for l in k+1:N*N
                if site == MA[l]
                    nB = nB + 1
                end
            end
        end
    end
    nNotB = nTot - nB
    #print("Total pairs in a configuration = ", nTot, "\n", "bounded = ", nB,"\n", "not bounded = ", nNotB, "\n\n")
    return nB, nNotB
end





function bond_pairs(N, M, CM, r, nC)
    #= 
    Main functions that find all
    bond pairs of V_AV of a configuration
    =#
    VV, AA, LL = find_vortex_simplevorticity(N, M, CM, r, nC)
    MV, MA = find_pairs(N, VV, AA, LL)
    nB, nNotB = count_pairs(MV, MA)
    return nB, nNotB
end

function ev_bond_pairs(N, G, M, β, tCalc, inter)
    #=
    Returns the expected value of the number of bonded pairs
    =#
    Nconf = tCalc/inter
    
    bpN = 0
    bpN2 = 0
    ebpN = 0
    
    NOTbpN = 0
    NOTbpN2 = 0
    eNOTbpN = 0
    
    percent = tCalc/10
    cont = 0
    print("    Percent of measure: ",cont, "%\n")
    
    for t in 1:tCalc
        M, sizes, CM, r, perimeters, nC = step_sw(N, G, M, β)
        if t%inter ==0
            nB, nNotB = bond_pairs(N, M, CM, r, nC)
            
            bpN = bpN + nB
            bpN2 = bpN2 + nB*nB
            
            NOTbpN = NOTbpN + nNotB
            NOTbpN2 = NOTbpN2 + nNotB*nNotB
            
        end
        if t%percent ==0
            cont = cont + 10
            print("    Percent of measure: ",cont, "%\n")
        end
    end
    bpN = bpN/Nconf
    bpN2 = bpN2/Nconf
    ebpN = sqrt((bpN2 - bpN*bpN)/(Nconf - 1))
    
    NOTbpN = NOTbpN/Nconf
    NOTbpN2 = NOTbpN2/Nconf
    eNOTbpN = sqrt((NOTbpN2 - NOTbpN*NOTbpN)/Nconf)
    return  bpN, ebpN, NOTbpN, eNOTbpN
end

        
include("./simplevorticity_MOD.jl")
include("./2indexto1index_MOD.jl")

#=
Este módulo encuentra los pares ligados de vortices. Es decir parejas V-AV 
que comparten los mismos dos clusters. La función principarl es bond_pairs()
y sus subrutinas son find_vortex_simplevorticity() y find_pairs()

CM : CLUSTER MATRIX
LL : list Of Label
LC : list of Cluster
RC : reference Configuration
RCp : reference Configuration partial 
nC : number of Clusters
=#

function write_vortex(N, r, LL, LC, RC)
    #=
    Escribe en cada entrada de una matriz corresponde a las plaquetas, las etiquetas 
    de los clusters a los que esta plaqueta contribuye. 
    =#
    MV = [[] for i in 1:N*N]
    l = 1
    
    for elements in LC
        RCp = partial_flip(RC, elements, r)
        V, A = find_vortex_dual(N, RCp)
    
        if length(V)==length(A)
            for i in 1:length(V)
                v = V[i]
                a = A[i]
                nv = ij_to_n(v[1], v[2], N)
                na = ij_to_n(a[1], a[2], N)
                    
                push!(MV[nv], LL[l])
                push!(MV[na], LL[l])
            end
        else
            print("Algo va mal")
        end
        l = l + 1
    end
    return MV
end

function count_vortex(N, MV)
    #=
    Check sites of vortex list and find  pairs vortex anti-vortex
    =#
    aux = zeros(Int64, N*N)
    nB = 0
    nTot = 0
    
    for k in 1:N*N
        site = MV[k]
        if length(site)!=0 
            nTot = nTot + 0.5
            if aux[k]==0
                match = 0
                for l in (k+1):N*N
                    if site==MV[l] && match==0 
                        nB = nB + 1.0
                        match = match + 1
                        aux[l] = 1
                    end
                end
            end
        end
    end
    nF = nTot - nB
    if nF<0
        print(nTot, "  ",nB, "\n")
        print(MV, "\n\n\n")
    end
    return nB, nF
end

function bounded_free_vortex(N, M, CM, r, nC)
    #=
    return number of bounded and free vortex
    =#
    LL = list_of_labels(N, CM)
    LC = list_of_clusters(N, CM, LL)
    RC = reference_conf(N, M, r)
    MV = write_vortex(N, r, LL, LC, RC)
    nB, nF = count_vortex(N, MV)
    return nB, nF
end

function ev_bond_pairs(N, G, M, β, tCalc, inter)
    #=
    Returns the expected value of the number of bonded and free vortex
    =#
    Nconf = tCalc/inter

    B = 0
    B2 = 0
    eB = 0
    
    F = 0
    F2 = 0
    eF = 0
    
    percent = tCalc/10
    cont = 0
    print("    Percentage of measurement: ",cont, "%\n")
    
    for t in 1:tCalc
        M, sizes, CM, r, perimeters, nC = step_sw(N, G, M, β)
        if t%inter == 0
            nB, nF = bounded_free_vortex(N, M, CM, r, nC)
            
            B = B + nB
            B2 = B2 + nB*nB
            
            F = F + nF
            F2 = F2 + nF*nF
        end
        if t%percent ==0
            cont = cont + 10
            print("    Percentage of measurement: ",cont, "%\n")
        end
    end
    
    B = B/Nconf # <B>
    B2 = B2/Nconf # <B^2>
    eB = sqrt((B2 - B*B)/Nconf)
    
    F = F/Nconf
    F2 = F2/Nconf
    eF = sqrt((F2 - F*F)/Nconf)
    
    return B, eB, F, eF
end





#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>New section 27ene20<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#Función de prueba para ver qué está pasando
function simple_vorticity_tot(N, M, CLUSTERS, r, nClusters)
    listLabel = list_of_labels(N, CLUSTERS)
    listCluster = list_of_clusters(N, CLUSTERS, listLabel)
    refConf = reference_conf(N, M, r)
    totV = 0
    
    for elements in listCluster
        refConfp = partial_flip(refConf, elements, r)
        #plot_vortex(N, refConfp, elements, "refConfp", true, false) #experiment
        V, A = find_vortex(N, refConfp)
        totV = totV + length(V) # + length(A)
    end
    return totV
end


#nueva moficación 4feb20>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
function diff_VcV(N, M, CM, r, nC)
    #=
    return diference between Vc-V for a configuration
    =#
    v, a = find_vortex(N, M) 
    V = length(v) # + length(a) #Estamos contando las plaquetas ocupadas, es lo importante
    
    Vc = simple_vorticity_tot(N, M, CM, r, nC)
    #LL = list_of_labels(N, CM)
    #LC = list_of_clusters(N, CM, LL)
    #RC = reference_conf(N, M, r)
    #MV = write_vortex(N, r, LL, LC, RC)
    #nB, nF = count_vortex(N, MV)
    #Vc = nB + nF
    print("Vc = ", Vc, "\n")
    print("V = ",V, "\n\n")
    return Vc, V
end


function ev_diff_VcV(N, G, M, β, tCalc, inter)
    #=
    Returns the expected value of the differece between Vc-V for a beta
    =#
    Nconf = tCalc/inter
    percent = tCalc/10
    cont = 0
    
    Vc = 0
    Vc2 = 0
    eVc = 0
    
    V = 0
    V2 = 0
    eV = 0
    
    print("    Percentage of measurement: ",cont, "%\n")
    for t in 1:tCalc
        M, sizes, CM, r, perimeters, nC = step_sw(N, G, M, β)
        if t%inter == 0
            vc, v = diff_VcV(N, M, CM, r, nC) #bounded_free_vortex(N, M, CM, r, nC)
            
            Vc = Vc + vc
            Vc2 = Vc2 + vc*vc
            
            V = V + v
            V2 = V2 + v*v
        end
        if t%percent ==0
            cont = cont + 10
            print("    Percentage of measurement: ",cont, "%\n")
        end
    end
    
    Vc = Vc/Nconf # <B>
    Vc2 = Vc2/Nconf # <B^2>
    eVc = sqrt((Vc2 - Vc*Vc)/Nconf)
    
    V = V/Nconf # <B>
    V2 = V2/Nconf # <B^2>
    eV = sqrt((V2 - V*V)/Nconf)
    
    print(Vc, "\n")
    print(V, "\n\n\n\n")
    
    return Vc, eVc, V, eV
end


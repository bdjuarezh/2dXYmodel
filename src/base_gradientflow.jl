#Gradient flow con el método Runge-Kutta de 4to orden RK4.
def gradient_flow_rk4(M, dt):
    N = len(M)
    MNew = np.zeros((N,N,2))
    
    for i in range(N):
        for j in range(N):
            x = (i, j, 0)
            y = (i, j, 1)
            xU, xD, xL, xR = neighborsX(i, j, N)
            yU, yD, yL, yR = neighborsY(i, j, N)
            
            a1 = const(i,j, M, N, M[x], M[y])
            f1x = M[xU] + M[xD] + M[xL] + M[xR] - M[x]*(4 + a1)
            f1y = M[yU] + M[yD] + M[yL] + M[yR] - M[y]*(4 + a1)
            k1x = f1x*dt
            k1y = f1y*dt
            
            a2 = const(i,j, M, N, M[x] + k1x/2, M[y] + k1y/2)
            f2x = M[xU] + M[xD] + M[xL] + M[xR] - (M[x] + k1x/2)*(4 + a2)
            f2y = M[yU] + M[yD] + M[yL] + M[yR] - (M[y] + k1y/2 )*(4 + a2)
            k2x = f2x*dt
            k2y = f2y*dt
            
            a3 = const(i,j, M, N, M[x] + k2x/2, M[y] + k2y/2)
            f3x = M[xU] + M[xD] + M[xL] + M[xR] - (M[x] + k2x/2)*(4 + a3)
            f3y = M[yU] + M[yD] + M[yL] + M[yR] - (M[y] + k2y/2 )*(4 + a3)
            k3x = f3x*dt
            k3y = f3y*dt
            
            a4 = const(i,j, M, N, M[x] + k3x, M[y] + k3y)
            f4x = M[xU] + M[xD] + M[xL] + M[xR] - (M[x] + k3x)*(4 + a4)
            f4y = M[yU] + M[yD] + M[yL] + M[yR] - (M[y] + k3y)*(4 + a4)
            k4x = f4x*dt
            k4y = f4y*dt
            
            
            Xnew = M[i,j,0] + k1x/6 + k2x/3 + k3x/3 + k4x/6
            Ynew = M[i,j,1] + k1y/6 + k2y/3 + k3y/3 + k4y/6
            
            norm = np.sqrt(Xnew*Xnew + Ynew*Ynew)
            
            MNew[i,j,0] = Xnew/norm
            MNew[i,j,1] = Ynew/norm
    M = MNew
    return M

    
    
def flow_rk4(M, N, dt, Nt):
    #files = []
    for k in range(Nt):
        M = gradient_flow_rk4(M, dt)
        #plot_vortex_dt(M,N,dt,k)
        #files.append('VA_lattice={}_beta={}_t={}.jpg'.format(N,beta,dt*k))
    return M#, files
    

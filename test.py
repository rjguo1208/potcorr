import read_eps
import matplotlib.pyplot as plt
import numpy as np
import math
import io_xsf
#--Read dieletric function files--#
#Eps0 =  read_eps.Epsmat('mos2_eps_uc/eps0mat.h5')
#Eps1 =  read_eps.Epsmat('mos2_eps_uc/epsmat.h5')
Eps1 =  read_eps.Epsmat('./mos2_chi/chimat.h5')
Eps0 =  read_eps.Epsmat('./mos2_chi/chi0mat.h5')
#----

# Crystal structure
A = 3.186*6/0.529
C = 3.186*8/0.529
#A = 2.51*12
#C = 25
#A = 5.43*4
#C = 5.43*4
# Define the grid for the potential
nx, ny, nz = 120, 120, 160 # number of grid points in each dimension
#N = nx*ny*nz
Rho = io_xsf.Rho('useful/rho_tot.xsf')
#Rho1 = read_eps.Rho('v.xsf')
#Rho1 = read_eps.Rho('Rho/WFC_d.dat_K001_B466WFC_d.xsf')
#Rho2 = read_eps.Rho('Rho/WFC_d.dat_K001_B467WFC_d.xsf')
nx = int(Rho.nx)
ny = int(Rho.ny)
nz = int(Rho.nz)
#nxyz = (nx-1)*(ny-1)*(nz-1)

dx, dy, dz = A/nx, A/ny, C/nz  # spacing between grid points in each dimension

x = np.linspace(0, (nx)*dx, nx) 
y = np.linspace(0, (ny)*dy, ny) 
z = np.linspace(0, (nz)*dz, nz)

xx, yy, zz = np.meshgrid(x, y, z, indexing='ij')

# Define the charge density
#sigma = 0.1
#rho = np.exp(-(((xx-nx*dx*0.5)-0.5*(yy-ny*dy*0.5))**2 + 3/4*(yy-ny*dy*0.5)**2 + (zz-nz*dz*(1-0.5))**2)/(2*sigma**2))
rho = Rho.rho
#rho = np.exp(-((xx-nx*dx*0.5)**2 + (yy-ny*dy*0.5)**2 + (zz-nz*dz*(0.5))**2)/(2*sigma**2))
#print(np.sum(rho))
#rho = rho/np.sum(rho)
#print(np.sum(rho))
#rho = Rho1.rho*Rho1.omega/nxyz + Rho2.rho*Rho2.omega/nxyz
#rho[18,18,72] = rho[18,18,72]-np.sum(rho)

rho_k0 = np.fft.fftn(rho)

#pot_k = np.fft.fftn(Rho1.rho)
kx = 2*np.pi*np.fft.fftfreq(nx, d=dx)
ky = 2*np.pi*np.fft.fftfreq(ny, d=dy)
kz = 2*np.pi*np.fft.fftfreq(nz, d=dz)
kxx, kyy, kzz = np.meshgrid(kx, ky, kz, indexing='ij')

k_squared = kxx**2 + (np.sqrt(3)/3*kxx+2*np.sqrt(3)/3*kyy)**2 + kzz**2
#k_squared = kxx**2 + kyy**2 + kzz**2

kxy = np.sqrt(kxx**2 + (np.sqrt(3)/3*kxx+2*np.sqrt(3)/3*kyy)**2)
#v_coul = 8*np.pi/k_squared*(1-np.exp(-0.5*kxy*C)*np.cos(0.5*kzz*C))
#v_coul[0,0,0] = 0
v_coul = 8*np.pi/(k_squared)
v_coul[0,0,0] = 0
#v_coul[1:-1, :,:] = 0
#v_coul[:,1:-1,:] = 0
#v_coul[:,:,1:-1] = 0

kx_tt = kx*3.186/(2*np.pi)/0.529
ky_tt = ky*3.186/(2*np.pi)/0.529
kz_tt = kz*3.186*8/(2*np.pi)/0.529
#kx_tt = kx*5.43/(2*np.pi)
#ky_tt = ky*5.43/(2*np.pi)
#kz_tt = kz*5.43/(2*np.pi)
#kx_tt = kx*2.51/(2*np.pi)
#ky_tt = ky*2.51/(2*np.pi)
#kz_tt = kz*25/(2*np.pi)

kxx_tt, kyy_tt, kzz_tt = np.meshgrid(kx_tt, ky_tt, kz_tt, indexing='ij')

rho_k = -np.ones((nx,ny,nz))*np.exp(1j*(kxx*nx*dx*0.5+kyy*ny*dy*0.5+kzz*nz*dz*0.5))*(nx)*(ny)*(nz)/Rho.omega
print(np.sum(np.fft.ifftn(rho_k)))
print(np.sum(np.fft.ifftn(rho_k0)))
#-Initialization-
g_tuple_tt = np.empty((nx,ny,nz), dtype = tuple)
g_tuple_floor_tt = np.empty((nx,ny,nz), dtype = tuple)
g_tuple_q_tt = np.empty((nx,ny,nz), dtype = tuple)
g_rhoind_tt = np.empty((nx,ny,nz), dtype = int)
q_ind_tt = np.empty((nx,ny,nz), dtype = int)
gg_eps = int(len(Eps1.G_vec2ind))
for i in range(nx):
    for j in range(ny):
        for k in range(nz):
            g_tuple_tt[i,j,k] = (float(format(kxx_tt[i,j,k],'.8f')), float(format(kyy_tt[i,j,k],'.8f')), float(format(kzz_tt[i,j,k],'.8f')))
            g_tuple_floor_tt[i,j,k] = (math.floor(g_tuple_tt[i,j,k][0]), math.floor(g_tuple_tt[i,j,k][1]), math.floor(g_tuple_tt[i,j,k][2]))
            g_tuple_q_tt[i,j,k] = (float(format(g_tuple_tt[i,j,k][0]-g_tuple_floor_tt[i,j,k][0],'.8f')), float(format(g_tuple_tt[i,j,k][1]-g_tuple_floor_tt[i,j,k][1],'.8f')),float(format(g_tuple_tt[i,j,k][2]-g_tuple_floor_tt[i,j,k][2],'.8f')))
            if g_tuple_floor_tt[i,j,k] in Eps1.G_vec2ind:
                g_rhoind_tt[i,j,k] = Eps1.G_vec2ind[g_tuple_floor_tt[i,j,k]]
            else:
                g_rhoind_tt[i,j,k] = gg_eps

            if g_tuple_q_tt[i,j,k] in Eps1.q_vec2ind:
                q_ind_tt[i,j,k] = Eps1.q_vec2ind[g_tuple_q_tt[i,j,k]]
            else:
                print("Error: No targeted q-points")
                print(g_tuple_tt[i,j,k])
                print(g_tuple_floor_tt[i,j,k])
                print(g_tuple_q_tt[i,j,k])
    print(i)
Eps1.gind_rho2eps = np.hstack((Eps1.gind_rho2eps,200000*np.ones([len(Eps1.qpts),1],dtype = int)))
Eps0.gind_rho2eps = np.hstack((Eps0.gind_rho2eps,200000*np.ones([len(Eps0.qpts),1],dtype = int)))
g_epsind_tt = Eps1.gind_rho2eps[q_ind_tt, g_rhoind_tt]
g_epsind0_tt = Eps0.gind_rho2eps[0, g_rhoind_tt]

lmax = len(g_epsind_tt[g_epsind_tt==200000])
g_epsind_tt[g_epsind_tt==200000] = 200000+np.arange(lmax)
g_epsind0_tt[g_epsind0_tt==200000] = 200000+np.arange(lmax)
#-End Initialization



G_ind_cut = 1300

eps_mat = np.zeros((len(Eps1.qpts),G_ind_cut, G_ind_cut), dtype=complex)
epsinv_mat = np.zeros((len(Eps1.qpts),G_ind_cut, G_ind_cut), dtype=complex)
for i in range(len(Eps1.qpts)):
    print("Constructing eps_mat of number %d/%d q-point" % (i+1, len(Eps1.qpts)))
    G_ind_list_tmp = g_epsind_tt[q_ind_tt == i]
    v_coul_list_tmp = v_coul[q_ind_tt == i]

    for j in range(G_ind_cut):
        v_coul_tmp = v_coul_list_tmp[G_ind_list_tmp == j+1]
        if (i == 0) :
            eps_mat[i,j,:] = -1.0*v_coul_tmp[0]*(Eps0.mat[0,0,0,j,:G_ind_cut,0]+Eps0.mat[0,0,0,j,:G_ind_cut,1]*1.0j)
        else:
            eps_mat[i,j,:] = -1.0*v_coul_tmp[0]*(Eps1.mat[i,0,0,j,:G_ind_cut,0]+Eps1.mat[i,0,0,j,:G_ind_cut,1]*1.0j)
        for jp in range(G_ind_cut):
            if (j == jp) :
                eps_mat[i,j,jp] = eps_mat[i,j,jp]+1.0
    
    print("Inversing eps_mat of number %d/%d q-point" % (i+1, len(Eps1.qpts)))
    epsinv_mat[i, :, :] = np.linalg.inv(eps_mat[i,:,:])
            
            







#-Calculation of screened/bare potential

phi_k = np.zeros((nx,ny,nz),dtype = complex)
phi_k[g_epsind_tt>1300] = 0#rho_k[g_epsind_tt>1300]*v_coul[g_epsind_tt>1300]

phik_list = []

g_epsind_reduced = g_epsind_tt[g_epsind_tt<=1300]
g_epsind0_reduced = g_epsind0_tt[g_epsind_tt<=1300]
q_ind_reduced = q_ind_tt[g_epsind_tt<=1300]
pot_reduced = v_coul[g_epsind_tt<=1300]*rho_k[g_epsind_tt<=1300]
v_coul_reduced = v_coul[g_epsind_tt<=1300]
for i in range(len(g_epsind_tt[g_epsind_tt<=1300])):
    G_index = g_epsind_reduced[i]
    q_index = q_ind_reduced[i]
    pot0 = pot_reduced[q_ind_reduced==q_index]
    eps_tmp=epsinv_mat[q_ind_reduced[q_ind_reduced==q_index],G_index-1,g_epsind_reduced[q_ind_reduced==q_index]-1]
#    phi_tmp =pot_reduced[i] - v_coul_reduced[i]*(eps_tmp@pot0.T)
    phi_tmp = eps_tmp@pot0.T
    phik_list.append(phi_tmp)
    print(i)

phi_k[g_epsind_tt<=1300] = phik_list

phi_k[0,0,0] = 0
#-End Calculation


#-FFT
phi_k_0 = v_coul*rho_k
phi_k_0[0,0,0] = 0
phi_k_1 = v_coul*rho_k0
phi_k_1[0,0,0] = 0
phi = np.fft.ifftn(phi_k)
phi_0 = np.fft.ifftn(phi_k_0)
phi_1 = np.fft.ifftn(phi_k_1)
phiR = np.real(phi)
phiR_0 = np.real(phi_0)
phiR_1 = np.real(phi_1)
#-End FFT

avg_z = []
avg_z_0 = []
avg_z_1 = []
avg_z_2 = []
for i in range(nz):
    phi_z = np.sum(phiR[:,:, i])
    phi_z_0 = np.sum(phiR_0[:,:,i]) 
    phi_z_1 = np.sum(phiR_1[:,:, i])
 #   phi_z_2 = np.sum(phiR_2[:,:, i])

    avg_z.append(phi_z)
    avg_z_0.append(phi_z_0)
    avg_z_1.append(phi_z_1)
  #  avg_z_2.append(phi_z_2)

fig,ax = plt.subplots(figsize = (8,6), dpi = 100)
ax.plot(range(nz), np.array(avg_z), label='screened potential')
ax.plot(range(nz), np.array(avg_z_0), label='bare potential')
ax.plot(range(nz), np.array(avg_z_1), label='DFT potential')
#ax.plot(range(nz), np.array(avg_z_2), label='bare potential')
ax.legend()
plt.show()

plot_xsf = True
plot_xsf = False
if plot_xsf:
    print("writing output file")
    with open('pot_point_charge_mos2.xsf','w') as v_file:
        print('Crystal\n PRIMVEC\n', end=' ',file=v_file)
        #print(' 3.186 0.00 0.00\n  -1.593  2.759 0.00\n 0.00 0.00 25.488\n', end=' ',file=v_file)
        #print(' 19.116 0.00 0.00\n  -9.558 16.555 0.00\n 0.00 0.00 25.488\n', end=' ',file=v_file)
        print("%16.10f %16.10f %16.10f\n%16.10f %16.10f %16.10f\n%16.10f %16.10f %16.10f\n" % (A, 0, 0, -0.5*A, np.sqrt(3)/2*A, 0, 0, 0, C),end='', file=v_file)
        #print("%16.10f %16.10f %16.10f\n%16.10f %16.10f %16.10f\n%16.10f %16.10f %16.10f\n" % (A, 0, 0, 0, A, 0, 0, 0, C),end='', file=v_file)
        #print(' 6.372 0.00 0.00\n  -3.186 5.518 0.00\n 0.00 0.00 25.488\n', end=' ',file=v_file)
        #print(' 10.86 0.00 0.00\n  0.00 10.86 0.00\n 0.00 0.00 10.86\n', end=' ',file=v_file)
        print(' PRIMCOORD\n 0    1\n', end=' ',file=v_file)
        print(' BEGIN_BLOCK_DATAGRID_3D\n 3D_PWSCF\n BEGIN_DATAGRID_3D_UNKNOWN\n', end=' ',file=v_file)
        print("%8d %8d %8d\n" % (nx, ny, nz), end=' ', file=v_file)
        print(' 0.00 0.00 0.00\n', end=' ',file=v_file)
        print("%16.10f %16.10f %16.10f\n%16.10f %16.10f %16.10f\n%16.10f %16.10f %16.10f\n" % (A, 0, 0, -0.5*A, np.sqrt(3)/2*A, 0, 0, 0, C),end='', file=v_file)
        #print("%16.10f %16.10f %16.10f\n%16.10f %16.10f %16.10f\n%16.10f %16.10f %16.10f\n" % (A, 0, 0, 0, A, 0, 0, 0, C),end='', file=v_file)
        #print(' 19.116 0.00 0.00\n  -9.558 16.555 0.00\n 0.00 0.00 25.488\n', end=' ',file=v_file)
        #print(' 6.372 0.00 0.00\n  -3.186 5.518 0.00\n 0.00 0.00 25.488\n', end=' ',file=v_file)
        #print(' 3.186 0.00 0.00\n  -1.593  2.759 0.00\n 0.00 0.00 25.488\n', end=' ',file=v_file)
        #print(' 10.86 0.00 0.00\n  0.00 10.86 0.00\n 0.00 0.00 10.86\n', end=' ',file=v_file)
        count = 0
        for k in range(nz):
            for j in range(ny):
                for i in range(nx):
                    print("%16.10f" %phiR[i,j,k], end = ' ', file=v_file)
                    count+=1
                    if count%5==0:
                        print(end='\n', file=v_file)
    
        print(' END_DATAGRID_3D\n', end=' ',file=v_file)
        print(' END_BLOCK_DATAGRID_3D', end=' ',file=v_file)
   
        v_file.close


import h5py
from scipy.interpolate import make_interp_spline
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import scipy
import math

#--
# Reading dielectric function hdf5 files produced by BerkeleyGW.
#--

def read_epsmat(name_):
    feps_ = h5py.File(name_, 'r')
    mat_diagonal_ = feps_['mats/matrix-diagonal']  # diagonal term of dielectric matrix
    qpts_ = feps_['eps_header/qpoints/qpts'] # qpoints coordinates
    nmtx_ = feps_['/eps_header/gspace/nmtx'] # Number of matrix elements BekerleyGW actually compute for each q-point.
    G_vec_ = feps_['/mf_header/gspace/components'] # G-vectors in RHO G-space
    gind_eps2rho_ = feps_['/eps_header/gspace/gind_eps2rho'] # convert Epsilon G-space to Rho G-space
    gind_rho2eps_ = feps_['/eps_header/gspace/gind_rho2eps'] # convert Rho G-space to Epsilon G-space
    return feps_, mat_diagonal_, qpts_, nmtx_, G_vec_, gind_eps2rho_, gind_rho2eps_

def G_converter(G_vec_): # Construct the converter between G-vector and G-index
    G_vec_tuple_list_ = []
    for i_ in G_vec_:
        i_ = tuple(i_)
        G_vec_tuple_list_.append(i_)
    G_ind2vec_ = dict(enumerate(G_vec_tuple_list_))
    G_vec2ind_ = {v:k for k,v in G_ind2vec_.items()}
    return G_ind2vec_, G_vec2ind_

def eps_model(Gz_): # interpolation of epsilon
    g0_ = [0,-1,Gz_]
    q0_ = qpts0[:10]
    q1_ = qpts1[:11]
    q0_abs_ = np.sqrt((q0_[:10,0])**2+ (np.sqrt(3)/3*q0_[:10,0]+2*np.sqrt(3)/3*q0_[:10,1])**2)
    q1_abs_ = np.sqrt((q1_[:11,0])**2+ (np.sqrt(3)/3*q1_[:11,0]+2*np.sqrt(3)/3*q1_[:11,1])**2)
    eps_qy_ = []
    for __ in range(3):
        g0_[1] = g0_[1]+1
        if __ == 0:
            for i in range(10):
                qind_ = i
                g_vec_ = tuple(g0_)
                gind_rho_ = G_vec2ind[g_vec_]
                gind_eps0_ = gind_rho2eps0[qind_, gind_rho_]
                mat0_ = mat_diagonal0[qind_, gind_eps0_-1, 0]
                eps_qy_.append(mat0_)
            q_abs_ = np.hstack((np.array([]),q0_abs_))

        for i in range(11):
            qind_ = i
            g_vec_ = tuple(g0_)
            gind_rho_ = G_vec2ind[g_vec_]
            gind_eps1_ = gind_rho2eps1[qind_, gind_rho_]
            mat1_ = mat_diagonal1[qind_, gind_eps1_-1,0]
            eps_qy_.append(mat1_)
        q_abs_ = np.hstack((q_abs_, q1_abs_+__*np.sqrt(3)*2/3))


    model_ = make_interp_spline(q_abs_,eps_qy_)
    return model_


#--Read dieletric function files--#
feps0, mat_diagonal0, qpts0, nmtx0, _, gind_eps2rho0, gind_rho2eps0 = read_epsmat('mos2_eps_gw/eps0mat.h5')
feps1, mat_diagonal1, qpts1, nmtx1, G_vec, gind_eps2rho1, gind_rho2eps1 = read_epsmat('mos2_eps_gw/epsmat.h5')
#----

#--Use G_converter function
G_ind2vec, G_vec2ind = G_converter(G_vec)
#----
'''
#--visulization--#
fig, ax = plt.subplots(figsize = (18,6), dpi = 100)
for i in range(20):
    print(i)
    model = eps_model(i)
    q_plot = np.linspace(0,5,10000)
    eps_plot = model(q_plot)

    ax.plot(q_plot, eps_plot, label = 'Gz=%d' % i)
#ax.scatter(q_abs, eps_qy)
ax.legend()
plt.show()
#----

epsz=[]
fig, ax = plt.subplots(figsize = (8,6), dpi = 100)
for i in range(-20, 21):
    model = eps_model(i)
    qz_0 = model([0])
    epsz.append(qz_0)
ax.scatter(range(-20, 21), epsz)
ax.plot(range(-20,21), epsz)
plt.show()
'''
#--FFT--
print("Starting FFT")
# Crystal structure
A = 3.186*6
C = 3.186*8

# Define the parameters of the Gaussian charge density
sigma = 0.5 # width of the Gaussian
center = [0, 0, 0] # center of the Gaussian

# Define the grid for the potential
nx, ny, nz = 30, 30, 120 # number of grid points in each dimension
N = nx*ny*nz
dx, dy, dz = A/nx, A/ny, C/nz  # spacing between grid points in each dimension
x = np.linspace(0, (nx)*dx, nx) - center[0]
y = np.linspace(0, (ny)*dy, ny) - center[1]
z = np.linspace(0, (nz)*dz, nz) - center[2]
#xx, yy, zz = np.meshgrid(x-0.5*y, np.sqrt(3)/2*y, z, indexing='ij')
xx, yy, zz = np.meshgrid(x, y, z, indexing='ij')

# Define the charge density
rho = np.exp(-((xx-nx*dx*0.5)**2 + (yy-ny*dy*0.5)**2 + (zz-nz*dz*0.5)**2)/(2*sigma**2))
v_r = np.zeros(xx.shape)
for i in range(nx):
    for j in range(ny):
        for k in range(nz):
            if (np.abs(xx[i,j,k]-0.5*nx*dx)<0.25*nx*dx and np.abs(yy[i,j,k]-0.5*ny*dy)<0.25*ny*dy and np.abs(zz[i,j,k]-0.5*nz*dz)<0.25*nz*dz):
                v_r[i,j,k] = 1
# Take the FFT of the charge density to get the potential in reciprocal space
rho_k = np.fft.fftn(rho)
v_k = np.fft.fftn(v_r)
# Multiply by the appropriate factors to get the long-range part of the potential
kx, ky, kz = 2*np.pi*np.fft.fftfreq(nx, d=dx), 2*np.pi*np.fft.fftfreq(ny, d=dy), 2*np.pi*np.fft.fftfreq(nz, d=dz)

kx_tt = kx*3.186*6/(2*np.pi)
ky_tt = ky*3.186*6/(2*np.pi)
kz_tt = kz*3.186*8/(2*np.pi)


print("Calculating eps_inv")
kxx_tt, kyy_tt, kzz_tt = np.meshgrid(kx_tt, ky_tt, kz_tt, indexing='ij')
eps_inv = -np.ones(kxx_tt.shape)
nx_, ny_, nz_ = eps_inv.shape
for i in range(nx_):
    for j in range(ny_):
        for k in range(nz_):
            Gz_tmp = int(np.round(kzz_tt[i,j,k]))
            if np.abs(Gz_tmp) < 22:
                model_tmp = eps_model(Gz_tmp)
                qxy_tmp = np.sqrt(kxx_tt[i,j,k]**2 + (np.sqrt(3)/3*kxx_tt[i,j,k]+2*np.sqrt(3)/3*kyy_tt[i,j,k])**2)
                mat_tmp = model_tmp(qxy_tmp)
                if mat_tmp > 1:
                    mat_tmp = 1.0
                eps_inv[i,j,k] = mat_tmp
            else:
                eps_inv[i,j,k] = 1.0

print("Starting iFFT")
kxx, kyy, kzz = np.meshgrid(kx, ky, kz, indexing='ij')

#kxx, kyy, kzz = np.meshgrid(kx, ky, kz, indexing='ij')
kxy = np.sqrt(kxx**2 + (np.sqrt(3)/3*kxx+2*np.sqrt(3)/3*kyy)**2)
k_squared = kxx**2 + (np.sqrt(3)/3*kxx+2*np.sqrt(3)/3*kyy)**2 + kzz**2
#k_squared = kxx**2 + kyy**2 + kzz**2
phi_k = 4*np.pi/(k_squared+0.01)*rho_k*(1-np.exp(-0.5*kxy*C)*np.cos(0.5*kzz*C))*eps_inv
phi_k_0 = 4*np.pi/(k_squared+0.01)*rho_k*(1-np.exp(-0.5*kxy*C)*np.cos(0.5*kzz*C))
phi_k_1 = 4*np.pi/(k_squared+0.01)*rho_k*eps_inv
phi_k_2 = 4*np.pi/(k_squared+0.01)*rho_k
#phi_k_1 = v_k*eps_inv
#phi_k_2 = v_k
#phi_k = 4*np.pi*rho_k/k_squared
#phi_k = 4*np.pi/k_squared*np.exp(1j*(kxx*nx*dx*0.5+kyy*ny*dy*0.5+kzz*nz*dz*0.5))*eps_inv
#phi_k = eps_inv
#phi_k[0, 0, 0] = 0
#phi_k_0[0, 0, 0] = 0
#phi_k_1[0, 0, 0] = 0
#phi_k_2[0, 0, 0] = 0
# Take the inverse FFT to get the potential in real space
phi = np.fft.ifftn(phi_k)
phi_0 = np.fft.ifftn(phi_k_0)
phi_1 = np.fft.ifftn(phi_k_1)
phi_2 = np.fft.ifftn(phi_k_2)

eps_real  = np.fft.ifftn(eps_inv)

phiR = np.real(phi)
phiR_0 = np.real(phi_0)
phiR_1 = np.real(phi_1)
phiR_2 = np.real(phi_2)

eps_realR = np.real(eps_real)

print("writing output file")
with open('v.xsf','w') as v_file:
    print('Crystal\n PRIMVEC\n', end=' ',file=v_file)
   # print(' 19.116 0.00 0.00\n  -9.558 16.555 0.00\n 0.00 0.00 25.488\n', end=' ',file=v_file)
    print(' 19.116 0.00 0.00\n  -9.558 16.555 0.00\n 0.00 0.00 25.488\n', end=' ',file=v_file)
    print(' PRIMCOORD\n 0    1\n', end=' ',file=v_file)
    print(' BEGIN_BLOCK_DATAGRID_3D\n 3D_PWSCF\n BEGIN_DATAGRID_3D_UNKNOWN\n', end=' ',file=v_file)
    print("%8d %8d %8d\n" % (nx, ny, nz), end=' ', file=v_file)
    print(' 0.00 0.00 0.00\n', end=' ',file=v_file)
    print(' 19.116 0.00 0.00\n  -9.558 16.555 0.00\n 0.00 0.00 25.488\n', end=' ',file=v_file)
#    print(' 10.00 0.00 0.00\n  0.00 10.00 0.00\n 0.00 0.00 20.00\n', end=' ',file=v_file)
    count = 0
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                print("%16.10f" %eps_realR[i,j,k], end = ' ', file=v_file)
                count+=1
                if count%5==0:
                    print(end='\n', file=v_file)
   
    v_file.close

avg_z = []
avg_z_0 = []
avg_z_1 = []
avg_z_2 = []
for i in range(nz):
    phi_z = np.sum(phiR[:,:, i])
    phi_z_0 = np.sum(phiR_0[:,:,i]) 
    phi_z_1 = np.sum(phiR_1[:,:, i])
    phi_z_2 = np.sum(phiR_2[:,:, i])

    avg_z.append(phi_z)
    avg_z_0.append(phi_z_0)
    avg_z_1.append(phi_z_1)
    avg_z_2.append(phi_z_2)

fig,ax = plt.subplots(figsize = (8,6), dpi = 100)
ax.plot(range(nz), np.array(avg_z), label='screened potential (2D)')
ax.plot(range(nz), np.array(avg_z_0), label='bare potential (2D)')
ax.plot(range(nz), np.array(avg_z_1), label='screened potential')
ax.plot(range(nz), np.array(avg_z_2), label='bare potential')
ax.legend()
plt.show()

fig, ax = plt.subplots(figsize = (8,6), dpi = 100)
ax.plot(range(nz), (np.array(avg_z_0))/(np.array(avg_z)))
ax.plot(range(nz), (np.array(avg_z_2))/(np.array(avg_z_1)))
plt.show()


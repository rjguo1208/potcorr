import read_eps
import io_xsf
import lab
import const
import numpy as np
import math
import psutil
import os


class PotCorr():
    def __init__(self, cell):
        self.prefix = cell['prefix']
        self.folder = cell['folder']
        self.Chi0_filename = cell['folder_chi']+'chi0mat.h5'
        self.Chi1_filename = cell['folder_chi']+'chimat.h5'
        #self.Chi0_doped_filename = cell['folder_chi_doped']+'chi0mat.h5'
        #self.Chi1_doped_filename = cell['folder_chi_doped']+'chi0mat.h5'
        self.is2d = cell['is2d']
        self.ibrav = cell['ibrav']
        self.lattpara_unit = np.array(cell['lattpara_unit'])
        self.sc = np.array(cell['sc'])
        self.fft_g = np.array(cell['fft_g'])
        self.lattpara = self.lattpara_unit*self.sc/const.Bohr_R
        if self.ibrav == 4:
            self.omega = self.lattpara_unit[0]*self.sc[0] * \
                         self.lattpara_unit[1]*self.sc[1] * \
                         self.lattpara_unit[2]*self.sc[2] * \
                         np.sqrt(3)/2 / (const.Bohr_R)**3
    

    def fft_init(self):
        print('Begin FFT grids Initialization')

        self.fft_nx=self.fft_g[0]
        self.fft_ny=self.fft_g[1]
        self.fft_nz=self.fft_g[2]

        self.fft_dx=self.lattpara[0]/self.fft_g[0]
        self.fft_dy=self.lattpara[1]/self.fft_g[1]
        self.fft_dz=self.lattpara[2]/self.fft_g[2]

        self.fft_x=np.linspace(0,self.fft_nx*self.fft_dx, self.fft_nx)
        self.fft_y=np.linspace(0,self.fft_ny*self.fft_dy, self.fft_ny)
        self.fft_z=np.linspace(0,self.fft_nz*self.fft_dz, self.fft_nz)

        self.fft_xx, self.fft_yy, self.fft_zz = np.meshgrid(self.fft_x, self.fft_y, self.fft_z, indexing='ij')

        self.fft_kx = 2*np.pi*np.fft.fftfreq(self.fft_nx, d=self.fft_dx)
        self.fft_ky = 2*np.pi*np.fft.fftfreq(self.fft_ny, d=self.fft_dy)
        self.fft_kz = 2*np.pi*np.fft.fftfreq(self.fft_nz, d=self.fft_dz)

        self.fft_kx_tt = self.fft_kx*self.lattpara_unit[0]/(2*np.pi)/const.Bohr_R
        self.fft_ky_tt = self.fft_ky*self.lattpara_unit[1]/(2*np.pi)/const.Bohr_R
        self.fft_kz_tt = self.fft_kz*self.lattpara_unit[2]/(2*np.pi)/const.Bohr_R
        
        self.fft_kxx, self.fft_kyy, self.fft_kzz = np.meshgrid(self.fft_kx, self.fft_ky, self.fft_kz, indexing='ij')
        self.fft_kxx_tt, self.fft_kyy_tt, self.fft_kzz_tt = np.meshgrid(self.fft_kx_tt, self.fft_ky_tt, self.fft_kz_tt, indexing='ij')
        print('FFT grids:%6d %6d %6d' % (self.fft_nx, self.fft_ny, self.fft_nz))


    def get_vcoul(self):
        print('Constructing Coulomb kernel...')
        if self.ibrav == 1:
            print("ibrav=", self.ibrav)
            self.k_squared = self.fft_kxx**2 + self.fft_kyy**2 + self.fft_kzz**2
            self.v_coul = 8*np.pi/(self.k_squared)
            self.v_coul[0,0,0] = 0
            if self.is2d:
                self.kxy = np.sqrt(self.fft_kxx**2 + self.fft_kyy**2)
                self.v_coul2d = 8*np.pi/(self.k_squared)*(1-np.exp(-0.5*self.kxy*self.lattpara[2])*np.cos(0.5*self.fft_kzz*self.lattpara[2]))
                self.v_coul2d[0,0,0] = 0
        elif self.ibrav == 4:
            print("ibrav=", self.ibrav)
            self.k_squared = self.fft_kxx**2 + (np.sqrt(3)/3*self.fft_kxx+2*np.sqrt(3)/3*self.fft_kyy)**2 + self.fft_kzz**2
            self.v_coul = 8*np.pi/(self.k_squared)
            self.v_coul[0,0,0] = 0
            if self.is2d:
                self.kxy = np.sqrt(self.fft_kxx**2 + (np.sqrt(3)/3*self.fft_kxx+2*np.sqrt(3)/3*self.fft_kyy)**2)
                self.v_coul2d = 8*np.pi/(self.k_squared)*(1-np.exp(-0.5*self.kxy*self.lattpara[2])*np.cos(0.5*self.fft_kzz*self.lattpara[2]))
                self.v_coul2d[0,0,0] = 0
        else:
            print('Lattice with this ibrav has not been implemented..')
        
    def get_greenfunc(self):
        print('Constructing Green Function...')
        x = np.linspace(-self.lattpara[0]/2, self.lattpara[0]/2, self.fft_nx)
        y = np.linspace(-self.lattpara[1]/2, self.lattpara[1]/2, self.fft_ny)
        z = np.linspace(-self.lattpara[2]/2, self.lattpara[2]/2, self.fft_nz)
        X, Y, Z = np.meshgrid(x, y, z)
        if self.ibrav == 1:
            print("ibrav=", self.ibrav)
            R = np.sqrt(X**2 + Y**2 + Z**2)

        elif self.ibrav == 4:
            print("ibrav=", self.ibrav)
            R = np.sqrt((X-0.5*Y)**2 + (np.sqrt(3)/2*Y)**2 + Z**2)
        epsilon = 1e-10  # to avoid divide by zero
        G = 1.0 / (4.0 * np.pi * (R + epsilon))
        V_r = np.roll(G, shift=-np.array(G.shape)//2, axis=(0, 1, 2))
        V_k = np.fft.fftn(V_r)
        self.greenfunc = G
        self.v_coulgfr = V_r
        self.v_coulgfk = V_k
    
    def read_chi(self):
        print('Reading Chi matrix hdf5 files')
        self.Chi0 = read_eps.Epsmat(self.Chi0_filename)
        self.Chi1 = read_eps.Epsmat(self.Chi1_filename)

    def read_chi_doped(self):
        print('Read doped Chi matrix hdf5 files')
        self.Chi0 = read_eps.Epsmat(self.Chi0_doped_filename)
        self.Chi1 = read_eps.Epsmat(self.Chi1_doped_filename)

    def read_rho_bare(self, rho_bare_file):
        Rho = io_xsf.Rho(rho_bare_file)
        self.rho_bare_r = Rho.rho#*Rho.omega/Rho.nx/Rho.ny/Rho.nz
        self.rho_bare_k = np.fft.fftn(self.rho_bare_r)

    def read_rho_tot(self, rho_tot_file):
        Rho = io_xsf.Rho(rho_tot_file)
        self.rho_tot_r = Rho.rho#*Rho.omega/Rho.nx/Rho.ny/Rho.nz
        self.rho_tot_k = np.fft.fftn(self.rho_tot_r)  

    def read_pot_bare(self, pot_bare_file):
        Pot = io_xsf.Pot(pot_bare_file)
        self.pot_bare_r = Pot.pot
        self.pot_bare_k = np.fft.ifftn(self.pot_bare_r)

    def read_pot_tot(self, pot_tot_file):
        Pot = io_xsf.Pot(pot_tot_file)
        self.pot_tot_r = Pot.pot
        self.pot_tot_k = np.fft.ifftn(self.pot_tot_r)



    def get_pointc_pot_bare(self, position='center', kernel='3d'):
        print('Setting point-charge potential as bare potential')
        if position == 'center':
            phase_p=np.exp(1j*(self.fft_kxx*self.fft_nx*self.fft_dx*0.5+
                               self.fft_kyy*self.fft_ny*self.fft_dy*0.5+
                               self.fft_kzz*self.fft_nz*self.fft_dz*0.5))
        else:
            phase_p = np.exp(1j*(self.fft_kxx*self.fft_nx*self.fft_dx*position[0]+
                                 self.fft_kyy*self.fft_ny*self.fft_dy*position[1]+
                                 self.fft_kzz*self.fft_nz*self.fft_dz*position[2]))
        
        if kernel == '3d':
            v_coul_ = self.v_coul
        elif kernel == '2d':
            v_coul_ = self.v_coul2d  
        elif kernel == 'gf':
            v_coul_ = self.v_coulgfk  
        
        self.pot_bare_k =  v_coul_*phase_p
        self.pot_bare_r = np.fft.ifftn(self.pot_bare_k)

    def epsmat_init(self):
        Eps1=self.Chi1
        Eps0=self.Chi0
        nx_, ny_, nz_ = self.fft_nx, self.fft_ny, self.fft_nz
        g_tuple_tt = np.empty((nx_,ny_,nz_), dtype = tuple)
        g_tuple_floor_tt = np.empty((nx_,ny_,nz_), dtype = tuple)
        g_tuple_q_tt = np.empty((nx_,ny_,nz_), dtype = tuple)
        g_rhoind_tt = np.empty((nx_,ny_,nz_), dtype = int)
        q_ind_tt = np.empty((nx_,ny_,nz_), dtype = int)

        gg_eps = int(len(Eps1.G_vec2ind))
        for i in range(nx_):
            for j in range(ny_):
                for k in range(nz_):
                    g_tuple_tt[i,j,k] = (float(format(self.fft_kxx_tt[i,j,k],'.8f')), 
                                         float(format(self.fft_kyy_tt[i,j,k],'.8f')), 
                                         float(format(self.fft_kzz_tt[i,j,k],'.8f')))
                    g_tuple_floor_tt[i,j,k] = (math.floor(g_tuple_tt[i,j,k][0]), 
                                               math.floor(g_tuple_tt[i,j,k][1]), 
                                               math.floor(g_tuple_tt[i,j,k][2]))
                    g_tuple_q_tt[i,j,k] = (float(format(g_tuple_tt[i,j,k][0]-g_tuple_floor_tt[i,j,k][0],'.8f')),
                                           float(format(g_tuple_tt[i,j,k][1]-g_tuple_floor_tt[i,j,k][1],'.8f')),
                                           float(format(g_tuple_tt[i,j,k][2]-g_tuple_floor_tt[i,j,k][2],'.8f')))
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
                        return
            print(i)

        Eps1.gind_rho2eps = np.hstack((Eps1.gind_rho2eps,200000*np.ones([len(Eps1.qpts),1],dtype = int)))
        Eps0.gind_rho2eps = np.hstack((Eps0.gind_rho2eps,200000*np.ones([len(Eps0.qpts),1],dtype = int)))
        g_epsind_tt = Eps1.gind_rho2eps[q_ind_tt, g_rhoind_tt]
        g_epsind0_tt = Eps0.gind_rho2eps[0, g_rhoind_tt]

        lmax = len(g_epsind_tt[g_epsind_tt==200000])
        g_epsind_tt[g_epsind_tt==200000] = 200000+np.arange(lmax)
        g_epsind0_tt[g_epsind0_tt==200000] = 200000+np.arange(lmax)

        self.q_ind_tt = q_ind_tt
        self.lmax = lmax
        self.g_epsind_tt = g_epsind_tt
        self.g_epsind0_tt = g_epsind0_tt 


    def epsmat_init_interp(self):
        Eps1=self.Chi1
        Eps0=self.Chi0
        nx_, ny_, nz_ = self.fft_nx, self.fft_ny, self.fft_nz
        g_tuple_tt = np.empty((nx_,ny_,nz_), dtype = tuple)
        g_tuple_floor_tt = np.empty((nx_,ny_,nz_), dtype = tuple)
        g_tuple_q_tt = np.empty((nx_,ny_,nz_), dtype = tuple)
        g_rhoind_tt = np.empty((nx_,ny_,nz_), dtype = int)
        q_ind_tt = np.empty((nx_,ny_,nz_), dtype = int)
        gg_eps = int(len(Eps1.G_vec2ind))

        for i in range(nx_):
            for j in range(ny_):
                for k in range(nz_):
                    g_tuple_tt[i,j,k] = (float(format(self.fft_kxx_tt[i,j,k],'.8f')), 
                                         float(format(self.fft_kyy_tt[i,j,k],'.8f')), 
                                         float(format(self.fft_kzz_tt[i,j,k],'.8f')))
                    g_tuple_floor_tt[i,j,k] = (math.floor(g_tuple_tt[i,j,k][0]), 
                                               math.floor(g_tuple_tt[i,j,k][1]), 
                                               math.floor(g_tuple_tt[i,j,k][2]))
                    g_tuple_q_tt[i,j,k] = (float(format(g_tuple_tt[i,j,k][0]-g_tuple_floor_tt[i,j,k][0],'.8f')),
                                           float(format(g_tuple_tt[i,j,k][1]-g_tuple_floor_tt[i,j,k][1],'.8f')),
                                           float(format(g_tuple_tt[i,j,k][2]-g_tuple_floor_tt[i,j,k][2],'.8f')))
                    if g_tuple_floor_tt[i,j,k] in Eps1.G_vec2ind:
                        g_rhoind_tt[i,j,k] = Eps1.G_vec2ind[g_tuple_floor_tt[i,j,k]]
                    else:
                        g_rhoind_tt[i,j,k] = gg_eps

        
            print(i)
        self.g_rhoind_tt = g_rhoind_tt
        self.g_tuple_q_tt = g_tuple_q_tt




    def get_chimat(self, G_ind_cut = 1300):
        Eps1=self.Chi1
        chi_mat = np.zeros((len(Eps1.qpts),G_ind_cut, G_ind_cut), dtype=complex)
        #epsinv_mat = np.zeros((len(Eps1.qpts),G_ind_cut, G_ind_cut), dtype=complex)
        for  i in range(len(Eps1.qpts)):
            print("Constructing chi_mat of number %d/%d q-point" % (i+1, len(Eps1.qpts)))
            #G_ind_list_tmp = self.g_epsind_tt[self.q_ind_tt == i]
            #print(G_ind_list_tmp)

            for j in range(G_ind_cut):
               
                #print(v_coul_tmp)
                if i == 0:
                    chi_mat[i,j,:] = self.Chi0.mat[0,0,0,j,:G_ind_cut,0]+self.Chi0.mat[0,0,0,j,:G_ind_cut,1]*1.0j
                else:
                    chi_mat[i,j,:] = self.Chi1.mat[i,0,0,j,:G_ind_cut,0]+self.Chi1.mat[i,0,0,j,:G_ind_cut,1]*1.0j
               # for jp in range(G_ind_cut):
               #     if j == jp:
               #         chi_mat[i,j,jp] = chi_mat[i,j,jp]+1.0

        self.chi_mat = chi_mat

    

    def get_epsmat(self, G_ind_cut = 1300, kernel='3d'):
        Eps1=self.Chi1
        eps_mat = np.zeros((len(Eps1.qpts),G_ind_cut, G_ind_cut), dtype=complex)
        #epsinv_mat = np.zeros((len(Eps1.qpts),G_ind_cut, G_ind_cut), dtype=complex)
        for  i in range(len(Eps1.qpts)):
            print("Constructing eps_mat of number %d/%d q-point" % (i+1, len(Eps1.qpts)))
            G_ind_list_tmp = self.g_epsind_tt[self.q_ind_tt == i]
            #print(G_ind_list_tmp)
            if kernel == '3d':
                v_coul_list_tmp = self.v_coul[self.q_ind_tt == i]
                #print(v_coul_list_tmp)
            elif kernel == '2d':
                v_coul_list_tmp = self.v_coul2d[self.q_ind_tt == i]
            elif kernel == 'gf':
                v_coul_list_tmp = self.v_coulgfk[self.q_ind_tt == i]
                
            else:
                print('Undefined kernel type!')
                return

            for j in range(G_ind_cut):
                v_coul_tmp = v_coul_list_tmp[G_ind_list_tmp == j+1]
                #print(v_coul_tmp)
                if i == 0:
                    eps_mat[i,j,:] = -1.0*v_coul_tmp[0]*(self.Chi0.mat[0,0,0,j,:G_ind_cut,0]+
                                                         self.Chi0.mat[0,0,0,j,:G_ind_cut,1]*1.0j)
                else:
                    eps_mat[i,j,:] = -1.0*v_coul_tmp[0]*(self.Chi1.mat[i,0,0,j,:G_ind_cut,0]+
                                                         self.Chi1.mat[i,0,0,j,:G_ind_cut,1]*1.0j)
                for jp in range(G_ind_cut):
                    if j == jp:
                        eps_mat[i,j,jp] = eps_mat[i,j,jp]+1.0

        self.eps_mat = eps_mat

    

    def epsmat_inv(self, G_ind_cut = 1300):
        Eps1 = self.Chi1
        epsinv_mat = np.zeros((len(Eps1.qpts), G_ind_cut, G_ind_cut), dtype=complex)
        for i in range(len(Eps1.qpts)):
            print("Inversing eps_mat of number %d/%d q-point" % (i+1, len(Eps1.qpts)))
            epsinv_mat[i,:,:] = np.linalg.inv(self.eps_mat[i,:,:])
        self.epsinv_mat = epsinv_mat


    def rho2pot_bare(self, kernel='3d'):
        if kernel == '3d':
            vcoul_tmp_ = self.v_coul
        elif kernel == '2d':
            vcoul_tmp_ = self.v_coul2d
        elif kernel == 'gf':
            vcoul_tmp_ = self.v_coulgfk
        else:
            print('Undefined kernel type!')
            return
        self.pot_bare_k = self.rho_bare_k*vcoul_tmp_
        self.pot_bare_r = np.fft.ifftn(self.pot_bare_k)
    
    def rho2pot_tot(self, kernel='3d'):
        if kernel == '3d':
            vcoul_tmp_ = self.v_coul
        elif kernel == '2d':
            vcoul_tmp_ = self.v_coul2d
        elif kernel == 'gf':
            vcoul_tmp_ = self.v_coulgfk
        else:
            print('Undefined kernel type!')
            return
        self.pot_tot_k = self.rho_tot_k*vcoul_tmp_
        self.pot_tot_r = np.fft.ifftn(self.pot_tot_k)

    def pot2rho_bare(self, kernel='3d', ncharge=0):
        if kernel == '3d':
            vcoul_tmp_ = self.v_coul
        elif kernel == '2d':
            vcoul_tmp_ = self.v_coul2d
        elif kernel == 'gf':
            vcoul_tmp_ = self.v_coulgfk
        else:
            print('Undefined kernel type!')
            return
        self.rho_bare_k = self.pot_bare_k/vcoul_tmp_
        self.rho_bare_k[0,0,0] = ncharge/self.omega*self.fft_nx*self.fft_ny*self.fft_nz
        self.rho_bare_r = np.fft.ifftn(self.rho_bare_k)


    def pot2rho_tot(self, kernel='3d', ncharge=0):
        if kernel == '3d':
            vcoul_tmp_ = self.v_coul
        elif kernel == '2d':
            vcoul_tmp_ = self.v_coul2d
        elif kernel == 'gf':
            vcoul_tmp_ = self.v_coulgfk
        else:
            print('Undefined kernel type!')
            return
        self.rho_tot_k = self.pot_tot_k/vcoul_tmp_
        self.rho_tot_k[0,0,0] = ncharge/self.omega*self.fft_nx*self.fft_ny*self.fft_nz
        self.rho_tot_r = np.fft.ifftn(self.rho_tot_k)

    def pot_bare2tot(self, kernel='3d', G_ind_cut=1300):
        if kernel == '3d':
            vcoul_tmp_ = self.v_coul
        elif kernel == '2d':
            vcoul_tmp_ = self.v_coul2d
        elif kernel == 'gf':
            vcoul_tmp_ = self.v_coulgfk
        else:
            print('Undefined kernel type!')
            return
        
        phi_k = np.zeros((self.fft_nx, self.fft_ny, self.fft_nz), dtype = complex)
        phi_k[self.g_epsind_tt > G_ind_cut] = self.pot_bare_k[self.g_epsind_tt > G_ind_cut]
        phik_list = []

        g_epsind_reduced = self.g_epsind_tt[self.g_epsind_tt<=G_ind_cut]
        g_epsind0_reduced = self.g_epsind_tt[self.g_epsind_tt<=G_ind_cut]
        q_ind_reduced = self.q_ind_tt[self.g_epsind_tt<=G_ind_cut]
        pot_reduced = self.pot_bare_k[self.g_epsind_tt<=G_ind_cut]
        v_coul_reduced = vcoul_tmp_[self.g_epsind_tt<=G_ind_cut]
        for i in range(len(self.g_epsind_tt[self.g_epsind_tt<=G_ind_cut])):
            G_index = g_epsind_reduced[i]
            q_index = q_ind_reduced[i]
            pot0 = pot_reduced[q_ind_reduced==q_index]
            eps_tmp=self.epsinv_mat[q_ind_reduced[q_ind_reduced==q_index],
                               G_index-1,
                               g_epsind_reduced[q_ind_reduced==q_index]-1]
            phi_tmp = eps_tmp@pot0.T
            phik_list.append(phi_tmp)
            print(i)
        phi_k[self.g_epsind_tt<=G_ind_cut] = phik_list
        #phi_k[0,0,0] = 0
        if kernel == '3d' or kernel == '2d':
            phi_k[0,0,0] = 0
        self.pot_tot_k = phi_k
        self.pot_tot_r = np.fft.ifftn(self.pot_tot_k)



    def pot_tot2bare(self, kernel='3d', G_ind_cut=1300):
        if kernel == '3d':
            vcoul_tmp_ = self.v_coul
        elif kernel == '2d':
            vcoul_tmp_ = self.v_coul2d
        elif kernel == 'gf':
            vcoul_tmp_ = self.v_coulgfk
        else:
            print('Undefined kernel type!')
            return
        
        phi_k = np.zeros((self.fft_nx, self.fft_ny, self.fft_nz), dtype = complex)
        phi_k[self.g_epsind_tt > G_ind_cut] = self.pot_tot_k[self.g_epsind_tt > G_ind_cut]
        phik_list = []

        g_epsind_reduced = self.g_epsind_tt[self.g_epsind_tt<=G_ind_cut]
        g_epsind0_reduced = self.g_epsind_tt[self.g_epsind_tt<=G_ind_cut]
        q_ind_reduced = self.q_ind_tt[self.g_epsind_tt<=G_ind_cut]
        pot_reduced = self.pot_tot_k[self.g_epsind_tt<=G_ind_cut]
        v_coul_reduced = vcoul_tmp_[self.g_epsind_tt<=G_ind_cut]
        for i in range(len(self.g_epsind_tt[self.g_epsind_tt<=G_ind_cut])):
            G_index = g_epsind_reduced[i]
            q_index = q_ind_reduced[i]
            pot0 = pot_reduced[q_ind_reduced==q_index]
            eps_tmp=self.eps_mat[q_ind_reduced[q_ind_reduced==q_index],
                               G_index-1,
                               g_epsind_reduced[q_ind_reduced==q_index]-1]
            phi_tmp = eps_tmp@pot0.T
            phik_list.append(phi_tmp)
            print(i)
        phi_k[self.g_epsind_tt<=G_ind_cut] = phik_list
        #phi_k[0,0,0] = 0
        self.pot_bare_k = phi_k
        self.pot_bare_r = np.fft.ifftn(self.pot_bare_k)

    def write_xsf(self, ftype=None,filedir='./outfile.xsf'):


        if ftype == None:
            print("Please define ftype")
            return
        elif ftype == 'pot_bare':
            print("Writing pot_bare to xsf file")
            out_file = self.pot_bare_r

        elif ftype == 'pot_tot':
            print("Writing pot_tot to xsf file")
            out_file = self.pot_tot_r

        elif ftype == 'rho_bare':
            print("Writing rho_bare to xsf file")
            out_file = self.rho_bare_r

        elif ftype == 'rho_tot':
            print("Writing rho_tot to xsf file")
            out_file = self.rho_tot_r
        
        else:
            print('Undefined ftype')
            return

        with open(filedir, 'w') as v_file:
            print('Crystal\n PRIMVEC\n', end=' ',file=v_file)
            A = self.lattpara_unit[0]*self.sc[0]
            B = self.lattpara_unit[1]*self.sc[1]
            C = self.lattpara_unit[2]*self.sc[2]
            nx, ny, nz = self.fft_nx, self.fft_ny, self.fft_nz
            if self.ibrav == 4:
                print("%16.10f %16.10f %16.10f\n%16.10f %16.10f %16.10f\n%16.10f %16.10f %16.10f\n" % (A, 0, 0, -0.5*B, np.sqrt(3)/2*B, 0, 0, 0, C),end='', file=v_file)

            elif self.ibrav == 1:
                print("%16.10f %16.10f %16.10f\n%16.10f %16.10f %16.10f\n%16.10f %16.10f %16.10f\n" % (A, 0, 0, 0, B, 0, 0, 0, C),end='', file=v_file)

            print(' PRIMCOORD\n 0    1\n', end=' ',file=v_file)
            print(' BEGIN_BLOCK_DATAGRID_3D\n 3D_PWSCF\n BEGIN_DATAGRID_3D_UNKNOWN\n', end=' ',file=v_file)
            print("%8d %8d %8d\n" % (nx, ny, nz), end=' ', file=v_file)
            print(' 0.00 0.00 0.00\n', end=' ',file=v_file)

            if self.ibrav == 4:
                print("%16.10f %16.10f %16.10f\n%16.10f %16.10f %16.10f\n%16.10f %16.10f %16.10f\n" % (A, 0, 0, -0.5*B, np.sqrt(3)/2*B, 0, 0, 0, C),end='', file=v_file)

            elif self.ibrav == 1:
                print("%16.10f %16.10f %16.10f\n%16.10f %16.10f %16.10f\n%16.10f %16.10f %16.10f\n" % (A, 0, 0, 0, B, 0, 0, 0, C),end='', file=v_file)
            
            count = 0

            for k in range(nz):
                for j in range(ny):
                    for i in range(nx):
                        print("%16.10f" % np.real(out_file[i,j,k]), end='', file=v_file)
                        count+=1
                        if count%5==0:
                            print(end='\n', file=v_file)
            if (nx*ny*nz)%5 != 0:
                print(end='\n', file=v_file)

            print(' END_DATAGRID_3D\n', end=' ',file=v_file)
            print(' END_BLOCK_DATAGRID_3D', end=' ',file=v_file)

            v_file.close


#    def eps_interp(self, inttype='dis'):
#        g_epsind_reduced = self.g_epsind_tt[self.g_epsind_tt<=G_ind_cut]
#        g_epsind0_reduced = self.g_epsind_tt[self.g_epsind_tt<=G_ind_cut]
#        q_ind_reduced = self.q_ind_tt[self.g_epsind_tt<=G_ind_cut]
#        pot_reduced = self.pot_bare_k[self.g_epsind_tt<=G_ind_cut]
#        v_coul_reduced = vcoul_tmp_[self.g_epsind_tt<=G_ind_cut]
        

        

        



if __name__ == '__main__':

    cell = lab.bn_12x12
   # cell = lab.mos2_tt
    #cell = lab.bn_6x6
    print('Starting job for '+cell['prefix']+'...')
    #Eps0 =  read_eps.Epsmat(cell['folder']+'chi0mat.h5')
    #Eps1 = read_eps.Epsmat(cell['folder']+'chimat.h5')
    #np.save(cell['folder']+'chi0mat.npy', Eps1.mat)

    ttkw=PotCorr(cell)
    print(ttkw.prefix)
    print(ttkw.lattpara_unit)
   # print(type(mos2.lattpara_unit))

    ttkw.fft_init()
    ttkw.get_vcoul()
    ttkw.get_greenfunc()
    ttkw.read_chi()
   # mos2.get_pointc_pot_bare(position='center', kernel='3d')
    #print(mos2.v_coul)
   # mos2.read_rho_bare(mos2.folder+'inp/rho_6pad12_coarse.xsf')
    ttkw.read_rho_tot(ttkw.folder+'drho.xsf')
    #print(np.shape(mos2.rho_bare_r)) 
    #mos2.rho_bare_r[:,:,:30] = 0
    #mos2.rho_bare_r[:,:,-30:] = 0
    #mos2.read_rho_bare(cell['rho_bare'])
    #mos2.read_rho_tot(cell['rho_tot'])
    
    ttkw.rho2pot_tot(kernel='3d')
    ttkw.epsmat_init()
    #print(mos2.q_ind_tt)
    ttkw.get_epsmat(G_ind_cut = 1000, kernel='3d')
    ttkw.epsmat_inv(G_ind_cut = 1000)
    del ttkw.Chi0
    del ttkw.Chi1
    ttkw.pot_tot2bare(kernel='3d', G_ind_cut=600)

    #print(mos2.pot_tot_r)
    ttkw.pot2rho_bare(kernel='3d', ncharge=1)

    ttkw.write_xsf(ftype='rho_bare', filedir=ttkw.folder+'rho_bare.xsf')
    #mos2.write_xsf(ftype='rho_bare', filedir=mos2.folder+'rho_bare.xsf')
    info = psutil.virtual_memory()
    print(u'Occupied memory：',psutil.Process(os.getpid()).memory_info().rss/1024/1024, 'MB')
    print(u'Total memory：',info.total/1024/1024, 'MB')


    
    print(u'Percent：',info.percent,'%')
    print(u'Number of Cup：',psutil.cpu_count())


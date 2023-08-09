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

class Epsmat(object):

    def __init__(self, name_):
        self.name = name_
        self.feps = h5py.File(self.name, 'r')
        self.mat_diagonal = np.array(self.feps['mats/matrix-diagonal'])  # diagonal term of dielectric matrix
        
        self.qpts = np.array(self.feps['eps_header/qpoints/qpts']) # qpoints coordinates
        self.nmtx = np.array(self.feps['/eps_header/gspace/nmtx']) # Number of matrix elements BekerleyGW actually compute for each q-point.
        #self.mat_diagonal[:,10000, ]
        self.G_vec = np.array(self.feps['/mf_header/gspace/components']) # G-vectors in RHO G-space
        self.gind_eps2rho = np.array(self.feps['/eps_header/gspace/gind_eps2rho']) # convert Epsilon G-space to Rho G-space
        self.gind_rho2eps = np.array(self.feps['/eps_header/gspace/gind_rho2eps']) # convert Rho G-space to Epsilon G-space
        self.mat = np.array(self.feps['mats/matrix']) # Matrix elements

        self.G_vec_tuple_list = []
        for i_ in self.G_vec:
            i_ = tuple(i_)
            self.G_vec_tuple_list.append(i_)
        self.G_ind2vec = dict(enumerate(self.G_vec_tuple_list))
        self.G_vec2ind = {v:k for k,v in self.G_ind2vec.items()}

        self.qlist = []
        for i_ in self.qpts:
            i_ = tuple(i_)
            self.qlist.append(i_)
        self.q_ind2vec = dict(enumerate(self.qlist))
        self.q_vec2ind = {v:k for k,v in self.q_ind2vec.items()}
        




    

    

if __name__ == '__main__':
    f = Epsmat('mos2_eps_gw/eps0mat.h5')

'''
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
'''

def eps_model(Eps0, Eps1, Gz_): # interpolation of epsilon
    g0_ = [0,-1,Gz_]
    q0_ = Eps0.qpts[:10]
    q1_ = Eps1.qpts[1:6]
    q0_abs_ = np.sqrt((q0_[:10,0])**2+ (np.sqrt(3)/3*q0_[:10,0]+2*np.sqrt(3)/3*q0_[:10,1])**2)
    q1_abs_ = np.sqrt((q1_[1:6,0])**2+ (np.sqrt(3)/3*q1_[1:6,0]+2*np.sqrt(3)/3*q1_[1:6,1])**2)
    eps_qy_ = []
    for __ in range(3):
        g0_[1] = g0_[1]+1
        if __ == 0:
            for i in range(10):
                qind_ = i
                g_vec_ = tuple(g0_)
                gind_rho_ = Eps0.G_vec2ind[g_vec_]
                gind_eps0_ = Eps0.gind_rho2eps[qind_, gind_rho_]
                mat0_ = Eps0.mat_diagonal[qind_, gind_eps0_-1, 0]
                eps_qy_.append(mat0_)
            q_abs_ = np.hstack((np.array([]),q0_abs_))

        for i in range(5):
            qind_ = i+1
            g_vec_ = tuple(g0_)
            gind_rho_ = Eps1.G_vec2ind[g_vec_]
            gind_eps1_ = Eps1.gind_rho2eps[qind_, gind_rho_]
            mat1_ = Eps1.mat_diagonal[qind_, gind_eps1_-1,0]
            eps_qy_.append(mat1_)
        q_abs_ = np.hstack((q_abs_, q1_abs_+__*np.sqrt(3)*2/3))


    model_ = make_interp_spline(q_abs_,eps_qy_)
    return model_





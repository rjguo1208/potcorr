import potcorr
import psutil
import os
import lab

cell = lab.bn_12x12
# cell = lab.mos2_tt
#cell = lab.bn_6x6
print('Starting job for '+cell['prefix']+'...')
#Eps0 =  read_eps.Epsmat(cell['folder']+'chi0mat.h5')
#Eps1 = read_eps.Epsmat(cell['folder']+'chimat.h5')
#np.save(cell['folder']+'chi0mat.npy', Eps1.mat)

ttkw=potcorr.PotCorr(cell)
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

ttkw.epsmat_init_interp()
print(ttkw.g_rhoind_tt)
print(ttkw.g_tuple_q_tt)

'''
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
'''
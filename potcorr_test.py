import potcorr
import psutil
import os
import lab
import numpy as np
import pickle
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression
from sklearn.pipeline import make_pipeline


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

ttkw.epsmat_init_interp(G_vec2ind_dict=cell['folder_G_info']+'G_vec2ind.pkl',
                        G_ind2vec_dict=cell['folder_G_info']+'G_ind2vec.pkl')
#print(ttkw.g_tuple_floor_tt)
#print(ttkw.g_rhoind_tt)
#print(ttkw.g_tuple_q_tt)

with open(cell['folder_model']+'model_all_re.pkl','rb') as file:
    model_all_re = pickle.load(file)
with open(cell['folder_model']+'model_all_im.pkl','rb') as file:
    model_all_im = pickle.load(file)
with open(cell['folder_G_info']+'Gind_list.pkl','rb') as file:
    Gind_list = pickle.load(file)


print(Gind_list)

#model_all_re



#for j in Gind_list:
#    print(j)
#    for jp in Gind_list:
#        if f'model_{j}_{jp}' not in model_all_re:
#            print(f'Missing model_{j}_{jp}')
#            break
#        



#unique_tuples = set()
#for matrix in ttkw.g_tuple_floor_tt:
#    for row in matrix:
#        for tup in row:
#            unique_tuples.add(tup)

#print(f"Number of unique tuples: {len(unique_tuples)}")
#for tup in unique_tuples:
#    print(tup)



info = psutil.virtual_memory()
print(u'Occupied memory：',psutil.Process(os.getpid()).memory_info().rss/1024/1024, 'MB')
print(u'Total memory：',info.total/1024/1024, 'MB')
print(u'Percent：',info.percent,'%')
print(u'Number of Cpus：',psutil.cpu_count())

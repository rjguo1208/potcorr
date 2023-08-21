import potcorr
import psutil
import os
import lab
import numpy as np

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

ttkw.epsmat_init_interp(G_vec2ind_dict='D:/data/charged_defect/bn_2/G_vec2ind.pkl',
                        G_ind2vec_dict='D:/data/charged_defect/bn_2/G_ind2vec.pkl')
print(ttkw.g_tuple_floor_tt)
print(ttkw.g_rhoind_tt)
print(ttkw.g_tuple_q_tt)

Gind_list = []
for i in range(72959):
    iix0, iiy0, iiz0 = ttkw.G_ind2vec[i]
    if np.abs(iix0)<=2 and np.abs(iiy0)<=2 and np.abs(iiz0)<=25:
        Gind_list.append(i)
print(len(Gind_list))

model_list=[]

for j in Gind_list:
    for jp in Gind_list:
        model_name = f"D:/data/charged_defect/bn_2/chi_interp/chi_interp/model_{j}_{jp}.pkl"
        if os.path.isfile(model_name):
            #print(f"'{model_name}' is a file and exists!")
            model_list.append(f"{j}_{jp}")
        #else:
            #print(f"'{model_name}' is not a file or does not exist!")

print(len(model_list))


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

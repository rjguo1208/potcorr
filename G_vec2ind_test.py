import potcorr
import numpy as np
import lab
import pickle 

cell = lab.bn_12x12
print('Starting job for '+cell['prefix']+'...')

ttkw=potcorr.PotCorr(cell)
print(ttkw.prefix)
print(ttkw.lattpara_unit)

ttkw.read_chi()


with open(f'D:/data/charged_defect/bn_2/G_vec2ind.pkl','wb') as file:
    pickle.dump(ttkw.Chi1.G_vec2ind, file)

with open(f'D:/data/charged_defect/bn_2/G_ind2vec.pkl','wb') as file:
    pickle.dump(ttkw.Chi1.G_ind2vec, file)
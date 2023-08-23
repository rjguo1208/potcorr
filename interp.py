import numpy as np
import pickle
import potcorr

model_dir = 'D:/data/charged_defect/bn_2/model_interp/'
G_info_dir = 'D:/data/charged_defect/bn_2/G_info/'

with open(G_info_dir+'Gind_list.pkl','rb') as file:
    Gind_list = pickle.load(file)


print(Gind_list)


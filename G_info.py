import potcorr
import numpy as np
import lab
import pickle 


def extract_G(ttkw, dir_name):
    print('Starting extracting G vector and index for ', ttkw.prefix)
    print(ttkw.lattpara_unit)
    ttkw.read_chi()

    with open(dir_name+'G_vec2ind.pkl','wb') as file:
        pickle.dump(ttkw.Chi1.G_vec2ind, file)

    with open(dir_name+'G_ind2vec.pkl','wb') as file:
        pickle.dump(ttkw.Chi1.G_ind2vec, file)

    print('Extracting complete')

def filter_G_by_Gekin(ttkw, dir_name, ecut=5):
    print('Staring filtering G vector and index for ', ttkw.prefix)
    print(ttkw.lattpara_unit)
    ttkw.fft_init()
    ttkw.get_vcoul()
    J_to_Ry = 1 / (2.17987e-18)
    a_0 = 0.52917721067e-10
    Ang= 1e-10
    h_bar = 1.0545718e-34
    m_electron = 9.10938356e-31

    a = float(ttkw.lattpara_unit[0])*Ang
    c = float(ttkw.lattpara_unit[2])*Ang

    b1 = (2 * np.pi / a) * np.array([1, 1/np.sqrt(3), 0])
    b2 = (2 * np.pi / a) * np.array([0, 2/np.sqrt(3), 0])
    b3 = (2 * np.pi / c) * np.array([0, 0, 1])


    E_k_limit = ecut
    filtered_G_vectors = []
    Gind_list = []
    for i in range(72959):
        coeff = ttkw.Chi0.G_ind2vec[i]
        G_vector = coeff[0]*b1 + coeff[1]*b2 + coeff[2]*b3
        E_k = (h_bar**2 * np.linalg.norm(G_vector)**2) / (2 * m_electron)
        E_k_Ry = E_k * J_to_Ry
        if E_k_Ry <= E_k_limit:
            filtered_G_vectors.append(coeff)
            Gind_list.append(i)

    with open(dir_name+'filtered_G_vec.pkl','wb') as file:
        pickle.dump(filtered_G_vectors, file)

    with open(dir_name+'Gind_list.pkl','wb') as file:
        pickle.dump(Gind_list, file)

    print('Filtering complete')
    print('Number of filtered G:', len(Gind_list))

def filter_G_by_Gqekin(ttkw,dirname, ecut = 5):
    print('Staring filtering G vector and index for ', ttkw.prefix)
    print(ttkw.lattpara_unit)
    ttkw.fft_init()
    ttkw.get_vcoul()
    J_to_Ry = 1 / (2.17987e-18)
    a_0 = 0.52917721067e-10
    Ang= 1e-10
    h_bar = 1.0545718e-34
    m_electron = 9.10938356e-31
    a = float(ttkw.lattpara_unit[0])*Ang
    c = float(ttkw.lattpara_unit[2])*Ang

    b1 = (2 * np.pi / a) * np.array([1, 1/np.sqrt(3), 0])
    b2 = (2 * np.pi / a) * np.array([0, 2/np.sqrt(3), 0])
    b3 = (2 * np.pi / c) * np.array([0, 0, 1])

    E_k_limit = ecut
    filtered_G_vectors = []
    Gind_list = []
    E_Gqvec_tt = np.empty((nx_, ny_, nz_),dtype=float)
    nx_, ny_, nz_ = ttkw.fft_nx, ttkw.fft_ny, ttkw.fft_nz
    for i in range(nx_):
        for j in range(ny_):
            for k in range(nz_):
                G_q_vector = ttkw.fft_kxx_tt[i,j,k]*b1 + ttkw.fft_kyy_tt[i,j,k]*b2 + ttkw.fft_kxx_tt[i,j,k]*b3
                E_k = (h_bar**2 * np.linalg.norm(G_q_vector)**2) / (2 * m_electron)
                E_k_Ry = E_k * J_to_Ry
                E_Gqvec_tt[i,j,k] = E_k_Ry
                


if __name__ == '__main__':
    cell = lab.bn_12x12
    ttkw=potcorr.PotCorr(cell)
    extract_G(ttkw, cell['folder_G_info'])
    filter_G_by_Gekin(ttkw,cell['folder_G_info'])




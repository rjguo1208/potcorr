import potcorr
import numpy as np
import lab
import matplotlib.pyplot as plt
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression
from sklearn.pipeline import make_pipeline
import pickle 

cell = lab.bn_12x12
print('Starting job for '+cell['prefix']+'...')

ttkw=potcorr.PotCorr(cell)
print(ttkw.prefix)
print(ttkw.lattpara_unit)

ttkw.fft_init()
ttkw.get_vcoul()
ttkw.read_chi()

ttkw.get_chimat()
#print(ttkw.chi_mat[:,0,0])
plt.figure(figsize=(12,6), dpi=100)

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


E_k_limit = 12
filtered_G_vectors = []
Gind_list = []
for i in range(72959):
            coeff = ttkw.Chi0.G_ind2vec[i]
            G_vector = coeff[0]*b1 + coeff[1]*b2 + coeff[2]*b3
            E_k = (h_bar**2 * np.linalg.norm(G_vector)**2) / (2 * m_electron)
            E_k_Ry = E_k * J_to_Ry
            if E_k_Ry <= E_k_limit:
                filtered_G_vectors.append(G_vector)
                Gind_list.append(i)




print(len(Gind_list))

iiz_list = []
epsz_re_list = []
epsz_im_list = []
model_list = []

for i0_ in Gind_list:
    print(i0_)
    for i1_ in Gind_list:
        
        iix0, iiy0, iiz0 = ttkw.Chi0.G_ind2vec[i0_]
        iix1, iiy1, iiz1 = ttkw.Chi0.G_ind2vec[i1_]
        g0_ = [iix0, iiy0, iiz0]
        g1_ = [iix1, iiy1, iiz1]
        eps_Re_ = []
        eps_Im_ = []
        q_list=[]
        q_ = ttkw.Chi1.qpts[:]
        q0_ = ttkw.Chi0.qpts[:]
        q0_abs_ = np.sqrt((q0_[:,0])**2+ (np.sqrt(3)/3*(q0_[:,0])+2*np.sqrt(3)/3*(q0_[:,1]))**2)
        q1_abs_ = np.sqrt((q_[:,0])**2+ (np.sqrt(3)/3*(q_[:,0])+2*np.sqrt(3)/3*(q_[:,1]))**2)

        for i in range(0):
            #print(i)
            qind_ = i
            g_vec0_ = tuple(g0_)
            g_vec1_ = tuple(g1_)
            gind_rho0_ = ttkw.Chi0.G_vec2ind[g_vec0_]
            gind_rho1_ = ttkw.Chi0.G_vec2ind[g_vec1_]
            gind_eps0_ = ttkw.Chi0.gind_rho2eps[qind_, gind_rho0_]
            gind_eps1_ = ttkw.Chi0.gind_rho2eps[qind_, gind_rho1_]
            #if (gind_eps0_ > 1400 or gind_eps1_ > 1400):
            #    continue     
            mat_Re_ = ttkw.Chi0.mat[qind_,0,0,gind_eps0_-1,gind_eps1_-1,0]
            mat_Im_ = ttkw.Chi0.mat[qind_,0,0,gind_eps0_-1,gind_eps1_-1,1]
            eps_Re_.append(mat_Re_)
            eps_Im_.append(mat_Im_)
        #    q_abs_ = np.hstack((np.array([]),q0_abs_))
        q_abs_ = []

        for i in range(0,144):
            #print(i)
            qind_ = i
            g_vec0_ = tuple(g0_)
            g_vec1_ = tuple(g1_)
            gind_rho0_ = ttkw.Chi1.G_vec2ind[g_vec0_]
            gind_rho1_ = ttkw.Chi1.G_vec2ind[g_vec1_]
            gind_eps0_ = ttkw.Chi1.gind_rho2eps[qind_, gind_rho0_]
            gind_eps1_ = ttkw.Chi1.gind_rho2eps[qind_, gind_rho1_]
            #if (gind_eps0_ > 1000 or gind_eps1_ > 1000):    
            #    continue     
            try:
                mat_Re_ = ttkw.Chi1.mat[qind_,0,0,gind_eps0_-1,gind_eps1_-1,0]
                mat_Im_ = ttkw.Chi1.mat[qind_,0,0,gind_eps0_-1,gind_eps1_-1,1]
            except (IndexError):
                continue
            eps_Re_.append(mat_Re_)
            eps_Im_.append(mat_Im_)
            q_list.append(qind_)
            
        q_abs_ = np.hstack((q_abs_, q1_abs_[0:] ))

        y_re_data = np.array(eps_Re_)[:,np.newaxis]
        y_im_data = np.array(eps_Im_)[:,np.newaxis]
        if len(y_re_data) < 100:
            
            #excout +=1
            print(len(y_re_data))
            print(g_vec0_)
            print(g_vec1_)
            continue

        qpt = q_[q_list,:2]
        #print(qpt)
        degree = 9
        model_2d = make_pipeline(PolynomialFeatures(degree), LinearRegression())
        model_2d.fit(qpt, y_re_data)
        model_list.append(f"{i0_}_{i1_}")

print(len(model_list))

    #    q_inp = np.random.rand(500,2)
    #    y_pred = model_2d.predict(q_inp)
    #    q_zero = np.array([[0,0]])
    #    y_zero = model_2d.predict(q_zero)
    #    params = model_2d.named_steps['linearregression'].coef_
    #    intercept = model_2d.named_steps['linearregression'].intercept_
        
        #print('Intercept: ', intercept)
        #print('Coefficients: ', params)

    #    fig = plt.figure()
    #    ax = fig.add_subplot(111, projection='3d')

        # Create a scatter plot
    #    ax.scatter(qpt[:, 0], qpt[:, 1], y_re_data)
    #    ax.scatter(q_inp[:, 0], q_inp[:, 1], y_pred)
    #    ax.scatter(q_zero[:,0], q_zero[:,1], y_zero)
        # Set labels
    #    ax.set_xlabel('X')
    #    ax.set_ylabel('Y')
    #    ax.set_zlabel('Z')

    #    plt.show()

    #    with open(f'D:/data/charged_defect/bn/chi_interp/model_{i0_}_{i1_}.pkl','wb') as f:
    #        pickle.dump(model_2d,f)



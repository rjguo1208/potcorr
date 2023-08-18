
import numpy as np
import utli
import io_xsf 
import lab
import  matplotlib.pyplot as plt
import matplotlib
import matplotlib.colors as mcolors
import potcorr

#cell_list = [lab.bn_6x6,lab.bn_7x7,lab.bn_8x8,lab.bn_9x9,lab.bn_10x10, 
 #   lab.bn_11x11,lab.bn_12x12,lab.bn_13x13,lab.bn_14x14,lab.bn_15x15]
#cell_list = [lab.bn_6x6_relaxed, lab.bn_9x9_relaxed]
cell_list = [lab.bn_6x6, lab.bn_12x12]
#fig,ax = plt.subplots(figsize = (8,6), dpi = 100)
for cell in cell_list:
    
    '''
    Pot_b1 = io_xsf.Pot(cell['folder']+'1-dn/Pot_b.xsf')
    Pot_bh1 = io_xsf.Pot(cell['folder']+'1-dn/Pot_bh.xsf')
    Pot_b2 = io_xsf.Pot(cell['folder']+'2-dc/Pot_b.xsf')
    Pot_bh2 = io_xsf.Pot(cell['folder']+'2-dc/Pot_bh.xsf')
    Rho1 = io_xsf.Rho(cell['folder']+'1-dn/Rho.xsf')
    Rho2 = io_xsf.Rho(cell['folder']+'2-dc/Rho.xsf')

    #Rho1_up = io_xsf.Rho(cell['folder']+'1-dn/Rho_up.xsf')
    #Rho1_down = io_xsf.Rho(cell['folder']+'1-dn/Rho_down.xsf')
    #Rho2_up = io_xsf.Rho(cell['folder']+'2-dc/Rho_up.xsf')
    #Rho2_down = io_xsf.Rho(cell['folder']+'2-dc/Rho_down.xsf')

    Pot_h1 = io_xsf.Pot(cell['folder']+'1-dn/Pot_bh.xsf')
    Pot_h2 = io_xsf.Pot(cell['folder']+'2-dc/Pot_bh.xsf')
    drho = io_xsf.Rho(cell['folder']+'2-dc/Rho.xsf')

    v11 = Pot_bh1.pot
    v21 = Pot_b1.pot

    v12 = Pot_bh2.pot
    v22 = Pot_b2.pot

    rho1 = Rho1.rho
    rho2 = Rho2.rho

    dv1= utli.delta_v(v11, v21, align='None')
    dv2= utli.delta_v(v12, v22, align='None')
    drhod = utli.delta_v(rho1,rho2,align='None')

    Pot_h1.pot = dv1
    Pot_h2.pot = dv2
    drho.rho = drhod

    print(np.sum(Rho1.rho)*Rho1.omega/Rho1.nx/Rho1.ny/Rho1.nz)
    print(np.sum(Rho2.rho)*Rho1.omega/Rho1.nx/Rho1.ny/Rho1.nz)
    print(np.sum(drho.rho)*Rho1.omega/Rho1.nx/Rho1.ny/Rho1.nz)

    #print(np.sum(Rho1_up.rho)*Rho1_up.omega/Rho1_up.nx/Rho1_up.ny/Rho1_up.nz)
    #print(np.sum(Rho2_up.rho)*Rho1_up.omega/Rho1_up.nx/Rho1_up.ny/Rho1_up.nz)
    #print(np.sum(Rho1_down.rho)*Rho1_up.omega/Rho1_up.nx/Rho1_up.ny/Rho1_up.nz)
    #print(np.sum(Rho2_down.rho)*Rho1_up.omega/Rho1_up.nx/Rho1_up.ny/Rho1_up.nz)


    #print(Rho1.omega)
    io_xsf.write_xsf(Pot_h1, xsf_type='pot', filedir=cell['folder']+'1-dn/Pot_h.xsf')
    io_xsf.write_xsf(Pot_h2, xsf_type='pot', filedir=cell['folder']+'2-dc/Pot_h.xsf')
    io_xsf.write_xsf(drho, xsf_type='rho', filedir=cell['folder']+'drho.xsf')

    del Pot_b1, Pot_h1
    del Pot_b2, Pot_h2
    del Pot_bh1
    del Pot_bh2
    del Rho1, rho1
    del Rho2, rho2, drho,drhod
    del v11, v12, v21, v22
    '''


    #Rho1=io_xsf.Rho(cell['folder']+'1-dn/Rho.xsf')
    #Pot1=io_xsf.Pot(cell['folder']+'1-dn/Pot_h.xsf')

    #Rho2=io_xsf.Rho(cell['folder']+'2-dc/Rho.xsf')
    #Pot2=io_xsf.Pot(cell['folder']+'2-dc/Pot_h.xsf')

    #print(np.sum(Pot1.pot*Rho1.rho)*Rho1.omega/(Rho1.nx)/(Rho1.ny)/(Rho1.nz)*0.5)
    #print(np.sum(Pot2.pot*Rho2.rho)*Rho2.omega/(Rho2.nx)/(Rho2.ny)/(Rho2.nz)*0.5)




    #'''
    #Pot_tot = io_xsf.Pot(cell['folder_out']+'pot_tot.xsf')
    #Pot_tot = io_xsf.Rho(cell['folder_dft']+'drho.xsf')
    #Pot_tot = io_xsf.Pot(cell['folder_dft']+'dv.xsf')
    Rho1 = io_xsf.Rho(cell['folder']+'rho_bare.xsf')
    chgden = Rho1.rho

    nx, ny, nz = Rho1.nx, Rho1.ny, Rho1.nz
    print(nz)

    #potential_2d = np.mean(potential[:,:,80], axis=2)
    chgden_2d = chgden[:,:,90]
    chgden_2d = np.sum(chgden[:,:,80:90], axis=2)
    x = np.arange(0, chgden_2d.shape[0])
    y = np.arange(0, chgden_2d.shape[1])

    #x = np.copy(a)
    #y = np.copy(b)


    X, Y = np.meshgrid(x, y)

    x = X.flatten()-0.5*Y.flatten()
    y = Y.flatten()*np.sqrt(3)*0.5 #- X.flatten()*0.5
    z = chgden_2d.flatten()

    norm= mcolors.TwoSlopeNorm(vmin=-0.3*z.max(), vmax=0.3*z.max(), vcenter=0)
    fig, ax = plt.subplots(dpi=150)
    hb = ax.scatter(x, y, c=z, 
                    #gridsize=100,
                    #vmin=-0.08,
                    #vmax=-0.0,
                    norm=norm,
                  #    cmap='coolwarm')
                    cmap='bwr')
    ax.set_aspect('equal', adjustable='box')
    #ax.set_xlim(0, 180)
    #ax.set_ylim(0, 180)
    cb = fig.colorbar(hb, ax=ax)
    cb.set_label('Potential')
    plt.show()
    #'''

    '''
    Pot_tot = io_xsf.Pot(cell['folder_out']+'pot_tot.xsf')
    Pot_dv = io_xsf.Pot(cell['folder_dft']+'dv.xsf')

    nx1, ny1, nz1 = Pot_tot.nx, Pot_tot.ny, Pot_tot.nz
    nx2, ny2, nz2 = Pot_dv.nx, Pot_dv.ny, Pot_dv.nz

    print(nx1, ny1, nz1)   
    print(nx2, ny2, nz2)

    potential_tot = Pot_tot.pot[:,:,int(0.45*nz1):int(0.55*nz1)]/1.01213
    potential_dv = Pot_dv.pot[:,:,int(0.45*nz2):int(0.55*nz2)]

    A = float(Pot_tot.a1[0])
    C = 0.2*float(Pot_tot.a3[2])

    print(A)
    print(nz1)
    print(nz2)


    dist_arr_tot = utli.distance_array(potential_tot.shape[0], potential_tot.shape[1], potential_tot.shape[2],A,A,C, defect_loc='center')
    dist_arr_dv = utli.distance_array(potential_dv.shape[0], potential_dv.shape[1], potential_dv.shape[2], A,A,C,defect_loc='center')

    plt.scatter(dist_arr_tot.flatten(), potential_tot.flatten(), label='model')
    plt.scatter(dist_arr_tot.flatten(), potential_tot.flatten()*1.05, label='model with interp')
    plt.scatter(dist_arr_dv.flatten(), potential_dv.flatten(), label='dft')

    plt.xlabel('Distance from Defect')
    plt.ylabel('Potential')
    plt.legend()
    plt.show()
    '''

    #plt.imshow(potential_2d, cmap='hot', interpolation='nearest')
    #plt.colorbar(label='Potential')
    #plt.xlabel('x')
    #plt.ylabel('y')
    #plt.title('Potential distribution in x-y plane')
    #plt.show()


    #Rho = io_xsf.Rho(cell['rho_bare'])
    #print(np.sum(Rho.rho)*Rho.omega/Rho.nx/Rho.ny/Rho.nz)

    '''
    ttkw=potcorr.PotCorr(cell)
    ttkw.fft_init()
    ttkw.get_vcoul()
    ttkw.read_rho_bare(ttkw.folder+'1-dn/Rho.xsf')
    ttkw.rho2pot_bare(kernel='3d')
    ttkw.write_xsf(ftype='pot_bare', filedir=ttkw.folder+'1-dn/pot_h_calc.xsf')

    ttkw2=potcorr.PotCorr(cell)
    ttkw2.fft_init()
    ttkw2.get_vcoul()
    ttkw2.read_rho_bare(ttkw2.folder+'2-dc/Rho.xsf')
    ttkw2.rho2pot_bare(kernel='3d')
    ttkw2.write_xsf(ftype='pot_bare', filedir=ttkw2.folder+'2-dc/pot_h_calc.xsf')
    '''

    '''
    Rho1=io_xsf.Rho(cell['folder']+'1-dn/Rho.xsf')
    Pot1=io_xsf.Pot(cell['folder']+'1-dn/Pot_h.xsf')
    Pot1_calc=io_xsf.Pot(cell['folder']+'1-dn/pot_h_calc.xsf')
    Rho2=io_xsf.Rho(cell['folder']+'2-dc/Rho.xsf')
    Pot2=io_xsf.Pot(cell['folder']+'2-dc/Pot_h.xsf')
    Pot2_calc=io_xsf.Pot(cell['folder']+'2-dc/pot_h_calc.xsf')

    print(np.sum(Pot1.pot*Rho1.rho)*Rho1.omega/(Rho1.nx)/(Rho1.ny)/(Rho1.nz)*0.5)
    print(np.sum(Pot2.pot*Rho2.rho)*Rho2.omega/(Rho2.nx)/(Rho2.ny)/(Rho2.nz)*0.5)

    print(np.sum(Pot1_calc.pot*Rho1.rho)*Rho1.omega/(Rho1.nx)/(Rho1.ny)/(Rho1.nz)*0.5)
    print(np.sum(Pot2_calc.pot*Rho2.rho)*Rho2.omega/(Rho2.nx)/(Rho2.ny)/(Rho2.nz)*0.5)
    '''
    '''
    ttkw=potcorr.PotCorr(cell)
    ttkw.fft_init()
    ttkw.get_vcoul()
    ttkw.read_rho_bare(ttkw.folder+'drho.xsf')
    ttkw.rho2pot_bare(kernel='3d')
    ttkw.write_xsf(ftype='pot_bare', filedir=ttkw.folder+'dpot_h_calc.xsf')
    '''
    
    #Rho1=io_xsf.Rho(cell['folder']+'drho.xsf')
    #print(float(Rho1.a1[0]))
    #Pot1=io_xsf.Pot(cell['folder']+'1-dn/Pot_b.xsf')
    print(np.sum(Rho1.rho)*Rho1.omega/(Rho1.nx)/(Rho1.ny)/(Rho1.nz))
    print(np.sum(Rho1.rho[int(Rho1.nx/2)-54:int(Rho1.nx/2)+54,int(Rho1.nx/2)-54:int(Rho1.nx/2)+54,: ])*Rho1.omega/(Rho1.nx)/(Rho1.ny)/(Rho1.nz))
    #print(np.average(Pot1.pot)/np.average(Rho1.rho))
    #bgc = np.average(Rho1.rho)
    #Rho1.rho = Rho1.rho-bgc
    #print(np.sum(Rho1.rho)*Rho1.omega/(Rho1.nx)/(Rho1.ny)/(Rho1.nz))
    #print(np.sum(Pot1.pot*bgc)*Rho1.omega/(Rho1.nx)/(Rho1.ny)/(Rho1.nz))
    #print(np.sum(Pot1.pot)*Rho1.omega/(Rho1.nx)/(Rho1.ny)/(Rho1.nz)*0.5)
    #ax.plot(np.linspace(-0.5*float(Rho1.a1[0]), 0.5*float(Rho1.a1[0]), int(Rho1.nx)), np.sum(Rho1.rho,axis=(1,2)))
    #plt.subplots(dpi = 200)  
    #im1=plt.imshow(Rho1.rho[:,:,30],extent=[0, 216,0, 216], cmap=matplotlib.cm.hot, interpolation='nearest', origin="lower")#pl.cm.jet   
    #plt.colorbar(im1)  
    #plt.show()
    #plt.subplots(dpi = 200)  
    #im1=plt.imshow(np.abs(pot_k)[:20,:20,0],extent=[0, 20,0, 20], cmap=matplotlib.cm.hot, interpolation='nearest', origin="lower")#pl.cm.jet  
    #plt.colorbar(im1) 
 
#plt.show()
import numpy as np

Bohr_R = 0.52917721067




def dat2xsf(dat_inp, xsf_outp):
    file = open(dat_inp)
    f = file.readlines()

    nx, ny, nz = int(f[1].split()[0]), int(f[1].split()[1]), int(f[1].split()[2])
    natom, ntype = int(f[1].split()[6]),int(f[1].split()[7])
    alat = float(f[2].split()[1])
    a1 = Bohr_R*alat*np.array([float(f[3].split()[0]),float(f[3].split()[1]),float(f[3].split()[2])])
    a2 = Bohr_R*alat*np.array([float(f[4].split()[0]),float(f[4].split()[1]),float(f[4].split()[2])])
    a3 = Bohr_R*alat*np.array([float(f[5].split()[0]),float(f[5].split()[1]),float(f[5].split()[2])])
    Nheader = 7+natom+ntype
    vl = ''.join(f[Nheader:])
    vn = [float(i) for i in vl.split()]
    __rho__ = np.array(vn).reshape([nz,ny,nx])
    rho = __rho__.transpose(2,1,0) 


    with open(xsf_outp, 'w') as v_file:
        print('Crystal\n PRIMVEC\n', end=' ',file=v_file)
        print("%16.10f %16.10f %16.10f\n%16.10f %16.10f %16.10f\n%16.10f %16.10f %16.10f\n" 
            % (float(a1[0]), float(a1[1]), float(a1[2]),
                float(a2[0]), float(a2[1]), float(a2[2]), 
                float(a3[0]), float(a3[1]), float(a3[2])),
                end='', file=v_file)
        print(' PRIMCOORD\n 0    1\n', end=' ',file=v_file)
        print(' BEGIN_BLOCK_DATAGRID_3D\n 3D_PWSCF\n BEGIN_DATAGRID_3D_UNKNOWN\n', end=' ',file=v_file)
        print("%8d %8d %8d\n" % (nx, ny, nz), end=' ', file=v_file)
        print(' 0.00 0.00 0.00\n', end=' ',file=v_file)
        print("%16.10f %16.10f %16.10f\n%16.10f %16.10f %16.10f\n%16.10f %16.10f %16.10f\n" 
            % (float(a1[0]), float(a1[1]), float(a1[2]),
                float(a2[0]), float(a2[1]), float(a2[2]), 
                float(a3[0]), float(a3[1]), float(a3[2])),
                end='', file=v_file)

        count = 0
        for k in range(nz):
            for j in range(ny):
                for i in range(nx):
                    print("%16.10f" % rho[i,j,k], end='', file=v_file)
                    count+=1
                    if count%5==0:
                        print(end='\n', file=v_file)
        if (nx*ny*nz)%5 != 0:
            print(end='\n', file=v_file)
        
        print(' END_DATAGRID_3D\n', end=' ',file=v_file)
        print(' END_BLOCK_DATAGRID_3D', end=' ',file=v_file)


fname = 'D:\Code\charged_defect\Rho.dat'
file_dir = 'D:\Code\charged_defect\Rho.xsf'
dat2xsf(fname, file_dir)
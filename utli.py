import numpy as np
import io_xsf 
import lab

def rho_expand(rho_array, target_fftn=np.array([241,241,241]), mode='constant'):
    pad_x = (target_fftn[0]-rho_array.shape[0])//2
    rest_x = (target_fftn[0]-rho_array.shape[0])%2

    pad_y = (target_fftn[1]-rho_array.shape[1])//2
    rest_y = (target_fftn[1]-rho_array.shape[1])%2

    print("pad_x: ", pad_x)
    print("rest_x: ", rest_x)
    print("pad_y: ", pad_y)
    print("rest_y: ", rest_y)

    extended_rho = np.pad(rho_array, 
                          pad_width=((pad_x, pad_x + rest_x), (pad_y, pad_y + rest_y), (0, 0)),
                            mode=mode, constant_values=0)
    print(rho_array.shape,'-->', extended_rho.shape)
    return extended_rho

def downsample_3d_array(arr, new_shape):
    s = [arr.shape[i] // new_shape[i] * new_shape[i] for i in range(3)]  # Find the nearest smaller shape divisible by new_shape
    print(s)
    arr = arr[:s[0], :s[1], :s[2]]  # Truncate the array to the new shape
    sh = (new_shape[0], arr.shape[0] // new_shape[0], 
          new_shape[1], arr.shape[1] // new_shape[1], 
          new_shape[2], arr.shape[2] // new_shape[2])  # Calculate the new shape
    print(sh)
    return arr.reshape(sh).mean(5).mean(1).mean(2)  # Take the mean along each dimension


def delta_v(arr1, arr2, align='None'):
    if align=='None':
       dv = arr2-arr1
    elif align=='vaccum':
       vacalign=np.average(arr2[:,:,0]-arr1[:,:,0])
       dv = arr2-arr1-vacalign
       
    return dv

def distance_array(nx, ny, nz, a,b,c, defect_loc='center'):
    x = np.linspace(0, a, nx)
    y = np.linspace(0, b, ny)
    z = np.linspace(0, c, nz)

    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    
    if defect_loc == 'center':
        defect_loc = [0.5*a, 0.5*b, 0.5*c]

    return np.sqrt((X - defect_loc[0])**2 + (Y - defect_loc[1])**2 + (Z - defect_loc[2])**2)
    

if __name__ == "__main__": 
  cell = lab.mos2_tt
  old_sc = np.array([6,6,1])
  new_sc = np.array([12,12,1])

  Rho = io_xsf.Rho(cell['folder']+'rho_test.xsf')
  a1_ = list(map(float, Rho.a1))
  a2_ = list(map(float, Rho.a2))
  a1 = [i*2 for i in a1_]
  a2 = [i*2 for i in a2_]
  Rho.a1 = a1
  Rho.a2 = a2



  rho = Rho.rho

  rho_12x12 = rho_expand(rho, np.array([361,361,241]),mode='constant')


  rho_12x12_coarse = downsample_3d_array(rho_12x12, new_shape=np.array([180,180,120]))

  print(rho_12x12_coarse.shape)
  Rho.nx, Rho.ny, Rho.nz = rho_12x12_coarse.shape
  Rho.rho = rho_12x12_coarse

  io_xsf.write_xsf(Rho, xsf_type='rho', filedir=cell['folder']+'rho_6pad12_coarse.xsf')



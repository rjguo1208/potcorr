import numpy as np

def equation_solver(a, b, c):
    delta = b*b - 4*a*c
    print(delta)
    if delta < 0:
        print('No real root')
    elif delta == 0:
        print('There are two equal roots')
        x = (-b+np.sqrt(b*b-4*a*c))/2*a
        print(x)
    else:
        print('There are two different roots')
        x1 = (-b+np.sqrt(b*b-4*a*c))/2*a
        x2 = (-b-np.sqrt(b*b-4*a*c))/2*a
        print('First root %6f, Second root %6f' % (x1, x2))


kzl = [0.6763888, 0.323611111, 0.3427253, 0.657274622]
kzl_=[0.323524, 0.3428363, 0.6571637, 0.6764760]
kzlarr = np.array(kzl)
kzl_arr = np.array(kzl_)
print(kzlarr[:, np.newaxis])
closest_indices = np.argmin(np.abs(kzlarr[:, np.newaxis] - kzl_arr), axis=1)
print(np.abs(kzlarr[:, np.newaxis] - kzl_arr))
print(closest_indices)
kzl_fin=kzl_arr[closest_indices]
#rmat_ = rmat_[closest_indices, :][:, closest_indices]
print("kzl",kzl)
print("kzl_",kzl_fin)
#print(kzl_[1])

 
        

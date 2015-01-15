import numpy as np
import pyfits as pf
import copy
import pylab as pl
import matplotlib.pyplot as plt

from ring_area_calc_many import get_circle_area_for_pix

def main3():
    px, py = np.meshgrid(
        np.arange(40) + 0.5,
        np.arange(40) + 0.5,
        )
    
    r = np.pi * 10
    A, intersec, intersected = get_circle_area_for_pix((px, py), r)
    
    print '\n'.join(map(
        lambda a: '{:.4f} {:.4f} {:.4f}'.format(*a),
        zip(px.flatten(), py.flatten(), A.flatten())  #, intersec #, intersected
        ))
    
    plt.close()
    plt.scatter(px, py, c=A)
    plt.show()
    
    print np.sum(A)

    
    
if __name__ == '__main__':
    #main()
    main3()
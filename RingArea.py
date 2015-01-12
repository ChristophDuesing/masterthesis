
# Program by Christoph Duesing
# This algorithm producess a the Powerspectrum of a random field 
import numpy as np
import pyfits as pf
import copy
import pylab as pl
import matplotlib.pyplot as plt



z = np.linspace(0.,99.,100.)

x = len(z)

a, b = np.meshgrid(np.linspace(0,x,x) , np.linspace(0,x,x))
c, d = np.meshgrid(np.random.ranf(x) , np.random.ranf(x))

#subsub=np.random.randint(10,size=(100,100))

#subsub = c*d
subsub = a*b
#subsub2 = np.cos(subsub)**2


r_min, r_max, a0,b0 =  1, x/2.,  x/2.,x/2. 
#r_grid = np.logspace(np.log10(r_min), np.log10(r_max),16)   
r_grid = np.linspace(r_min, r_max,64)
r_grid_mean = 0.5 * (r_grid[1:] + r_grid[:-1])



dists_sq = (a-a0)**2 + (b-b0)**2

radius_mask_func = lambda dists_sq, r: dists_sq <= r**2

masks_ringshape = []
Areas_True = []
for idx, r in enumerate(r_grid[:-1]):
    subsub_2 = subsub.copy()
    check_pix =[]
    Area_True = np.pi*(r_grid[idx+1]**2-r_grid[idx]**2)
    
    masks_ringshape.append(np.logical_and(
        radius_mask_func(dists_sq,  r_grid[idx+1]),
        np.logical_not(radius_mask_func(dists_sq,  r_grid[idx]))
        ))
    
    Areas_True.append(Area_True)
    
    from scipy.integrate import quad
    
    idx=15
    r = 15.
    r_max = r_grid[19]
    print r_max
    r_min = r_grid[18]
    b= r_max-r_min
    print r_min
    
    Check_Range = np.logical_and(
    radius_mask_func(dists_sq,  r_max + np.sqrt(2)),
    np.logical_not(radius_mask_func(dists_sq,  r_min-np.sqrt(2))
    ))
    
    subsub_2[Check_Range] =1
    xval, yval = np.where(subsub_2==1) 
    coordinates = zip(xval, yval)
    for i in coordinates:
        
        x,y = i
        x,y =x-49.5,y-49.5
        print x,y
        print 
        
        x_min,x_max,y_min,y_max = x-0.5,x+0.5,y-0.5,y+0.5
        
        print x_min,x_max,y_min,y_max
        
        #Check for interceptionpoints between cicle and pixel
        
        x_lower_max = np.sqrt(abs(r_max**2 - y_min**2))
        x_upper_max = np.sqrt(abs(r_max**2 - y_max**2))
        y_left_max = np.sqrt(abs(r_max**2 - x_min**2))
        y_right_max = np.sqrt(abs(r_max**2 - x_max**2))
        
        
        #x_lower_min = np.sqrt(abs(r_min**2 - y_min**2))
        #x_upper_min = np.sqrt(abs(r_min**2 - y_max**2))
        
        print  x_lower_max, x_upper_max, y_left_max, y_right_max
        print
        print r_min, r_max
        
        #Integrate Area 
        def integrand_rmax(x,r_max):
            return 2*np.sqrt(r_max**2-x**2)
        def integrand_rmin(x,r_min):
            return 2*np.sqrt(r_min**2-x**2)
            
        #if x_lower_min >= x_min and x_lower_min <= x_max:
            #x_min = x_lower_min
        
        #if x_lower_max <= x_max and x_lower_max >= x_min:
            #x_max = x_lower_max
        
        '''CASE 1 just r_max hits the pixel'''
        if x_min < x_lower_min and x_min < x_upper_min:
            if x_lower_max <= x_max and x_lower_max >= x_min:
                x_max = x_lower_max
            I = quad(integrand_rmax, x_min, x_max,args=(r_max))
        
        '''CASE 2 just r_min hits the pixel'''
        if x_max > x_lower_max and x_max > x_upper_min:   
            if x_lower_min >= x_min and x_lower_min <= x_max:
                x_min = x_lower_min
            I = quad(integrand_rmin, x_min, x_max,args=(r_min))
        '''CASE 3 both r_min and r_max hit the pixel'''
        if x_lower_min > x_min and x_lower_max < x_max:
            I_1 = quad(integrand_rmax, x_min, x_max,args=(r_max))
            I_2 = quad(integrand_rmin, x_min, x_max,args=(r_min))
            I = I_1[0]-I_2[0]
        print I
        
        
        
        
        
        
        
        
        
        
        
        ##print coordinates
    ##check_pix.append(subsub[Check_Range])
        
    #print check_pix

     
    #x_min,x_max,y_min,y_max = 
    
#for mask_part in masks_ringshape:
    #pl.close()
    #pl.figure()
    #pl.rc('text', usetex=True)
    #pl.rc('font', family='serif')
    #pl.xlabel(r" $\log_{10}$ Radius $|k|$")
    #pl.ylabel(r"Velocities     $\displaystyle\frac{\text{km}}{\text{s}}$")
    #pl.imshow((mask_part), interpolation='nearest', cmap='cubehelix')
    #pl.colorbar()
    #pl.show()
#final = np.mean(subsub[masks_ringshape][0])  

Areas_Estimated =[]   
means_of_cicular_grids_spek = [ ]   
for mask_ringshape in masks_ringshape:
    means_of_cicular_grids_spek.append(np.mean(subsub[mask_ringshape])) 
    Areas_Estimated.append(np.count_nonzero(mask_ringshape))
        #means_of_cicular_grids_spek.append(np.count_nonzero(aperture_sample))
#plt.close()          
#plt.loglog(r_grid_mean,means_of_cicular_grids_spek,'xb')
##plt.savefig('/export/data1/cduesing/Figures/Problems/CheckPlot_'+Zusatz+datei.split('.')[0]+'_'+str(z)+'.png', dpi=300)   
#plt.show()
#plt.clf()



t = zip(Areas_True,Areas_Estimated)

diff = np.subtract(Areas_True,Areas_Estimated)
diff
div = np.divide(Areas_True,Areas_Estimated)
div






    

    
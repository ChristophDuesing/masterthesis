# Program by Christoph Duesing
# This algorithm producess a the Powerspectrum of a random field 
import numpy as np
import pyfits as pf
    
    
def open_fits(fit):
    hdu = pf.open(fit)  
    primhdu = hdu[0]
   
    return primhdu
    
    
def convert_pixel_to_world(header):
    from kapteyn import wcs
    xdim, ydim, zdim = proj.naxis
    proj = wcs.Projection(header)
    proj_spec_1d = proj.sub(axes=[proj.specaxnum])
    proj_imag_2d = proj.sub(axes=[proj.lonaxnum, proj.lataxnum])
    
    
def convert_spec_channel_to_velocity():
    from kapteyn import wcs
    print proj_spec_1d.toworld((z,))[0] / 1000. 
    

def get_shape(fit):
    print fit.shape  
    
    
def quadratic_shape(fit):
    z,y,x = fit.shape
        # For galactic coordinates in 20 times 20 degrees every cube is based on :
    # (945, 373, 237), thus the y value needs a reshape
    a = 0.
    if y > x:
        a = y - x 
        max_y, min_y = (y-a/2), (a/2)
        fit_new = fit[:,min_y:max_y,:]
        print "y was reshaped" 
    elif x > y:
        a = x - y  
        max_x, min_x = (x-a/2), (a/2)
        fit_new = fit[:,:,min_y:max_y]# [z,y,x]!!    
        print "x was reshaped"
    if y==x:
        fit_new = fit 
    return fit_new
    

def restrict_velocities(datacube):
    z,y,x = datacube.shape
    datacube = datacube[222.5:722.5,: ,: ]
    return datacube
    
    
def inspect(datacube):    
    l = []
    for i in range(len(datacube)):
        planemean = np.mean(datacube[i])
        l.append(planemean)
    return l 
    #return np.mean(datacube, axis=0)

    
def get_power_spektra(datacube):
    power_spektra = [] 
    for i in datacube:
        fft_spek = np.fft.fft2(i)
        fft_spek = np.fft.fftshift(fft_spek)
        pow_spek = (np.abs(fft_spek)**2)/len(fft_spek)
        power_spektra.append(pow_spek)
    return power_spektra 

def get_cube(imgs, datei, cube_header, name):
    fits = pf.PrimaryHDU(np.array(imgs), header=cube_header)
    fits.writeto(datei.split('.')[0]+"_"+name+".fits",clobber=True)
    return 
    

def get_sub_cube(cube, datei, header):
    get_shape(cube)
    z,y,x = cube.shape
    i,j = np.arange(0,4), np.arange(0,4)
    Small_Cubes = []  
    for elem_i in i:
        print elem_i
        for elem_j in j:
            print elem_j
            datacube = cube[:,y*elem_i/4:y*(elem_i+1)/4 ,x*elem_j/4:x*(elem_j+1)/4 ]    
            get_cube(datacube, datei, header, "SUB_SAMPLE_"+str(elem_i+1)+"_x_"+str(elem_j+1) )
            Small_Cubes.append(datacube)
            
    return Small_Cubes
       
    
def vis2d(a):
    import pylab 
    import matplotlib.pyplot as plt
    
    
    plt.close()
    fig = plt.figure(figsize=(15, 5))
    ax1 = fig.add_subplot(1,4,1)
    #ax2 = fig.add_subplot(1,4,2)
    #ax3 = fig.add_subplot(1,4,3)
    #ax4 = fig.add_subplot(1,4,4)
    imdict = dict(interpolation='nearest', origin='lower')
    datas = [a]
    axes = [ax1]

    for d, ax in zip(datas, axes):
        ax.set_aspect('equal')
        ax.imshow(d, **imdict)
        im = ax.imshow(d, **imdict)
        pylab.colorbar(im, ax=ax , orientation='horizontal')
    pylab.savefig('/export/data1/cduesing/Figures/Spektra/Test_Picture.png')   
    plt.show()
          
    
def get_power_spek_grid(power_spektra_list):
    import pylab as pl
    import copy
    import matplotlib.pyplot as plt
    
    
    means_of_cicular_grids_spek_all = []
    
    x = len(power_spektra_list[0])
    r_min, r_max, a0,b0 =  1, x/2.,  x/2.,x/2. 
    r_grid = np.logspace(np.log10(r_min), np.log10(r_max),64)
    #find out gridding here ! Check with the results of 1024 from yesterday mate ;) 
    #r_grid = np.linspace(r_min, r_max,256)
    
    
    a, b = np.meshgrid(np.linspace(0,x,x) , np.linspace(0,x,x))
    
    #aperture = np.zeros_like(a)
    
    dists_sq = (a-a0)**2 + (b-b0)**2
    radius_mask_func = lambda dists_sq, r: dists_sq <= r**2

    masks_ringshape = []
    for idx, r in enumerate(r_grid[:-1]):  
        masks_ringshape.append(np.logical_and(
            radius_mask_func(dists_sq,  r_grid[idx+1]),
            np.logical_not(radius_mask_func(dists_sq,  r_grid[idx]))
            ))
    
    r_grid_mean = 0.5 * (r_grid[1:] + r_grid[:-1])
    
    for z, i in enumerate(power_spektra_list):
        means_of_cicular_grids_spek = []
        #pl.close()
        #pl.figure()
        #pl.imshow(i)
        #pl.show()
    
        print z
        #need quadratic input thus x=y

        for mask_ringshape in masks_ringshape:
        #for idx, r in enumerate(r_grid[:-1]):  
            #mask_ringshape = np.logical_and(
                #radius_mask_func(dists_sq,  r_grid[idx+1]),
                #np.logical_not(radius_mask_func(dists_sq,  r_grid[idx]))
                #)
            #aperture_sample  = copy.deepcopy(aperture)
            #aperture_sample[mask_ringshape]=1
            
            ##pl.close()
            ##pl.figure()
            ##pl.imshow(aperture_sample)
            ##pl.show()
            
            ##HOW to Mask
            ##import numpy as np
            ##a = np.arange(100)
            ##a
            ##m = (10 < a) & (a < 20)
            ##a[m]
            
            #pow_spek_grid = i * aperture_sample
            
            ##pl.close()
            ##pl.figure()
            ##pl.imshow(pow_spek_grid)
            ##pl.show()
            
            ##print pow_spek_grid
            ##print np.count_nonzero(aperture_sample)
            means_of_cicular_grids_spek.append(np.sum(i[mask_ringshape])/np.count_nonzero(i[mask_ringshape]))
            
            #means_of_cicular_grids_spek.append(np.mean(i[mask_ringshape])) #Warum klappt es nicht ????????
            #means_of_cicular_grids_spek.append(np.count_nonzero(aperture_sample))
        
       
        if z >290 and z <321:
                
            plt.close()   
            ##plt.plot(r_grid_mean,means_of_cicular_grids_spek,'xr')
            ##print means_of_cicular_grids_spek
            
            plt.loglog(10**(r_grid_mean),means_of_cicular_grids_spek,'xb')
            plt.savefig('/export/data1/cduesing/Figures/Spektra/PS_'+datei.split('.')[0]+'_'+str(z)+'.png', dpi=300)   
            #plt.show()
            #plt.clf()
        means_of_cicular_grids_spek_all.append(means_of_cicular_grids_spek)
        
    
    return means_of_cicular_grids_spek_all, r_grid_mean
    
    
def waterfall_plot(f, rgridding, datei, attachment):
    import pylab as pl
    
    pl.close()
    pl.figure()
    pl.rc('text', usetex=True)
    pl.rc('font', family='serif')
    pl.xlabel(r"Radii $|k|$")
    pl.ylabel(r"Velocities     $\displaystyle\frac{\text{km}}{\text{s}}$")
    pl.imshow(np.log10(f), extent=(np.log10(rgridding[0]), np.log10(rgridding[1]), -400, 400), interpolation='nearest', origin='lower' ,cmap='cubehelix')
    pl.colorbar()
    pl.gca().set_aspect('auto')
    pl.savefig('/export/data1/cduesing/Figures/Spektra/WP_new'+attachment+datei.split('.')[0]+'.png', dpi=300)
    #pl.show()
    
    
    
    
if __name__ == '__main__':
    
#for datei in Cubes_4x4:   
    
    #datei = input("Need some input for fitting \n")
    datei = ("CAR_NorSky30-50_80.fits")
    
    cube = open_fits(datei)
    cube_header = cube.header
    cube_data = cube.data
    get_shape(cube_data)
    
    #cube_data = restrict_velocities(cube_data)
    print "save the the cube"
    get_shape(cube_data)

    
    #print "Data need a quadratic shape"
    data = quadratic_shape(cube_data)
    #data = cube_data

    get_shape(data)
    get_cube(cube_data, datei, cube_header, "RESHAPED")
    mean_velocities = inspect(data)
    #vis(mean_velocities)
    
    ''' Build some subcubes '''
    Cubes_4x4 = get_sub_cube(data, datei, cube_header)
    
    
    FFT_Power_spektra = get_power_spektra(data)
    #get_cube(FFT_Power_spektra, datei, cube_header, "POWERSPEK")
    
    powerspektrum, rgridding = get_power_spek_grid(FFT_Power_spektra)
    powerspektrum = np.asarray(powerspektrum)
    #powerspektrum = powerspektrum.reshape(500.,127.)
    
    waterfall_plot(powerspektrum, rgridding, datei,'Basic_')
    #vis(powerspektrum, "Powerspectra", r"Velocities     $\displaystyle (\frac{km}{s})$", "Waterfallplot_",datei) 
    
    
    #for index, Cubus in enumerate(Cubes_4x4):
        #FFT_Power_spektra = get_power_spektra(Cubus)
        #powerspektrum, rgridding = get_power_spek_grid(FFT_Power_spektra)
        #powerspektrum = np.asarray(powerspektrum)
        #waterfall_plot(powerspektrum, rgridding, datei,"Sub_Sample_"+str(index)+"_" )
    
 
 
 
###############################################################################################
 
 
 

    
    
    
    
    
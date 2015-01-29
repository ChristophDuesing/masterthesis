# Program by Christoph Duesing
# This algorithm producess a the Powerspectrum of a random field 
import numpy as np
import pyfits as pf
import copy
from ring_area_calc_many import get_circle_area_for_pix
    
    
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
        max_y, min_y = (y-a/2.), (a/2.)
        fit_new = fit[:,min_y:max_y,:]
        print "y was reshaped" 
    elif x > y:
        a = x - y  
        max_x, min_x = (x-a/2.), (a/2.)
        fit_new = fit[:,:,min_x:max_x]# [z,y,x]!!    
        print "x was reshaped"
    else:
        fit_new = fit
    return fit_new

    
def get_power_spektra(datacube):
    power_spektra = [] 
    for i in datacube:
        fft_spek = np.fft.fft2(i, axes=(-2, -1))
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
        #print elem_i
        for elem_j in j:
            #print elem_j
            datacube = cube[:,y*elem_i/4:y*(elem_i+1)/4 ,x*elem_j/4:x*(elem_j+1)/4 ]    
            #print datacube
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
    #pylab.savefig('/export/data1/cduesing/Figures/Problems/TEST/Test_Picture.png')   
    plt.show()

    
def get_radius_gridding(length_of_frame, radius_spaceing):
    x = length_of_frame    
    r_min, r_max, a0,b0 =  2*np.sqrt(2)+1, x/2.,  x/2.,x/2. 
    
    #r_grid = np.logspace(np.log10(r_min), np.log10(r_max),radius_spaceing)       
    r_grid = np.linspace(r_min, r_max, radius_spaceing)

    r_grid_mean = []
    r_grid_mean = 0.5 * (r_grid[1:] + r_grid[:-1])
    #print r_grid_mean
    r_grid_mean = np.insert(r_grid_mean,0,np.sqrt(2) )
    #r_grid_mean = np.insert(r_grid_mean,0,np.log10(2*np.sqrt(2)+1) )
    #print r_grid_mean
    
    return r_grid, r_grid_mean, a0,b0 
    
    
def get_power_spek_grid(power_spektra_list, radius_spaceing):

    import pylab as pl
    import copy
    import matplotlib.pyplot as plt
    
    means_of_cicular_grids_spek_all = []
    
    x = len(power_spektra_list[0])
    
    r_grid, r_grid_mean, a0,b0 = get_radius_gridding(x,radius_spaceing)

    Area_ringshape = []

    px, py = np.meshgrid(
        np.arange(x/2.) + 0.5,
        np.arange(x/2.) + 0.5,
        )

    
    for idx, r in enumerate(r_grid[:-1]): 
        if idx == 0:        
            Area_ringshape.append(get_circle_area_for_pix((px, py),r_grid_mean[idx])[0])

        else :
            A_up, intersec_up, intersected_up = get_circle_area_for_pix((px, py), r_grid_mean[idx+1])
            A_low, intersec_low, intersected_low = get_circle_area_for_pix((px, py), r_grid_mean[idx])

            delta_A = A_up - A_low
        
            Area_ringshape.append(abs(delta_A))
            '''
            plt.close()
            plt.scatter(px, py, c=A_up)
            plt.show()
        
            plt.close()
            plt.scatter(px, py, c=A_low)
            plt.show()

            plt.close()
            plt.scatter(px, py, c=delta_A)
            plt.show()
        '''
    
            
    for z, i in enumerate(power_spektra_list):
        means_of_cicular_grids_spek = []
        # Warning first y then x and y goes from up to down and x from left to right 
        lo = i[:x/2,  :x/2.]
        lo = np.fliplr(np.flipud(lo))
        
        lu = i[x/2.: , :x/2.]
        lu = np.fliplr(lu)

        ro = i[ :x/2.,  x/2.:]
        ro = np.flipud(ro)
        
        ru = i[ x/2.:, x/2.:]
        ru = np.flipud(np.flipud(ru))
        
        for f, Area in enumerate(Area_ringshape):
            means_of_cicular_grids_spek.append(np.sum(np.multiply(ro,Area))/np.sum(Area)+ np.sum(ru*Area)/np.sum(Area) +np.sum(lo*Area)/np.sum(Area) + np.sum(lu*Area)/np.sum(Area))
        
        means_of_cicular_grids_spek_all.append(means_of_cicular_grids_spek)

    return means_of_cicular_grids_spek_all, r_grid_mean, Area_ringshape

    
#def get_masks_saved(masks, rgrid_mean, word):
    
    #import pylab as pl
    #for mask,rgrid in zip(masks,rgrid_mean):
        #pl.close()          
        #pl.imshow(mask,interpolation='nearest', cmap='cubehelix')
        #pl.savefig('/export/data1/cduesing/Figures/Problems/TEST/Masks_'+str(rgrid)+word+'.png', dpi=300)   
        ###plt.show()
        #pl.clf() 

    
def waterfall_plot(f, r_grid_mean, datei, attachment):
    import pylab as pl
    
    pl.close()
    pl.figure()
    pl.rc('text', usetex=True)
    pl.rc('font', family='serif')
    pl.xlabel(r" $\log_{10}$ Radius $|k|$")
    pl.ylabel(r"Velocities     $\displaystyle\frac{\text{km}}{\text{s}}$")
    #pl.imshow(f, interpolation='nearest',  cmap='cubehelix')
    #pl.imshow(np.log10(f), extent=(0.1, r_grid_mean[-1], 311, -311),interpolation='nearest', cmap='cubehelix')
    pl.imshow(np.log10(f), extent=(np.log10(r_grid_mean[0]), np.log10(r_grid_mean[-1]), -311, 311),origin='lower',interpolation='nearest', cmap='cubehelix')#TODO Watch out if is valid to just start from [1]
  
    cbar = pl.colorbar()
    cbar.set_label(r'$\log_{10}$ Power $P_k$', labelpad=-40, y=0.5)
    pl.gca().set_aspect('auto')
    #pl.savefig('/export/data1/cduesing/Figures/Problems/TEST/Waterfallplot_'+attachment+datei.split('.')[0]+'.png', dpi=300)
    pl.show()
    
 


    
if __name__ == '__main__':

    datei = ("CAR_NorSky30-50_80.fits")
    
    cube = open_fits(datei)
    cube_header = cube.header
    cube_data = cube.data
    get_shape(cube_data)

    
    print "intro"
    
    data = cube_data.copy()
    #Cubes_4x4 = get_sub_cube(data, datei, cube_header)
    
    print "FFT"
    FFT_Power_spektra = get_power_spektra(data)
    
    print "powerspektrum"
    powerspektrum, rgrid_mean, masks = get_power_spek_grid(FFT_Power_spektra, 128)
    #powerspektrum = np.asarray(powerspektrum)
    
    print "waterfall_plot"
    waterfall_plot(powerspektrum, rgrid_mean, datei,'_')
    print  rgrid_mean
    
    #def Wert_durch_d(plane):
import pylab as pl
px, py = np.meshgrid(
    np.arange(-185,186),
    np.arange(-185,186),
    )
z = 185.
d = [] 

d = np.sqrt(px**2+py**2)
tupel = list()
for i,j in zip(plane,d):    
    tupel.extend(zip(np.log10(j),np.log10(i)))
tupel2=zip(*tupel)    

pl.close()
#pl.rc('text', usetex=True)
#pl.rc('font', family='serif')
pl.xlabel(r"$\log_{10}(|d|=\sqrt{p_x^2+p_y^2})$ Distance to center ")
pl.ylabel(r"$\log_{10}$ Value of Pixel")

#pl.scatter(*zip(*tupel))
pl.plot(*zip(*tupel),color='red', linestyle='.', marker='.')
#pl.loglog(tupel,"xr")

pl.show()
    
    #f = 1
    #for index, Cubus in enumerate(Cubes_4x4):
        #print Cubus.shape
        #Cubus = quadratic_shape(Cubus)
        #Cubus = get_odd_shape(Cubus)
        #FFT_Power_spektra_Cubus = get_power_spektra(Cubus)
        #Cubus_powerspektrum, Cubus_rgridding, Cubus_masks_ringshape = get_power_spek_grid(FFT_Power_spektra_Cubus,"Subcube_1_"+str(f), 16.)
        #Cubus_powerspektrum = np.asarray(Cubus_powerspektrum)
        #waterfall_plot(Cubus_powerspektrum, Cubus_rgridding, datei,"Sub_Sample_"+str(index)+"_" )
        #if f==1:
            #get_masks_saved(Cubus_masks_ringshape,Cubus_rgridding,"subcube_")
        #f+=1
    
 
 
 
 
 

    
    
    
    
    
########################################################################
########################################################################        
########################################################################
# Program by Christoph Duesing
# This algorithm producess a the Powerspectrum of a random field 
import numpy as np
import copy
import matplotlib.pyplot as plt
import pyfits as pf
import scipy.optimize
import scipy.signal 



###############################################################
###############################################################
###############################################################

#REAL DATA ANALYSIS

###############################################################
###############################################################
###############################################################


hdu = pf.open('NorSky70-90_40.fits')
primhdu = hdu[0]
header = primhdu.header
data = primhdu.data


print data.shape 

spectrum = data[407,:,:]# [z,y,x]!!
spectrum_whole = spectrum[68:305,:]
spectrum_void = spectrum[150:250,:100]
spectrum_cloud = spectrum[200:300,135:235]
spectrum_cloud2 = spectrum[0:100,0:100]

plt.close()
fig = plt.figure(figsize=(15, 5))
ax1 = fig.add_subplot(1,4,1)
ax2 = fig.add_subplot(1,4,2)
ax3 = fig.add_subplot(1,4,3)
ax4 = fig.add_subplot(1,4,4)


imdict = dict(interpolation='nearest', origin='lower')
datas = [spectrum_whole, spectrum_void, spectrum_cloud, spectrum_cloud2]
axes = [ax1, ax2, ax3, ax4]

for d, ax in zip(datas, axes):
    ax.set_aspect('equal')
    ax.imshow(d, **imdict)
    
plt.show()


fft_spek = np.fft.fft2(spectrum_whole)
fft_spek = np.fft.fftshift(fft_spek)

fft_spek_void = np.fft.fft2(spectrum_void)
fft_spek_void = np.fft.fftshift(fft_spek_void)

fft_spek_cloud = np.fft.fft2(spectrum_cloud)
fft_spek_cloud = np.fft.fftshift(fft_spek_cloud)

fft_spek_cloud2 = np.fft.fft2(spectrum_cloud2)
fft_spek_cloud2 = np.fft.fftshift(fft_spek_cloud2)

#freq = np.fft.fftfreq(spectrum.shape[1])
pow_spek = (np.abs(fft_spek)**2)/len(fft_spek)

#freq_void = np.fft.fftfreq(spectrum_void.shape[1])
pow_spek_void = (np.abs(fft_spek_void)**2)/len(fft_spek_void)

#freq_cloud = np.fft.fftfreq(spectrum_cloud.shape[1])
pow_spek_cloud = (np.abs(fft_spek_cloud)**2)/len(fft_spek_cloud)

pow_spek_cloud2 = (np.abs(fft_spek_cloud2)**2)/len(fft_spek_cloud2)


plt.close()
fig = plt.figure(figsize=(15, 5))
ax1 = fig.add_subplot(1,4,1)
ax2 = fig.add_subplot(1,4,2)
ax3 = fig.add_subplot(1,4,3)
ax4 = fig.add_subplot(1,4,4)
imdict = dict(interpolation='nearest', origin='lower')
datas = [pow_spek, pow_spek_void, pow_spek_cloud, pow_spek_cloud2]
axes = [ax1, ax2, ax3, ax4]

for d, ax in zip(datas, axes):
    ax.set_aspect('equal')
    ax.imshow(d, **imdict)

plt.show()


###############################################################

#Circular masks:
np.max(spectrum)
np.max(spectrum_void)
np.max(spectrum_cloud)
np.max(spectrum_cloud2)


a0,b0 = len(spectrum)/2.,len(spectrum[1])/2.
c0,d0 = len(spectrum_void)/2.,len(spectrum_void)/2.
e0,f0 = len(spectrum_cloud)/2.,len(spectrum_cloud)/2.
g0,h0 = len(spectrum_cloud2)/2.,len(spectrum_cloud2)/2.

#spectrum = data[300:600,256,:]# [z,y,x]!!
#spectrum_void = spectrum[:100,:100]
#spectrum_cloud = spectrum[80:130,200:235]



a, b = np.meshgrid(np.linspace(0,len(spectrum[1]), len(spectrum[1])), np.linspace(0,len(spectrum[0]), len(spectrum[0])))
c, d = np.meshgrid(np.linspace(0,len(spectrum_void), len(spectrum_void)), np.linspace(0,len(spectrum_void), len(spectrum_void)))
e, f = np.meshgrid(np.linspace(0,len(spectrum_cloud), len(spectrum_cloud)), np.linspace(0,len(spectrum_cloud), len(spectrum_cloud)))
g, h = np.meshgrid(np.linspace(0,len(spectrum_cloud2), len(spectrum_cloud2)), np.linspace(0,len(spectrum_cloud2), len(spectrum_cloud2)))

r_min =  1 
r_max = 50
r_grid = np.linspace(r_min,r_max,256)

# From 256 grid poits maximum stays constant. Also for more grid points there are no variations !! 

#aperture = np.zeros_like(x)


aperture = np.zeros_like(a)
aperture_void = np.zeros_like(c)
aperture_cloud = np.zeros_like(e)
aperture_cloud2 = np.zeros_like(g)


radius_mask_func = lambda a, b, a0, b0, r: (a-a0)**2 + (b-b0)**2 <= r**2

means_of_cicular_grids_spek = []
means_of_cicular_grids_spek_void = []
means_of_cicular_grids_spek_cloud = []
means_of_cicular_grids_spek_cloud2 = []


for idx, r in enumerate(r_grid[:-1]):  
    mask_ringshape = np.logical_and(
        radius_mask_func(a, b, a0, b0,  r_grid[idx+1]),
        np.logical_not(radius_mask_func(a, b, a0, b0,  r_grid[idx]))
        )
    aperture_sample  = copy.deepcopy(aperture)
    aperture_sample[mask_ringshape]=1

 
    pow_spek_grid = pow_spek * aperture_sample
    means_of_cicular_grids_spek.append(np.sum(pow_spek_grid)/np.count_nonzero(aperture_sample))
    print r
    print len(means_of_cicular_grids_spek)
    
    #######################################################################

    mask_ringshape = np.logical_and(
        radius_mask_func(c, d, c0, d0,  r_grid[idx+1]),
        np.logical_not(radius_mask_func(c, d, c0, d0,  r_grid[idx]))
        )
    aperture_void_sample  = copy.deepcopy(aperture_void)
    aperture_void_sample[mask_ringshape]=1

 
    pow_spek_grid_void = pow_spek_void * aperture_void_sample
    means_of_cicular_grids_spek_void.append(np.sum(pow_spek_grid_void)/np.count_nonzero(aperture_void_sample))
    print r
    print len(means_of_cicular_grids_spek_void)
    
    #######################################################################
    mask_ringshape = np.logical_and(
        radius_mask_func(e, f, e0, f0,  r_grid[idx+1]),
        np.logical_not(radius_mask_func(e, f, e0, f0,  r_grid[idx]))
        )
    aperture_cloud_sample = copy.deepcopy(aperture_cloud)
    aperture_cloud_sample[mask_ringshape]=1

 
    pow_spek_grid_cloud = pow_spek_cloud * aperture_cloud_sample
    
    
    means_of_cicular_grids_spek_cloud.append(np.sum(pow_spek_grid_cloud)/np.count_nonzero(aperture_cloud_sample))
    
    
    
    #######################################################################    
    mask_ringshape = np.logical_and(
        radius_mask_func(g, h, g0, h0,  r_grid[idx+1]),
        np.logical_not(radius_mask_func(g, h, g0, h0,  r_grid[idx]))
        )
    aperture_cloud2_sample = copy.deepcopy(aperture_cloud)
    aperture_cloud2_sample[mask_ringshape]=1

 
    pow_spek_grid_cloud2 = pow_spek_cloud2 * aperture_cloud2_sample
    means_of_cicular_grids_spek_cloud2.append(np.sum(pow_spek_grid_cloud2)/np.count_nonzero(aperture_cloud2_sample))
    
    
       
r_grid_mean = 0.5 * (r_grid[1:] + r_grid[:-1])



plt.close()   
plt.plot(r_grid_mean,np.log10(means_of_cicular_grids_spek),'x')
plt.plot(r_grid_mean,np.log10(means_of_cicular_grids_spek_void),':g')
plt.plot(r_grid_mean,np.log10(means_of_cicular_grids_spek_cloud),':r')
plt.plot(r_grid_mean,np.log10(means_of_cicular_grids_spek_cloud2),':b')
plt.show()    
    























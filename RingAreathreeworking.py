
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


r_min_start, r_max_start, a0,b0 =  1, x/2.,  x/2.,x/2. 
#r_grid = np.logspace(np.log10(r_min), np.log10(r_max),16)   
r_grid = np.linspace(r_min_start, r_max_start,64)
r_grid_mean = 0.5 * (r_grid[1:] + r_grid[:-1])



dists_sq = (a-a0)**2 + (b-b0)**2

radius_mask_func = lambda dists_sq, r: dists_sq <= r**2

masks_ringshape = []
Areas_True = []
for idx, r in enumerate(r_grid[:-1]):
    
    check_pix =[]
    Area_True = np.pi*(r_grid[idx+1]**2-r_grid[idx]**2)
    
    masks_ringshape.append(np.logical_and(
        radius_mask_func(dists_sq,  r_grid[idx+1]),
        np.logical_not(radius_mask_func(dists_sq,  r_grid[idx]))
        ))
    
    Areas_True.append(Area_True)
    

    
    ##########################################################################
    subsub_2 = subsub.copy()
    def Antiderivative(a_min, a_max ,b):
    #b is the radius 
        A = ((a_max/2 * np.sqrt(b**2 - a_max**2) + b**2/2 * np.arcsin(a_max/b)) - (a_min/2 * np.sqrt(b**2-a_min**2) + b**2/2 * np.arcsin(a_min/b)))
        return A
    
    
    def RadialGleichung(Radius, Koordinate):
        return np.sqrt(abs(Radius**2 - Koordinate**2))
    
    def get_intercections(a_low, a_up, a_min, a_max, b_left, b_right, b_min, b_max ):
        
        Schnittpunkt_x_low, Schnittpunkt_x_up, Schnittpunkt_y_left, Schnittpunkt_y_right = 0, 0, 0, 0
        
        if a_low > a_min and a_low < a_max:
            Schnittpunkt_x_low +=1
            
        if a_up > a_min and a_up < a_max:
            Schnittpunkt_x_up +=1
        
        if b_left > b_min and b_left < b_max:
            Schnittpunkt_y_left +=1
        
        if y_right_max > y_min and y_right_max < y_max: 
            Schnittpunkt_y_right +=1
    
        return Schnittpunkt_x_low, Schnittpunkt_x_up, Schnittpunkt_y_left, Schnittpunkt_y_right
    
        
       
    
    r_max = r_grid[35]
    print r_max
    r_min = r_grid[34]
    b= r_max-r_min
    print r_min
    
    Check_Range = np.logical_and(
    radius_mask_func(dists_sq,  r_max + np.sqrt(2)),
    np.logical_not(radius_mask_func(dists_sq,  r_min-np.sqrt(2))
    ))
    
    Ring_true = np.logical_and(
    radius_mask_func(dists_sq,  r_max ),
    np.logical_not(radius_mask_func(dists_sq,  r_min)
    ))
    
    subsub_2[Check_Range] =1
    subsub_2[Ring_true] = 9999
    xval, yval = np.where(subsub_2==1) 
    coordinates = zip(xval, yval)
    
    Whole_pix_hit = []
    

    Cicular_area = []
    for i in coordinates:
        
        x,y = i  
        x,y =x-49.5,y-49.5
        #print x,y
        
        #First Quadrant
        
        if x>=0 and y>=0:
            #print
            x_min,x_max,y_min,y_max = x-0.5,x+0.5,y-0.5,y+0.5
            #print x_min, "x_max", "y_min", "y_max"
            #print x_min, x_max, y_min, y_max
            #print
            
            #Check for interceptionpoints between cicle and pixel
            
            x_lower_max = RadialGleichung(r_max, y_min)
            x_upper_max = RadialGleichung(r_max, y_max)
            y_left_max = RadialGleichung(r_max, x_min)
            y_right_max = RadialGleichung(r_max, x_max)
            
            x_lower_min = RadialGleichung(r_min, y_min)
            x_upper_min = RadialGleichung(r_min, y_max)
            y_left_min = RadialGleichung(r_min, x_min)
            y_right_min = RadialGleichung(r_min, x_max)
            
            #print  "x_lower_max", "x_upper_max", "y_left_max", "y_right_max"
            #print  x_lower_max, x_upper_max, y_left_max, y_right_max
            #print 
            
            '''General dis mission'''
            
            max_Schnittpunkt_x_low, max_Schnittpunkt_x_up, max_Schnittpunkt_y_left, max_Schnittpunkt_y_right = get_intercections(x_lower_max, x_upper_max, x_min, x_max, y_left_max, y_right_max, y_min, y_max )
            min_Schnittpunkt_x_low, min_Schnittpunkt_x_up, min_Schnittpunkt_y_left, min_Schnittpunkt_y_right = get_intercections(x_lower_min, x_upper_min, x_min, x_max, y_left_min, y_right_min, y_min, y_max )
            
            print min_Schnittpunkt_x_low, min_Schnittpunkt_x_up, min_Schnittpunkt_y_left, min_Schnittpunkt_y_right
            
            #CASE 1 WHOLE PIXEL FILLED 
            
            if max_Schnittpunkt_x_low == 0 and max_Schnittpunkt_x_up == 0 and max_Schnittpunkt_y_left == 0 and max_Schnittpunkt_y_right == 0 and \
               min_Schnittpunkt_x_low ==0 and min_Schnittpunkt_x_up ==0 and min_Schnittpunkt_y_left == 0 and min_Schnittpunkt_y_right ==0 :
                if np.sqrt(x_max**2 +y_max**2) <= r_max and np.sqrt(x_min**2 +y_min**2) >= r_min:
                     
                    i_new = (x,y,1,"CASE1")     
                    Cicular_area.append(i_new)
            
            #CASE 2 BOTTOM LEFT CORNER 
            
            if max_Schnittpunkt_x_low == 1 and max_Schnittpunkt_y_left == 1:
                I_max = Antiderivative(x_min, x_lower_max ,r_max)
                Kasten_max = (x_lower_max-x_min)*y_min
                Stamm_max = (I_max - Kasten_max)
                I_min = 0
                Stamm_min = 0
                Kasten_min = 0

                if min_Schnittpunkt_x_low == 1 and min_Schnittpunkt_y_left == 1:    
                    I_min = Antiderivative(x_min, x_lower_min ,r_max)
                    Kasten_min = (x_lower_min-x_min)*y_min
                    Stamm_min = (I_min - Kasten_min)
                

                Stamm = Stamm_max-Stamm_min
                
                i_new = (x,y,Stamm,"CASE2")     
                Cicular_area.append(i_new)
                
            #CASE 3 TOP RIGHT CORNER 
            
            if min_Schnittpunkt_x_up == 1 and min_Schnittpunkt_y_right == 1:
                
                I_min = Antiderivative(x_upper_min, x_max ,r_min)
                Kasten_min = (x_max-x_upper_min)*y_min
                Stamm_min = (I_min - Kasten_min)
               
                I_max = 0
                Stamm_max = 0
                Kasten_max = 0
                Stamm_kleine_Ecke = 0
                
                if max_Schnittpunkt_x_up == 1 and max_Schnittpunkt_y_right == 1:
                    I_max = Antiderivative(x_upper_max, x_max ,r_max)
                    Kasten_max = (x_max-x_upper_min)*y_min
                    Stamm_max = (I_max - Kasten_min)
                    Stamm_kleine_Ecke = 1-(x_upper_max-x_min)-Stamm_min
                
                Stamm = 1-(x_upper_min-x_min)-Stamm_min-Stamm_kleine_Ecke
                i_new = (x,y,Stamm,"CASE3")
                Cicular_area.append(i_new)                

            #CASE 4 BOTH Y BOARDERS ARE INTERCEPTED
            
            if max_Schnittpunkt_y_left == 1 and max_Schnittpunkt_y_right == 1:
                I_max = Antiderivative(x_min, x_max ,r_max)
                Kasten_min = (x_max-x_min)*y_min
                Stamm_min = (I_max - Kasten_min)
                Stamm = 1-Stamm_min
                i_new = (x,y,Stamm,"CASE4")
                Cicular_area.append(i_new)      
                
                
            #CASE 5 BOTH X BOARDERS ARE INTERCEPTED
            
            if max_Schnittpunkt_x_low == 1 and max_Schnittpunkt_x_up == 1:
                Kasten_Im_Kreis = (x_upper_max - x_min) 
                print "Kasten_Im_Kreis",Kasten_Im_Kreis
                I_max = Antiderivative(x_upper_max, x_lower_max ,r_max)
                print "I_max",I_max
                Stamm_max = (x_lower_max - x_upper_max ) * y_min
                print "Stamm_max", Stamm_max
                Stamm = 1 - (Kasten_Im_Kreis + (I_max -Stamm_max))
                #print "Stamm = 1 - (Kasten_Im_Kreis + (I_max -Stamm_max))", Stamm
                more = 0
                if min_Schnittpunkt_x_low == 1 and min_Schnittpunkt_x_up == 1:
                    more = 1
                if min_Schnittpunkt_x_low == 1 and min_Schnittpunkt_y_left == 1:    
                    more = 2
                i_new = (x,y,Stamm,"CASE5", more)
                
                
                print "max_Schnittpunkt_x_low", "max_Schnittpunkt_x_up", "max_Schnittpunkt_y_left", "max_Schnittpunkt_y_right"
                print max_Schnittpunkt_x_low, max_Schnittpunkt_x_up, max_Schnittpunkt_y_left, max_Schnittpunkt_y_right
                print "min_Schnittpunkt_x_low", "min_Schnittpunkt_x_up", "min_Schnittpunkt_y_left", "min_Schnittpunkt_y_right"
                print min_Schnittpunkt_x_low, min_Schnittpunkt_x_up, min_Schnittpunkt_y_left, min_Schnittpunkt_y_right
                
                Cicular_area.append(i_new)  

                
pl.figure()
pl.imshow(subsub_2,interpolation='nearest',)
pl.show()

        
        ##Integrate Area 
        #def integrand_rmax(x,r_max):
            #return 2*np.sqrt(r_max**2-x**2)
        #def integrand_rmin(x,r_min):
            #return 2*np.sqrt(r_min**2-x**2)
            
        ##if x_lower_min >= x_min and x_lower_min <= x_max:
            ##x_min = x_lower_min
        
        ##if x_lower_max <= x_max and x_lower_max >= x_min:
            ##x_max = x_lower_max
        
        #'''CASE 1 just r_max hits the pixel'''
        #if x_min < x_lower_min and x_min < x_upper_min:
            #if x_lower_max <= x_max and x_lower_max >= x_min:
                #x_max = x_lower_max
            #I = quad(integrand_rmax, x_min, x_max,args=(r_max))
        
        #'''CASE 2 just r_min hits the pixel'''
        #if x_max > x_lower_max and x_max > x_upper_min:   
            #if x_lower_min >= x_min and x_lower_min <= x_max:
                #x_min = x_lower_min
            #I = quad(integrand_rmin, x_min, x_max,args=(r_min))
        #'''CASE 3 both r_min and r_max hit the pixel'''
        #if x_lower_min > x_min and x_lower_max < x_max:
            #I_1 = quad(integrand_rmax, x_min, x_max,args=(r_max))
            #I_2 = quad(integrand_rmin, x_min, x_max,args=(r_min))
            #I = I_1[0]-I_2[0]
        #print I
        

        
        
Areas_Estimated =[]   
means_of_cicular_grids_spek = [ ]   
for mask_ringshape in masks_ringshape:
    means_of_cicular_grids_spek.append(np.mean(subsub[mask_ringshape])) 
    Areas_Estimated.append(np.count_nonzero(mask_ringshape))




t = zip(Areas_True,Areas_Estimated)

diff = np.subtract(Areas_True,Areas_Estimated)
diff
div = np.divide(Areas_True,Areas_Estimated)
div






    

    
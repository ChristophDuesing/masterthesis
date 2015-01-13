
# Program by Christoph Duesing
# This algorithm producess a the Powerspectrum of a random field 
import numpy as np
import pyfits as pf
import copy
import pylab as pl
import matplotlib.pyplot as plt



z = np.linspace(0.,99.,100.)

x = len(z)

a, b = np.meshgrid(np.linspace(0,x-1,x)+0.5 , np.linspace(0,x-1,x)+0.5)
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
    
    
    def Antiderivative(a_min, a_max, Radius):
        A = ((a_max/2 * np.sqrt(Radius**2 - a_max**2) + Radius**2/2 * np.arcsin(a_max/Radius)) - (a_min/2 * np.sqrt(Radius**2-a_min**2) + Radius**2/2 * np.arcsin(a_min/Radius)))
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
    
    
    #def CASE1_whole_pixel(Schnittpunkte, x, y, borders, r_min, r_max ):
        #if all(S == 0 for S in Schnittpunkte) == True:
            #if np.sqrt(borders[1]**2 +borders[3]**2) <= r_max and np.sqrt(borders[0]**2 +borders[2]**2) >= r_min:
                #print "CASE1"
                #return 1
            
    def CASE2_bottom_left_corner(borders, lower_x, r):
        I = Antiderivative(borders[0], lower_x ,r)
        Kasten = (lower_x-borders[0])*borders[2]
        Circ_Area = (I - Kasten)
        print "CASE2"
        return Circ_Area     
    
    
    def CASE3_upper_right_corner(borders, upper_x, r):    
        I = Antiderivative(upper_x, borders[1] ,r)
        Kasten = (borders[1]-upper_x)*borders[2]
        Circ_Area = (I - Kasten)
        print "CASE3"
        return Circ_Area         
    
    
    def CASE4_both_y_borders(borders, y_left, y_right, r):
        I = Antiderivative(borders[0], borders[1] ,r)
        Kasten = (borders[1] - borders[0]) * borders[2]
        #Circ_Area = 1 - (I - Kasten)
        Circ_Area = (I - Kasten)
        print "CASE4"
        return Circ_Area           

       
    
    
    
    r_max = r_grid[35]
    print r_max
    r_min = r_grid[34]
    b= r_max-r_min
    print r_min
    
    Check_Range = np.logical_and(
    radius_mask_func(dists_sq,  r_max + np.sqrt(2)),
    np.logical_not(radius_mask_func(dists_sq,  r_min-np.sqrt(2))
    ))
    

    subsub_2[Check_Range] =1
    
    Ring_true = np.logical_and(
    radius_mask_func(dists_sq,  r_max ),
    np.logical_not(radius_mask_func(dists_sq,  r_min)
    ))
    subsub_2[Ring_true] = 9999
    
    
    
    
    xval, yval = np.where(subsub_2==1) 
    coordinates = zip(xval, yval)
    
    Whole_pix_hit = []
    

    Cicular_area = []
    for i in coordinates:
        
        x,y = i  
        x,y =x-49.5,y-49.5
        #print x,y
        
        '''
            First Quadrant is used only 
        '''
        
        if x>=0 and y>=0:
            print x,y 
            x_min,x_max,y_min,y_max = x-0.5,x+0.5,y-0.5,y+0.5
            
            borders = (x_min,x_max,y_min,y_max)
            #Check for interceptionpoints between cicle and pixel
            
            x_lower_max = RadialGleichung(r_max, y_min)
            x_upper_max = RadialGleichung(r_max, y_max)
            y_left_max = RadialGleichung(r_max, x_min)
            y_right_max = RadialGleichung(r_max, x_max)
            
            x_lower_min = RadialGleichung(r_min, y_min)
            x_upper_min = RadialGleichung(r_min, y_max)
            y_left_min = RadialGleichung(r_min, x_min)
            y_right_min = RadialGleichung(r_min, x_max)
            
            
            max_Schnittpunkt_x_low, max_Schnittpunkt_x_up, max_Schnittpunkt_y_left, max_Schnittpunkt_y_right = get_intercections(x_lower_max, x_upper_max, x_min, x_max, y_left_max, y_right_max, y_min, y_max )
            min_Schnittpunkt_x_low, min_Schnittpunkt_x_up, min_Schnittpunkt_y_left, min_Schnittpunkt_y_right = get_intercections(x_lower_min, x_upper_min, x_min, x_max, y_left_min, y_right_min, y_min, y_max )
            
            Schnittpunkte_max = (max_Schnittpunkt_x_low, max_Schnittpunkt_x_up, max_Schnittpunkt_y_left, max_Schnittpunkt_y_right)
            Schnittpunkte_min = (min_Schnittpunkt_x_low, min_Schnittpunkt_x_up, min_Schnittpunkt_y_left, min_Schnittpunkt_y_right)
            
            
            ''' 
                First Check the outer radius and calculate the Integral
            
                If Done continue with the lower radius 
                
                and calculate the Integral. Then Substract the lower from the Top
            
            '''
            
            Area = 0
            
            #CASE 1 WHOLE PIXEL FILLED 
            
            if all(S == 0 for S in Schnittpunkte_max+Schnittpunkte_min) == True:
                if np.sqrt(borders[1]**2 +borders[3]**2) <= r_max and np.sqrt(borders[0]**2 +borders[2]**2) >= r_min:
                    print "CASE1"
                    Area =1
                

            #CASE 2 BOTTOM LEFT CORNER 
            if Schnittpunkte_max[0]== 1 and Schnittpunkte_max[2] == 1:
                Area_r_max = CASE2_bottom_left_corner(borders, x_lower_max, r_max)
            if Schnittpunkte_min[0]== 1 and Schnittpunkte_min[2] == 1:
                Area_r_min = CASE2_bottom_left_corner(borders, x_lower_min, r_min)

                
            #CASE 3 TOP RIGHT CORNER 
         
            if Schnittpunkte_max[1] == 1 and Schnittpunkte_max[3] == 1:
                Area_r_max = CASE3_upper_right_corner(borders, x_upper_max, r_max)   
            if Schnittpunkte_min[1] == 1 and Schnittpunkte_min[3] == 1:
                Area_r_min = CASE3_upper_right_corner(borders, x_upper_min, r_min)  
       
       
            #CASE 4 BOTH Y BOARDERS ARE INTERCEPTED

            if Schnittpunkte_max[2] == 1 and Schnittpunkte_max[3] == 1:
                Area_r_max = CASE4_both_y_borders(borders, y_left_max, y_right_max, r_max) 
            if Schnittpunkte_min[2] == 1 and Schnittpunkte_min[3] == 1:    
                Area_r_min = CASE4_both_y_borders(borders, y_left_min, y_right_min, r_min)
                
                
            #CASE 5 BOTH X BOARDERS ARE INTERCEPTED
            
            if Schnittpunkte_max[0] == 1 and Schnittpunkte_max[1] == 1:
                Kasten_Im_Kreis = (x_upper_max - x_min) 
                
                I_max = Antiderivative(x_upper_max, x_lower_max ,r_max)
               
                Stamm_max = (x_lower_max - x_upper_max ) * y_min
               
                Stamm = 1 - (Kasten_Im_Kreis + (I_max -Stamm_max))
                #print "Stamm = 1 - (Kasten_Im_Kreis + (I_max -Stamm_max))", Stamm
                more = 0
                if min_Schnittpunkt_x_low == 1 and min_Schnittpunkt_x_up == 1:
                    more = 1
                if min_Schnittpunkt_x_low == 1 and min_Schnittpunkt_y_left == 1:    
                    more = 2
                i_new = (x,y,Stamm,"CASE5", more)
                
               
                
                Cicular_area.append(i_new)  
            #in case if Area_r_max == 0 take the absolut  
            Area = abs(Area_r_max - Area_r_min)
            print Area 
                
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






    

    
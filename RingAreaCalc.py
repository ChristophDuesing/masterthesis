# Program by Christoph Duesing
# This algorithm producess a the Powerspectrum of a random field 
import numpy as np
import pyfits as pf
import copy
import pylab as pl
import matplotlib.pyplot as plt

def get_pixels(r_max, r_min, dists_sq):
    Check_Range = np.logical_and(
    radius_mask_func(dists_sq,  r_max + np.sqrt(2)),
    np.logical_not(radius_mask_func(dists_sq,  r_min-np.sqrt(2))
    ))
    
    Test_Plane_copy[Check_Range] =1
    xval, yval = np.where(Test_Plane_copy==1) 
    return zip(xval, yval)    
    

def Antiderivative(a_min, a_max, Radius):
    A = ((a_max/2 * np.sqrt(Radius**2 - a_max**2) + Radius**2/2 * np.arcsin(a_max/Radius)) - (a_min/2 * np.sqrt(Radius**2-a_min**2) + Radius**2/2 * np.arcsin(a_min/Radius))) 
    return A


def RadialEqu(Radius, coord):
    return np.sqrt(abs(Radius**2 - coord**2))


def get_intercections(borders,a_low, a_up, b_left, b_right,):
    
    intersec_x_low, intersec_x_up, intersec_y_left, intersec_y_right = 0, 0, 0, 0
    
    if a_low > borders[0] and a_low < borders[1]:
        intersec_x_low +=1
        
    if a_up > borders[0] and a_up < borders[1]:
        intersec_x_up +=1
    
    if b_left > borders[2] and b_left < borders[3]:
        intersec_y_left +=1
    
    if b_right > borders[2] and b_right < borders[3]: 
        intersec_y_right +=1

    return intersec_x_low, intersec_x_up, intersec_y_left, intersec_y_right


        
def CASE2_bottom_left_corner(borders, lower_x, r):
    I = Antiderivative(borders[0], lower_x ,r)
    Box = (lower_x-borders[0])*borders[2]
    Circ_Area = (I - Box)
    print "CASE2"
    return Circ_Area     


def CASE3_upper_right_corner(borders, upper_x, r):    
    I = Antiderivative(upper_x, borders[1] ,r)
    Box = (borders[1]-upper_x)*borders[2]
    Circ_Area = (I - Box)
    print "CASE3"
    return Circ_Area         


def CASE4_both_y_borders(borders, y_left, y_right, r):
    I = Antiderivative(borders[0], borders[1] ,r)
    Box = (borders[1] - borders[0]) * borders[2]
    #Circ_Area = 1 - (I - Box)
    Circ_Area = (I - Box)
    print "CASE4"
    return Circ_Area    
    
    
def CASE5_both_x_borders(borders, x_lower, x_upper):
    Box_in_circ = (x_upper - borders[0]) 
    Box_under_int = (x_lower - x_upper) * borders[2]
    I = Antiderivative(x_upper, x_lower ,r_max)
    Circ_Area = (Box_in_circ + (I - Box_under_int))
    print "CASE5"
    return Circ_Area 
    
    
def get_circle_area_for_pix(coordinates, r):
            x_min,x_max,y_min,y_max = x-0.5,x+0.5,y-0.5,y+0.5
            borders = (x_min,x_max,y_min,y_max)
            
            #Check for interceptionpoints between cicle and axis
            x_lower = RadialEqu(r, y_min)
            x_upper = RadialEqu(r, y_max)
            y_left = RadialEqu(r, x_min)
            y_right = RadialEqu(r, x_max)
            
            #Check for actual interceptionpoints between cicle and pixel
            intersec_x_low, intersec_x_up, intersec_y_left, intersec_y_right = get_intercections(borders, x_lower, x_upper, y_left, y_right)
            intersec = (intersec_x_low, intersec_x_up, intersec_y_left, intersec_y_right)
            #print intersec
            #CASE 1 WHOLE PIXEL FILLED 
            if np.sqrt(x_min**2+y_max**2) <= r:
                Area = 1
            else:
                Area = 0

            #CASE 2 BOTTOM LEFT CORNER 
            if intersec[0]== 1 and intersec[2] == 1:
                Area = CASE2_bottom_left_corner(borders, x_lower, r)

            #CASE 3 TOP RIGHT CORNER 
            
            if intersec[1] == 1 and intersec[3] == 1:
                Area = CASE3_upper_right_corner(borders, x_upper, r)   
          
        
            #CASE 4 BOTH Y BORDERS ARE INTERCEPTED

            if intersec[2] == 1 and intersec[3] == 1:
                Area = CASE4_both_y_borders(borders, y_left, y_right, r) 
            
                
            #CASE 5 BOTH X BORDERS ARE INTERCEPTED
            
            if intersec[0] == 1 and intersec[1] == 1:
                Area = CASE5_both_x_borders(borders, x_lower, x_upper)
            

            return Area 
                
    
if __name__ == '__main__':

    z = np.linspace(0.,99.,100.)
    x = len(z)
    a, b = np.meshgrid(np.linspace(0,x-1,x)+0.5 , np.linspace(0,x-1,x)+0.5)

    Test_Plane = a*b
    Test_Plane_copy = Test_Plane.copy()

    r_min_start, r_max_start, a0,b0 =  1, x/2.,  x/2.,x/2. 
    r_grid = np.linspace(r_min_start, r_max_start,64)
    r_grid_mean = 0.5 * (r_grid[1:] + r_grid[:-1])
    
    dists_sq = (a-a0)**2 + (b-b0)**2

    radius_mask_func = lambda dists_sq, r: dists_sq <= r**2
    
    r_min, r_max = r_grid[34], r_grid[35]
    b = r_max-r_min

    coordinates = get_pixels(r_max, r_min, dists_sq)
            
    for i in coordinates:
        x,y = i  
        x,y =x-49.5,y-49.5
        #First Quadrant is used only 
        if x>=0 and y>=0:
            print x,y 
   
            Area = abs(get_circle_area_for_pix(i, r_max) - get_circle_area_for_pix(i, r_min))
            print Area
            
     
    pl.figure()
    pl.imshow(Test_Plane_copy,interpolation='nearest',)
    pl.show()

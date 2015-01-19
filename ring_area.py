#!/usr/bin/python
# -*- coding: utf-8 -*-

# Program by Christoph Duesing
# This algorithm producess a the Powerspectrum of a random field
import numpy as np
import pyfits as pf
import copy
import pylab as pl
import matplotlib.pyplot as plt

'''
Get the Integrated Area of the circle in cartesian coordinates
'''
def Antiderivative(a_min, a_max, Radius):
    aR_max = a_max / Radius
    aR_min = a_min / Radius
    A = (
        a_max * np.sqrt(1. - aR_max ** 2) +
        Radius * np.arcsin(aR_max) -
        a_min * np.sqrt(1. - aR_min ** 2) -
        Radius * np.arcsin(aR_min)
        ) * Radius / 2.
    return A


def RadialEqu(Radius, coord):
    return np.sqrt(abs(Radius ** 2 - coord ** 2))

'''
Find the intersection for a pixel with the circle.
'''
def get_intercections(borders, a_low, a_up, b_left, b_right):

    intersec_x_low = np.zeros(a_low.shape, dtype=np.bool)
    intersec_x_up = np.zeros(a_low.shape, dtype=np.bool)
    intersec_y_left = np.zeros(a_low.shape, dtype=np.bool)
    intersec_y_right = np.zeros(a_low.shape, dtype=np.bool)
    
    mask = (a_low >= borders[0]) & (a_low <= borders[1])
    intersec_x_low[mask] = True
    
    mask = (a_up >= borders[0]) & (a_up <= borders[1])
    intersec_x_up[mask] = True

    mask = (b_left >= borders[2]) & (b_left <= borders[3])
    intersec_y_left[mask] = True

    mask = (b_right >= borders[2]) & (b_right <= borders[3])
    intersec_y_right[mask] = True

    return intersec_x_low, intersec_x_up, intersec_y_left, intersec_y_right

    
def CASE2_bottom_left_corner(borders, lower_x, r):
    I = Antiderivative(borders[0], lower_x, r)
    Box = (lower_x - borders[0]) * borders[2]
    Circ_Area = I - Box
    print 'CASE2'
    return Circ_Area


def CASE3_upper_right_corner(borders, upper_x, r):
    Circ_box = (upper_x - borders[0])
    I = Antiderivative(upper_x, borders[1], r)
    Box = (borders[1] - upper_x) * borders[2]
    Circ_Area = (I + Circ_box )- Box
    print 'CASE3'
    return Circ_Area


def CASE4_both_y_borders(borders, y_left, y_right, r):
    I = Antiderivative(borders[0], borders[1], r)
    Box = (borders[1] - borders[0]) * borders[2]
    # Circ_Area = 1 - (I - Box)
    Circ_Area = I - Box
    print 'CASE4'
    return Circ_Area


def CASE5_both_x_borders(borders, x_lower, x_upper, r):
    Box_in_circ = x_upper - borders[0]
    Box_under_int = (x_lower - x_upper) * borders[2]
    I = Antiderivative(x_upper, x_lower, r)
    Circ_Area = Box_in_circ + I - Box_under_int
    print 'CASE5'
    return Circ_Area    
    
    
def get_circle_area_for_pix(coordinates, r):
    x, y = coordinates
    x_min, x_max, y_min, y_max = x - 0.5, x + 0.5, y - 0.5, y + 0.5
    borders = x_min, x_max, y_min, y_max

    # Check for interceptionpoints between circle and axis
    x_lower = RadialEqu(r, y_min)
    x_upper = RadialEqu(r, y_max)
    y_left = RadialEqu(r, x_min)
    y_right = RadialEqu(r, x_max)
    intersec = x_lower, x_upper, y_left, y_right
    print intersec
    # Check for actual interceptionpoints between cicle and pixel
    #intersected == intersec_x_low, intersec_x_up, intersec_y_left, intersec_y_right
    intersected = get_intercections(borders, x_lower, x_upper, y_left, y_right)
    # print intersec, intersected

    '''
    CASE 1  
    The whole pixel lays within the circle    
    '''
    Area = np.zeros_like(x)
    max_mask = np.sqrt(x_max ** 2 + y_max ** 2) <= r
    Area[max_mask] = 1
    intersec_mask = (np.sqrt(x_max ** 2 + y_max ** 2) > r) & (np.sqrt(x_min ** 2 + y_min ** 2) < r)
    
    '''
    CASE 2  
    The circle intersects with the left and bottom border of the pixel    
    '''
    def masked_borders(borders, mask):
        return tuple(a[mask] for a in borders)
        
    imask = intersected[0] & intersected[2] & intersec_mask
    Area[imask] = CASE2_bottom_left_corner(
        masked_borders(borders, imask), x_lower[imask], r
        )

    '''
    CASE 3 
    The circle intersects with the right and top border of the pixel    
    '''
    imask = intersected[1] & intersected[3] & intersec_mask
    Area[imask] = CASE3_upper_right_corner(
        masked_borders(borders, imask), x_upper[imask], r
        )

    '''
    CASE 4 
    The circle intersects with left and right border of the pixel    
    '''
    imask = intersected[2] & intersected[3] & intersec_mask
    Area[imask] = CASE4_both_y_borders(
        masked_borders(borders, imask), y_left[imask], y_right[imask], r
        )

    '''
    CASE 5  
    The circle intersects with the top and bottom border of the pixel    
    '''
    imask = intersected[0] & intersected[1] & intersec_mask
    Area[imask] = CASE5_both_x_borders(
        masked_borders(borders, imask), x_lower[imask], x_upper[imask], r
        )

    return Area, intersec, intersected

def main():
    px, py = np.meshgrid(
        np.arange(185) + 0.5,
        np.arange(185) + 0.5,
        )
    
    r = 1.06
    A, intersec, intersected = get_circle_area_for_pix((px, py), r)
    
    print '\n'.join(map(
        lambda a: '{:.4f} {:.4f} {:.4f}'.format(*a),
        zip(px.flatten(), py.flatten(), A.flatten())  #, intersec #, intersected
        ))
    
    plt.close()
    plt.scatter(px, py, c=A)
    plt.show()
    
    
if __name__ == '__main__':
    main()
    
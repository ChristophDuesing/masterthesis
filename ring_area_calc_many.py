#!/usr/bin/python
# -*- coding: utf-8 -*-

# Program by Christoph Duesing
# This algorithm producess a the Powerspectrum of a random field
import numpy as np
import pyfits as pf
import copy
import pylab as pl
import matplotlib.pyplot as plt


def get_pixels(r_max, r_min, dists_sq):
    Check_Range = np.logical_and(
        radius_mask_func(dists_sq, r_max + np.sqrt(2)),
        np.logical_not(radius_mask_func(dists_sq, r_min - np.sqrt(2)))
        )

    Test_Plane_copy[Check_Range] = 1
    xval, yval = np.where(Test_Plane_copy == 1)

    return zip(xval, yval)


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


def RadialEqu(Radius , coord):
    return np.sqrt(abs(Radius ** 2 - coord ** 2))


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
    if np.mean(Area)==0:
            
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
    start, stop = 0., 99.
    steps = 100.
    z = np.linspace(start, stop, steps)
    x = len(z)
    a, b = np.meshgrid(
        np.linspace(0, x - 1, x) + 0.5,
        np.linspace(0, x - 1, x) + 0.5
        )

    Test_Plane = a * b
    Test_Plane_copy = Test_Plane.copy()

    r_min_start, r_max_start, a0, b0 = 1, x / 2., x / 2., x / 2.
    r_grid = np.linspace(r_min_start, r_max_start, 64)
    r_grid_mean = 0.5 * (r_grid[1:] + r_grid[:-1])

    dists_sq = (a - a0) ** 2 + (b - b0) ** 2

    radius_mask_func = lambda dists_sq, r: dists_sq <= r ** 2

    r_min, r_max = r_grid[34], r_grid[35]
    b = r_max - r_min

    coordinates = get_pixels(r_max, r_min, dists_sq)

    Circ_Area = []

    for i in coordinates:
        x, y = i
        x, y = x - (stop - start) / 2, y - (stop - start) / 2
        # First Quadrant is used only
        if x >= 0 and y >= 0:
            Area = (
                get_circle_area_for_pix(i, r_max)[0] -
                get_circle_area_for_pix(i, r_min)[0]
                )

        x_old, y_old = x + (stop - start) / 2, y + (stop - start) / 2

        Areas = (x_old, y_old, Area, Test_Plane[x_old, y_old], Area
                 * Test_Plane[x_old, y_old]), (x_old, -y_old, Area,
                Test_Plane[x_old, -y_old], Area * Test_Plane[x_old, -y_old]), \
            (-x_old, y_old, Area, Test_Plane[-x_old, y_old], Area
             * Test_Plane[-x_old, y_old]), (-x_old, -y_old, Area,
                Test_Plane[-x_old, -y_old], Area * Test_Plane[-x_old, -y_old])

        Circ_Area.append(Areas)

        print Areas

    # pl.figure()
    # pl.imshow(Test_Plane_copy,interpolation='nearest',)
    # pl.show()


def main2():
    '''Test the 5 different cases and plot'''

    def plot_with_area(px, py, r):
        rect_stroke = (
            [px - 0.5, px + 0.5, px + 0.5, px - 0.5, px - 0.5],
            [py - 0.5, py - 0.5, py + 0.5, py + 0.5, py - 0.5],
            )

        circle_stroke = (
            r * np.sin(np.linspace(0, np.pi / 2, 100)),
            r * np.cos(np.linspace(0, np.pi / 2, 100)),
            )

        A, intersec, intersected = get_circle_area_for_pix(
            (np.atleast_1d(px), np.atleast_1d(py)), r
            )
        print intersec

        plt.close()
        plt.plot(*rect_stroke, c='k', ls='-')
        plt.plot(*circle_stroke, c='r', ls='-')
        if intersected[0]:
            plt.plot(intersec[0], py - 0.5, 'go')
        if intersected[1]:
            plt.plot(intersec[1], py + 0.5, 'go')
        if intersected[2]:
            plt.plot(px - 0.5, intersec[2], 'go')
        if intersected[3]:
            plt.plot(px + 0.5, intersec[3], 'go')

        plt.gca().set_aspect('equal')
        plt.title('px, py = ({}, {}); r = {}; A = {}'.format(px, py, r, A))
        plt.show()
    
    # circle smaller than lower left edge
    px, py = 0.5, 2.5
    #r = np.sqrt((px - 0.5) ** 2 + (py - 0.5) ** 2) - 0.2
    r= 2*np.sqrt(2)
    plot_with_area(px, py, r)

    # circle hits lower left edge
    px, py = 0.5, 2.5
    r = 2*np.sqrt(2)
    plot_with_area(px, py, r)  # <- wrong as well ! 

    ## circle larger than lower left edge
    #px, py = 3.5, 5.5
    #r = np.sqrt((px - 0.5) ** 2 + (py - 0.5) ** 2) + 0.2
    #plot_with_area(px, py, r)

    ## circle hits lower right edge
    #px, py = 3.5, 5.5
    #r = np.sqrt((px + 0.5) ** 2 + (py - 0.5) ** 2)
    #plot_with_area(px, py, r)

    ## circle larger than lower right edge
    #px, py = 3.5, 5.5
    #r = np.sqrt((px + 0.5) ** 2 + (py - 0.5) ** 2) + 0.2
    #plot_with_area(px, py, r)

    ## circle hits upper left edge
    #px, py = 3.5, 5.5
    #r = np.sqrt((px - 0.5) ** 2 + (py + 0.5) ** 2)
    #plot_with_area(px, py, r)

    ## circle larger than upper left edge
    #px, py = 3.5, 5.5
    #r = np.sqrt((px - 0.5) ** 2 + (py + 0.5) ** 2) + 0.2
    #plot_with_area(px, py, r)  # <-- wrong CORRECTED ! 

    ## circle hits upper right edge
    #px, py = 3.5, 5.5
    #r = np.sqrt((px + 0.5) ** 2 + (py + 0.5) ** 2)
    #plot_with_area(px, py, r)  # <-- wrong SEEMS RIGHT ! 

    ## circle larger than upper right edge
    #px, py = 3.5, 5.5
    #r = np.sqrt((px + 0.5) ** 2 + (py + 0.5) ** 2) + 0.2
    #plot_with_area(px, py, r)




    ## circle smaller than lower left edge
    #px, py = 6.5, 5.5
    #r = np.sqrt((px - 0.5) ** 2 + (py - 0.5) ** 2) - 0.2
    #plot_with_area(px, py, r)

    ## circle hits lower left edge
    #px, py = 6.5, 5.5
    #r = np.sqrt((px - 0.5) ** 2 + (py - 0.5) ** 2)
    #plot_with_area(px, py, r)

    ## circle larger than lower left edge
    #px, py = 6.5, 5.5
    #r = np.sqrt((px - 0.5) ** 2 + (py - 0.5) ** 2) + 0.2
    #plot_with_area(px, py, r)

    ## circle hits lower right edge
    #px, py = 6.5, 5.5
    #r = np.sqrt((px + 0.5) ** 2 + (py - 0.5) ** 2)
    #plot_with_area(px, py, r)  

    ## circle larger than lower right edge
    #px, py = 6.5, 5.5
    #r = np.sqrt((px + 0.5) ** 2 + (py - 0.5) ** 2) + 0.2
    #plot_with_area(px, py, r)  # <-- wrong SEEMS ALSO RIGHT 

    ## circle hits upper left edge
    #px, py = 6.5, 5.5
    #r = np.sqrt((px - 0.5) ** 2 + (py + 0.5) ** 2)
    #plot_with_area(px, py, r)

    ## circle larger than upper left edge
    #px, py = 6.5, 5.5
    #r = np.sqrt((px - 0.5) ** 2 + (py + 0.5) ** 2) + 0.2
    #plot_with_area(px, py, r)  # <-- wrong CORRECTED

    ## circle hits upper right edge
    #px, py = 6.5, 5.5
    #r = np.sqrt((px + 0.5) ** 2 + (py + 0.5) ** 2)
    #plot_with_area(px, py, r)  # <-- wrong

    ## circle larger than upper right edge
    #px, py = 6.5, 5.5
    #r = np.sqrt((px + 0.5) ** 2 + (py + 0.5) ** 2) + 0.2
    #plot_with_area(px, py, r)


def main3():
    px, py = np.meshgrid(
        np.arange(185) + 0.5,
        np.arange(185) + 0.5,
        )
    
    r1 = 1.79586770
    r2 = 2.10161708
    
    A1, intersec1, intersected1 = get_circle_area_for_pix((px, py), r1)
    A2, intersec2, intersected2 = get_circle_area_for_pix((px, py), r2)
    
    A =A2-A1 
    print '\n'.join(map(
        lambda a: '{:.4f} {:.4f} {:.4f}'.format(*a),
        zip(px.flatten(), py.flatten(), A.flatten())  #, intersec #, intersected
        ))
    
    plt.close()
    plt.scatter(px, py, c=A)
    plt.show()
    
    print np.sum(A)
    print A.shape
    print A
    print np.max(A)
    
    # circle larger than upper right edge
    px, py = 75.5, 8.5
    r = np.sqrt((px + 0.5) ** 2 + (py + 0.5) ** 2)
    plot_with_area(px, py, r)

    
if __name__ == '__main__':
    #main()
    main2()

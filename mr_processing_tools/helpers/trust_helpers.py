#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Sky Jones
"""
import sys

from copy import copy

import matplotlib.pyplot as plt
from scipy import ndimage
import numpy as np
import numpy.ma as ma
import matplotlib.colors as colors




def manual_sss(trust_diff, r=3):
    
    
    circle_pixels, x, y = auto_sss(trust_diff, r) # initial guess
    
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8,8))
    ax.imshow(trust_diff, cmap='gist_gray')
    ax.set_title('Superior sagittal sinus isolation\nClick whitespace to confirm, or click image again to make new mask')
    
    palette = copy(plt.cm.Reds)
    palette.set_bad(color='white', alpha=0)
    
    confirmed = False
        
    def on_click(event):
        global coords
        coords = (event.xdata, event.ydata)
        #print(f'you clicked {coords}')
    
    def on_close(event):
        print('Closed figure')
        raise Exception('SSS isolation closed before confirmation')
        sys.exit()
        
    cid = fig.canvas.mpl_connect('button_press_event', on_click)
    fid = fig.canvas.mpl_connect('close_event', on_close)
    while not confirmed:
        
        ex, why = round(x), round(y)
        
        circle_pixels = circle_mask(trust_diff, why, ex, r)
        
        circ_ma = ma.array(circle_pixels, mask=~circle_pixels)
        
        imask = ax.imshow(circ_ma, alpha=0.5, cmap=palette, vmin=0, vmax=1)
        fig.canvas.draw()
        fig.waitforbuttonpress() # button press will change value of coords
        if None not in coords:
            imask.set_visible(False)
            x, y = coords
        else:
            confirmed = True
        
    #print(f'Coords selected: {ex, why}')
    
    fig.canvas.mpl_disconnect(cid)
    fig.canvas.mpl_disconnect(fid)
        
    plt.close(fig)
    
    return circle_pixels, ex, why


def auto_sss(trust_diff, r=3, sigma=2):
    '''
    takes a 2d TRUST difference image (label-control) and returns a mask
    of the superior sagittal sinus
    
    Does this by:
        1) Blurring the image to dampen spurious bright spots
        2) Making a circle of radius r around the brightest spot in blurred image
    
    '''
    blurred = ndimage.gaussian_filter(trust_diff, sigma)
    max_gray_level = blurred.max()
    rows, columns = np.where(blurred == max_gray_level) # this could return multiple locations if there is a tie for max intensity
    
    why = rows[0]
    ex = columns[0]
    
    
    #print(f'Coords selected: {ex, why}')
    
    circle_pixels = circle_mask(trust_diff, why, ex, r)
    
    return circle_pixels, ex, why
    

def circle_mask(arr, py, px, r):
    # given an array and the central coordinates of a circle, return
    # a circular mask of radius
    
    dims = arr.shape
    
    columns_in_image, rows_in_image = np.meshgrid(np.arange(0,dims[0]), np.arange(0,dims[1]))
    circle_radius = r
    circle_pixels = ((rows_in_image - py)**2 + (columns_in_image - px)**2) <= circle_radius**2
    return circle_pixels

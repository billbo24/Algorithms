# -*- coding: utf-8 -*-
"""
Created on Mon Feb 22 21:44:59 2021

@author: dingl
Will be attempting to implement error diffusion dithering. 

If we import a png I believe we have RGB and opacity channels
Where 255,255,255,255 is white with full opacity

Now, to implement floyd steinberg dithering, my gut tells me we need to:
    1) Map each of our pixels from input space to output space
    2) Determine the error in doing step one

Note that this approach implies a sort of "1a" step, which is scoring 
the input pixels.  I think the best approach here is using the RGB color model,
which are really just points in three space, determing which color 
in our goal palette is closest to our current color, and then noting that error.
This method will allow us to use an arbitrary number of colors in our output
palette.  I think if we really wanted to we could then apply another function
to the output RGB values to convert them to other color models if we were 
so inclined (I think .jpg uses YBrBc or whatever it is)


    
"""

from matplotlib.pyplot import imshow
import numpy as np
from PIL import Image
from IPython.display import display



#image_path = r'C:\\Users\\dingl\\Pictures\\2_20\\Marshmallow.png'
image_path = r'C:\\Users\\dingl\\Pictures\\2_20\\Winnie.PNG'
output_filename = 'Winnie_dither.png'

img = Image.open(image_path)
my_palette = [[0,0,0],[0,0,255],[255,255,255]] #black and then white
weight = 0.3
img_info = np.array(img)


def vector_add(a,b):
    return [min(max(a[i]+b[i],0),255) for i in range(3)] #restrains us to the 0,255 space

def score_pixel(input_pixel,color_palette):
    euclid = lambda a,b: [(a[i]-b[i])**2 for i in range(3)]
    min_error = 1000 #I think this is high enough, max distance in 3-space is...255rad(3)
    ans_color = []
    for i in color_palette:
        this_error = sum(euclid(input_pixel,i))**.5 + .01
        if this_error < min_error:
            min_error,ans_color = this_error,i
    #Note this gives us the vector which when added to our
    #input pixel gives the closest color
    unit_error_vec = [(ans_color[i]-input_pixel[i])*(1/min_error) for i in range(3)]
    return unit_error_vec,ans_color,min_error

def subtract_error(pixel,error_vec,magnitude):
    #Because we had to add [a,b,c] to our pixel to move it to the 
    #destination color in the palette, we need to move 
    #neighbor pixels away by the same amount
    new_error = [min(max(magnitude*error_vec[i]*-1,0),255) for i in range(3)]
    return vector_add(pixel,new_error)

def diffuse_error(row,col,error_mag,error_vec,source_image):
    #floyd steinberg adds 7/16 of the error on pixel to the right
    #3/16 southwest, 5/16 south, and 1/16 southeast
    #no clue what happens on the edges, just gonna do a try, except
    
    
    try:
        source_image[row,col+1] = subtract_error(source_image[row,col+1],error_vec,(7/16)*error_mag)
    except:
        pass

    try:
        source_image[row+1,col-1] = subtract_error(source_image[row+1,col-1],error_vec,(3/16)*error_mag)
    except:
        pass
    
    try:
        source_image[row+1,col] = subtract_error(source_image[row+1,col],error_vec,(5/16)*error_mag)
    except:
        pass
    
    try:
        source_image[row+1,col+1] = subtract_error(source_image[row+1,col+1],error_vec,(1/16)*error_mag)
    except:
        pass
    return source_image


def floyd_steinberg_dithering(input_image,color_palette):
    img = input_image
    display(img)
    img_info = np.array(img)[:,:,:3]
    
    rows = len(img_info)
    cols = len(img_info[0])
    
    mapped_img = np.zeros((rows,cols,3)) #this is an empty matrix which #will become our augmented image
    errors = np.zeros((rows,cols,3)) #this will store our error vectors
    
    
    for row in range(rows):
        for col in range(cols):
            
            t = score_pixel(img_info[row,col],color_palette)
            new_color = t[1]
            error_vector = t[0]
            error_mag = t[2]
            
            mapped_img[row,col] = new_color #gives us an image mapped to the space
            
            #now we must diffuse the error
            img_info = diffuse_error(row,col,error_mag,error_vector,img_info)
            
            
            #errors[row,col] = unit_error_vec
            
    
    return mapped_img
    

'''
y = score_pixel(test_pixel  ,my_palette)
x = subtract_error(other_test_pixel,y[0],20)
print(y[0],x)
'''

my_img = Image.open(image_path)

#things have to be in the format uint8 to convert back to an image.  
a = floyd_steinberg_dithering(my_img,my_palette).astype('uint8')
resultim = Image.fromarray(a)
resultim.save(r'C:\\Users\\dingl\\Pictures\\2_20\\' + output_filename)


#now my understanding is that for floyd dithering we need 




#This is how we go from array to an image














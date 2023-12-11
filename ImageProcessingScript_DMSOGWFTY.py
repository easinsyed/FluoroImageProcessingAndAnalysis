#!/usr/bin/env python
# coding: utf-8

# In[16]:


import os
import glob
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import re
import time
import bisect
from numba import jit
import cv2
import skimage 
import pandas as pd
from skimage import exposure
import glob
from skimage.restoration import rolling_ball


# In[6]:


#Function to rescale images
def rescale_to_bits(im: np.ndarray, n_bits: int) -> np.ndarray:
    """
    Helper method to rescale an image to have values supported with a specific bit-depth (unsigned integer).
    """
    max_value = 2**n_bits - 1
    im = im - im.min()
    im = im * (max_value / im.max())
    im = np.round(im)
    return im


# In[7]:


def colorize(im, color, clip_percentile=0.1):
    """
    Helper function to create an RGB image from a single-channel image using a 
    specific color.
    """
    # Check that we do just have a 2D image
    if im.ndim > 2 and im.shape[2] != 1:
        raise ValueError('This function expects a single-channel image!')
        
    # Rescale the image according to how we want to display it
    im_scaled = im.astype(np.float32) - np.percentile(im, clip_percentile)
    im_scaled = im_scaled / np.percentile(im_scaled, 100 - clip_percentile)
    im_scaled = np.clip(im_scaled, 0, 1)
    
    # Need to make sure we have a channels dimension for the multiplication to work
    im_scaled = np.atleast_3d(im_scaled)
    
    # Reshape the color (here, we assume channels last)
    color = np.asarray(color).reshape((1, 1, -1))
    return im_scaled * color


# In[8]:


def gray2color(I,channel):
    """
    Compute color image from intensity in fluorescence in a given channel.
    Arguments:
    -----------
        I: np.ndarray
            Input fluorescence image (2D).
        channel: int
            Channel to code the image in (0: Red, 1: Green, 2: Blue).
    Returns:
    ----------- 
        I_color: np.ndarray
            The computed output image in color.
    """
    
    I_color = np.dstack((
        skimage.exposure.rescale_intensity(I if channel==0 else np.zeros_like(I)),
        skimage.exposure.rescale_intensity(I if channel==1 else np.zeros_like(I)),
        skimage.exposure.rescale_intensity(I if channel==2 else np.zeros_like(I)),
        ))
    return I_color


# In[9]:


def subtract_background(image, radius=25, light_bg=False):
    # import libraries
    from skimage.morphology import white_tophat, black_tophat, disk 
    
    # generate structuring element
    str_el = disk(radius)
     
    # use appropriate filter depending on the background colour
    if light_bg:
        return black_tophat(image, str_el)
    else:
        return white_tophat(image, str_el)


# In[ ]:


import os
import numpy as np
from PIL import Image

# Define the folder containing the images
folder_path = "D:/Easin_Microscopy_Images/IncellAnalyzer/Easin/20230309_ZF3/Easin_20x_ZF3/Easin_20x_ZF3_1"

# Define the well and site numbers
well_numbers = ["A-01", "A-02", ..., "H-12"]
site_numbers = ["fld 1", "fld 2", ..., "fld N"]

# Define the channels and their corresponding file name keywords
channels = ["Blue - FITC", "Red - Cy5", "UV - DAPI"]
channel_keywords = ["wv Blue - FITC", "wv Red - Cy5", "wv UV - DAPI"]

# Loop through all the wells and sites
for well in well_numbers:
    for site in site_numbers:
        try:
            # Loop through all the channels
            channels_data = []
            for i, channel in enumerate(channels):
                # Define the file name of the image for this channel, well, and site
                file_name = "{}({} {} {})".format(well, site, channel_keywords[i], channel)
                file_path = os.path.join(folder_path, file_name)

                # Load the image for this channel, well, and site
                image = np.asarray(Image.open(file_path))
                channels_data.append(image)

            # Merge the images for all channels into a single color image
            color_image = np.dstack(channels_data)

            # Save the color image with the same file name as the first channel's image
            save_file_name = "{}({} color)".format(well, site)
            save_file_path = os.path.join(folder_path, save_file_name)
            Image.fromarray(color_image).save(save_file_path)

        except FileNotFoundError:
            # Skip this well or site if the image file is not found
            continue


# In[9]:


import os
import re
import os.path 
t1 = time.perf_counter()

# Define the folder containing the images
pathIn = r"D:/Easin_Microscopy_Images/IncellAnalyzer/Easin/Easin_ZF3_20230317_2"
pathOut = r"D:/Easin_processed_images"+ '/'+os.path.basename(pathIn)


isExist = os.path.exists(pathOut) # Check whether the specified path exists or not
if not isExist:
    os.makedirs(pathOut)
    print("The new directory is created!")    

def prepend(list, str):
    str += '{0}'
    list = [str.format(i) for i in list]
    return(list)

# Define the well and site numbers
well_numbers_list = []
site_number_list = []

# Loop through all the files in the folder
for file_name in os.listdir(pathIn):
    # Check if the file name matches the expected pattern
    if "(" in file_name and ")" in file_name:
        # Extract the well and site numbers from the file name
        well, site = file_name.split("(")[:2]
        well_numbers_list.append(well.strip())
        well_names_unique = sorted(set(well_numbers_list))
        site_number_list.append(re.search('\d+', site).group())
        min_site_number = int(min(site_number_list))
        max_site_number = int(max(site_number_list))
        
# Remove duplicates (set) and sort(sorted) the well and site numbers
well_numbers = sorted(list(set(well_numbers_list)))
site_numbers = sorted(set(prepend(site_number_list, 'fld ')))

image_type = '.tif'

# Define the channels and their corresponding file name keywords
channels = ["green", "red", "blue"]
channel_keywords = ["wv Blue - FITC", "wv Red - Cy5", "wv UV - DAPI"]

# Loop through all the wells and sites
for well in well_numbers:
    for site in site_numbers:
        im = []
        for i, channel in enumerate(channels):
            # Define the file name of the image for this channel, well, and site
            file_name = "{}({} {}){}".format(well, site, channel_keywords[i], image_type)
            file_path = pathIn +'/'+ file_name
               
            try:
                im = plt.imread(file_path)
                assert im.dtype == np.uint16
            except FileNotFoundError:
                print(f"Skipping well {well}, site {site} for channel {channel}")
                continue
                
            im = rescale_to_bits(im, 8)
            #check which channel it is and colorize accordingly
            save_file_name = "{} ({}, {}).tiff".format(well, site, channel_keywords[i])
            save_file_name_path = pathOut + '/' + save_file_name
            if channel == 'green':
                im_g = colorize(im[...,], (0, 1, 0),clip_percentile=0.01)
                plt.imsave(save_file_name_path, im_g,format='TIFF')
            elif channel == 'red':
                im_r = colorize(im[...,], (1, 0, 0),clip_percentile=0.01)
                plt.imsave(save_file_name_path, im_r,format='TIFF')
            else :
                im_b = colorize(im[...,], (0, 0, 1),clip_percentile=0.01)
                plt.imsave(save_file_name_path, im_b,format='TIFF')

        save_file_name_rg = "{} ({}) rg.tiff".format(well, site)
        save_file_name_path_rg = pathOut + '/' + save_file_name_rg
        merge1 = np.clip(im_r + im_g , 0, 1)
        plt.imsave(save_file_name_path_rg, merge1,format='TIFF') #Saving the 2 colors merged file
        
        save_file_name_rgb = "{} ({}) rgb.tiff".format(well, site)
        save_file_name_path_rgb = pathOut + '/' + save_file_name_rgb
        merge2 = np.clip(im_r + im_g + im_b , 0, 1)
        plt.imsave(save_file_name_path_rgb, merge2,format='TIFF') #Saving the 3 colors merged file

t2 = time.perf_counter()
print(t2-t1) #Get an average on the run time


# In[141]:


#folder_path = r"V:/Laboratoires/Lécuyer/Juan-Carlos/HCS/20230802_MCF7_Ab_DMSOGWFTY-2/20230801_MCF7_Ab_DMSO/TimePoint_1/DMSO_F08_s4_w3.TIF"
#print(folder_path)
#well = 'B - 02'
#site = 'fld 1'
#channel_keywords = ["wv Blue - FITC", "wv Red - Cy5", "wv UV - DAPI"]
#channel = "Red - Cy5"
#image_type = '.tif'
#file_name ="{}({} {}){}".format(well, site, channel_keywords[1], image_type)
#file_path = os.path.join(folder_path, file_name)
#print(file_path)

im = plt.imread( r"V:/Laboratoires/Lécuyer/Juan-Carlos/HCS/20230802_MCF7_Ab_DMSOGWFTY-2/20230801_MCF7_Ab_DMSO/TimePoint_1/DMSO_F08_s4_w1.TIF")
assert im.dtype == np.uint16
# Show the image
plt.imshow(im, cmap='gray')
plt.show()


im_8bits = rescale_to_bits(im, 8)
plt.imshow(im_8bits, cmap='gray')
plt.show()


img_g = subtract_background(im_8bits, radius=50, light_bg=False)
#background = rolling_ball(im_8bits)
#im_8bits = im_8bits - background

# The color we provide gives RGB values in order, between 0 and 1
im_red = colorize(img_g[...,], (0, 1, 0))
#im_red = skimage.exposure.equalize_adapthist(im_red)
plt.imshow(im_red)
plt.axis(False)
plt.title('Magenta')
plt.show()

plt.imsave(r'C:\Users\PadillJ\Desktop\test3.tif', im_red,format='TIFF')


# In[46]:


pathOut = r"D:/Easin_processed_images"+ '/'+basename(pathIn)


isExist = os.path.exists(pathOut) # Check whether the specified path exists or not
if not isExist:
    os.makedirs(pathOut)
    print("The new directory is created!")    
#OutFolder = 

#isExist = os.path.exists(pathOut+) # Check whether the specified path exists or not
#if not isExist:
#    os.makedirs(geneFolder)
#    print("The new directory is created!")    


# In[81]:


import os
import re
import os.path 

# Define input and output directories
pathIn = r"S:/20230802_MCF7_Ab_DMSOGWFTY-2/20230801_MCF7_Ab_DMSO/TimePoint_1"
pathOut = r"S:/20230918_MCF7_Ab_DMSOGWFTY-2_ImageEnhancement/DMSO"+ '/'+os.path.basename(pathIn)

isExist = os.path.exists(pathOut) # Check whether the specified path exists or not

if not isExist:
    os.makedirs(pathOut)
    print("The new directory is created!")    

def prepend(list, str):
    str += '{0}'
    list = [str.format(i) for i in list]
    return(list)

# Define the well and site numbers
well_numbers_list = []
site_number_list = []

# Loop through all the files in the folder
for file_name in os.listdir(pathIn):
    condition, well, site, wavelength = file_name.split("_")
    well_numbers_list.append(well.strip())
    well_names_unique = sorted(set(well_numbers_list))
    site_number_list.append(re.search('\d+', site).group())
    min_site_number = int(min(site_number_list))
    max_site_number = int(max(site_number_list))
    #print(condition, well, site, wavelength)

# Remove duplicates (set) and sort(sorted) the well and site numbers
well_numbers = sorted(list(set(well_numbers_list)))
site_numbers = sorted(set(prepend(site_number_list, 's'))) 
image_type = '.TIF'

# Define the channels and their corresponding file name keywords
channels = ["blue", "green", "red"]
channel_keywords = ["w1", "w2", "w3"]

# Loop through all the wells and sites
for well in well_numbers:
    for site in site_numbers:
        im = []
        for i, channel in enumerate(channels):
            # Define the file name of the image for this channel, well, and site
            file_name = "{}_{}_{}_{}{}".format(condition, well, site, channel_keywords[i], image_type)
            print(file_name)
            file_path = pathIn +'/'+ file_name
            print(file_path)   
            try:
                im = plt.imread(file_path)
                assert im.dtype == np.uint16
            except FileNotFoundError:
                print(f"Skipping well {well}, site {site} for channel {channel}")
                continue
                
            im = rescale_to_bits(im, 8)
            #check which channel it is and colorize accordingly
            save_file_name = "{}_Pr_{}_{}_{}.tiff".format(condition, well, site, channel_keywords[i])
            save_file_name_path = pathOut + '/' + save_file_name
            if channel == 'green':
                background_g = rolling_ball(im)
                im_g_corr = im - background_g
                im_g = colorize(im_g_corr[...,], (0, 1, 0),clip_percentile=0.01)
                plt.imsave(save_file_name_path, im_g,format='TIFF')
            elif channel == 'red':
                background_r = rolling_ball(im)
                im_r_corr = im - background_r
                im_r = colorize(im_r_corr[...,], (1, 0, 0),clip_percentile=0.01)
                plt.imsave(save_file_name_path, im_r,format='TIFF')
            else :
                im_b = subtract_background(im, radius=50, light_bg=False)
                im_b = skimage.exposure.equalize_adapthist(im_b)
                im_b = colorize(im_b[...,], (0, 0, 1),clip_percentile=0.01)
                plt.imsave(save_file_name_path, im_b,format='TIFF')
                
        save_file_name_rg = "{}_Pr_{}_{}_mergeRG.tiff".format(condition, well, site)
        save_file_name_path_rgb = pathOut + '/' + save_file_name_rgb
        merge2 = np.clip(im_r + im_g  , 0, 1)
        plt.imsave(save_file_name_path_rgb, merge2,format='TIFF') #Saving the 3 colors merged file

        save_file_name_rb = "{}_Pr_{}_{}_mergeRB.tiff".format(condition, well, site)
        save_file_name_path_rgb = pathOut + '/' + save_file_name_rgb
        merge2 = np.clip(im_r + im_b , 0, 1)
        plt.imsave(save_file_name_path_rgb, merge2,format='TIFF') #Saving the 3 colors merged file

        save_file_name_gb = "{}_Pr_{}_{}_mergeGB.tiff".format(condition, well, site)
        save_file_name_path_rgb = pathOut + '/' + save_file_name_rgb
        merge2 = np.clip(im_g + im_b , 0, 1)
        plt.imsave(save_file_name_path_rgb, merge2,format='TIFF') #Saving the 3 colors merged file
                
        save_file_name_rgb = "{}_Pr_{}_{}_mergeRGB.tiff".format(condition, well, site)
        save_file_name_path_rgb = pathOut + '/' + save_file_name_rgb
        merge2 = np.clip(im_r + im_g + im_b , 0, 1)
        plt.imsave(save_file_name_path_rgb, merge2,format='TIFF') #Saving the 3 colors merged file
        


# In[147]:


inputFolder = 'S:/20230802_MCF7_Ab_DMSOGWFTY-2/'

listOfDirectories = []

for inputFolder, dirs, files in os.walk(inputFolder):
    for subdir in dirs:
        print(os.path.join(inputFolder, subdir))
        listOfDirectories.append(os.path.join(inputFolder, subdir))
        
def partial(lst, query):
    return [s for s in lst if query in s]

print(partial(listOfDirectories, 'TimePoint'))

listOfDirectoriesWithImages = partial(listOfDirectories, 'TimePoint')

outputFolder = r"S:/20230918_MCF7_Ab_DMSOGWFTY-2_ImageEnhancement"

isExist = os.path.exists(outputFolder) # Check whether the specified path exists or not

if not isExist:
    os.makedirs(outputFolder)
    print("The new directory is created!") 

for x in listOfDirectoriesWithImages:
    pathIn = x
    pathOut = outputFolder + '/' + os.path.basename(os.path.dirname(pathIn))
    
    isExist = os.path.exists(pathOut) # Check whether the specified path exists or not

    if not isExist:
        os.makedirs(pathOut)
        print("The new directory is created!")    

    def prepend(list, str):
        str += '{0}'
        list = [str.format(i) for i in list]
        return(list)

    # Define the well and site numbers
    well_numbers_list = []
    site_number_list = []

    # Loop through all the files in the folder
    for file_name in os.listdir(pathIn):
        condition, well, site, wavelength = file_name.split("_")
        well_numbers_list.append(well.strip())
        well_names_unique = sorted(set(well_numbers_list))
        site_number_list.append(re.search('\d+', site).group())
        min_site_number = int(min(site_number_list))
        max_site_number = int(max(site_number_list))
        #print(condition, well, site, wavelength)

    # Remove duplicates (set) and sort(sorted) the well and site numbers
    well_numbers = sorted(list(set(well_numbers_list)))
    site_numbers = sorted(set(prepend(site_number_list, 's'))) 
    image_type = '.TIF'

    # Define the channels and their corresponding file name keywords
    channels = ["blue", "green", "red"]
    channel_keywords = ["w1", "w2", "w3"]

    # Loop through all the wells and sites
    for well in well_numbers:
        for site in site_numbers:
            im = []
            for i, channel in enumerate(channels):
                # Define the file name of the image for this channel, well, and site
                file_name = "{}_{}_{}_{}{}".format(condition, well, site, channel_keywords[i], image_type)
                print(file_name)
                file_path = pathIn +'/'+ file_name
                print(file_path)   
                try:
                    im = plt.imread(file_path)
                    assert im.dtype == np.uint16
                except FileNotFoundError:
                    print(f"Skipping well {well}, site {site} for channel {channel}")
                    continue

                im = rescale_to_bits(im, 8)
                #check which channel it is and colorize accordingly
                save_file_name = "{}_Pr_{}_{}_{}.tiff".format(condition, well, site, channel_keywords[i])
                save_file_name_path = pathOut + '/' + save_file_name
                if channel == 'green':
                    
                    background_g = rolling_ball(im)
                    im_g_corr = im - background_g
                    im_g = colorize(im_g_corr[...,], (0, 1, 0),clip_percentile=0.01)
                    plt.imsave(save_file_name_path, im_g,format='TIFF')
                elif channel == 'red':
                    background_r = rolling_ball(im)
                    im_r_corr = im - background_r
                    im_r = colorize(im_r_corr[...,], (1, 0, 0),clip_percentile=0.01)
                    plt.imsave(save_file_name_path, im_r,format='TIFF')
                else :
                    im_b = subtract_background(im, radius=50, light_bg=False)
                    im_b = colorize(im_b[...,], (0, 0, 1),clip_percentile=0.01)
                    im_b = skimage.exposure.equalize_adapthist(im_b)
                    plt.imsave(save_file_name_path, im_b,format='TIFF')

            save_file_name_rg = "{}_Pr_{}_{}_mergeRG.tiff".format(condition, well, site)
            save_file_name_path_rg = pathOut + '/' + save_file_name_rg
            merge1 = np.clip(im_r + im_g  , 0, 1)
            plt.imsave(save_file_name_path_rg, merge1,format='TIFF') #Saving the rg colors merged file

            save_file_name_rb = "{}_Pr_{}_{}_mergeRB.tiff".format(condition, well, site)
            save_file_name_path_rb = pathOut + '/' + save_file_name_rb
            merge2 = np.clip(im_r + im_b , 0, 1)
            plt.imsave(save_file_name_path_rb, merge2,format='TIFF') #Saving the rb colors merged file

            save_file_name_gb = "{}_Pr_{}_{}_mergeGB.tiff".format(condition, well, site)
            save_file_name_path_gb = pathOut + '/' + save_file_name_gb
            merge3 = np.clip(im_g + im_b , 0, 1)
            plt.imsave(save_file_name_path_gb, merge3,format='TIFF') #Saving the gb colors merged file

            save_file_name_rgb = "{}_Pr_{}_{}_mergeRGB.tiff".format(condition, well, site)
            save_file_name_path_rgb = pathOut + '/' + save_file_name_rgb
            merge4 = np.clip(im_r + im_g + im_b , 0, 1)
            plt.imsave(save_file_name_path_rgb, merge4,format='TIFF') #Saving the rgb colors merged file



# In[65]:


import os
import re
import os.path 

# Define input and output directories
pathIn = r"S:/20230802_MCF7_Ab_DMSOGWFTY-2/20230803_MCF7_Ab_FTY720/TimePoint_1"
pathOut = r"S:/20230918_MCF7_Ab_DMSOGWFTY-2_ImageEnhancement/FTY"+ '/'+os.path.basename(pathIn)

def prepend(list, str):
    str += '{0}'
    list = [str.format(i) for i in list]
    return(list)

# Define the well and site numbers
well_numbers_list = []
site_number_list = []

# Loop through all the files in the folder
for file_name in os.listdir(pathIn):
    condition, well, site, wavelength = file_name.split("_")
    well_numbers_list.append(well.strip())
    well_names_unique = sorted(set(well_numbers_list))
    site_number_list.append(re.search('\d+', site).group())
    min_site_number = int(min(site_number_list))
    max_site_number = int(max(site_number_list))
    print(condition, well, site, wavelength)

# Remove duplicates (set) and sort(sorted) the well and site numbers
well_numbers = sorted(list(set(well_numbers_list)))
site_numbers = sorted(set(prepend(site_number_list, 's'))) 
        


# In[66]:


print(site_numbers)


# In[55]:


# import the required module
import os
folder_path='S:/20230802_MCF7_Ab_DMSOGWFTY-2/20230801_MCF7_Ab_DMSO/TimePoint_1'
for root, _, files in os.walk(folder_path):
    for f in files:
        complete_path=os.path.join(root,f)
        try:
            # set the size of the files you want to delete
            if os.path.getsize(complete_path) < 15 * 1024:
                print(complete_path)
                # function to delete the files
                os.remove(complete_path)
        except FileNotFound:
            print("This file does not exist"+ complete_path)
            


# In[54]:


# import the required module
import os
folder_path='S:/20230802_MCF7_Ab_DMSOGWFTY-2/20230803_MCF7_Ab_FTY720/MCF7 AB Screen- Double IF_Plate_21826/TimePoint_1'
for root, _, files in os.walk(folder_path):
    for f in files:
        complete_path=os.path.join(root,f)
        try:
            # set the size of the files you want to delete
            if os.path.getsize(complete_path) < 15 * 1024:
                print(complete_path)
                # function to delete the files
                os.remove(complete_path)
        except FileNotFound:
            print("This file does not exist"+ complete_path)
            
print(done)


# In[137]:


inputFolder = 'S:/20230802_MCF7_Ab_DMSOGWFTY-3/'

listOfDirectories = []

for inputFolder, dirs, files in os.walk(inputFolder):
    for subdir in dirs:
        print(os.path.join(inputFolder, subdir))
        listOfDirectories.append(os.path.join(inputFolder, subdir))
        
def partial(lst, query):
    return [s for s in lst if query in s]

print(listOfDirectories)
print(partial(listOfDirectories, 'TimePoint'))

rootdir = 'S:/20230802_MCF7_Ab_DMSOGWFTY-3/'

listOfDirectories = []

for rootdir, dirs, files in os.walk(rootdir):
    for subdir in dirs:
        print(os.path.join(rootdir, subdir))
        listOfDirectories.append(os.path.join(rootdir, subdir))
        
def partial(lst, query):
    return [s for s in lst if query in s]

print(listOfDirectories)
print(partial(listOfDirectories, 'TimePoint'))

listOfDirectoriesWithImages = partial(listOfDirectories, 'TimePoint')

OutputFolder = r"S:/20230918_MCF7_Ab_DMSOGWFTY-2_ImageEnhancement"
print(pathOut)

for x in listOfDirectoriesWithImages:
    pathIn = x
    pathOut = OutputFolder + '/' + os.path.basename(os.path.dirname(pathIn))



# In[ ]:





# In[ ]:





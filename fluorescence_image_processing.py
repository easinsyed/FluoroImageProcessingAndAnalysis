#!/usr/bin/env python
# coding: utf-8

"""
Fluorescence Microscopy Image Processing Script
This script processes fluorescence microscopy images, enhancing and merging them based on their channels.
"""

# Load dependecies
import os
import numpy as np
import re
import skimage 
from skimage import exposure
import glob
from skimage.restoration import rolling_ball
from skimage.morphology import white_tophat, black_tophat, disk 


def rescale_to_bits(im: np.ndarray, n_bits: int) -> np.ndarray:
    """
    Fuction that rescales an image to have values supported with a specific bit-depth (unsigned integer).
    Args:
        im (np.ndarray): Input image.
        n_bits (int): Target bit depth.
    Returns:
        np.ndarray: Rescaled image.
    """
    max_value = 2**n_bits - 1
    im = im - im.min()
    im = im * (max_value / im.max())
    im = np.round(im)
    return im


def colorize(im, color, clip_percentile=0.1):
    """
    Fuction that creates an RGB image from a single-channel image using a specific color.
    Args:
        im (np.ndarray): Input single-channel image.
        color (tuple): RGB color for tinting the image.
        clip_percentile (float): Percentile for clipping image intensities.
    Returns:
        np.ndarray: Colorized image.
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
    
    # Reshape the color
    color = np.asarray(color).reshape((1, 1, -1))
    return im_scaled * color


def subtract_background(image, radius=25, light_bg=False):
    """
    Fuction that subtracts the background from an image using morphological operations.
    Args:
        image (np.ndarray): Input image.
        radius (int): Radius of the structuring element.
        light_bg (bool): Indicates if the background is lighter than the signal.
    Returns:
        np.ndarray: Image with background subtracted.
    """
    # generate structuring element
    str_el = disk(radius)
     
    # use appropriate filter depending on the background colour
    if light_bg:
        return black_tophat(image, str_el)
    else:
        return white_tophat(image, str_el)



    inputFolder = 'S:/20230802_MCF7_Ab_DMSOGWFTY-2/'
    outputFolder = 'S:/20230918_MCF7_Ab_DMSOGWFTY-2_ImageEnhancement'

def main(inputFolder, outputFolder):
    """
    Main function to process images from input_folder and save them to output_folder.
    Args:
        input_folder (str): Path to the input directory containing images.
        output_folder (str): Path to the output directory where processed images will be saved.
    """


    listOfDirectories = []
    for inputFolder, dirs, files in os.walk(inputFolder):
        for subdir in dirs:
            print(os.path.join(inputFolder, subdir))
            listOfDirectories.append(os.path.join(inputFolder, subdir))
            
    def partial(lst, query):
        return [s for s in lst if query in s]

    print(partial(listOfDirectories, 'TimePoint'))

    listOfDirectoriesWithImages = partial(listOfDirectories, 'TimePoint')

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
                plt.imsave(save_file_name_path_rg, merge1,format='TIFF') #Saving the Red and Green colors merged file

                save_file_name_rb = "{}_Pr_{}_{}_mergeRB.tiff".format(condition, well, site)
                save_file_name_path_rb = pathOut + '/' + save_file_name_rb
                merge2 = np.clip(im_r + im_b , 0, 1)
                plt.imsave(save_file_name_path_rb, merge2,format='TIFF') #Saving the Red and Blue colors merged file

                save_file_name_gb = "{}_Pr_{}_{}_mergeGB.tiff".format(condition, well, site)
                save_file_name_path_gb = pathOut + '/' + save_file_name_gb
                merge3 = np.clip(im_g + im_b , 0, 1)
                plt.imsave(save_file_name_path_gb, merge3,format='TIFF') #Saving the Green and Blue colors merged file

                save_file_name_rgb = "{}_Pr_{}_{}_mergeRGB.tiff".format(condition, well, site)
                save_file_name_path_rgb = pathOut + '/' + save_file_name_rgb
                merge4 = np.clip(im_r + im_g + im_b , 0, 1)
                plt.imsave(save_file_name_path_rgb, merge4,format='TIFF') #Saving the Red, Green, and Blue colors merged file

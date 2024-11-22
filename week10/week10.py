#!/usr/bin/env python

import numpy
import scipy
import matplotlib.pyplot as plt
import imageio.v3
import os

#Set working directory
os.chdir("/Users/cmdb/qbb2024-answers/week10")

#Define gene names, fields, and channels
genes = ["APEX1", "PIM2", "POLR2B", "SRSF1"]
fields = ["field0", "field1"]
channels = ["DAPI", "PCNA", "nascentRNA"] 



#Part 1

# Initialize a list to hold all the images for each gene-field combination
images = []

#Test image to get image dimensions
img = imageio.v3.imread("APEX1_field0_DAPI.tif").astype(numpy.uint16)

for gene in genes:
    for field in fields:
        #Initialize an array to hold the channels for gene/field combo
        rgbimg = numpy.zeros((img.shape[0], img.shape[1], 3), numpy.uint16)  
        #Loop through each channel
        for i, channel in enumerate(channels):
            #Construct the filename for the current channel image
            filename = f"{gene}_{field}_{channel}.tif"
            #Load the image
            img = imageio.v3.imread(filename)
            #Store the normalized image in the rgbimg array at the correct channel position
            rgbimg[:, :, i] = img
        # Normalize the combined RGB image across all channels
        rgbimg = rgbimg.astype(numpy.float32)
        #Append the image to the list of images
        images.append(rgbimg)



#Part 2.1 and 2.2

#From live code:
def find_labels(mask):
    # Set initial label
    l = 0
    # Create array to hold labels
    labels = numpy.zeros(mask.shape, numpy.int32)
    # Create list to keep track of label associations
    equivalence = [0]
    # Check upper-left corner
    if mask[0, 0]:
        l += 1
        equivalence.append(l)
        labels[0, 0] = l
    # For each non-zero column in row 0, check back pixel label
    for y in range(1, mask.shape[1]):
        if mask[0, y]:
            if mask[0, y - 1]:
                # If back pixel has a label, use same label
                labels[0, y] = equivalence[labels[0, y - 1]]
            else:
                # Add new label
                l += 1
                equivalence.append(l)
                labels[0, y] = l
    # For each non-zero row
    for x in range(1, mask.shape[0]):
        # Check left-most column, up  and up-right pixels
        if mask[x, 0]:
            if mask[x - 1, 0]:
                # If up pixel has label, use that label
                labels[x, 0] = equivalence[labels[x - 1, 0]]
            elif mask[x - 1, 1]:
                # If up-right pixel has label, use that label
                labels[x, 0] = equivalence[labels[x - 1, 1]]
            else:
                # Add new label
                l += 1
                equivalence.append(l)
                labels[x, 0] = l
        # For each non-zero column except last in nonzero rows, check up, up-right, up-right, up-left, left pixels
        for y in range(1, mask.shape[1] - 1):
            if mask[x, y]:
                if mask[x - 1, y]:
                    # If up pixel has label, use that label
                    labels[x, y] = equivalence[labels[x - 1, y]]
                elif mask[x - 1, y + 1]:
                    # If not up but up-right pixel has label, need to update equivalence table
                    if mask[x - 1, y - 1]:
                        # If up-left pixel has label, relabel up-right equivalence, up-left equivalence, and self with smallest label
                        labels[x, y] = min(equivalence[labels[x - 1, y - 1]], equivalence[labels[x - 1, y + 1]])
                        equivalence[labels[x - 1, y - 1]] = labels[x, y]
                        equivalence[labels[x - 1, y + 1]] = labels[x, y]
                    elif mask[x, y - 1]:
                        # If left pixel has label, relabel up-right equivalence, left equivalence, and self with smallest label
                        labels[x, y] = min(equivalence[labels[x, y - 1]], equivalence[labels[x - 1, y + 1]])
                        equivalence[labels[x, y - 1]] = labels[x, y]
                        equivalence[labels[x - 1, y + 1]] = labels[x, y]
                    else:
                        # If neither up-left or left pixels are labeled, use up-right equivalence label
                        labels[x, y] = equivalence[labels[x - 1, y + 1]]
                elif mask[x - 1, y - 1]:
                    # If not up, or up-right pixels have labels but up-left does, use that equivalence label
                    labels[x, y] = equivalence[labels[x - 1, y - 1]]
                elif mask[x, y - 1]:
                    # If not up, up-right, or up-left pixels have labels but left does, use that equivalence label
                    labels[x, y] = equivalence[labels[x, y - 1]]
                else:
                    # Otherwise, add new label
                    l += 1
                    equivalence.append(l)
                    labels[x, y] = l
        # Check last pixel in row
        if mask[x, -1]:
            if mask[x - 1, -1]:
                # if up pixel is labeled use that equivalence label 
                labels[x, -1] = equivalence[labels[x - 1, -1]]
            elif mask[x - 1, -2]:
                # if not up but up-left pixel is labeled use that equivalence label 
                labels[x, -1] = equivalence[labels[x - 1, -2]]
            elif mask[x, -2]:
                # if not up or up-left but left pixel is labeled use that equivalence label 
                labels[x, -1] = equivalence[labels[x, -2]]
            else:
                # Otherwise, add new label
                l += 1
                equivalence.append(l)
                labels[x, -1] = l
    equivalence = numpy.array(equivalence)
    # Go backwards through all labels
    for i in range(1, len(equivalence))[::-1]:
        # Convert labels to the lowest value in the set associated with a single object
        labels[numpy.where(labels == i)] = equivalence[i]
    # Get set of unique labels
    ulabels = numpy.unique(labels)
    for i, j in enumerate(ulabels):
        # Relabel so labels span 1 to # of labels
        labels[numpy.where(labels == j)] = i
    return labels


#Part 2.3

def filter_by_size(labels, minsize, maxsize):
    #Find label sizes
    sizes = numpy.bincount(labels.ravel())
    #Iterate through labels, skipping background
    for i in range(1, sizes.shape[0]):
        #If # of pixels falls outsize the cutoff range, relabel as background
        if sizes[i] < minsize or sizes[i] > maxsize:
            # Find all pixels for label
            where = numpy.where(labels == i)
            labels[where] = 0
    #Get set of unique labels
    ulabels = numpy.unique(labels)
    for i, j in enumerate(ulabels):
        #Relabel so labels span 1 to # of labels
        labels[numpy.where(labels == j)] = i
    return labels



#Part 3
for i in range(8): 
    current_img = images[i]
    DAPI_img = current_img[:,:,0]
    nascent_img = current_img[:,:,1]
    PCNA_img = current_img[:,:,2]

    mask = DAPI_img >= numpy.mean(DAPI_img)

    label_array = find_labels(mask)
    label_array = filter_by_size(label_array, 100, 100000000)

    sizes = numpy.bincount(label_array.ravel())
    sizes = sizes[1:]
    mean_size = numpy.mean(sizes)
    sd_size = numpy.std(sizes)

    label_array = filter_by_size(label_array, mean_size - sd_size, mean_size + sd_size)
    num_nuclei = numpy.amax(label_array)
    num_nuclei = num_nuclei + 1

    for j in range(1, num_nuclei):
        where = numpy.where(label_array == j)
        nascent_signal = numpy.mean(nascent_img[where])
        PCNA_signal = numpy.mean(PCNA_img[where])
        log2ratio = numpy.log2(nascent_signal / PCNA_signal)
        if i in [0, 1]: 
            Gene = "APEX1"
        if i in [2, 3]:
            Gene = "PIM2"
        if i in [4,5]:
            Gene = "POLR2B"
        if i in [6,7]: 
            Gene = "SRSF1"
        print(Gene, nascent_signal, PCNA_signal, log2ratio, sep = "\t")

#To save as .txt, wrote the following in terminal:
#python week10.py > nuclei_signal.txt
#!/usr/bin/env python

#import cv2
import os
from itertools import repeat
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import skimage
from datetime import datetime
from scipy import ndimage as ndi
from skimage import io, viewer, color, feature, filters, measure, morphology, segmentation, util
#change this later to be a loop
RFPimg = '/Users/claytonwandishin/December 14 RT glow run/Image Segmentation custom/BCD_images/B2_02_1_1_RFP_001.tif'

#GFPimg = '/Users/claytonwandishin/December 14 RT glow run/Image Segmentation custom/BCD_images/B2_02_2_1_GFP_001.tif'
momentum_log = pd.read_csv('/Users/claytonwandishin/December 14 RT glow run/Image Segmentation custom/MomentumLog_CycloDRCandHES1ko_CMW_Quaranta_20210611.csv')
read_log = momentum_log.loc[momentum_log['Get'] == 'Read']
read_log = read_log.rename(columns={'1':'plate','Unnamed: 13':'barcode'})

#change this to be by barcode later not mcf
mcf_reads = read_log.loc[read_log['barcode'] == 'VU60024H']
mcf_img_reads = mcf_reads.loc[mcf_reads['Stack 1:Nest 1']=='CycloDRCandHES1ko_100msRFP_10msGFP_BF_CMW_Quaranta_20210611.prt']
mcf_img_DT_dic = mcf_img_reads.iloc[:, [11,13]]
mcf_img_DT_dic = mcf_img_DT_dic.rename({mcf_img_DT_dic.columns[0]: 'DateTime_Image'}, axis='columns')
mcf_img_DT_dic = mcf_img_DT_dic.reset_index()
mcf_img_DT_dic = mcf_img_DT_dic.drop(columns=['index'])

mcf_img_DT_dic['DateTime']=pd.to_datetime(mcf_img_DT_dic['DateTime_Image'])
well_list=[]
well_columns=[*range(2,12)]
well_rows=['B','C','D']

for x in well_rows:
    for y in well_columns:
        well_list.append(x+str(y))
mcf10a_reads=[]
mcf10a_read_freq = [*range(2,47,3)]

#this bit extracts the timestamp from the filename
file_list = os.listdir('/Volumes/Hackathon1/20210611_images')
fl_ts = pd.DataFrame(file_list)
timestamp_list =[]
cyclofilesonly=[]
for k in mcf10a_read_freq:
    cyclofilesonly.append(file_list[k])
cyclosubfiles =[]
mcf10a_explicit_file_list=[]
cyclo_img_file_list=[]
# all the files that need to be segmented contain 'RFP' and each individual folder is under the format "/Volumes/Hackathon1/20210611_images/210612024031_Experiment12/!PLATE_ID!_210612024031/" and the variable "well_list" contains all the mcf10a only reads
for x in cyclofilesonly:
    dtslice=x[0:12]
    cyclopath="/Volumes/Hackathon1/20210611_images/"+x+"/!PLATE_ID!_"+dtslice+"/"
    cyclosubfiles.append(cyclopath)
    
#this next loop selects only the files that contain RFP
for j in cyclosubfiles:
    img_sub_files=os.listdir(j)
    for m in img_sub_files:
        subfilepath = j+m
        if 'RFP' in subfilepath:
            cyclo_img_file_list.append(subfilepath)
#this last loop should select only the relevant wells for mcf10a
for x in cyclo_img_file_list:
    for m in well_list:
        if m in x:
            mcf10a_explicit_file_list.append(x)
excludedwells = ['B20','B21','B22','B23','C20','C21','C22','C23','D20','D21','D22','D23']
notwant=[]
for j in excludedwells:
    for m in mcf10a_explicit_file_list:
        if j in m:
            notwant.append(m)
for x in notwant:
    mcf10a_explicit_file_list.remove(x)

        

    
#mcf10a_explicit_file_list now contains the full path to every image we want to segment
                    
    
    

for i, row in fl_ts.iterrows():
    timestamp_name = fl_ts.at[i,0]
    timestamp = timestamp_name[0:10]
    timestamp_list.append(timestamp)

fl_ts['timestamp'] = timestamp_list
dt_list = []
for x in timestamp_list:
    datetime = pd.to_datetime(x, format='%y%m%d%H%M')
    dt_list.append(datetime)
fl_ts['DateTime'] = dt_list


mcf10a_cyclo=pd.merge(mcf_img_DT_dic,fl_ts)

#this is the list of relevant experiment files for the cyclohexamide curves (all cell lines)
mcf10exp_files = mcf10a_cyclo[0].copy()
mcf10_df = pd.DataFrame()
mcf10_df['file_name']=mcf10a_explicit_file_list.copy()
CellCountList =[]
for i in mcf10a_explicit_file_list:
    #loads the image into skimage
    imageRFP = skimage.io.imread(fname=i)
    #viewer = skimage.viewer.ImageViewer(imageRFP)
    #viewer.show()
    #imageGFP = skimage.io.imread(fname=GFPimg)
    #viewer = skimage.viewer.ImageViewer(imageGFP)
    #viewer.show()

    #blurs the image before binary so it's a smoother mask
    gblur = skimage.filters.gaussian(imageRFP, sigma=3)

    image = gblur

    thresholds = filters.threshold_multiotsu(image, classes=2)
    regions = np.digitize(image, bins=thresholds)


    '''
    fig, ax = plt.subplots(ncols=2, figsize=(10, 5))
    ax[0].imshow(image)
    ax[0].set_title('Original')
    ax[0].axis('off')
    ax[1].imshow(regions)
    ax[1].set_title('Multi-Otsu thresholding')
    ax[1].axis('off')
    plt.show()
    '''
    cells = image > thresholds[0]
    nuclei = measure.label(cells)

    distance = ndi.distance_transform_edt(cells)

    local_max_coords = feature.peak_local_max(distance, min_distance=16)
    local_max_mask = np.zeros(distance.shape, dtype=bool)
    local_max_mask[tuple(local_max_coords.T)] = True
    markers = measure.label(local_max_mask)

    segmented_cells = segmentation.watershed(-distance, markers, mask=cells)

    fig, ax = plt.subplots(ncols=2, figsize=(10, 5))
    ax[0].imshow(cells, cmap='gray')
    ax[0].set_title('Overlapping nuclei')
    ax[0].axis('off')
    ax[1].imshow(color.label2rgb(segmented_cells, bg_label=0))
    ax[1].set_title(i)
    ax[1].axis('off')
    plt.show()
    nuclei_count = markers.max()-1
    CellCountList.append(nuclei_count)
    print(nuclei_count)
mcf10_df['cell_count']=CellCountList.copy()
mcf10_df.to_csv('/Users/claytonwandishin/December 14 RT glow run/Image Segmentation custom/mcf10a_df.csv')
'''
#runs an OTSU adaptive thresholding algo and converts to binary by comparing blurred pixels values to determined threshold and assigning 0 or 1 based on the gate
thresholdimg = skimage.filters.threshold_otsu(gblur)
mask = gblur > thresholdimg

# Now we want to separate the two objects in image
# Generate the markers as local maxima of the distance to the background
distance = ndi.distance_transform_edt(mask)
coords = peak_local_max(distance, footprint=np.ones((3, 3)), labels=mask)
mask2 = np.zeros(distance.shape, dtype=bool)
mask2[tuple(coords.T)] = True
markers, _ = ndi.label(mask2)
labels = watershed(-distance, markers, mask=mask)


# ...but then you have to clean up the tiny intersections between coins
regions = regionprops(labels)
regions = [r for r in regions if r.area > 140]

print('Nuclei:'+str(len(regions) - 1))





fig, axes = plt.subplots(ncols=3, figsize=(9, 3), sharex=True, sharey=True)
ax = axes.ravel()

ax[0].imshow(mask, cmap=plt.cm.gray, interpolation='nearest')
ax[0].set_title('Overlapping objects')
ax[1].imshow(-distance, cmap=plt.cm.gray, interpolation='nearest')
ax[1].set_title('Distances')
ax[2].imshow(labels, cmap=plt.cm.nipy_spectral, interpolation='nearest')
ax[2].set_title('Separated objects')


for a in ax:
    a.set_axis_off()

fig.tight_layout()
plt.show()


#print(datetime.now() - startTime)
'''
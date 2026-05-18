#### Script for Visium HD
## bin2cell with DAPI

#### Preare DAPI mask
# 1. Open tiff (sourece image from cytassist) in napari and save it as tiff.
# 2. Open napari_tiff in GIMP.
# 3. Load DAPI tiles as mask 
# 4. Rotate 90º, mirror vertically and rotate -1.2ª
# 5. Align with the H&E, remove border areas from one of the tiles. 
# 4. Remove H&E mask and save the alineated DAPI as tiff keeping the size of source tiff. Add a black bounground image. 


### Load python
##conda activate /datos_2/Spatial/VisiumHD_bin2cell_environment/ 
##python

## Import packages
import matplotlib.pyplot as plt
import scanpy as sc
import numpy as np
import pandas as pd
import os
import tifffile
import bin2cell as b2c
import imageio
from scipy.spatial import cKDTree
import json
from skimage.transform import rescale
from skimage.transform import resize
import imagecodecs
from matplotlib.backends.backend_pdf import PdfPages


### Set the path to file
os.chdir("/datos_2/Spatial/Bone_GP/Analysis_DAPI/")

path_D1 = "/datos_2/Spatial/Bone_GP/data/Counts_D1/outs/binned_outputs/square_002um/"
spaceranger_image_path_D1 = "/datos_2/Spatial/Bone_GP/data/Counts_D1/outs/binned_outputs/square_002um/spatial"

path_A1 = "/datos_2/Spatial/Bone_GP/data/Counts_A1/outs/binned_outputs/square_002um/"
spaceranger_image_path_A1 = "/datos_2/Spatial/Bone_GP/data/Counts_A1/outs/binned_outputs/square_002um/spatial"


#the image you used for --image of spaceranger, that's the one the spatial coordinates are based on
source_image_HE_D1 = "/datos_2/Spatial/Bone_GP/data/HumanBone.250304_VisiumHD_C16_Ir16_DQ12_LA-Spot000001.tif"
source_image_DAPI_D1 = "/datos_2/Spatial/Bone_GP/data/Visium_HD_DQ1295_napari_DQ1295_opacity50.tif"
source_image_path_D1 = "/datos_2/Spatial/Bone_GP/data/VisiumHD_HE_DAPI_DQ9512.tif"

source_image_HE_A1 = "/datos_2/Spatial/Bone_GP/data/HumanBone.250304_VisiumHD_C15_Ir15_CC5215_rep_LA-Spot000001.tif"
source_image_DAPI_A1 = "/datos_2/Spatial/Bone_GP/data/Visium_HD_A1_napari_opacity50.tif"
source_image_path_A1 = "/datos_2/Spatial/Bone_GP/data/VisiumHD_HE_DAPI.tif"

sc.settings.figdir = "/datos_2/Spatial/Bone_GP/Analysis_DAPI"


############################
####### Run Bin2Cell #######
############################

## Generate tiff image with the coordinates from expression

mpp = 0.3
b2c.scaled_if_image(adata_DQ9512, channel = 3, mpp=mpp, save_path="stardist/dapi_DQ9512.tiff")
b2c.scaled_if_image(adata_CC5215, channel = 3, mpp=mpp, save_path="stardist/dapi_CC5215.tiff")

b2c.scaled_he_image(adata_he_DQ9512, mpp=mpp, save_path="stardist/he_DQ9512.tiff",
    spatial_cropped_key = "spatial_cropped_150he_buffer", img_key = "0.3_mpp_150he_buffer")
b2c.scaled_he_image(adata_he_CC5215, mpp=mpp, save_path="stardist/he_CC5215.tiff",
    spatial_cropped_key = "spatial_cropped_150he_buffer", img_key = "0.3_mpp_150he_buffer")

# Remove stripes from figure, mostly due to scanning
b2c.destripe(adata_DQ9512)
b2c.destripe(adata_CC5215)

## Stardist nuclei segmentation of full image
b2c.stardist(image_path="stardist/dapi_DQ9512.tiff", labels_npz_path="stardist/dapi_DQ9512.npz", 
    stardist_model="2D_versatile_fluo", prob_thresh=0.05)
b2c.stardist(image_path="stardist/dapi_CC5215.tiff", labels_npz_path="stardist/dapi_CC5215.npz", 
    stardist_model="2D_versatile_fluo", prob_thresh=0.05)


## Add nuclei segmentation to adata
b2c.insert_labels(adata_DQ9512, labels_npz_path="stardist/dapi_DQ9512.npz", basis="spatial", 
    spatial_key="spatial_cropped_150_buffer", mpp=mpp, labels_key="labels_DAPI")
b2c.insert_labels(adata_CC5215, labels_npz_path="stardist/dapi_CC5215.npz", basis="spatial", 
    spatial_key="spatial_cropped_150_buffer", mpp=mpp, labels_key="labels_DAPI")

### Extend nuclei labels to cells
b2c.expand_labels(adata_DQ9512, labels_key='labels_DAPI', expanded_labels_key="labels_DAPI_expanded", algorithm="volume_ratio")
b2c.expand_labels(adata_CC5215, labels_key='labels_DAPI', expanded_labels_key="labels_DAPI_expanded", algorithm="volume_ratio")

#the labels obs are integers, 0 means unassigned
adata_DQ9512_noNuc = adata_DQ9512[adata_DQ9512.obs['labels_DAPI_expanded']>0]
adata_DQ9512_noNuc.obs['labels_DAPI_expanded'] = adata_DQ9512_noNuc.obs['labels_DAPI_expanded'].astype(str)

sc.pl.spatial(adata_DQ9512_noNuc, color=[None, "labels_DAPI_expanded", "n_counts_adjusted"], legend_loc=None, save="_cellSegmentation_DQ9512.pdf")

adata_CC5215_noNuc = adata_CC5215[adata_CC5215.obs['labels_DAPI_expanded']>0]
adata_CC5215_noNuc.obs['labels_DAPI_expanded'] = adata_CC5215_noNuc.obs['labels_DAPI_expanded'].astype(str)

sc.pl.spatial(adata_CC5215_noNuc, color=[None, "labels_DAPI_expanded", "n_counts_adjusted"], legend_loc=None, save="_cellSegmentation_CC5215.pdf")

#the label viewer wants a crop of the processed image
#get the corresponding coordinates spanning the subset object
crop_DQ9512 = b2c.get_crop(adata_DQ9512, basis="spatial", spatial_key="spatial_cropped_150_buffer", mpp=mpp)
crop_CC5215 = b2c.get_crop(adata_CC5215, basis="spatial", spatial_key="spatial_cropped_150_buffer", mpp=mpp)

rendered_DQ9512 = b2c.view_labels(image_path="stardist/dapi_DQ9512.tiff", labels_npz_path="stardist/dapi_DQ9512.npz", crop=crop_DQ9512)
rendered_CC5215 = b2c.view_labels(image_path="stardist/dapi_CC5215.tiff", labels_npz_path="stardist/dapi_CC5215.npz", crop=crop_CC5215)

plt.figure(figsize=(200, 200))
plt.imshow(rendered_DQ9512)
plt.savefig("NucleiSegment_DQ9512_dapi.pdf", format="pdf", bbox_inches='tight')  
plt.show()
plt.figure(figsize=(200, 200))
plt.imshow(rendered_CC5215)
plt.savefig("NucleiSegment_CC5215_dapi.pdf", format="pdf", bbox_inches='tight')  
plt.show()

rendered_he_DQ9512 = b2c.view_labels(image_path="stardist/he_DQ9512.tiff", labels_npz_path="stardist/dapi_DQ9512.npz", crop=crop_DQ9512)
rendered_he_CC5215 = b2c.view_labels(image_path="stardist/he_DQ9512.tiff", labels_npz_path="stardist/dapi_DQ9512.npz", crop=crop_CC5215)

plt.figure(figsize=(200, 200))
plt.imshow(rendered_he_DQ9512)
plt.savefig("NucleiSegment_DQ9512_he.pdf", format="pdf", bbox_inches='tight')  
plt.show()
plt.figure(figsize=(200, 200))
plt.imshow(rendered_he_CC5215)
plt.savefig("NucleiSegment_CC5215_he.pdf", format="pdf", bbox_inches='tight')  
plt.show()


### Group bins into cells
cellData_DQ9512 = b2c.bin_to_cell(adata_DQ9512, labels_key="labels_DAPI_expanded", spatial_keys=["spatial", "spatial_cropped_150_buffer"])
cellData_CC5215 = b2c.bin_to_cell(adata_CC5215, labels_key="labels_DAPI_expanded", spatial_keys=["spatial", "spatial_cropped_150_buffer"])

#### Add labels into he objects
adata_he_DQ9512.obs = adata_he_DQ9512.obs.join(adata_DQ9512.obs['labels_DAPI'])
adata_he_DQ9512.obs = adata_he_DQ9512.obs.join(adata_DQ9512.obs['labels_DAPI_expanded']).fillna({'labels_DAPI_expanded': 0})
adata_he_DQ9512.obs['labels_DAPI_expanded'] = adata_he_DQ9512.obs['labels_DAPI_expanded'].astype(int)
cellData_he_DQ9512 = b2c.bin_to_cell(adata_he_DQ9512, labels_key="labels_DAPI_expanded", spatial_keys=["spatial"])

adata_he_CC5215.obs = adata_he_CC5215.obs.join(adata_CC5215.obs['labels_DAPI'])
adata_he_CC5215.obs = adata_he_CC5215.obs.join(adata_CC5215.obs['labels_DAPI_expanded']).fillna({'labels_DAPI_expanded': 0})
adata_he_CC5215.obs['labels_DAPI_expanded'] = adata_he_CC5215.obs['labels_DAPI_expanded'].astype(int)
cellData_he_CC5215 = b2c.bin_to_cell(adata_he_CC5215, labels_key="labels_DAPI_expanded", spatial_keys=["spatial"])

# Total UMI counts per cell
cellData_DQ9512.obs["n_counts"] = cellData_DQ9512.X.sum(axis=1).A1 
cellData_he_DQ9512.obs["n_counts"] = cellData_he_DQ9512.X.sum(axis=1).A1 
cellData_CC5215.obs["n_counts"] = cellData_CC5215.X.sum(axis=1).A1 
cellData_he_CC5215.obs["n_counts"] = cellData_he_CC5215.X.sum(axis=1).A1 

# Total Features counts per cell
cellData_DQ9512.obs["n_features"] = (cellData_DQ9512.X > 0).sum(axis=1).A1
cellData_he_DQ9512.obs["n_features"] = (cellData_he_DQ9512.X > 0).sum(axis=1).A1
cellData_CC5215.obs["n_features"] = (cellData_CC5215.X > 0).sum(axis=1).A1
cellData_he_CC5215.obs["n_features"] = (cellData_he_CC5215.X > 0).sum(axis=1).A1


### Save adata with single-cell resolution
cellData_DQ9512.write("cells_DQ9512_dapi.h5ad")
cellData_he_DQ9512.write("cells_DQ9512_he.h5ad")
cellData_CC5215.write("cells_CC5215_dapi.h5ad")
cellData_he_CC5215.write("cells_CC5215_he.h5ad")


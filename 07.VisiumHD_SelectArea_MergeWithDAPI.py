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


########################################################
####### Create image with HE & DAPI for Bin2Cell #######
########################################################

# Load images
he_img_D1 = tifffile.imread(source_image_HE_D1)
with tifffile.TiffFile(source_image_DAPI_D1) as tif:
    dapi_img_D1 = tif.pages[0].asarray()

he_img_A1 = tifffile.imread(source_image_HE_A1)
with tifffile.TiffFile(source_image_DAPI_A1) as tif:
    dapi_img_A1 = tif.pages[0].asarray()

### take only the DAPI channel
dapi_img_D1 = dapi_img_D1[..., 2]
dapi_img_A1 = dapi_img_A1[..., 2]

## Ensuer DAPI and HE are same size
if he_img_D1.shape[:2] != dapi_img_D1.shape:
    raise ValueError("H&E and DAPI dimensions do not match")

if he_img_A1.shape[:2] != dapi_img_A1.shape:
    raise ValueError("H&E and DAPI dimensions do not match")

# Convert H&E RGB into 3 channels
he_channels_A1 = np.moveaxis(he_img_A1, -1, 0)  # shape: (3, H, W)
he_channels_D1 = np.moveaxis(he_img_D1, -1, 0)  # shape: (3, H, W)

# Add DAPI as channel 4
combined_A1 = np.concatenate([he_channels_A1, dapi_img_A1[None, ...]], axis=0)  # shape: (4, H, W)
combined_D1 = np.concatenate([he_channels_D1, dapi_img_D1[None, ...]], axis=0)  # shape: (4, H, W)


# Open original metadata
with tifffile.TiffFile(source_image_HE_A1) as tif:
    metadata_A1 = tif.pages[0].tags  # recoge tags TIFF estándar
    ome_metadata_A1 = tif.ome_metadata  # si es OME-TIFF

with tifffile.TiffFile(source_image_HE_D1) as tif:
    metadata_D1 = tif.pages[0].tags  # recoge tags TIFF estándar
    ome_metadata_D1 = tif.ome_metadata  # si es OME-TIFF

# Save each channel as a separate TIFF page (multi-page TIFF)
tifffile.imwrite(
    source_image_path_D1,
    combined_D1,
    photometric='minisblack',  # Para multicanal no RGB
    metadata={'axes': 'CYX'},
    description=ome_metadata_D1 if ome_metadata_D1 else None
)
# Save each channel as a separate TIFF page (multi-page TIFF)
tifffile.imwrite(
    source_image_path_A1,
    combined_A1,
    photometric='minisblack',  # Para multicanal no RGB
    metadata={'axes': 'CYX'},
    description=ome_metadata_A1 if ome_metadata_A1 else None
)




#####################################
####### Load data into scanpy #######
#####################################


## Load visium data into adata
adata_D1 = b2c.read_visium(path_D1, source_image_path = source_image_path_D1,
    spaceranger_image_path = spaceranger_image_path_D1)
adata_D1.var_names_make_unique()
adata_D1

mask_D1 = ((adata_D1.obs['array_row'] >= 1500) & 
        (adata_D1.obs['array_row'] <= 3000) & 
        (adata_D1.obs['array_col'] >= 1500) & 
        (adata_D1.obs['array_col'] <= 4000))
### Subset to take spots within samples
adata_DQ9512 = adata_D1[mask_D1].copy()


## Load visium data into adata
adata_A1 = b2c.read_visium(path_A1, source_image_path = source_image_path_A1,
    spaceranger_image_path = spaceranger_image_path_A1)
adata_A1.var_names_make_unique()
adata_A1

mask_paths_Sample = {
    "Ctl": "/datos_2/Spatial/Bone_GP/data/C15_Ctl.png",
    "Irr": "/datos_2/Spatial/Bone_GP/data/C15_Irr.png",
    "CC5215": "/datos_2/Spatial/Bone_GP/data/C15_CC5215.png"
}

# Get coordiantes from adata
coords_A1 = adata_A1.obsm['spatial'].astype(int)
adata_A1_tree = cKDTree(coords_A1)
radius = 2
area_A1 = np.full(len(adata_A1), "notDefined", dtype=object)
sf = adata_A1.uns['spatial']['H1-PYFPXX3_A1']['scalefactors']['tissue_hires_scalef']

# Extract black coordinates and label
sample_A1 = np.full(len(adata_A1), "notDefined", dtype=object)
for region, path in mask_paths_Sample.items():
    img = imageio.v3.imread(path, mode="L")
    
    img_rescaled = rescale(img, scale=(1/sf, 1/sf), anti_aliasing=False,
                           preserve_range=True, order=0).astype(np.uint8)
    
    black_pixels = np.argwhere(img_rescaled < 200)
    region_tree = cKDTree(black_pixels[:, [1, 0]])  # invertir orden (x, y)
    
    matches = adata_A1_tree.query_ball_tree(region_tree, r=radius)
    matched_indices = set(i for i, m in enumerate(matches) if len(m) > 0)
    
    for i in matched_indices:
        sample_A1[i] = region


# Añadir etiquetas al objeto
adata_A1.obs['sample'] = sample_A1
mask_A1 = ((adata_A1.obs['array_row'] >= 1000) & 
        (adata_A1.obs['array_row'] <= 4000) & 
        (adata_A1.obs['array_col'] >= 000) & 
        (adata_A1.obs['array_col'] <= 1600))

### Subset to take spots within samples
adata_CC5215 = adata_A1[mask_A1].copy()
adata_CC5215 = adata_CC5215[adata_CC5215.obs['sample'] == "CC5215"].copy()

# Filter lowly expressed genes, remove areas without counts and leave spots within samples
sc.pp.filter_genes(adata_DQ9512, min_cells=3)
sc.pp.filter_cells(adata_DQ9512, min_counts=1)

sc.pp.filter_genes(adata_CC5215, min_cells=3)
sc.pp.filter_cells(adata_CC5215, min_counts=1)

# Normalize and log-transform ----
adata_DQ9512.layers["counts"] = adata_DQ9512.X.copy()
sc.pp.normalize_total(adata_DQ9512, inplace=True)
sc.pp.log1p(adata_DQ9512)
adata_DQ9512.obs['normalized_counts'] = adata_DQ9512.X.sum(axis=1)

adata_CC5215.layers["counts"] = adata_CC5215.X.copy()
sc.pp.normalize_total(adata_CC5215, inplace=True)
sc.pp.log1p(adata_CC5215)
adata_CC5215.obs['normalized_counts'] = adata_CC5215.X.sum(axis=1)

### Object to obtain the H&E image
adata_he_D1 = b2c.read_visium(path_D1, source_image_path = source_image_HE_D1,
    spaceranger_image_path = spaceranger_image_path_D1)
adata_he_D1.var_names_make_unique()
adata_he_D1
### Subset to take spots within samples
adata_he_DQ9512 = adata_he_D1[mask_D1].copy()


### Object to obtain the H&E image
adata_he_A1 = b2c.read_visium(path_A1, source_image_path = source_image_HE_A1,
    spaceranger_image_path = spaceranger_image_path_A1)
adata_he_A1.var_names_make_unique()
adata_he_A1
### Subset to take spots within samples
adata_he_CC5215 = adata_he_A1[mask_A1].copy()

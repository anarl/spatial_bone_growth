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


######################################################
####### Add area labels based on mask (Mahtab) #######
######################################################


## Import packages
import imageio
from scipy.ndimage import zoom
from skimage.transform import rescale
import matplotlib.pyplot as plt
import scanpy as sc
import numpy as np
import pandas as pd
import os
from scipy.spatial import cKDTree


### Set the path to file
os.chdir("/datos_2/Spatial/Bone_GP/Analysis_DAPI/")
sc.settings.figdir = "/datos_2/Spatial/Bone_GP/Analysis_DAPI"

cellData_D1 = sc.read_h5ad("cells_DQ9512_dapi.h5ad")
cellData_A1 = sc.read_h5ad("cells_CC5215_dapi.h5ad")


# Load binary image (black = ROI)
mask_paths_Area_A1 = {
    "PZ": "/datos_2/Spatial/Bone_GP/data/C15_PZ.png",
    "HZ": "/datos_2/Spatial/Bone_GP/data/C15_HZ.png",
    "RZ": "/datos_2/Spatial/Bone_GP/data/C15_RZ.png",
    "SOC": "/datos_2/Spatial/Bone_GP/data/C15_SOC.png",
    "POC": "/datos_2/Spatial/Bone_GP/data//C15_POC.png",
    "BM": "/datos_2/Spatial/Bone_GP/data/C15_BM.png"
}
# Load binary image (black = ROI)
mask_paths_Area_D1 = {
    "PZ": "/datos_2/Spatial/Bone_GP/data/C16_PZ.png",
    "HZ": "/datos_2/Spatial/Bone_GP/data/C16_HZ.png",
    "RZ": "/datos_2/Spatial/Bone_GP/data/C16_RZ.png",
    "SOC": "/datos_2/Spatial/Bone_GP/data/C16_SOC.png",
    "POC": "/datos_2/Spatial/Bone_GP/data//C16_POC.png",
    "BM": "/datos_2/Spatial/Bone_GP/data/C16_BM.png"
}

# Get coordiantes from adata
coords_D1 = cellData_D1.obsm['spatial'].astype(int)
adata_D1_tree = cKDTree(coords_D1)
radius = 2
area_D1 = np.full(len(cellData_D1), "notDefined", dtype=object)

# Extract black coordinates and label
for region, path in mask_paths_Area_D1.items():
    img = imageio.v3.imread(path, mode="L")
    
    img_rescaled = rescale(img, scale=(1/sf, 1/sf), anti_aliasing=False,
                           preserve_range=True, order=0).astype(np.uint8)
    
    black_pixels = np.argwhere(img_rescaled < 200)
    region_tree = cKDTree(black_pixels[:, [1, 0]])  
    
    matches = adata_D1_tree.query_ball_tree(region_tree, r=radius)
    matched_indices = set(i for i, m in enumerate(matches) if len(m) > 0)
    
    for i in matched_indices:
        area_D1[i] = region

# Añadir etiquetas al objeto
cellData_D1.obs['area'] = area_D1


# Get coordiantes from adata
coords_A1 = cellData_A1.obsm['spatial'].astype(int)
adata_A1_tree = cKDTree(coords_A1)
area_A1 = np.full(len(cellData_A1), "notDefined", dtype=object)

# Extract black coordinates and label
for region, path in mask_paths_Area_A1.items():
    img = imageio.v3.imread(path, mode="L")
    
    img_rescaled = rescale(img, scale=(1/sf, 1/sf), anti_aliasing=False,
                           preserve_range=True, order=0).astype(np.uint8)
    
    black_pixels = np.argwhere(img_rescaled < 200)
    region_tree = cKDTree(black_pixels[:, [1, 0]])  
    
    matches = adata_A1_tree.query_ball_tree(region_tree, r=radius)
    matched_indices = set(i for i, m in enumerate(matches) if len(m) > 0)
    
    for i in matched_indices:
        area_A1[i] = region

# Añadir etiquetas al objeto
cellData_A1.obs['area'] = area_A1

# Visualization
custom_palette = {"RZ":"#A72326","PZ":"#F1DC2D","HZ":"#2DA297","SOC":"#2F75B0","BM":"#7A3803","POC":"#AEF359","notDefined":"#311432"}
sc.pl.spatial(cellData_D1, color=['area'], palette=custom_palette, save="_area_CellDQ9512.pdf")
sc.pl.spatial(cellData_A1, color=['area'], palette=custom_palette, save="_area_CellCC5215.pdf")



#####################################################
#### Clean objects, keep only bins within labels ####
#####################################################

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Asegúrate de que los factores están correctamente tipados
cellData_DQ9512_clean = cellData_D1.copy()
cellData_DQ9512_clean.obs['area'] = cellData_DQ9512_clean.obs['area'].astype('category')
cellData_DQ9512_clean = cellData_DQ9512_clean[cellData_DQ9512_clean.obs["area"].isin(["RZ","PZ","HZ","SOC"])].copy()
cellData_DQ9512_clean = cellData_DQ9512_clean[(cellData_DQ9512_clean.obs["bin_count"] > 4) & (cellData_DQ9512_clean.obs["bin_count"] < 200)].copy()
sc.pp.filter_genes(cellData_DQ9512_clean, min_cells=1)
sc.pp.filter_cells(cellData_DQ9512_clean, min_counts=1)

cellData_CC5215_clean = cellData_A1.copy()
cellData_CC5215_clean.obs['area'] = cellData_CC5215_clean.obs['area'].astype('category')
cellData_CC5215_clean = cellData_CC5215_clean[cellData_CC5215_clean.obs["area"].isin(["RZ","PZ","HZ","SOC"])].copy()
cellData_CC5215_clean = cellData_CC5215_clean[(cellData_CC5215_clean.obs["bin_count"] > 4) & (cellData_CC5215_clean.obs["bin_count"] < 200)].copy()
sc.pp.filter_genes(cellData_CC5215_clean, min_cells=1)
sc.pp.filter_cells(cellData_CC5215_clean, min_counts=1)

### Save h5ad objects
cellData_DQ9512_clean.write("cells_DQ9512_regionSubset.h5ad")
sc.pl.spatial(cellData_DQ9512_clean, palette=custom_palette, color=['area'], save="_areaSubset_CellDQ9512.pdf")

cellData_CC5215_clean.write("cells_CC5215_regionSubset.h5ad")
sc.pl.spatial(cellData_CC5215_clean, palette=custom_palette, color=['area'], save="_areaSubset_CellCC5215.pdf")


###########################################
#### Plot interesting genes from paper ####
###########################################

# Genes to plot
resting_markers = ["CHRDL2","SFRP5","UCMA","ZNF550"]
proliferating_markers = ["ANGPTL2","CHPF","CNMD","COL16A1","CXCL14","FMOD","HIST1H1B","NKX3-2","PAMR1","SOX5","SYT8","THBS3"]
hypertrophic_markers = ["ALPL","ARSI","CDKN1C","CLIC3","COL10A1","DYSF","EZR","F13A1","FGFBP2","FN1","IBSP","LOX","MAN1C1","LOXL4","MARCKSL1","MGP","PCSK6","PDLIM4","PMP22","PNPLA7","PRKG1","SCUBE1","SGMS2","SLC13A5","SLC44A2","TSPAN13","VEGFA","WNK4"]

hypertrophic_markers1 = ["ALPL","ARSI","CDKN1C","CLIC3","COL10A1","DYSF","EZR","F13A1","FGFBP2","FN1","IBSP","LOX","MAN1C1","LOXL4"]
hypertrophic_markers2 = ["MARCKSL1","MGP","PCSK6","PDLIM4","PMP22","PNPLA7","PRKG1","SCUBE1","SGMS2","SLC13A5","SLC44A2","TSPAN13","VEGFA","WNK4"]

markers = ["CLU", "COL2A1", "COL9A1", "ACAN", "MATN3", "COL1A1", "CXCL12", "C1QTNF3", "MKI67"]

# Normalized each obejct
cellData_DQ9512_clean.layers["counts"] = cellData_DQ9512_clean.X.copy()
sc.pp.normalize_total(cellData_DQ9512_clean)
sc.pp.log1p(cellData_DQ9512_clean)

sc.pp.highly_variable_genes(cellData_DQ9512_clean, flavor="seurat", n_top_genes=2000)
sc.pp.scale(cellData_DQ9512_clean, max_value=10)

# Normalized each obejct
cellData_CC5215_clean.layers["counts"] = cellData_CC5215_clean.X.copy()
sc.pp.normalize_total(cellData_CC5215_clean)
sc.pp.log1p(cellData_CC5215_clean)

sc.pp.highly_variable_genes(cellData_CC5215_clean, flavor="seurat", n_top_genes=2000)
sc.pp.scale(cellData_CC5215_clean, max_value=10)

# Spatial plots
sc.pl.spatial(cellData_DQ9512_clean, color=resting_markers, cmap="Blues", save="_RZnewMarkers_DQ9512_newColours.pdf")
sc.pl.spatial(cellData_DQ9512_clean, color=proliferating_markers, cmap="Greens", save="_PZnewMarkers_DQ9512_newColours.pdf")
sc.pl.spatial(cellData_DQ9512_clean, color=hypertrophic_markers1, cmap="Blues", save="_HZnewMarkers1_DQ9512_newColours.pdf")
sc.pl.spatial(cellData_DQ9512_clean, color=hypertrophic_markers2, cmap="Blues", save="_HZnewMarkers2_DQ9512_newColours.pdf")
sc.pl.spatial(cellData_DQ9512_clean, color=markers, cmap="Blues", save="_knownMarkers_DQ9512_newColours.pdf")

sc.pl.spatial(cellData_DQ9512_clean, color=proliferating_markers, cmap="white_to_blue", save="_PZnewMarkers_DQ9512_newColours.pdf")

# Genes to plot
resting_markers = ["CHRDL2","ZNF550"]
proliferating_markers = ["ANGPTL2","CHPF","CNMD","COL16A1","CXCL14","FMOD","HIST1H1B","NKX3-2","PAMR1","SOX5","THBS3"]
hypertrophic_markers1 = ["ALPL","ARSI","CDKN1C","CLIC3","COL10A1","DYSF","EZR","F13A1","FGFBP2","FN1","IBSP","LOX","MAN1C1","LOXL4"]
hypertrophic_markers2 = ["MARCKSL1","MGP","PCSK6","PMP22","SCUBE1","SGMS2","SLC13A5","SLC44A2","TSPAN13","VEGFA","WNK4"]

markers = ["CLU", "COL2A1", "COL9A1", "ACAN", "MATN3", "COL1A1", "CXCL12", "C1QTNF3"]

# Spatial plots
sc.pl.spatial(cellData_CC5215_clean, color=resting_markers, cmap="YlGn", save="_RZnewMarkers_CC5215_newColours.pdf")
sc.pl.spatial(cellData_CC5215_clean, color=proliferating_markers, cmap="YlGn", save="_PZnewMarkers_CC5215_newColours.pdf")
sc.pl.spatial(cellData_CC5215_clean, color=hypertrophic_markers1, cmap="YlGn", save="_HZnewMarkers1_CC5215_newColours.pdf")
sc.pl.spatial(cellData_CC5215_clean, color=hypertrophic_markers2, cmap="YlGn", save="_HZnewMarkers2_CC5215_newColours.pdf")
sc.pl.spatial(cellData_CC5215_clean, color=markers, cmap="YlGn", save="_knownMarkers_CC5215_newColours.pdf")

### Other gens to plot
cellData_DQ9512_clean = sc.read_h5ad("cells_DQ9512_regionSubset.h5ad")
sc.pl.spatial(cellData_DQ9512_clean, color=["NT5E", "CLU", "APOE", "FOXA2", "IHH"], save="_MarkersForPaper_DQ9512.pdf")

sc.pl.spatial(cellData_DQ9512_clean, color=["LRP5", "SOD2", "IDS"], save="_LRP5_SOD2_IDS_DQ9512.pdf")


genes = ["AEBP1","BOC","BRWD1","CD81","CERCAM","CHRDL2","COL15A1","DOCK6","HYI","IDS","LRP5","MT1X","NUMA1","OBSCN","PDGFA",
         "PIK3CD","PMEPA1","RGL2","RILP","SFRP5","SOD2","SOD3","TAMM41","TPST2","UCMA","VSTM4","ZNF394","ZNF550","ZNF606","ZSWIM9"]
sc.pl.spatial(cellData_DQ9512_clean, color=genes, save="_upRZvsHZ_PZ_DQ9512.pdf")


# This is the repository of the article:

## _Amazon biodiversity at risk from metal contamination due to mining activity_
## resubmitted on the 19th of November 2025

### Appendix material is available in here
* Appendix 1: List of the biological species considered, with the associated proportion of their geographic distribution areas exposed to metal contamination.
* Appendix 2: List of the biological species affected both by metals (this study) and fire (Feng et al. 2021).
* Appendix 3: Percentage of overlap of the species richness in 10th, 25th and 50th percentile (i.e. hotspot classes) with mining contaminated areas in each biogeographic region. 
* Appendix 4: Maximum range exposure (%) within Indigenous territories.

### Example code to calculate the mining exposure in fishes:
* getClassificationFishes3.R

#### envidata.zip 
* This file contains the environmental rasters and shapefiles necessary to run the code _getClassificationFishes3.R_. It incudes:
  * AllDownstream.tif - a rasterfile with the river network of the Amazon
  *  amazon_basin_minewatchUTM20S.shp - a shapefile that defines the limits of the Amazon basin
  *  wetlands2_am_wgsUTM20S.tif - a raster file containing the distribution of Amazonian wetlands, as obtained from the Hydroshed project (see the paper for description)
  *  As.tif - raster file with the modelled As contamination
  *  Cu.tif - raster file with the modelled Cu contamination
  *  Hg.tif - raster file with the modelled Hg contamination
  *  Pb.tif - raster file with the modelled Pb contamination
  *  Zn.tif - raster file with the modelled Zn contamination

#### fish_sample.zip
* This file contains the range area of ~40 random fish species that can be used to test the code _getClassificationFishes3.R_

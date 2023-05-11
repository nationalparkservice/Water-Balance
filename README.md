# nationalparkservice/Water-Balance
The GitHub 'nationalparkservice/Water-Balance' repository contains tools for deriving Point and Raster/Gridded Modified Thornthwaite water balance model products.  The Excel Point Water Balance Tools derived daily water balance products, while the Python GIS Desktop tools derive a monthly water balance products.

# Excel Point Daily Water Balance Tools
Excel Point Water Balance Tool for deriving daily water balance products. This includes the NPS water balance model version 3.0 processing logic. Tools and files dervied by David Thoma Ecohydrologist, Intentory and Monitoring Division, National Park Service.

**Excel_xlsm_water_balance_model_v3.xlsm and shareable_Excel_xlsm_water_balance_model_v3** - Excel based Water Balance Tool file(s). Files with a .xlsm extension are blocked by IT security to address this matter the shareable_Excel_xlsm_water_balance_model_v3 file without the .xlsm extension is included and is downloadable to NPS/DOI computers. 

**Simple Water Balance Model Instructions ver 3.docx** - Excel Water Balance Tool Instructions document.  

**Water balance inputs outputs_2_2.docx** -  User guide for the point based 'Excel_xlsm_water_balance_model_v3.xlsm' excel water balance tool.

**Soil_Data_Development_User_Guide_v5.pdf** - USDA Soil Data Management Toolbox - user guide. This tool can be used to extract Available Water Supply GIS data which is used as input to the Water Balance Model. 

**extract_jennings_coeff.R** - R Script to extract the Jennings et. al. 2016 Snow Melt Coefficients for defiend sites in the raster file 'merged_jennings.tif'.

**merged_jennings.tif** - Jennings et. al. 2016 Snow Melt Coefficent raster file.

*sites.xlsx* - Example sites file to be used in the Excel_xlsm_water_balance_model_v3.xlsm excel tool.

*Water balance graphic.pptx* - Conceptual Water Balance Graphic.

# Python GIS Desktop Gridded Monthly Water-Balance Script
Python code for deriving gridded monthly water balance products following logic defined in Lutz et. al 2010 and Dilts et. al 2015 Journal of Biogreography articles. Script has been modified to process successive years of data (e.g. 2000, 2001, 2002, ...).  There is the option to derive water balance using  either using a 'Penman-Monteith' (Physically Based) or 'Hamon' (Emperically Based) evapotranspiration calculations. Script derived by: Kirk Sherrill (Data Manager, Rocky Mountain Network, Inventory and Monitoring Division, NPS.

**SWB_Hamon_and_Penman_ET_Python_3.x.py** - Python script to dervie monthly gridded water balance mdoel products.  

Script dependences: Python 3.x, Numpy and Geospatial Data Abstraction Library (GDAL).

Lutz, J. A., van Wagtendonk, J.W and Franklin, J.F. (2010) Climatic water deficit, tree species ranges and climate change in Yosemite Nationa Park. Journal of Biogeopgrahpy, 37 (5), 936- 950. https://doi.org/10.1111/j.1365-2699.2009.02268.x

Dilts, T.E., P.J. Weisberg, C.M. Dencker., and J.C. Chambers. 2015. Functionally relevant climate variables for arid lands: a climatic water deficit approach for modelling desert shrub distributions. Journal of Biogeography. https://doi.org/10.1111/jbi.12561


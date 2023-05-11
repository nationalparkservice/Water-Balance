# nationalparkservice/Water-Balance
The GitHub 'nationalparkservice/Water-Balance' repository contains tools for deriving Point and Raster/Gridded water balance products.  Te Excel Point Water Balance Tools derived daily water balance products, while the Python GIS Desktop tools derive a montyhly water balance products.


# Excel Point Water Balance Tools
Files for the Excel Point Water Balance Tool for deriving daily water balance products. This includes the NPS water balance model version 3.0 processing logic.

Tools and files dervied by David Thoma Ecohydrologist, Intentory and Monitoring Division, National Park Service.

## Excel_xlsm_water_balance_model_v3.xlsm 
Is the excel based tool file. Files with a .xlsm extension are blocked by IT security to address this matter the shareable_Excel_xlsm_water_balance_model_v3 file without the .xlsm extension is included and is downloadable to NPS/DOI computers. 

**Water balance inputs outputs_2_2.docx**: User guide for the point based 'Excel_xlsm_water_balance_model_v3.xlsm' excel water balance tool

**Soil_Data_Development_User_Guide_v5.pdf**: USDA Soil Data Management Toolbox - user guide. This tool can be used to extract Available Water Supply GIS data which is used as input to the Water Balance Model. 

# Python GIS Desktop Monthly Water-Balance Scripts

Python code for deriving gridded monthly water balance products following logic defined in Lutz et. al 2010 and Dilts et. al 2015 Journal of Biogreography articles. Script has been modified to process successive years of data (e.g. 2000, 2001, 2002, ...).  There is the option to derive water balance using  either using a 'Penman-Monteith' (Physically Based) or 'Hamon' (Emperically Based) evapotranspiration calculations.

Script dependences: Python 3.x, Numpy and Geospatial Data Abstraction Library (GDAL).

Scripts derived By: Kirk Sherrill (Data Manager, Rocky Mountain Network, Inventory and Monitoring Division, NPS

Lutz, J. A., van Wagtendonk, J.W and Franklin, J.F. (2010) Climatic water deficit, tree species ranges and climate change in Yosemite Nationa Park. Journal of Biogeopgrahpy, 37 (5), 936- 950. https://doi.org/10.1111/j.1365-2699.2009.02268.x

Dilts, T.E., P.J. Weisberg, C.M. Dencker., and J.C. Chambers. 2015. Functionally relevant climate variables for arid lands: a climatic water deficit approach for modelling desert shrub distributions. Journal of Biogeography. https://doi.org/10.1111/jbi.12561

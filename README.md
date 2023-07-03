## README file for the isotopic box model from Skeie et al. (2023)
******************************************************************

Author:  Oivind Hodnebrog (oivinho@cicero.oslo.no)
Version: 2023-07-03


Reference
---------
The isotopic box model has been used and described in the paper by:
Skeie, R.B., Hodnebrog, O., Myhre, G.: "Methane trend over the last decades driven and modified by anthropogenic emissions", to be published in the journal "Communications Earth & Environment", 2023.


Contents
--------
- main_BoxModel.m: Main Matlab script for running the isotopic box model
- make_Figure6.m: Matlab script used to generate Figure 6 in the paper (see below)
- functions/: Folder with Matlab functions
- inputdata/: Folder with input data needed to run the box model
- matfiles/: Folder with box model output stored as binary "mat" files (Matlab format)
             (contains results from the 87 different box model runs used in Skeie et al.)
- plots/: Folder to store plots from box model simulations


Input data
----------
- Emissions of CH4 from CEDS:
  - Download from DOI: 10.5281/zenodo.4741285
  - Filename: CH4_global_CEDS_emissions_by_sector_2021_04_21.csv (put file in inputdata/)
- Emissions of CH4 from EDGAR:
  - Download from: https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/EDGAR/datasets/v70_FT2021_GHG/EDGAR_CH4_1970-2021.zip
  - After downloading, extract the zip file, open excel workbook, and save the sheet "IPCC 2006" as EDGAR_CH4_1970-2021.csv
  - Remove the first 10 (header) lines and put the csv file in inputdata/
- Emissions of biomass burning from IPCC AR6+GFED:
  - Run the script https://github.com/ragnhibs/Methane-trend/blob/main/figures/plot_figure_2.py
  - Put file in inputdata/ and set correct path in variable file_GFED in main_BoxModel.m
- Natural emissions of CH4:
  - Exists in inputdata/natemis.csv
- Table of global annual mean mole fractions:
  - Download from: https://gml.noaa.gov/aggi/NOAA_Annual_Mean_MoleFractions_2023.csv
  - After downloading, set correct path in variable file_CH4_NOAA in main_BoxModel.m
- Annual CH4 mixing ratios from Meinshausen et al. (2020, 10.5194/gmd-13-3571-2020):
  - Exists in inputdata/WMGHG_vmr_SSP2-4.5
- Time evolution of OH from OsloCTM3:
  - Exists in inputdata/OsloCTM3_OH_histO3_ceds2021.txt
- Time evolution of OH from AerChemMIP and CCMI:
  - See separate CH4 box model here: https://github.com/ragnhibs/Methane-trend/


How to reproduce figures in Skeie et al. (2023)
-----------------------------------------------
Figure 6a:
- Run make_Figure6.m with the preset input options except set INP_idealized_OH=1

Figure 6b:
- Run make_Figure6.m with preset options except INP_anthro_const2007=1 and INP_inc_dC13emis_unc=1
- Run make_Figure6.m with preset options except INP_anthro_const2007=2 and INP_inc_dC13emis_unc=2
- Run make_Figure6.m with preset options except INP_anthro_const2007=3 and INP_inc_dC13emis_unc=3
- Run make_Figure6.m with preset options except INP_anthro_const2007=4 and INP_inc_dC13emis_unc=4

Figure 6c:
- Run make_Figure6.m with preset options except INP_nat_emis=1
- Run make_Figure6.m with preset options except INP_nat_emis=2
- Run make_Figure6.m with preset options except INP_nat_emis=3 and INP_inc_dC13emis_unc=6
- Run make_Figure6.m with preset options except INP_nat_emis=10
- Run make_Figure6.m with preset options except INP_nat_emis=11 and INP_inc_dC13emis_unc=6

Figure S10:
- Run main_BoxModel.m with the preset input options
- Set idealized_OH=1 and rerun the script to get the simulation with constant OH

Figure S11:
- Run main_BoxModel.m with the preset options to get plot with EDGARv7 emissions
- Set anthro_emis=1, eyear=2019 and eyear_plot=2019 to get plot with CEDS-2021 emissions

Figure S12:
- Same procedure as for Figure 6b (see above) except set INP_anthro_emis=1

Figure S13a:
- Run main_BoxModel.m with the preset input options except change the following variables:
  emisBB_clim=5, plot_aerchemmip=0, plot_ccmi=0, plot_dC13CH4_twoyaxes=0, plot_emis=0, init_OH=1302918

Figure S13b-d:
- Same procedure as for Figure 6a-c (see above) except set INP_emisBB_clim=5 and INP_init_OH=1302918

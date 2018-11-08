#!/usr/bin/env python
''' 
Created on Apr 2, 2018

@author: Shugong Wang
@email: shugong.wang@nasa.gov

2018 Apr 06, a index bug has been fixed. Shugong Wang 
2018 Apr 10, a doule rotation bug has been fixed by Shugong Wang and Eric Kemp 
'''

import os 
import sys
import mule 
import ConfigParser
import numpy as np 
import numpy.ma as ma
import netCDF4 as nc4 
import shutil as shell

class dump2lis(object):
  '''
  dump2lis: taking JULES state variables from UM dump file
  and insert into LIS restart file 
  JULES STATE VARIABLES


  '''
  def __init__(self, config_file_name):
    self.config_file = config_file_name
    if os.path.isfile(self.config_file)==False:
      print("The configuration file {0} not exist".format(self.config_file))
      sys.exit(0) 

    self.dump_file     = ""
    self.lis_rst_file0 = ""
    self.lis_rst_file1 = ""
    self.stashmaster   = ""
    self.nrow = 1920
    self.ncol = 2560
    self.dump_var = np.zeros((self.nrow, self.ncol)) 

  def read_config_file(self):
    config = ConfigParser.ConfigParser()
    config.read(self.config_file)
    self.dump_file     = config.get("UM", "UM dump file")
    self.stashmaster   = config.get("UM", "STASH Master file")
    self.lis_rst_file0 = config.get("LIS", "LIS restart file (0)")
    self.lis_rst_file1 = config.get("LIS", "LIS restart file (1)")
    self.lis_ldt_file  = config.get("LIS", "LIS domain and parameter file")

  def prepare(self):
    print("reading STASH master file {0}".format(self.stashmaster))
    sm=mule.stashmaster.STASHmaster.from_file(self.stashmaster)
    print("reading UM dump file {0}".format(self.dump_file))
    self.dump=mule.DumpFile.from_file(self.dump_file)
    print("attaching the STASH master file {0} to the dump file {1}".format(self.stashmaster, self.dump_file))
    self.dump.attach_stashmaster_info(sm)
    print("creating new LIS restart file {0}".format(self.lis_rst_file1)) 
    shell.copy(self.lis_rst_file0, self.lis_rst_file1) 
    # open the new NetCDF restart file for updating 
    self.nc_rst    = nc4.Dataset(self.lis_rst_file1, "r+")
    self.lis_ldt   = nc4.Dataset(self.lis_ldt_file, "r")  
    self.land_mask = self.lis_ldt["LANDMASK"][:,:]
    self.lis_var   = np.zeros((int(self.land_mask.sum()), 1))

  def read_dump_variable(self, section, item, lblev=-1, lbuser5=-1):
    if lblev>=0 and lbuser5==-1:
      for field in self.dump.fields:
        if field.stash.section==section and field.stash.item==item and field.lblev==lblev:
          info = "section=%03d item=%03d lblev=%02d %s"  %(field.stash.section, field.stash.item, field.lblev, field.stash.name)
          print(info)
          self.dump_var[:,:] = field.get_data()[:,:]
          break
    elif lblev==-1 and lbuser5>=0:
      for field in self.dump.fields:
        if field.stash.section==section and field.stash.item==item and field.lbuser5==lbuser5:
          info = "section=%03d item=%03d lbuser5=%02d %s"  %(field.stash.section, field.stash.item, field.lbuser5, field.stash.name)
          print(info)
          self.dump_var[:,:] = field.get_data()[:,:]
          break
    elif lblev==-1 and lbuser5==-1:
      for field in self.dump.fields:
        if field.stash.section==section and field.stash.item==item:
          info = "section=%03d item=%03d %s"  %(field.stash.section, field.stash.item,  field.stash.name)
          print(info)
          self.dump_var[:,:] = field.get_data()[:,:]
          break
  
  def shift_um2lis(self):
    data = self.dump_var
    col1 = self.ncol/2
    col2 = self.ncol
    tmp = ma.zeros((self.nrow, self.ncol))
    tmp[:,:] = data[:,:]

    data[:,0   :col1] = tmp[:,col1:col2]
    data[:,col1:col2] = tmp[:,0   :col1]

  def convert_map2tile(self):
    data2d = self.dump_var
    data1d = self.lis_var 
    k = 0
    for row in range(0, self.nrow):
      for col in range(0, self.ncol):
        if self.land_mask[row,col] == 1:
          data1d[k] = data2d[row,col]
          k = k + 1 
  
  def update_lis_variable(self, var_name, lev=-1):
    self.shift_um2lis()    
    self.convert_map2tile()
    if lev==-1:
      self.nc_rst[var_name][:]     = self.lis_var[:]
    else:
      self.nc_rst[var_name][lev,:] = self.lis_var.transpose()


  def finalize(self):
    self.nc_rst.close() 

if __name__=="__main__":
  # using the command line argument 
  # dump2lis config_file_name 
  if len(sys.argv) <> 2:
    print("Usage: dump2lis config_file_name")
    sys.exit() 
  if os.path.isfile(sys.argv[1]) <> True:
    print("The configuation file {0} does not exist. Please double check.".format(sys.argv[1]))
  d2l = dump2lis(sys.argv[1]) 
  d2l.read_config_file() 
  d2l.prepare()
  # 1. NSNOW: all 0, keep LIS values
  # 2. TSURF_ELEV_SURFT: all 0, keep LIS values
  # 3. Z0MSEA: all 3.9999999e-05, keep LIS values 
  # 4. FRAC_PAST_PREV: all 0, keep LIS values
  # 5. N_INORG_SOILT_LYRS, all 0, keep LIS values 
  # 6. N_INORG_AVAIL_PFT, all 0, keep LIS values 
  # 7. TRIFFID_CO2_GB, all 0, keep LIS values 
  # 8. CANOPY_GB, read section=000 item=022 lblev=00 CANOPY WATER AFTER TIMESTEP    KG/M2
  d2l.read_dump_variable(0, 22, 0)
  d2l.update_lis_variable("CANOPY_GB") 

  #18 CANOPY, the same as CANOPY_GB in aggregated surface mode
  d2l.read_dump_variable(0, 22, 0)
  d2l.update_lis_variable("CANOPY") 
  
  # 9 SNOW_MASS_IJ, section=000 item=023 lblev=9999 lbuser4=0023 SNOW AMOUNT OVER LAND AFT TSTP KG/M2
  d2l.read_dump_variable(0, 23)
  d2l.update_lis_variable("SNOW_MASS_IJ") 
  
  #27 SNOW_TILE, the same value as SNOW_MASS_IJ in aggregated surface mode
  d2l.read_dump_variable(0, 23)
  d2l.update_lis_variable("SNOW_TILE") 
  
  #10 TSOIL_DEEP, all 0 
  #11 WOOD_PROD_FAST, all 0
  #12 WOOD_PROD_MED, all 0
  #13 WOOD_PROD_SLOW, all 0 
  #14 FRAC_AGR_PREV, all 0
  #15 FRAC_AGR, all 0 
  #16 N_INORG
  #17 CANHT_FT
  # section=000 item=218 lblev=9999 lbuser5=0001 CANOPY HEIGHT OF PLANT FUNC TYPES M
  for k in range(0,5):
    d2l.read_dump_variable(0, 218, lbuser5=k+1)
    d2l.update_lis_variable("CANHT_FT", k)

  #19 NS, organic soil nitrogen, all 0.001
  #20 CS, grid box soil carbon in each pool, all 12.1 
  #21 GC, Tile surface conductance to evaporation for land tiles
  # section=000 item=213 lblev=00 lbuser4=0213 CANOPY CONDUCTANCE AFTER TIMESTEP
  d2l.read_dump_variable(0, 213, 0)
  d2l.update_lis_variable("GC") 
  
  #22 GS, Gridbox surface conductance to evaporation
  # It is weird that the GS is in the dump file module instead GC. 
  # We should double check with Met Office. 
  # Assign the same value as GC to GS as a kluge solution 
  d2l.read_dump_variable(0, 213, 0)
  d2l.update_lis_variable("GS")

  #23 LAI, jules_vegetation.nml:l_phenol=.false., it is not used

  #24 RGRAIN, all 0 
  #   jules_radiation.nml:l_snow_albedo=.false.
  #   jules_radiation.nml:l_embedded_snow=.false.

  #25 SMC, should not be a state variable 
  #26 SMCL 
  # section=000 item=009 lblev=01 lbuser4=0009 SOIL MOISTURE CONTENT IN A LAYER
  for k in range(0, 4):
    d2l.read_dump_variable(0, 9, k+1)
    d2l.update_lis_variable("SMCL", k) 
  
  #27
  #28 SNOW_GRND, all 0 
  
  #29 SOOT, all 0 
  
  #30 T_SOIL
  # section=000 item=020 lblev=01 lbuser4=0020 DEEP SOIL TEMP AFTER TIMESTEP
  for k in range(0, 4):
    d2l.read_dump_variable(0, 20, k+1)
    d2l.update_lis_variable("T_SOIL", k) 
    
  #31 TSTAR_TILE
  #section=000 item=024 lblev=9999 lbuser4=0024 SURFACE TEMPERATURE AFTER TIMESTEP 
  d2l.read_dump_variable(0, 24) 
  d2l.update_lis_variable("TSTAR_TILE") 

  #32 ASTEPS_SINCE_TRIFFID, none 0 number

  #33 LAI_PHEN, all 0
  
  #34 C_VEG, all 0

  #35 CV, all 0
  
  #36 FEXP, parameter 
  
  #37 TI_MEAN, parameter

  #38 TI_SIG, parameter

  #39 ZW
  # section=000 item=278 lblev=9999 lbuser4=0278 MEAN WATER TABLE DEPTH            M
  d2l.read_dump_variable(0, 278)
  d2l.update_lis_variable("ZW") 

  #40 STHZW 
  # section=000 item=281 lblev=9999 lbuser4=0281 SATURATION FRAC IN DEEP LAYER
  d2l.read_dump_variable(0, 281)
  d2l.update_lis_variable("STHZW") 

  #41 FWETL
  # section=000 item=280 lblev=9999 lbuser4=0280 SURFACE WETLAND FRACTION
  d2l.read_dump_variable(0, 280)
  d2l.update_lis_variable("FWETL") 
  
  #42 FSAT
  # section=000 item=279 lblev=9999 lbuser4=0279 SURFACE SATURATION FRACTION
  d2l.read_dump_variable(0, 279)
  d2l.update_lis_variable("FSAT") 


  #43 STHU
  # section=000 item=214 lblev=01 lbuser4=0214 UNFROZEN SOIL MOISTURE FRAC AFTER TS
  for k in range(0, 4):
    d2l.read_dump_variable(0, 214, k+1)
    d2l.update_lis_variable("STHU", k) 

  #44 STHF 
  # section=000 item=215 lblev=01 lbuser4=0215 FROZEN SOIL MOISTURE FRAC AFTER TS
  for k in range(0, 4):
    d2l.read_dump_variable(0, 215, k+1)
    d2l.update_lis_variable("STHF", k) 

  #45 SICE, all 0

  #46 SLIQ, all 0 

  #47 SNOWDEPTH
  # section=000 item=376 lblev=9999  SNOW DEPTH ON GROUND ON TILES (M)
  d2l.read_dump_variable(0, 376)
  d2l.update_lis_variable("SNOWDEPTH") 
  
  #48: TSNOW, all 0

  #49: RGRAINL, all 0
  #50: RHO_SNOW_GRND, all 250
  #51: RHO_SNOW, all 0 
  #52: DS, all 0 
  
  d2l.finalize() 

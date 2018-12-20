#!/usr/bin/env python 
import netCDF4 as nc4
import json
import h5py
import numpy as np 
import os, errno
#global constants, defs etc:
 #Prolly wrong check later
def getLatIdx(self, latitude):
  idx = (latitude-self.llat)/self.resolution
  return int(idx)

def getLonIdx(self, longitude):
  idx = (longitude-self.llon)/self.resolution
  return int(idx)

class awral():
  def __init__(self):
    self.tile_fractions=[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    self.llat = -44.00
    self.llon = 112.00
    self.latitude=-33.00
    self.longitude=145.00
    self.num_cellsx=841
    self.num_cellsy=681
    self.resolution=0.05
    self.domsize=2
    self.halosize=6
    self.num_hrus=2
    self.num_hypsobins=20
    # read this in from json file:
    with open('DefaultParameters.json') as f:
      data = json.load(f)
    for member in data:
      if(member['MemberName'] == 'K_rout_scale'):
        self.k_rout_scale=member['Value'] #0.04855041
      if(member['MemberName'] == 'K_rout_int'):
        self.k_rout_int=member['Value'] #0.16542047
      if(member['MemberName'] == 'K_gw_scale'):
        self.k_gw_scale=member['Value'] #0.93908135
      if(member['MemberName'] == 'S0max_scale'):
        self.s0max_scale=member['Value'] #2.80487914
      if(member['MemberName'] == 'Ssmax_scale'):
        self.ssmax_scale=member['Value'] #1.99322065
      if(member['MemberName'] == 'Sdmax_scale'):
        self.sdmax_scale=member['Value'] #0.88436116
      if(member['MemberName'] == 'K0sat_scale'):
        self.k0sat_scale=member['Value'] #3.89200702
      if(member['MemberName'] == 'Kssat_scale'):
        self.kssat_scale=member['Value'] #0.05234894
      if(member['MemberName'] == 'Kdsat_scale'):
        self.kdsat_scale=member['Value'] #0.0122021
      if(member['MemberName'] == 'Pref_gridscale'):
        self.pref_gridscale=member['Value'] #2.56420087
      if(member['MemberName'] == 'ne_scale'):
        self.ne_scale=member['Value'] #0.04253378

  def mapping(self, SPATIAL_FILE):
        ds = h5py.File(SPATIAL_FILE,mode='r')
        SPATIAL_GRIDS = list(ds['parameters'])

        self.HEIGHT = ds['parameters']['height'][:]
        self.SLOPE = ds['parameters']['slope'][:]
        self.meanpet_grid = ds['parameters']['meanPET'][:]
        self.K_ROUT = self.k_rout_scale* self.meanpet_grid + self.k_rout_int
        self.k_gw_grid = ds['parameters']['k_gw'][:]
        self.K_GW = np.multiply(self.k_gw_scale,self.k_gw_grid)
        self.s0fracawc_grid = ds['parameters']['s0fracAWC'][:]
        self.ssfracawc_grid = ds['parameters']['ssfracAWC'][:]
        self.sdfracawc_grid = ds['parameters']['sdfracAWC'][:]
        self.S0MAX = 100.*np.multiply(self.s0max_scale,self.s0fracawc_grid)
        self.SSMAX = 900.*np.multiply(self.ssmax_scale,self.ssfracawc_grid)
        self.SDMAX = 5000.*np.multiply(self.sdmax_scale,self.sdfracawc_grid)
        self.k0sat_grid = ds['parameters']['k0sat'][:]
        self.kssat_grid = ds['parameters']['kssat'][:]
        self.kdsat_grid = ds['parameters']['kdsat'][:]
        self.K0SAT = np.multiply(self.k0sat_scale,self.k0sat_grid)
        self.KSSAT = np.multiply(self.kssat_scale,self.kssat_grid)
        self.KDSAT = np.multiply(self.kdsat_scale,self.kdsat_grid)
        self.KR_0S = np.minimum(150.0, np.maximum(1.0, self.K0SAT/self.KSSAT))
        self.KR_SD = np.minimum(150.0, np.maximum(1.0, self.KSSAT/self.KDSAT))
        self.pref_grid = ds['parameters']['pref'][:]
        self.PREFR = np.multiply(self.pref_gridscale,self.pref_grid)
        self.ne_grid = ds['parameters']['ne']
        self.NE = np.multiply(self.ne_scale,self.ne_grid)

        #HRU stuff:
        self.f_tree_grid = ds['parameters']['f_tree'][:]
        self.sub_f_tree = np.subtract(1.0,self.f_tree_grid)
        self.FHRU = (self.f_tree_grid, self.sub_f_tree)
        self.FHRU = np.reshape(self.FHRU,(2,self.num_cellsy,self.num_cellsx))
        self.hveg_grid = ds['parameters']['hveg_dr'][:]
        zeros = np.zeros((self.num_cellsy, self.num_cellsx))
        self.HVEG = (self.hveg_grid, zeros)
        self.HVEG = np.reshape(self.HVEG,(2,self.num_cellsy,self.num_cellsx))
        self.laimax_grid = ds['parameters']['lai_max'][:]
        self.LAIMAX = (self.laimax_grid, self.laimax_grid)
        self.LAIMAX = np.reshape(self.LAIMAX,(2,self.num_cellsy,self.num_cellsx))

        self.hru_params={"AWRAL600_FHRU":self.FHRU, "AWRAL600_HVEG":self.HVEG, "AWRAL600_LAIMAX":self.LAIMAX}
        self.grid_params={"AWRAL600_K_ROUT":self.K_ROUT, "AWRAL600_KSSAT":self.KSSAT, "AWRAL600_PREFR":self.PREFR, "AWRAL600_S0MAX":self.S0MAX, "AWRAL600_SLOPE":self.SLOPE, "AWRAL600_SSMAX":self.SSMAX, "AWRAL600_KDSAT":self.KDSAT, "AWRAL600_NE":self.NE} 
        self.hypso_params={"AWRAL600_HEIGHT":self.HEIGHT}

  def write_nc(self, nc_file):
    try:
      os.remove(nc_file)
    except OSError:
      pass
    tmp = np.zeros((2,2),dtype=float)
    rootgrp = nc4.Dataset(nc_file, 'w')
    lat = rootgrp.createDimension("north_south", self.domsize)
    lon = rootgrp.createDimension("east_west", self.domsize)
    lat_b = rootgrp.createDimension("north_south_b", self.halosize)
    lon_b = rootgrp.createDimension("east_west_b", self.halosize)
    sfctypes = rootgrp.createDimension("sfctypes", 9)
    hrutypes = rootgrp.createDimension("hrutypes", self.num_hrus)
    hypsobins = rootgrp.createDimension("hypsobins", self.num_hypsobins)
    self.x_idx = getLonIdx(self,self.longitude)
    self.y_idx = getLatIdx(self,self.latitude)
    # domain mask
    rootgrp.createVariable("DOMAINMASK", "f4", ("north_south", "east_west"))
    rootgrp["DOMAINMASK"][:,:] = 1
    #rootgrp["DOMAINMASK"][0,0] = 1 
    rootgrp.createVariable("LANDMASK", "f4", ("north_south", "east_west"))
    rootgrp["LANDMASK"][:,:] = 1
    #rootgrp["LANDMASK"][0,0] = 1 
    #HRU params 
    for k,v in self.hru_params.items():
      var_name = "%s" %(k)
      slice = v[:,self.y_idx:self.y_idx+self.domsize,self.x_idx:self.x_idx+self.domsize]
      rootgrp.createVariable(var_name, "f4", ("hrutypes", "north_south", "east_west"))
      rootgrp[var_name][:,:,:] = slice    
    #Grid params
    for k,v in self.grid_params.items():
      var_name = "%s" % (k)
      slice = v[self.y_idx:self.y_idx+self.domsize,self.x_idx:self.x_idx+self.domsize]
      rootgrp.createVariable(var_name, "f4", ("north_south", "east_west"))
      rootgrp[var_name][:,:] = slice   
    #Hypso params
    for k,v in self.hypso_params.items():
      var_name = "%s" % (k)
      slice = v[:,self.y_idx:self.y_idx+self.domsize,self.x_idx:self.x_idx+self.domsize]
      rootgrp.createVariable(var_name, "f4", ("hypsobins", "north_south", "east_west"))
      rootgrp[var_name][:,:,:] = slice

    rootgrp.createVariable("SURFACETYPE", "f4", ("sfctypes", "north_south", "east_west"))
    rootgrp.createVariable("LANDCOVER", "f4", ("sfctypes", "north_south", "east_west"))
    for n in range(0, 1):
      if self.tile_fractions[n]>0:
        rootgrp["SURFACETYPE"][n, :,:] = 1.0
        rootgrp["LANDCOVER"][n, :,:] = self.tile_fractions[n]
      else:
        rootgrp["SURFACETYPE"][n, :,:] = 0.0
        rootgrp["LANDCOVER"][n, :,:] = 0.0 
        
    
    rootgrp.createVariable("lat", "f4", ("north_south", "east_west"))
    rootgrp.createVariable("lon", "f4", ("north_south", "east_west"))
    rootgrp.createVariable("lat_b", "f4", ("north_south_b", "east_west_b"))
    rootgrp.createVariable("lon_b", "f4", ("north_south_b", "east_west_b"))    
    rootgrp["lat"][0,:] = self.latitude
    rootgrp["lon"][:,0] = self.longitude
    rootgrp["lat_b"][0:int(self.halosize/2),:] = self.latitude
    rootgrp["lon_b"][:,0:int(self.halosize/2)] = self.longitude
    for i in range(1,self.domsize):
      rootgrp["lat"][i,:] = self.latitude + i*self.resolution
      rootgrp["lon"][:,i] = self.longitude + i*self.resolution
      rootgrp["lat_b"][i*int(self.halosize/2):i*self.halosize,:] = self.latitude + i*self.resolution
      rootgrp["lon_b"][:,i*int(self.halosize/2):i*self.halosize] = self.longitude + i*self.resolution 
    # global attributes
    rootgrp.MAP_PROJECTION = "EQUIDISTANT CYLINDRICAL" 
    rootgrp.SOUTH_WEST_CORNER_LAT = self.latitude 
    rootgrp.SOUTH_WEST_CORNER_LON = self.longitude 
    rootgrp.DX = 0.05
    rootgrp.DY = 0.05 
    rootgrp.INC_WATER_PTS = "false" 
    rootgrp.LANDCOVER_SCHEME = "" 
    rootgrp.WATERCLASS = 0 
    rootgrp.LANDCLASS = 1 
    rootgrp.NUMVEGTYPES = 1 
    rootgrp.LANDMASK_SOURCE = "AWRAL_LANDMASK" 
    rootgrp.title = "Land Data Toolkit (LDT) parameter-processed output" 
    rootgrp.institution = "NASA GSFC Hydrological Sciences Laboratory" 
    rootgrp.history = "created on date: 2017-08-22T14:34:26.543" 
    rootgrp.references = "Kumar_etal_EMS_2006, Peters-Lidard_etal_ISSE_2007" 
    rootgrp.comment = "website: http://lis.gsfc.nasa.gov/" 
    rootgrp.close()
    
def main():
  awra = awral()
  awra.mapping("spatial_parameters.h5")
  awra.write_nc("awral_lis_input.nc")
  
if __name__=="__main__":
  main()

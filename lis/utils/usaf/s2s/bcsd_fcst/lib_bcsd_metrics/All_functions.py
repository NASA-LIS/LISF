def CALC_EMPIRICAL_PCTL(TARGET_VAL, OBS_CLIM, MEAN_TYPE):
	from scipy.stats import percentileofscore as pscore
	PCTL = pscore(OBS_CLIM, TARGET_VAL, kind=MEAN_TYPE)
	return PCTL

def Sel_var (SEL_CIM_DATA, VAR_NAME, MODEL):
	import numpy as np
	# Now Selecting the climatology further for the given variable
	if VAR_NAME == 'RootZone-SM':
		if (MODEL=='CLSM'):
			VAR_SEL_CLIM_DATA = SEL_CIM_DATA.SoilMoist_tavg.isel(SoilMoist_profiles=1) # for clsm the layer-2 is rootzone soil moisture
		elif (MODEL=='NOAHMP') or (MODEL=='NoahMP'):
			VAR_SEL_CLIM_DATA = SEL_CIM_DATA.SoilMoist_tavg.isel(SoilMoist_profiles=np.arange(3)).sum(dim='SoilMoist_profiles') ## Summing soil moisture from top 3 layers
	if VAR_NAME == 'Total-SM':
		if (MODEL=='CLSM'):
			VAR_SEL_CLIM_DATA = SEL_CIM_DATA.SoilMoist_tavg.isel(SoilMoist_profiles=2) # for clsm the total soil moisture is in the third layer
		else:
			VAR_SEL_CLIM_DATA = SEL_CIM_DATA.SoilMoist_tavg.sum(dim='SoilMoist_profiles') ## Summing soil moisture from all layers
	if VAR_NAME == 'Surface-SM':
		VAR_SEL_CLIM_DATA = SEL_CIM_DATA.SoilMoist_tavg.isel(SoilMoist_profiles=0)

	elif VAR_NAME == 'Total-Runoff':
		VAR_SEL_CLIM_DATA = SEL_CIM_DATA.Qs_tavg+SEL_CIM_DATA.Qsb_tavg ## Adding total surface runoff with sub-surface runoff
	elif VAR_NAME == 'TWS':
		VAR_SEL_CLIM_DATA = SEL_CIM_DATA.TWS_tavg
	elif VAR_NAME == 'Precip':
		VAR_SEL_CLIM_DATA = SEL_CIM_DATA.TotalPrecip_tavg
	elif (VAR_NAME=='Air-T'):
		VAR_SEL_CLIM_DATA = SEL_CIM_DATA.Tair_f_tavg
	elif (VAR_NAME=='ET'):
		VAR_SEL_CLIM_DATA = SEL_CIM_DATA.Evap_tavg
	return VAR_SEL_CLIM_DATA


def GET_BOUNDARY(REGION_NAME):
	BOUNDARY_EA = (22, 52, -12, 23); BOUNDARY_WA = (-19, 26, -5, 25)
	BOUNDARY_SA = (8, 52, -37, 0); BOUNDARY_FAME = (-20, 60, -40, 40)
	if (REGION_NAME == 'EA'):
		Boundary = BOUNDARY_EA
	elif (REGION_NAME == 'WA'):
		Boundary = BOUNDARY_WA
	elif (REGION_NAME == 'SA'):
		Boundary = BOUNDARY_SA
	else: 
		Boundary = BOUNDARY_FAME
	return Boundary

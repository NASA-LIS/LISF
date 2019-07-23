mv module_sf_noahlsm.F module_sf_noah39lsm.F90
mv module_model_constants.F module_model_constants_39.F90
mv module_sf_noahlsm_glacial_only.F module_sf_noah39lsm_glacial.F90

sed -i 's/module_sf_noahlsm/module_sf_noah39lsm/g' module_sf_noah39lsm.F90
sed -i 's/module_sf_noahlsm/module_sf_noah39lsm/g' module_sf_noah39lsm_glacial.F90
sed -i 's/module_model_constants/module_model_constants_39/g' module_model_constants_39.F90
sed -i 's/module_sf_noahlsm_glacial_only/module_sf_noah39lsm_glacial/g' module_sf_noah39lsm_glacial.F90
sed -i 's/module_model_constants/module_model_constants_39/g' module_sf_noah39lsm.F90
sed -i "s/module_model_constants/module_model_constants_39/g" module_sf_noah39lsm_glacial.F90

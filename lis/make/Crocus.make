FLAGS_Crocus := -O0 -g -autodouble -noerror_limit -FR -auto -WB -traceback
dummy_hook.o:               dummy_hook.F90
	$(FC) $(FFLAGS) $(FLAGS_Crocus) $(HEADER_DIRS) $<
ini_csts.o:                 ini_csts.F90
	$(FC) $(FFLAGS) $(FLAGS_Crocus) $(HEADER_DIRS) $<
modd_const_atm.o:           modd_const_atm.F90
	$(FC) $(FFLAGS) $(FLAGS_Crocus) $(HEADER_DIRS) $<
modd_const_tartes.o:        modd_const_tartes.F90
	$(FC) $(FFLAGS) $(FLAGS_Crocus) $(HEADER_DIRS) $<
modd_csts.o:                modd_csts.F90
	$(FC) $(FFLAGS) $(FLAGS_Crocus) $(HEADER_DIRS) $<
modd_flood_par.o:           modd_flood_par.F90
	$(FC) $(FFLAGS) $(FLAGS_Crocus) $(HEADER_DIRS) $<
modd_meb_par.o:             modd_meb_par.F90
	$(FC) $(FFLAGS) $(FLAGS_Crocus) $(HEADER_DIRS) $<
modd_prep_snow.o:           modd_prep_snow.F90
	$(FC) $(FFLAGS) $(FLAGS_Crocus) $(HEADER_DIRS) $<
modd_reprod_oper.o:         modd_reprod_oper.F90
	$(FC) $(FFLAGS) $(FLAGS_Crocus) $(HEADER_DIRS) $<
modd_snow_metamo.o:         modd_snow_metamo.F90
	$(FC) $(FFLAGS) $(FLAGS_Crocus) $(HEADER_DIRS) $<
modd_snow_par.o:            modd_snow_par.F90
	$(FC) $(FFLAGS) $(FLAGS_Crocus) $(HEADER_DIRS) $<
modd_surf_atm.o:            modd_surf_atm.F90
	$(FC) $(FFLAGS) $(FLAGS_Crocus) $(HEADER_DIRS) $<
modd_surf_conf.o:           modd_surf_conf.F90
	$(FC) $(FFLAGS) $(FLAGS_Crocus) $(HEADER_DIRS) $<
modd_surfex_mpi.o:          modd_surfex_mpi.F90
	$(FC) $(FFLAGS) $(FLAGS_Crocus) $(HEADER_DIRS) $<
modd_surfex_omp.o:          modd_surfex_omp.F90
	$(FC) $(FFLAGS) $(FLAGS_Crocus) $(HEADER_DIRS) $<
modd_surf_par.o:            modd_surf_par.F90
	$(FC) $(FFLAGS) $(FLAGS_Crocus) $(HEADER_DIRS) $<
modd_type_date_surf.o:      modd_type_date_surf.F90
	$(FC) $(FFLAGS) $(FLAGS_Crocus) $(HEADER_DIRS) $<
modd_water_par.o:           modd_water_par.F90
	$(FC) $(FFLAGS) $(FLAGS_Crocus) $(HEADER_DIRS) $<
mode_crodebug.o:            mode_crodebug.F90
	$(FC) $(FFLAGS) $(FLAGS_Crocus) $(HEADER_DIRS) $<
mode_snow3l.o:              mode_snow3l.F90
	$(FC) $(FFLAGS) $(FLAGS_Crocus) $(HEADER_DIRS) $<
mode_tartes.o:              mode_tartes.F90
	$(FC) $(FFLAGS) $(FLAGS_Crocus) $(HEADER_DIRS) $<
mode_thermos.o:             mode_thermos.F90
	$(FC) $(FFLAGS) $(FLAGS_Crocus) $(HEADER_DIRS) $<
modi_abor1_sfx.o:           modi_abor1_sfx.F90
	$(FC) $(FFLAGS) $(FLAGS_Crocus) $(HEADER_DIRS) $<
modi_close_file.o:          modi_close_file.F90
	$(FC) $(FFLAGS) $(FLAGS_Crocus) $(HEADER_DIRS) $<
modi_get_luout.o:           modi_get_luout.F90
	$(FC) $(FFLAGS) $(FLAGS_Crocus) $(HEADER_DIRS) $<
modi_ini_surf_csts.o:       modi_ini_surf_csts.F90
	$(FC) $(FFLAGS) $(FLAGS_Crocus) $(HEADER_DIRS) $<
modi_parkind1.o:            modi_parkind1.F90
	$(FC) $(FFLAGS) $(FLAGS_Crocus) $(HEADER_DIRS) $<
modi_surface_aero_cond.o:   modi_surface_aero_cond.F90
	$(FC) $(FFLAGS) $(FLAGS_Crocus) $(HEADER_DIRS) $<
modi_surface_cd.o:          modi_surface_cd.F90
	$(FC) $(FFLAGS) $(FLAGS_Crocus) $(HEADER_DIRS) $<
modi_surface_ri.o:          modi_surface_ri.F90
	$(FC) $(FFLAGS) $(FLAGS_Crocus) $(HEADER_DIRS) $<
modi_wind_threshold.o:      modi_wind_threshold.F90
	$(FC) $(FFLAGS) $(FLAGS_Crocus) $(HEADER_DIRS) $<
modn_io_offline.o:          modn_io_offline.F90
	$(FC) $(FFLAGS) $(FLAGS_Crocus) $(HEADER_DIRS) $<
parkind1.o:                 parkind1.F90
	$(FC) $(FFLAGS) $(FLAGS_Crocus) $(HEADER_DIRS) $<
snowcro.o:                  snowcro.F90
	$(FC) $(FFLAGS) $(FLAGS_Crocus) $(HEADER_DIRS) $<
sunpos.o:                   sunpos.F90
	$(FC) $(FFLAGS) $(FLAGS_Crocus) $(HEADER_DIRS) $<
surface_ri.o:               surface_ri.F90
	$(FC) $(FFLAGS) $(FLAGS_Crocus) $(HEADER_DIRS) $<
tridiag_ground_snowcro.o:   tridiag_ground_snowcro.F90
	$(FC) $(FFLAGS) $(FLAGS_Crocus) $(HEADER_DIRS) $<
yomhook.o:                  yomhook.F90
	$(FC) $(FFLAGS) $(FLAGS_Crocus) $(HEADER_DIRS) $<
#Crocus81_coldstart.o:    Crocus81_coldstart.F90
#	$(FC) $(FFLAGS) $(FLAGS_Crocus) $(HEADER_DIRS) $<
#Crocus81_dynsetup.o:     Crocus81_dynsetup.F90
#	$(FC) $(FFLAGS) $(FLAGS_Crocus) $(HEADER_DIRS) $<
#Crocus81_f2t.o:          Crocus81_f2t.F90
#	$(FC) $(FFLAGS) $(FLAGS_Crocus) $(HEADER_DIRS) $<
#Crocus81_finalize.o:     Crocus81_finalize.F90
#	$(FC) $(FFLAGS) $(FLAGS_Crocus) $(HEADER_DIRS) $<
#Crocus81_lsmMod.o:       Crocus81_lsmMod.F90
#	$(FC) $(FFLAGS) $(FLAGS_Crocus) $(HEADER_DIRS) $<
#Crocus81_main.o:         Crocus81_main.F90
#	$(FC) $(FFLAGS) $(FLAGS_Crocus) $(HEADER_DIRS) $<
#Crocus81_module.o:       Crocus81_module.F90
#	$(FC) $(FFLAGS) $(FLAGS_Crocus) $(HEADER_DIRS) $<
#Crocus81_readcrd.o:      Crocus81_readcrd.F90
#	$(FC) $(FFLAGS) $(FLAGS_Crocus) $(HEADER_DIRS) $<
#Crocus81_readrst.o:      Crocus81_readrst.F90
#	$(FC) $(FFLAGS) $(FLAGS_Crocus) $(HEADER_DIRS) $<
#Crocus81_setup.o:        Crocus81_setup.F90
#	$(FC) $(FFLAGS) $(FLAGS_Crocus) $(HEADER_DIRS) $<
#Crocus81_writerst.o:     Crocus81_writerst.F90
#	$(FC) $(FFLAGS) $(FLAGS_Crocus) $(HEADER_DIRS) $<
#crocus_driver.o:         crocus_driver.F90
#	$(FC) $(FFLAGS) $(FLAGS_Crocus) $(HEADER_DIRS) $<


!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! 
! !ROUTINE: LIS_mapSurfaceModelType
! \label{LIS_mapSurfaceModelType}
!
! !REVISION HISTORY: 
!
! 1 July 2012: Sujay Kumar, initial specification
!
! !INTERFACE: 
subroutine LIS_mapSurfaceModelType(surface_model, sf_index)

  implicit none
! !ARGUMENTS: 
  character(len=*)        :: surface_model
  integer                 :: sf_index
! 
! !DESCRIPTION: 
!  This subroutine maps the surface model type name to the
!  index used in LIS. The following convention is adopted: 
!  
!  Land Surface Model (LSM) - 1 \newline
!  Lake Model               - 2 \newline
!  Glacier Model            - 3 \newline
!  Wetland Model            - 4 \newline
!  Openwater Model          - 5 \newline
!
!EOP
  if(surface_model.eq."LSM") then 
     sf_index = 1
  elseif(surface_model.eq."Lake") then 
     sf_index = 2
  elseif(surface_model.eq."Glacier") then 
     sf_index = 3
  elseif(surface_model.eq."Wetland") then 
     sf_index = 4
  elseif(surface_model.eq."Openwater") then 
     sf_index = 5
  endif
end subroutine LIS_mapSurfaceModelType

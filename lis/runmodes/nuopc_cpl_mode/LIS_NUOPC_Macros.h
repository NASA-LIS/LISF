!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

!-------------------------------------------------------------------------------
! LIS NUOPC CPP Macros
!-------------------------------------------------------------------------------
#ifndef FILENAME
#define FILENAME __FILE__
#endif
#define CONTEXT  line=__LINE__,file=FILENAME
#define PASSTHRU msg=ESMF_LOGERR_PASSTHRU,CONTEXT
#define ESMF_STDERRORCHECK(rc) ESMF_LogFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)

!-------------------------------------------------------------------------------
! Define ESMF real kind to match Appplications single/double precision
!-------------------------------------------------------------------------------
#if defined(REAL4)
#define ESMF_KIND_FIELD ESMF_KIND_R4
#define ESMF_KIND_COORD ESMF_KIND_R4
#define ESMF_TYPEKIND_FIELD ESMF_TYPEKIND_R4
#define ESMF_TYPEKIND_COORD ESMF_TYPEKIND_R4
#elif defined(REAL8)
#define ESMF_KIND_FIELD ESMF_KIND_R8
#define ESMF_KIND_COORD ESMF_KIND_R8
#define ESMF_TYPEKIND_FIELD ESMF_TYPEKIND_R8
#define ESMF_TYPEKIND_COORD ESMF_TYPEKIND_R8
#else
#define ESMF_KIND_FIELD ESMF_KIND_R4
#define ESMF_KIND_COORD ESMF_KIND_R8
#define ESMF_TYPEKIND_FIELD ESMF_TYPEKIND_R4
#define ESMF_TYPEKIND_COORD ESMF_TYPEKIND_R8
#endif

!-------------------------------------------------------------------------------
! Define Missing Value
!-------------------------------------------------------------------------------

#define MISSINGVALUE 9.99e20_ESMF_KIND_R8


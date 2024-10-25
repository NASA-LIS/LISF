//-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
// NASA Goddard Space Flight Center
// Land Information System Framework (LISF)
// Version 7.5
//
// Copyright (c) 2024 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
//-------------------------END NOTICE -- DO NOT EDIT-----------------------
// These values define the indices used for creating the function pointer
// tables.  E.g., FT_NUM_DOMAIN is the number of domains supported by LIS.
#define FT_NUM_RUNMODE   ( 50 )
#define FT_NUM_PERTURB   ( 50 )
#define FT_NUM_RTM       ( 50 )
#define FT_NUM_LSM       ( 50 )
#define FT_NUM_LANDSLIDEMODEL       ( 50 )
#define FT_NUM_DOMAIN    ( 50 )
#define FT_NUM_DAALG     ( 50 )
#define FT_NUM_DAVAROBS  ( 50 )
#define FT_NUM_OPTUEALG  ( 50 )
#define FT_NUM_OPTUESET  ( 50 )
#define FT_NUM_OPTUETYPE ( 50 )
#define FT_NUM_OBJFUNC   ( 50 )
#define FT_NUM_UEALG     ( 50 )
#define FT_NUM_UEOBJFUNC ( 50 )
#define FT_NUM_UESET     ( 50 )
#define FT_NUM_BIASALG   ( 50 )
#define FT_NUM_METFORCE  ( 50 )
#define FT_NUM_SNOWSRC   ( 50 )
#define FT_NUM_ALBSRC    ( 50 )
#define FT_NUM_GFRACSRC  ( 50 )
#define FT_NUM_VEGSRC    ( 50 )
#define FT_NUM_LAISRC    ( 50 )
#define FT_NUM_SOILSRC   ( 50 )
#define FT_NUM_TOPOSRC   ( 50 )
#define FT_NUM_TBOTSRC   ( 50 )
#define FT_NUM_SHDMAXSRC ( 50 )
#define FT_NUM_SHDMINSRC ( 50 )
#define FT_NUM_SLOPETYPESRC ( 50 )
#define FT_NUM_ROUTING (50)

// Since C counts from 0 and LIS counts plugins from 1, add 1 to the 
// FT_NUM_* values.   These values must be used when allocating memory for
// the function pointer tables.
#define FT_MAX_RUNMODE    ( FT_NUM_RUNMODE    + 1 )
#define FT_MAX_PERTURB    ( FT_NUM_PERTURB    + 1 )
#define FT_MAX_RTM        ( FT_NUM_RTM        + 1 )
#define FT_MAX_LSM        ( FT_NUM_LSM        + 1 )
#define FT_MAX_LANDSLIDEMODEL        ( FT_NUM_LANDSLIDEMODEL        + 1 )
#define FT_MAX_DOMAIN     ( FT_NUM_DOMAIN     + 1 )
#define FT_MAX_DAALG      ( FT_NUM_DAALG      + 1 )
#define FT_MAX_OPTUEALG   ( FT_NUM_OPTUEALG   + 1 ) 
#define FT_MAX_OPTUESET   ( FT_NUM_OPTUESET   + 1 ) 
#define FT_MAX_OPTUETYPE  ( FT_NUM_OPTUETYPE  + 1 ) 
#define FT_MAX_OBJFUNC    ( FT_NUM_OBJFUNC    + 1 ) 
#define FT_MAX_UEALG      ( FT_NUM_UEALG      + 1 )
#define FT_MAX_UEOBJFUNC  ( FT_NUM_UEOBJFUNC  + 1 )
#define FT_MAX_UESET      ( FT_NUM_UESET      + 1 ) 
#define FT_MAX_BIASALG    ( FT_NUM_BIASALG    + 1 )
#define FT_MAX_DAVAROBS   ( FT_NUM_DAVAROBS   + 1 )
#define FT_MAX_METFORCE   ( FT_NUM_METFORCE  + 1 )
#define FT_MAX_SNOWSRC    ( FT_NUM_SNOWSRC    + 1 )
#define FT_MAX_ALBSRC     ( FT_NUM_ALBSRC     + 1 )
#define FT_MAX_GFRACSRC   ( FT_NUM_GFRACSRC   + 1 )
#define FT_MAX_VEGSRC     ( FT_NUM_VEGSRC     + 1 )
#define FT_MAX_LAISRC     ( FT_NUM_LAISRC     + 1 )
#define FT_MAX_SOILSRC    ( FT_NUM_SOILSRC    + 1 )
#define FT_MAX_TOPOSRC    ( FT_NUM_TOPOSRC    + 1 )
#define FT_MAX_TBOTSRC    ( FT_NUM_TBOTSRC    + 1 )
#define FT_MAX_SHDMAXSRC  ( FT_NUM_SHDMAXSRC  + 1 )
#define FT_MAX_SHDMINSRC  ( FT_NUM_SHDMINSRC  + 1 )
#define FT_MAX_SLOPETYPESRC  ( FT_NUM_SLOPETYPESRC  + 1 )
#define FT_MAX_ROUTING    (FT_NUM_ROUTING + 1)

void ft_check_index(int, int, char *);

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
#define FT_NUM_METRIC    ( 100 )

// Since C counts from 0 and LIS counts plugins from 1, add 1 to the 
// FT_NUM_* values.   These values must be used when allocating memory for
// the function pointer tables.
#define FT_MAX_METRIC     ( FT_NUM_METRIC     + 1 )
void ft_check_index(int, int, char *);

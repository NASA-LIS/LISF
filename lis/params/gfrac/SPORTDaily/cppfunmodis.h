//-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
// NASA Goddard Space Flight Center
// Land Information System Framework (LISF)
// Version 7.3
//
// Copyright (c) 2020 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
//-------------------------END NOTICE -- DO NOT EDIT-----------------------

#include "ftn_drv.h"
  extern "C" {
  int FTN(readgvfmodis)(char * infile, int * nx, int * ny, char * flag, float data_1D[], int len);
  }


//-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
// NASA Goddard Space Flight Center
// Land Information System Framework (LISF)
// Version 7.5
//
// Copyright (c) 2024 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
//-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "gpt.h"
char *t_pclstr (int code)
{

#if ( defined DISABLE_TIMERS )
  return "";
#endif

#ifdef HAVE_PCL
  switch (code) {

  case PCL_SUCCESS: 
    return "Success";
    
  case PCL_NOT_SUPPORTED:
    return "Event not supported";
    
  case PCL_TOO_MANY_EVENTS:
    return "Too many events";
    
  case PCL_TOO_MANY_NESTINGS:
    return "More nesting levels than allowed";
    
  case PCL_ILL_NESTING:
    return "Bad nesting";
    
  case PCL_ILL_EVENT:
    return "Illegal event identifier";
    
  case PCL_MODE_NOT_SUPPORTED:
    return "Mode not supported";
    
  case PCL_FAILURE:
    return "Failure for unspecified reason";
    
  default:
    return "Unknown error code";
    
  }
#endif
}

      
  

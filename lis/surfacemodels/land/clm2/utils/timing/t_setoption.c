//-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
// NASA Goddard Space Flight Center
// Land Information System Framework (LISF)
// Version 7.5
//
// Copyright (c) 2024 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
//-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include <stdio.h>
#include "gpt.h"

/*
** t_setoption: set option value to true or false.
**
** Input arguments:
**   option: option to be set
**   val:    value to which option should be set
**
** Return value: 0 (success) or -1 (failure)
*/

int t_setoption (OptionName option, Boolean val)
{
  int n;

#if ( defined DISABLE_TIMERS )
  return 0;
#endif

  if (t_initialized)
    return (t_error ("t_setoption: Options must be set BEFORE t_initialize\n"));

  for (n = 0; n < npossible; n++) {
    if (possible_event[n].name == option) {
      possible_event[n].enabled = val;

      if (val)
	printf ("t_setoption: option enabled:  %s\n", possible_event[n].string);
      else
	printf ("t_setoption: option disabled: %s\n", possible_event[n].string);

      return 0;
    }
  }

  return (t_error ("t_setoption: Option with enum index %d not available\n",
		     option));
}

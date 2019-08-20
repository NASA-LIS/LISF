//-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
// NASA Goddard Space Flight Center Land Information System (LIS) V5.0 BETA
// Released January 2008
//
// See SOFTWARE DISTRIBUTION POLICY for software distribution policies
//
// The LIS source code and documentation are in the public domain,
// available without fee for educational, research, non-commercial and
// commercial purposes.  Users may distribute the binary or source
// code to third parties provided this statement appears on all copies and
// that no charge is made for such copies.
//
// NASA GSFC MAKES NO REPRESENTATIONS ABOUT THE SUITABILITY OF THE
// SOFTWARE FOR ANY PURPOSE.  IT IS PROVIDED AS IS WITHOUT EXPRESS OR
// IMPLIED WARRANTY.  NEITHER NASA GSFC NOR THE US GOVERNMENT SHALL BE
// LIABLE FOR ANY DAMAGES SUFFERED BY THE USER OF THIS SOFTWARE.
//
// See COPYRIGHT.TXT for copyright details.
//
//-------------------------END NOTICE -- DO NOT EDIT-----------------------
// These values define the indices used for creating the function pointer
// tables.  E.g., FT_NUM_DOMAIN is the number of domains supported by LIS.
#define FT_NUM_METRIC    ( 100 )

// Since C counts from 0 and LIS counts plugins from 1, add 1 to the 
// FT_NUM_* values.   These values must be used when allocating memory for
// the function pointer tables.
#define FT_MAX_METRIC     ( FT_NUM_METRIC     + 1 )
void ft_check_index(int, int, char *);

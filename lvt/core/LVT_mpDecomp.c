//-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
// NASA Goddard Space Flight Center
// Land Information System Framework (LISF)
// Version 7.5
//
// Copyright (c) 2024 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
//-------------------------END NOTICE -- DO NOT EDIT-----------------------
//BOP
//
// !ROUTINE: LVT_mpDecomp
// \label{LVT_mpDecomp}
//
// !DESCRIPTION:
//  This file provides an algorithm for the decomposition of a given 
//  domain. The method is adopted from the Weather Research and Forecasting
//  (WRF) source code. 
//  
// !REVISION HISTORY: 
//  Feb 2006  Sujay Kumar Adopted in LIS
//EOP

#include <stdio.h>
#include "ftn_drv.h"

/* updated 20051021, new algorithm distributes the remainder, if any, at either ends of the dimension
   rather than the first remainder number of processors in the dimension. Idea is that the processes
   on the ends have less work because they're boundary processes.  New alg works like this:
                     a                         b
         + + + + + + o o o o o o o o o o o o o + + + + + +

   + represents a process with an extra point (npoints is n/p+1), o processors that don't (n/p)
   a and b are the starting process indices in the dimension of the new section of o or x.
   JM
*/

void FTN(lvt_mpdecomp) ( int *i_p , int *j_p , int *ids_p, int *ide_p , int *jds_p, int *jde_p , int *npx_p , int *npy_p , int *Px_p, int *Py_p, int *P_p )
{
  int i , j , ids, ide, jds, jde, npx, npy ;  /* inputs */
  int P, Px, Py ;                             /* output */
  int idim, jdim ;
  int rem, a, b ;

  i = *i_p - 1 ;
  j = *j_p - 1 ;
  npx = *npx_p ;
  npy = *npy_p ;
  ids = *ids_p - 1 ; ide = *ide_p - 1 ;
  jds = *jds_p - 1 ; jde = *jde_p - 1 ;
  idim = ide - ids + 1 ;
  jdim = jde - jds + 1 ;

  i = i >= ids ? i : ids ; i = i <= ide ? i : ide ;
  rem = idim % npx ;
  a = ( rem / 2 ) * ( (idim / npx) + 1 ) ;
  b = a + ( npx - rem ) * ( idim / npx ) ;
  if ( i-ids < a ) {
    Px = (i-ids) / ( (idim / npx) + 1 ) ;
  }
  else if ( i-ids < b ) {
    Px = ( a / ( (idim / npx) + 1 ) ) + (i-a-ids) / ( ( b - a ) / ( npx - rem ) )     ;
  }
  else {
    Px = ( a / ( (idim / npx) + 1 ) ) + (b-a-ids) / ( ( b - a ) / ( npx - rem ) ) +
                                        (i-b-ids) / ( ( idim / npx ) + 1 )  ;
  }

  j = j >= jds ? j : jds ; j = j <= jde ? j : jde ;
  rem = jdim % npy ;
  a = ( rem / 2 ) * ( (jdim / npy) + 1 ) ;
  b = a + ( npy - rem ) * ( jdim / npy ) ;
  if ( j-jds < a ) {
    Py = (j-jds) / ( (jdim / npy) + 1 ) ;
  }
  else if ( j-jds < b ) {
    Py = ( a / ( (jdim / npy) + 1 ) ) + (j-a-jds) / ( ( b - a ) / ( npy - rem ) )     ;
  }
  else {
    Py = ( a / ( (jdim / npy) + 1 ) ) + (b-a-jds) / ( ( b - a ) / ( npy - rem ) ) +
                                        (j-b-jds) / ( ( jdim / npy ) + 1 )  ;
  }

  *Px_p = Px ;
  *Py_p = Py ;
  *P_p = Px + Py * npx ;
}

#if 0
main()
{
  int ips[100], ipe[100] ;
  int jps[100], jpe[100] ;
  int shw, i , j , ids, ide, jds, jde, npx, npy ;  /* inputs */
  int Px, Py, P ;                             /* output */
  printf("i, j, ids, ide, jds, jde, npx, npy\n") ;
  scanf("%d %d %d %d %d %d %d %d",&i, &j, &ids,&ide,&jds,&jde,&npx,&npy ) ;
  shw =0 ;
  for ( i = 0 ; i < 100 ; i++ ) { ips[i] = 9999999 ; ipe[i] = -99999999 ; }
  for ( i = 0 ; i < 100 ; i++ ) { jps[i] = 9999999 ; jpe[i] = -99999999 ; }
#if 1
  for ( j = jds-shw ; j <= jde+shw ; j++ )
  {
  for ( i = ids-shw ; i <= ide+shw ; i++ )
  {
#endif
  TASK_FOR_POINT ( &i , &j ,
                   &ids, &ide, &jds, &jde , &npx , &npy ,
                   &Px, &Py, &P ) ;
  if ( i < ips[P] ) ips[P] = i ;
  if ( j < jps[P] ) jps[P] = j ;
  if ( i > ipe[P] ) ipe[P] = i ;
  if ( j > jpe[P] ) jpe[P] = j ;
/*  printf("%3d",P) ; */
#if 1
  }
/*  printf("\n") ; */
  }
for ( i = 0 ; i < npx*npy ; i++ ) {
  fprintf(stderr,"%3d. ips %d ipe %d (%d) jps %d jpe %d (%d)\n", i, ips[i], ipe[i], ipe[i]-ips[i]+1, jps[i], jpe[i], jpe[i]-jps[i]+1 ) ;
}
#endif
}
#endif


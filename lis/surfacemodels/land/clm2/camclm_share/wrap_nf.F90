!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!
! Purpose:
!
! Wrapper routines for the netCDF library for input and output data.
!
! Author: Jim Rosinski
!
! $Id: wrap_nf.F90,v 1.7 2004/05/07 22:18:33 jim Exp $
!
!-------------------------------------------------------------------------------
#include "LIS_misc.h"

!===============================================================================

   subroutine wrap_create (path, cmode, ncid)
   implicit none
!-------------------------------------------------------------------------------
!
! Purpose:
!
! Create a netCDF file for reading and/or writing
!
!-------------------------------------------------------------------------------
   
!#include <netcdf.inc>
   character*(*), intent(in):: path
   integer, intent(in):: cmode
   integer, intent(out):: ncid

   integer ret      ! NetCDF return code
#if 0 
   ret = nf_create (path, cmode, ncid)
   if (ret/=NF_NOERR) call handle_error (ret)
#endif
   end subroutine wrap_create

!===============================================================================

   subroutine wrap_inq_dim (nfid, dimid, dimname, dimlen)
   implicit none
!-------------------------------------------------------------------------------
!
! Purpose:
!
! Gets dimension name for a given dimension id
!
!-------------------------------------------------------------------------------
!#include <netcdf.inc>
   integer, intent(in)::  nfid
   integer, intent(in):: dimid
   integer, intent(out):: dimlen
#if ( defined ABSOFT )
   character*(*), intent(inout):: dimname
#else
   character*(*), intent(out):: dimname
#endif

   integer ret      ! NetCDF return code
#if 0 
   ret = nf_inq_dim (nfid, dimid, dimname, dimlen)
   if (ret/=NF_NOERR) call handle_error (ret)
#endif
   end subroutine wrap_inq_dim

!===============================================================================

   subroutine wrap_inq_dimid (nfid, dimname, dimid)
   implicit none
!-------------------------------------------------------------------------------
!
! Purpose:
!
! Gets the dimension id
!
!-------------------------------------------------------------------------------
!#include <netcdf.inc>
   
   integer, intent(in):: nfid
   integer, intent(out):: dimid
   character*(*), intent(in):: dimname

   integer ret      ! NetCDF return code
#if 0 
   ret = nf_inq_dimid (nfid, dimname, dimid)
   if (ret/=NF_NOERR) call handle_error (ret)
#endif
   end subroutine wrap_inq_dimid

!===============================================================================

   subroutine wrap_inq_dimlen (nfid, dimid, dimlen)
   implicit none
!-------------------------------------------------------------------------------
!
! Purpose:
!
! Gets the dimension length for a given dimension
!
!-------------------------------------------------------------------------------
!#include <netcdf.inc>
   
   integer, intent(in)::  nfid
   integer, intent(in)::  dimid 
   integer, intent(out):: dimlen
   
   integer ret      ! NetCDF return code
#if 0 
   ret = nf_inq_dimlen (nfid, dimid, dimlen)
   if (ret/=NF_NOERR) call handle_error (ret)
#endif
   end subroutine wrap_inq_dimlen

!===============================================================================

   subroutine wrap_inq_vardimid (nfid, varid, dimids)
!-------------------------------------------------------------------------------
!
! Purpose:
!
! Returns the dimension Id's from a variable
!
!-------------------------------------------------------------------------------
   implicit none
!#include <netcdf.inc>
   
   integer, intent(in):: nfid
   integer, intent(in):: varid
   integer, intent(out)::  dimids(*)
   
   integer ret      ! NetCDF return code
#if 0 
   ret = nf_inq_vardimid (nfid, varid, dimids)
   if (ret/=NF_NOERR) call handle_error (ret)
#endif
   end subroutine wrap_inq_vardimid

!===============================================================================

   subroutine wrap_inq_varid (nfid, varname, varid)
!-------------------------------------------------------------------------------
!
! Purpose:
!
! Returns the variable ID
!
!-------------------------------------------------------------------------------
   implicit none
!#include <netcdf.inc>
   
   integer, intent(in):: nfid
   integer, intent(out):: varid
   character*(*), intent(in):: varname
   
   integer ret      ! NetCDF return code
#if 0 
   ret = nf_inq_varid (nfid, varname, varid)
   if (ret/=NF_NOERR) then
     write(6,*)'wrap_inq_varid: id for ',trim(varname),' not found'
     call handle_error (ret)
   end if
#endif
   end subroutine wrap_inq_varid

!===============================================================================

   subroutine wrap_inq_var (nfid, varid, varname, xtype, ndims, &
                            dimids, natts)
!-------------------------------------------------------------------------------
!
! Purpose:
!
! Returns the variable name, type, number of dimensions, dimension ID's, and number of attributes
!
!-------------------------------------------------------------------------------
   implicit none
!#include <netcdf.inc>

   integer, intent(in):: nfid
   integer, intent(in):: varid
   integer, intent(out):: xtype
   integer, intent(out):: ndims
   integer, intent(out):: dimids(*)
   integer, intent(out):: natts
#if ( defined ABSOFT )
   character*(*), intent(inout):: varname
#else
   character*(*), intent(out):: varname
#endif

   integer ret      ! NetCDF return code
#if 0 
   ret = nf_inq_var (nfid, varid, varname, xtype, ndims, dimids, &
                     natts)
   if (ret/=NF_NOERR) call handle_error (ret)
#endif
   end subroutine wrap_inq_var

!===============================================================================

   subroutine wrap_inq_varname (nfid, varid, varname)
!-------------------------------------------------------------------------------
!
! Purpose:
!
! Returns the variable name from the dimension ID
!
!-------------------------------------------------------------------------------
   implicit none
!#include <netcdf.inc>

   integer ret      ! NetCDF return code

   integer, intent(in):: nfid
   integer, intent(in):: varid
#if ( defined ABSOFT )
   character*(*), intent(inout):: varname
#else
   character*(*), intent(out):: varname
#endif
#if 0 
   ret = nf_inq_varname (nfid, varid, varname)
   if (ret/=NF_NOERR) call handle_error (ret)
#endif
   end subroutine wrap_inq_varname

!===============================================================================

   subroutine wrap_get_att_text (nfid, varid, attname, atttext)
!-------------------------------------------------------------------------------
!
! Purpose:
!
! Returns the attribute text from the given variable ID and attribute name
!
!-------------------------------------------------------------------------------
   implicit none
!#include <netcdf.inc>
   
   integer, intent(in):: nfid
   integer, intent(in):: varid
   character*(*), intent(in):: attname
#if ( defined ABSOFT )
   character*(*), intent(inout):: atttext
#else
   character*(*), intent(out):: atttext
#endif

   integer ret      ! NetCDF return code
#if 0 
   ret = nf_get_att_text (nfid, varid, attname, atttext)
   if (ret/=NF_NOERR) call handle_error (ret)
#endif
   end subroutine wrap_get_att_text

!===============================================================================

   subroutine wrap_put_att_text (nfid, varid, attname, atttext)
!-------------------------------------------------------------------------------
!
! Purpose:
!
! Puts the given attribute text to variable ID.
!
! This routine violates the convetion that the wrapper codes take an identical
! set of arguments as the netcdf library code.  The length of the character
! argument is computed inside the wrapper.
!
!-------------------------------------------------------------------------------
   implicit none
!#include <netcdf.inc>
   
   integer, intent(in):: nfid
   integer, intent(in):: varid
   character*(*), intent(in):: attname
   character*(*), intent(in):: atttext

   integer ret      ! NetCDF return code
   integer siz
#if 0 
   siz = len_trim(atttext)
   ret = nf_put_att_text (nfid, varid, attname, siz, atttext)
   if (ret/=NF_NOERR) call handle_error (ret)
#endif
   end subroutine wrap_put_att_text

!===============================================================================

   subroutine wrap_put_att_realx (nfid, varid, attname, xtype, len, &
                                  attval)
!-------------------------------------------------------------------------------
!
! Purpose:
!
! Puts the given real attribute to the variable id
!
!-------------------------------------------------------------------------------
   use LIS_precisionMod
   implicit none
!#include <netcdf.inc>
   
   integer , intent(in):: nfid
   integer , intent(in):: varid
   integer , intent(in):: xtype
   integer , intent(in):: len
   character*(*) , intent(in):: attname
   real(r8) , intent(in):: attval

   integer ret      ! NetCDF return code
#if 0 
#ifdef CRAY
   ret = nf_put_att_real (nfid, varid, attname, xtype, len, attval)
#else
   ret = nf_put_att_real (nfid, varid, attname, xtype, len, attval)
#endif
   if (ret/=NF_NOERR) call handle_error (ret)
#endif
   end subroutine wrap_put_att_realx
!===============================================================================

   subroutine wrap_def_dim (nfid, dimname, len, dimid)
!-------------------------------------------------------------------------------
!
! Purpose:
!
! Defines the input dimension
!
!-------------------------------------------------------------------------------
   implicit none
!#include <netcdf.inc>
   integer, intent(in):: nfid
   integer, intent(in):: len
   integer, intent(out):: dimid
   character*(*), intent(in):: dimname
   
   integer ret      ! NetCDF return code
#if 0 
   ret = nf_def_dim (nfid, dimname, len, dimid)
   if (ret/=NF_NOERR) call handle_error (ret)
#endif
   end subroutine wrap_def_dim

!===============================================================================

   subroutine wrap_def_var (nfid, name, xtype, nvdims, vdims, varid)
!-------------------------------------------------------------------------------
!
! Purpose:
!
! Defines the given variable
!
!-------------------------------------------------------------------------------
   implicit none
!#include <netcdf.inc>

   integer, intent(in):: nfid
   integer, intent(in)::xtype
   integer, intent(in)::nvdims
   integer, intent(out)::varid
   integer, intent(in):: vdims(nvdims+1)
   character*(*), intent(in):: name
   
   integer ret      ! NetCDF return code
#if 0 
   ret = nf_def_var (nfid, name, xtype, nvdims, vdims, varid)
   if (ret/=NF_NOERR) call handle_error (ret)
#endif
   end subroutine wrap_def_var

!===============================================================================

   subroutine wrap_get_var_realx (nfid, varid, arr)
!-------------------------------------------------------------------------------
!
! Purpose:
!
! Gets the given real variable from a input file
!
!-------------------------------------------------------------------------------
   use LIS_precisionMod
   implicit none
!#include <netcdf.inc>

   integer, intent(in):: nfid
   integer, intent(in):: varid
   real(r8), intent(out):: arr(*)

   integer ret      ! NetCDF return code
#if 0 
#ifdef CRAY
   ret = nf_get_var_real (nfid, varid, arr)
#else
   ret = nf_get_var_real (nfid, varid, arr)
#endif
   if (ret/=NF_NOERR) then
     write(6,*)'WRAP_GET_VAR_REALX: error reading varid =', varid
     call handle_error (ret)
   end if
#endif
   end subroutine wrap_get_var_realx

!===============================================================================

   subroutine wrap_get_var_int (nfid, varid, arr)
!-------------------------------------------------------------------------------
!
! Purpose:
!
! Gets the given integer variable from a input file
!
!-------------------------------------------------------------------------------
   implicit none
!#include <netcdf.inc>

   integer, intent(in):: nfid
   integer, intent(in):: varid
   integer, intent(out):: arr(*)

   integer ret      ! NetCDF return code
#if 0 
   ret = nf_get_var_int (nfid, varid, arr)
   if (ret/=NF_NOERR) then
     write(6,*)'WRAP_GET_VAR_INT: error reading varid =', varid
     call handle_error (ret)
   end if
#endif
   end subroutine wrap_get_var_int

!===============================================================================

   subroutine wrap_get_vara_realx (nfid, varid, start, count, arr)
!-------------------------------------------------------------------------------
!
! Purpose:
!
! Gets a range of the given real variable from a input file
!
!-------------------------------------------------------------------------------
   use LIS_precisionMod
   implicit none
!#include <netcdf.inc>

   integer, intent(in):: nfid
   integer, intent(in)::varid
   integer, intent(in)::start(*)
   integer, intent(in)::count(*)
   real(r8), intent(out):: arr(*)

   integer ret      ! NetCDF return code
#if 0 
#ifdef CRAY
   ret = nf_get_vara_real (nfid, varid, start, count, arr)
#else
   ret = nf_get_vara_real (nfid, varid, start, count, arr)
#endif
   if (ret/=NF_NOERR) then
     write(6,*)'WRAP_GET_VARA_REALX: error reading varid =', varid
     call handle_error (ret)
   end if
#endif
   end subroutine wrap_get_vara_realx

!===============================================================================

   subroutine wrap_get_vara_int (nfid, varid, start, count, arr)
!-------------------------------------------------------------------------------
!
! Purpose:
!
! Gets a range of the given integer variable from a input file.
!
!-------------------------------------------------------------------------------
   implicit none
!#include <netcdf.inc>

   integer, intent(in):: nfid
   integer, intent(in):: varid
   integer, intent(in):: start(*)
   integer, intent(in):: count(*)
   integer, intent(out):: arr(*)

   integer ret      ! NetCDF return code
#if 0 
   ret = nf_get_vara_int (nfid, varid, start, count, arr)
   if (ret/=NF_NOERR) then
     write(6,*)'WRAP_GET_VARA_INT: error reading varid =', varid
     call handle_error (ret)
   end if
#endif
   end subroutine wrap_get_vara_int

!===============================================================================

   subroutine wrap_open (path, omode, ncid)
!-------------------------------------------------------------------------------
!
! Purpose:
!
! Open a netCDF file
!
!-------------------------------------------------------------------------------
   implicit none
!#include <netcdf.inc>

   character*(*), intent(in):: path
   integer, intent(in):: omode
   integer, intent(out):: ncid

   integer ret      ! NetCDF return code
#if 0 
   ret = nf_open (path, omode, ncid)
   if (ret/=NF_NOERR) then
     write(6,*)'WRAP_OPEN: nf_open failed for file ',path
     call handle_error (ret)
   end if
#endif
   end subroutine wrap_open

!===============================================================================

   subroutine wrap_close (ncid)
!-------------------------------------------------------------------------------
!
! Purpose:
!
! Close netCDF file
!
!-------------------------------------------------------------------------------
   implicit none
!#include <netcdf.inc>

   integer, intent(in):: ncid

   integer ret      ! NetCDF return code
#if 0 
   ret = nf_close (ncid)
   if (ret/=NF_NOERR) then
     write(6,*)'WRAP_CLOSE: nf_close failed for id ',ncid
     call handle_error (ret)
   end if
#endif
   end subroutine wrap_close

!===============================================================================

   subroutine wrap_put_var_int (nfid, varid, arr)
!-------------------------------------------------------------------------------
!
! Purpose:
!
! Put a integer variable on output file.
!
!-------------------------------------------------------------------------------
   implicit none
!#include <netcdf.inc>

   integer, intent(in):: nfid
   integer, intent(in):: varid
   integer, intent(in):: arr(*)

   integer ret      ! NetCDF return code
#if 0 
   ret = nf_put_var_int (nfid, varid, arr)
   if (ret/=NF_NOERR) call handle_error (ret)
#endif
   end subroutine wrap_put_var_int

!===============================================================================

   subroutine wrap_put_var1_int (nfid, varid, index, ival)
!-------------------------------------------------------------------------------
!
! Purpose:
!
! Put a variable on output file at a given index.
!
!-------------------------------------------------------------------------------
   implicit none
!#include <netcdf.inc>

   integer, intent(in):: nfid
   integer, intent(in):: varid
   integer, intent(in):: index(*)
   integer, intent(in):: ival

   integer ret      ! NetCDF return code
#if 0 
   ret = nf_put_var1_int (nfid, varid, index, ival)
   if (ret/=NF_NOERR) call handle_error (ret)
#endif
   end subroutine wrap_put_var1_int

!===============================================================================

   subroutine wrap_put_vara_int (nfid, varid, start, count, arr)
!-------------------------------------------------------------------------------
!
! Purpose:
!
! Put a range of a integer variable on a output file.
!
!-------------------------------------------------------------------------------
   implicit none
!#include <netcdf.inc>

   integer, intent(in):: nfid
   integer, intent(in):: varid
   integer, intent(in):: start(*)
   integer, intent(in):: count(*)
   integer, intent(in):: arr(*)

   integer ret      ! NetCDF return code
#if 0 
   ret = nf_put_vara_int (nfid, varid, start, count, arr)
   if (ret/=NF_NOERR) call handle_error (ret)
#endif
   end subroutine wrap_put_vara_int

!===============================================================================

   subroutine wrap_put_vara_text (nfid, varid, start, count, text)
!-------------------------------------------------------------------------------
!
! Purpose:
!
! Put a range of the given text variable to output file.
!
!-------------------------------------------------------------------------------
   implicit none
!#include <netcdf.inc>

   integer, intent(in):: nfid
   integer, intent(in):: varid
   integer, intent(in):: start(*)
   integer, intent(in):: count(*)
   character*(*), intent(in):: text(*)

   integer ret      ! NetCDF return code
#if 0 
   ret = nf_put_vara_text (nfid, varid, start, count, text)
   if (ret/=NF_NOERR) call handle_error (ret)
#endif
   end subroutine wrap_put_vara_text

!===============================================================================

   subroutine wrap_put_var1_realx (nfid, varid, index, val)
!-------------------------------------------------------------------------------
!
! Purpose:
!
! Put the given real variable to output file at given index.
!
!-------------------------------------------------------------------------------
   use LIS_precisionMod
   implicit none
!#include <netcdf.inc>

   integer, intent(in):: nfid
   integer, intent(in):: varid
   integer, intent(in):: index(*)
   real(r8), intent(in):: val

   integer ret      ! NetCDF return code
#if 0 
#ifdef CRAY
   ret = nf_put_var1_real (nfid, varid, index, val)
#else
   ret = nf_put_var1_real (nfid, varid, index, val)
#endif
   if (ret/=NF_NOERR) call handle_error (ret)
#endif
   end subroutine wrap_put_var1_realx

!===============================================================================

   subroutine wrap_put_vara_realx (nfid, varid, start, count, arr)
!-------------------------------------------------------------------------------
!
! Purpose:
!
! Output the given portion of the real array.
!
!-------------------------------------------------------------------------------
   use LIS_precisionMod
   implicit none
!#include <netcdf.inc>

   integer, intent(in):: nfid
   integer, intent(in):: varid
   integer, intent(in):: start(*)
   integer, intent(in):: count(*)
   real(r8), intent(in):: arr(*)

   integer ret      ! NetCDF return code
#if 0 
#ifdef CRAY
   ret = nf_put_vara_real (nfid, varid, start, count, arr)
#else
   ret = nf_put_vara_real (nfid, varid, start, count, arr)
#endif
   if (ret/=NF_NOERR) call handle_error (ret)
#endif
   end subroutine wrap_put_vara_realx

!===============================================================================

   subroutine wrap_put_vara_real (nfid, varid, start, count, arr)
!-------------------------------------------------------------------------------
!
! Purpose:
!
! Output the given portion of the real array.
!
!-------------------------------------------------------------------------------
   use LIS_precisionMod
   implicit none
!#include <netcdf.inc>

   integer, intent(in):: nfid
   integer, intent(in):: varid
   integer, intent(in):: start(*)
   integer, intent(in):: count(*)
   real(r4), intent(in):: arr(*)

   integer ret      ! NetCDF return code
#if 0 
   ret = nf_put_vara_real (nfid, varid, start, count, arr)
   if (ret/=NF_NOERR) call handle_error (ret)
#endif
   end subroutine wrap_put_vara_real

!===============================================================================

   subroutine wrap_put_var_realx (nfid, varid, arr)
!-------------------------------------------------------------------------------
!
! Purpose:
!
! Put the given real variable to output file.
!
!-------------------------------------------------------------------------------
   use LIS_precisionMod
   implicit none
!#include <netcdf.inc>

   integer, intent(in):: nfid
   integer, intent(in):: varid
   real(r8), intent(in):: arr(*)

   integer ret      ! NetCDF return code
#if 0 
#ifdef CRAY
   ret = nf_put_var_real (nfid, varid, arr)
#else
   ret = nf_put_var_real (nfid, varid, arr)
#endif
   if (ret/=NF_NOERR) call handle_error (ret)
#endif
   end subroutine wrap_put_var_realx

!===============================================================================

   subroutine handle_error(ret)
!-------------------------------------------------------------------------------
!
! Purpose:
!
! Handle netCDF errors.
!
!-------------------------------------------------------------------------------
   implicit none
!#include <netcdf.inc>

   integer, intent(in):: ret
#if 0    
   write(6,*)nf_strerror(ret)
   call endrun
#endif 
   end subroutine handle_error

!===============================================================================



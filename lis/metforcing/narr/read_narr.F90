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
! !ROUTINE: read_narr
! \label{read_narr}
!
! !REVISION HISTORY:
!  30 APR 2009: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine read_narr(n, m, order, narrfile)
! !USES:  
  use LIS_coreMod, only : LIS_rc, LIS_domain
  use LIS_logMod,  only : LIS_logunit, LIS_endrun
  use narr_forcingMod, only : narr_struc

  implicit none

  
! !ARGUMENTS:
!
  integer,   intent(in) :: n 
  integer,   intent(in) :: m
  integer,   intent(in) :: order
  character(len=*)      :: narrfile

! !DESCRIPTION:
!
!EOP
  
  integer               :: ftn
  logical               :: file_exists
  integer               :: i, iv


!Indexed by param identifier (prefix p), layer type (t), and then actual layer(z)
  integer, parameter    :: pPRES=1
  integer, parameter    :: pTMP=11
  integer, parameter    :: pDSWRF=204
  integer, parameter    :: pDLWRF=205
  integer, parameter    :: pPRATE=59
  integer, parameter    :: pACPCP=63
  integer, parameter    :: pUGRD=33
  integer, parameter    :: pVGRD=34
  integer, parameter    :: pSPFH=51
  integer, parameter    :: pTCDC=71

  integer, parameter    :: tSURFACE=1
  integer, parameter    :: tPROFILE=100
  integer, parameter    :: tMETERSABOVEGROUND=105
  integer, parameter    :: tHYBRIDLEVEL=109
  integer, parameter    :: tATMCOLUMN=200

  integer, parameter    :: zSURFACE=0
  integer, parameter    :: z2METER=2
  integer, parameter    :: z10METER=10
  integer, parameter    :: zHYBRIDLEVEL1=1
  integer, parameter    :: zPLEVEL(29)= (/1000,975,950,925,900,875,&
       850,825,800,775,750,725,700,650,600, & 
       550,500,450,400,350,300,275,250,225, &
       200,175,150,125,100/)

  
  inquire(file=trim(narrfile),exist=file_exists)

  if(.not.file_exists) then 
     write(LIS_logunit,*) 'NARR file',trim(narrfile), ' does not exist'
     write(LIS_logunit,*) 'Program stopping ..'
     call LIS_endrun()
  endif

  iv = 1
!surface air temp, then profile (surface to top)
  if(order.eq.1) then
     call interp_narrfield(n,trim(narrfile),iv,pTMP,&
          tMETERSABOVEGROUND,z2METER,narr_struc(n)%metdata1(iv,:))
!read vertical profiles surface to top. 
     do i=1,narr_struc(n)%nlevels
        iv = iv+1
        call interp_narrfield(n,trim(narrfile),iv,pTMP,&
             tPROFILE,zPLEVEL(i),narr_struc(n)%metdata1(iv,:))
     enddo
  elseif(order.eq.2) then 
     call interp_narrfield(n,trim(narrfile),iv,pTMP,&
          tMETERSABOVEGROUND,z2METER,narr_struc(n)%metdata2(iv,:))
!read vertical profiles surface to top. 
     do i=1,narr_struc(n)%nlevels
        iv = iv+1
        call interp_narrfield(n,trim(narrfile),iv,pTMP,&
             tPROFILE,zPLEVEL(i),narr_struc(n)%metdata2(iv,:))     
     enddo
  endif

  iv = iv+1
!surface specific humidity, then profile
  if(order.eq.1) then 
     call interp_narrfield(n,trim(narrfile),iv,pSPFH,&
          tMETERSABOVEGROUND,z2METER,narr_struc(n)%metdata1(iv,:))
     do i=1,narr_struc(n)%nlevels
        iv = iv+1
        call interp_narrfield(n,trim(narrfile),iv,pSPFH,&
             tPROFILE,zPLEVEL(i),narr_struc(n)%metdata1(iv,:))
     enddo
  elseif(order.eq.2) then 
     call interp_narrfield(n,trim(narrfile),iv,pSPFH,&
          tMETERSABOVEGROUND,z2METER,narr_struc(n)%metdata2(iv,:))
     do i=1,narr_struc(n)%nlevels
        iv = iv+1
        call interp_narrfield(n,trim(narrfile),iv,pSPFH,&
             tPROFILE,zPLEVEL(i),narr_struc(n)%metdata2(iv,:))     
     enddo
  endif

!surface pressure, then pressure levels
  iv = iv+1
  if(order.eq.1) then 
     call interp_narrfield(n,trim(narrfile),iv,pPRES,&
          tMETERSABOVEGROUND,z2METER,narr_struc(n)%metdata1(iv,:))
     do i=1,narr_struc(n)%nlevels
        iv = iv+1
        narr_struc(n)%metdata1(iv,:) = zPLEVEL(i)*100.0 !millibars to Pa
     enddo
  elseif(order.eq.2) then 
     call interp_narrfield(n,trim(narrfile),iv,pPRES,&
          tMETERSABOVEGROUND,z2METER,narr_struc(n)%metdata2(iv,:))
     do i=1,narr_struc(n)%nlevels
        iv = iv+1
        narr_struc(n)%metdata2(iv,:) = zPLEVEL(i)*100.0 !millibars to Pa
     enddo     
  endif
!swdown   
  iv = iv+1
  if(order.eq.1) then 
     call interp_narrfield(n,trim(narrfile),iv,pDSWRF,&
          tSURFACE,zSURFACE,narr_struc(n)%metdata1(iv,:))
  elseif(order.eq.2) then 
     call interp_narrfield(n,trim(narrfile),iv,pDSWRF,&
          tSURFACE,zSURFACE,narr_struc(n)%metdata2(iv,:))
  endif
!lwdown   
  iv = iv+1
  if(order.eq.1) then 
     call interp_narrfield(n,trim(narrfile),iv,pDLWRF,&
          tSURFACE,zSURFACE,narr_struc(n)%metdata1(iv,:))
  elseif(order.eq.2) then 
     call interp_narrfield(n,trim(narrfile),iv,pDLWRF,&
          tSURFACE,zSURFACE,narr_struc(n)%metdata2(iv,:))
  endif

!10m u wind
  iv = iv+1
  if(order.eq.1) then 
     call interp_narrfield(n,trim(narrfile),iv,pUGRD,&
          tMETERSABOVEGROUND,z10METER,narr_struc(n)%metdata1(iv,:))
  elseif(order.eq.2) then 
     call interp_narrfield(n,trim(narrfile),iv,pUGRD,&
          tMETERSABOVEGROUND,z10METER,narr_struc(n)%metdata2(iv,:))
  endif
!10m v wind
  iv = iv+1
  if(order.eq.1) then 
     call interp_narrfield(n,trim(narrfile),iv,pVGRD,&
          tMETERSABOVEGROUND,z10METER,narr_struc(n)%metdata1(iv,:))
  elseif(order.eq.2) then 
     call interp_narrfield(n,trim(narrfile),iv,pVGRD,&
          tMETERSABOVEGROUND,z10METER,narr_struc(n)%metdata2(iv,:))
  endif

!pcp
  iv = iv+1
  if(order.eq.1) then 
     call interp_narrfield(n,trim(narrfile),iv,pPRATE,&
          tSURFACE,zSURFACE,narr_struc(n)%metdata1(iv,:))
  elseif(order.eq.2) then 
     call interp_narrfield(n,trim(narrfile),iv,pPRATE,&
          tSURFACE,zSURFACE,narr_struc(n)%metdata2(iv,:))
  endif
!convective precip 
  iv = iv+1
  if(order.eq.1) then 
     call interp_narrfield(n,trim(narrfile),iv,pACPCP,&
          tSURFACE,zSURFACE,narr_struc(n)%metdata1(iv,:))
  elseif(order.eq.2) then 
     call interp_narrfield(n,trim(narrfile),iv,pACPCP,&
          tSURFACE,zSURFACE,narr_struc(n)%metdata2(iv,:))
  endif
  write(LIS_logunit,*) 'Read ',trim(narrfile), ' successfully '

end subroutine read_narr


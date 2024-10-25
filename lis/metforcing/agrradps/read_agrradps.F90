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
!
! !ROUTINE: read_agrradps
! \label{read_agrradps}
!
! !REVISION HISTORY:
! 29Jul2005; Sujay Kumar, Initial Code
!
! !INTERFACE:    
subroutine read_agrradps(n,m,order,yr,mo,da,hr)
! !USES: 
  use ESMF
  use LIS_coreMod,         only : LIS_rc,LIS_domain
  use agrradps_forcingMod, only : agrradps_struc
  use LIS_metforcingMod,   only : LIS_forc
  use LIS_logMod,          only : LIS_logunit, LIS_verify
  use LIS_FORC_AttributesMod 
  use LIS_metforcingMod,  only : LIS_FORC_Base_State
  use LIS_constantsMod,   only : LIS_CONST_PATH_LEN

  implicit none
! !ARGUMENTS: 
  integer,   intent(IN)    :: n 
  integer,   intent(IN)    :: m
  integer,   intent(IN)    :: order
  integer,   intent(IN)    :: yr,mo,da,hr
! 
! !DESCRIPTION: 
!  This routine opens the AGRMET files and reads the radiation 
!  fields. The fields are then spatially interpolated and 
!  gaps in the interpolation due to mismatches in landmask 
!  (between LIS and AGRMET) are filled. 
! 
!  The arguments are: 
!  \begin{description}
!   \item[n]  index of the nest
!  \item[order]
!    flag indicating which data to be read (order=1, read the previous 
!    hourly instance, order=2, read the next hourly instance)
!   \item[error]
!    output error code to indicate if the data was read and interpolated
!    successfully. 
!   \item[yr,mo,da,hr]
!    current year, month, day and hour. 
!  \end{description}
! 
!   The routines invoked are: 
!  \begin{description}
!   \item[agrradps\_filename\_sw](\ref{agrradps_filename_sw}) \newline
!    generates the name of the AGRMET file to be read
!   \item[agrradps\_filename\_cloud](\ref{agrradps_filename_cloud}) \newline
!    generates the name of the AGRMET file to be read
!   \item[agrmet2latlon1](\ref{agrmet2latlon1}) \newline
!    spatially interpolate the forcing data using bilinear interpolation
!  \end{description}
!EOP

  character(len=LIS_CONST_PATH_LEN) :: agrradpsfileNH,  agrradpsfileSH
  logical                  :: exists1,exists2,exists3,exists4
  integer                  :: c,r
  integer                  :: t
  integer                  :: error
  integer                  :: ferror
  real                     :: varfield(LIS_rc%lnc(n),LIS_rc%lnr(n))
  integer                  :: index,status
  integer                  :: fvalid
  real                     :: tair,vaporP,rldown
  type(ESMF_Field)         :: tmpField,q2Field,psurfField
  real,pointer             :: tmp(:),q2(:),psurf(:)
  real :: t2(LIS_rc%lnc(n),LIS_rc%lnr(n)),qair(LIS_rc%lnc(n),LIS_rc%lnr(n)),ps(LIS_rc%lnc(n),LIS_rc%lnr(n))     
  real :: cldamtH(LIS_rc%lnc(n), LIS_rc%lnr(n))     !high layer cloud amount [%]
  real :: cldamtM(LIS_rc%lnc(n), LIS_rc%lnr(n))     !mid  layer cloud amount [%]
  real :: cldamtL(LIS_rc%lnc(n), LIS_rc%lnr(n))     !low  layer cloud amount [%]
  real :: cldamt1d( 3 )                             !single column 3 layet cldamt [%]

  ! Get tair, qair, and psurf from the base forcing
  call ESMF_StateGet(LIS_FORC_Base_State(n,1),trim(LIS_FORC_Tair%varname(1)),tmpField,&
       rc=status)
  call LIS_verify(status)

  call ESMF_StateGet(LIS_FORC_Base_State(n,1),trim(LIS_FORC_Qair%varname(1)),q2Field,&
       rc=status)
  call LIS_verify(status)

  call ESMF_StateGet(LIS_FORC_Base_State(n,1),trim(LIS_FORC_Psurf%varname(1)),psurfField,&
       rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(tmpField,localDE=0,farrayPtr=tmp,rc=status)
  call LIS_verify(status)
 
  call ESMF_FieldGet(q2Field,localDE=0,farrayPtr=q2,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(psurfField,localDE=0,farrayPtr=psurf,rc=status)
  call LIS_verify(status)

! Convert forcing from tile to grid space
  do t=1,LIS_rc%ntiles(n)
    t2(LIS_domain(n)%tile(t)%col,LIS_domain(n)%tile(t)%row)=tmp(t)
    qair(LIS_domain(n)%tile(t)%col,LIS_domain(n)%tile(t)%row)=q2(t)
    ps(LIS_domain(n)%tile(t)%col,LIS_domain(n)%tile(t)%row)=psurf(t)
  enddo
!
! Get Shortwave radiation data
!
  error = 0
  call agrradps_filename_sw(agrradpsfileNH,agrradpsfileSH, &
       agrradps_struc(n)%agrpsdir,yr,mo,da,hr)
  inquire(file=trim(agrradpsfileNH),exist=exists1)
  inquire(file=trim(agrradpsfileSH),exist=exists2)

  if(exists1 .and. exists2) then 
  call agrmet2latlon1(n,trim(agrradpsfileNH),trim(agrradpsfileSH),varfield,ferror)
  else
     error = 1
     write(LIS_logunit,*) 'Missing AGRMET SW file NH: ',trim(agrradpsfileNH)
     write(LIS_logunit,*) 'Missing AGRMET SW file SH: ',trim(agrradpsfileSH)
  endif
  if(ferror.eq.0 .and. error.eq.0)  then 
     do r=1,LIS_rc%lnr(n)
        do c=1,LIS_rc%lnc(n)
           index = LIS_domain(n)%gindex(c,r)
           if(index.ne.-1) then 
              if(order.eq.1) then 
                 agrradps_struc(n)%metdata1(1,index) = varfield(c,r)
              elseif(order.eq.2) then 
                 agrradps_struc(n)%metdata2(1,index) = varfield(c,r)
              endif
           endif
        enddo
     enddo
  else
     if(order.eq.1) then 
        agrradps_struc(n)%metdata1(1,:) = LIS_rc%udef
     elseif(order.eq.2) then 
        agrradps_struc(n)%metdata2(1,:) = LIS_rc%udef
     endif
  endif
     
!
! Get Cloud amount data in 3 layers
!  
  error = 0
  call agrradps_filename_cloud(agrradpsfileNH,agrradpsfileSH, &
       agrradps_struc(n)%agrpsdir,yr,mo,da,hr,"H")
  
  inquire(file=trim(agrradpsfileNH),exist=exists3)
  inquire(file=trim(agrradpsfileSH),exist=exists4)

  if(exists3 .and. exists4) then 
  call agrmet2latlon1(n,trim(agrradpsfileNH),trim(agrradpsfileSH),&
        cldamtH,ferror)
  else
     error = 1
     write(LIS_logunit,*) 'Missing AGRMET CloudH file NH: ',trim(agrradpsfileNH)
     write(LIS_logunit,*) 'Missing AGRMET CloudH file SH: ',trim(agrradpsfileSH)
  endif

  call agrradps_filename_cloud(agrradpsfileNH,agrradpsfileSH, &
       agrradps_struc(n)%agrpsdir,yr,mo,da,hr,"M")
  
  inquire(file=trim(agrradpsfileNH),exist=exists3)
  inquire(file=trim(agrradpsfileSH),exist=exists4)

  if(exists3 .and. exists4) then 
  call agrmet2latlon1(n,trim(agrradpsfileNH),trim(agrradpsfileSH), &
       cldamtM,ferror)
  else
     error = error + 1
     write(LIS_logunit,*) 'Missing AGRMET CloudM file NH: ',trim(agrradpsfileNH)
     write(LIS_logunit,*) 'Missing AGRMET CloudM file SH: ',trim(agrradpsfileSH)
  endif

  call agrradps_filename_cloud(agrradpsfileNH,agrradpsfileSH, &
       agrradps_struc(n)%agrpsdir,yr,mo,da,hr,"L")
  
  inquire(file=trim(agrradpsfileNH),exist=exists3)
  inquire(file=trim(agrradpsfileSH),exist=exists4)

  if(exists3 .and. exists4) then 
  call agrmet2latlon1(n,trim(agrradpsfileNH),trim(agrradpsfileSH),&
       cldamtL,ferror)
  else
     error = error + 1
     write(LIS_logunit,*) 'Missing AGRMET CloudL file NH: ',trim(agrradpsfileNH)
     write(LIS_logunit,*) 'Missing AGRMET CloudL file SH: ',trim(agrradpsfileSH)
  endif

! initialize to udef 
  if(order.eq.1) then 
     agrradps_struc(n)%metdata1(2,:) = LIS_rc%udef
  elseif(order.eq.2) then 
     agrradps_struc(n)%metdata2(2,:) = LIS_rc%udef
  endif

  if (ferror.eq.0 .and. error.eq.0)  then 

    fvalid = 0
    do r=1, LIS_rc%lnr(n)
       do c=1, LIS_rc%lnc(n)
          if (cldamtH(c,r) >= 0.0) fvalid = 1
          if (cldamtM(c,r) >= 0.0) fvalid = 1
          if (cldamtL(c,r) >= 0.0) fvalid = 1
       end do
    end do

    if (fvalid.ne.0) then
     do r=1,LIS_rc%lnr(n)
        do c=1,LIS_rc%lnc(n)
          index = LIS_domain(n)%gindex(c,r)
          if(index.ne.-1) then 
           rldown = 0.
           tair = t2(c,r)

!=== CONVERT 2 METER SPECIFIC HUMIDITY TO VAPOR PRESSURE ==============

           vaporP =  qair(c,r)* ps(c,r)/ &
                    ( 0.622 + qair(c,r) * (1-0.622) )

!=== If tair, vaporP, and ALL 3 AGRMET LAYERS cldamt are defined, 
!===    calculate rldown; transfer to proper output array

           fvalid = 1
           if (tair         .LT. 0.0) fvalid = 0
           if (vaporP       .LT. 0.0) fvalid = 0
           if (cldamtH(c,r) .LT. 0.0) fvalid = 0
           if (cldamtM(c,r) .LT. 0.0) fvalid = 0
           if (cldamtL(c,r) .LT. 0.0) fvalid = 0
           if (fvalid == 1) then
        
!=== Extract the 3 layers cldamt at a single grid

              cldamt1d(1) = cldamtL(c,r)
              cldamt1d(2) = cldamtM(c,r)
              cldamt1d(3) = cldamtH(c,r)

              call agrlwdn( tair, vaporP, cldamt1d, rldown )
               
              if(order.eq.1) then 
                 agrradps_struc(n)%metdata1(2,index) = rldown
              elseif(order.eq.2) then 
                 agrradps_struc(n)%metdata2(2,index) = rldown
              endif

           else !(fvalid /= 1)
              if(order.eq.1) then 
                 agrradps_struc(n)%metdata1(2,index) = LIS_rc%udef
              elseif(order.eq.2) then 
                 agrradps_struc(n)%metdata2(2,index) = LIS_rc%udef
              endif
           endif !(fvalid =1);tair/vaporP/cldamt defined

          endif
        enddo
     enddo

    endif ! end if fvalid ne 0
  endif ! end if error=0

 end subroutine read_agrradps

!BOP
! !ROUTINE: fillgaps_agrradps
!  \label{fillgaps_agrradps}
! 
! !INTERFACE:
subroutine fillgaps_agrradps(n,varfield)
! !USES:
  use LIS_coreMod,         only : LIS_rc,LIS_domain
  use agrradps_forcingMod, only : agrradps_struc
  use LIS_logMod,          only : LIS_logunit, LIS_endrun

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  real, intent(inout)    :: varfield(LIS_rc%lnc(n),LIS_rc%lnr(n)) 
! !DESCRIPTION:
!   This subroutine fills in invalid grid points introduced due to 
!   reprojection from PS to lat/lon. This routine assumes that the undef
!   or invalid value is the LIS undefined value. 
! 
!  The arguments are: 
!  \begin{description}
!  \item[n]
!    index of the nest
!  \item[ip]
!    interpolation option
!  \item[varfield]
!    updated output field
!  \end{description}
!
!EOP
  integer      :: c,r
  logical      :: foundPt
  integer      :: i,j,str,enr,stc,enc,kk
  integer      :: try

  try = 0 
  if(agrradps_struc(n)%fillflag1) then !This will be done once 
     agrradps_struc(n)%smask = 0
     do r=1,LIS_rc%lnr(n)
        do c=1,LIS_rc%lnc(n)
           if((LIS_domain(n)%gindex(c,r).ne.-1).and.&
                varfield(c,r).eq.LIS_rc%udef) then !agrmet mismatch
              agrradps_struc(n)%smask(c,r) = 1
           endif
        enddo
     enddo
     agrradps_struc(n)%fillflag1 = .false. 
  endif
     
  do r=1,LIS_rc%lnr(n)
     do c=1,LIS_rc%lnc(n)
        if(agrradps_struc(n)%smask(c,r).eq.1) then 
           foundPt = .false.
           kk = 1
           try = 0 
           do while(.not.foundPt) 
              try = try +1
              str = max(r-kk,1)
              enr = min(r+kk,LIS_rc%lnr(n))
              stc = max(c-kk,1)
              enc = min(c+kk,LIS_rc%lnc(n))
              do j=str,enr
                 do i=stc,enc
                    if(LIS_domain(n)%gindex(i,j).ne.-1&
                         .and.agrradps_struc(n)%smask(i,j).ne.1) then 
                       varfield(c,r) = varfield(i,j)
                       foundPt = .true.
                       exit
                    endif
                 enddo
              enddo
              kk = kk+1
              if(try.gt.100) then 
                 write(LIS_logunit,*) 'AGRMET fillgaps failed, stopping..',try,kk,c,r
                 call LIS_endrun()
              endif
           enddo
        endif
     enddo
  enddo
     
end subroutine fillgaps_agrradps

!BOP
! !ROUTINE: agrmet2latlon1
!  \label{agrmet2latlon1}
! 
! !INTERFACE:
subroutine agrmet2latlon1 (n, nameNH, nameSH, varfield, ferror)
! !USES:
  use LIS_coreMod,         only : LIS_rc,LIS_domain
  use LIS_logMod,          only : LIS_logunit
  use agrradps_forcingMod, only : agrradps_struc

   implicit none
! !ARGUMENTS:
   integer,intent(in)          :: n
   character(len=*),intent(in) :: nameNH, nameSH
   real,intent(out)            :: varfield(LIS_rc%lnc(n),LIS_rc%lnr(n))
   integer,intent(out)         :: ferror    !set to 1 if error found
   real, parameter             :: udef = -999.  ! input data missing value
! !DESCRIPTION:
! Reads in NH and SH data, and interpolates to LIS lat/lon domain.
! Currently set to proceed only if both files exist even for processing
! one hemisphere.
!EOP
   real    :: pdata(2,agrradps_struc(n)%imax,agrradps_struc(n)%jmax) !2D input data
   integer :: i,j,c,r,iret,cnt

   integer :: openerrN=0, openerrS=0 !set to non-zero if error found
   integer :: readerrN=0, readerrS=0 !set to non-zero if error found
   real                        :: scann, scans

    write(LIS_logunit,*) 'Reading AGRMET file ',trim(nameNH)
    open(11, file=nameNH, form="unformatted", access="direct", &
         recl=4, status="old", iostat=openerrN)
     cnt = 0
     do j=1,agrradps_struc(n)%jmax
      do i=1,agrradps_struc(n)%imax
       cnt = cnt + 1
       read(11,rec=cnt,iostat=readerrN) pdata(1,i,j)
      enddo
     enddo
    close(11)

    write(LIS_logunit,*) 'Reading AGRMET file ',trim(nameSH)
    open(12, file=nameSH, form="unformatted", access="direct", &
        recl=4, status="old", iostat=openerrS)
     cnt = 0
     do j=1,agrradps_struc(n)%jmax
      do i=1,agrradps_struc(n)%imax
       cnt = cnt + 1
       read(12,rec=cnt,iostat=readerrS) pdata(2,i,j)
      enddo
     enddo
    close(12)

    if ((openerrN+openerrS+readerrN+readerrS) > 0) then
      ferror = 1
      write(LIS_logunit,*) 'AGRMET file i/o problem:', &
                 openerrN,openerrS,readerrN,readerrS
      !call LIS_endrun
    else

! scan the fields to make sure data is valid
     scann = maxval(pdata(1,:,:))
     scans = maxval(pdata(2,:,:))
     if ( scann .ne. udef .and. scans .ne. udef ) then

      call interp_agrradpsvar(n,1,pdata,udef,varfield)
      call fillgaps_agrradps(n,varfield)
      ferror = 0

     else

      write(LIS_logunit,*) 'AGRMET bad file: ',scann,scans
      ferror = 1

     endif

    endif

end subroutine agrmet2latlon1

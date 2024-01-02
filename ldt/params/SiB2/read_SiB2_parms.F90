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
! !ROUTINE: read_SiB2_parms
! \label{read_SiB2_parms}
!
! !REVISION HISTORY:
!  02 Dec 2013: KR Arsenault; Initial Specification
!
! !INTERFACE:
subroutine read_SiB2_parms(n, nvegtypes, &
                           z2,     &  ! 1
                           z1,     &  ! 2
                           vcover, &  ! 3
                           chil,   &  ! 4
                           sodep,  &  ! 5
                           rootd,  &  ! 6
                           phc,    &  ! 7
                           tran1,  &  ! 8
                           tran2,  &  ! 9
                           tran3,  &  ! 10
                           tran4,  &  ! 11
                           ref1,   &  ! 12
                           ref2,   &  ! 13
                           ref3,   &  ! 14
                           ref4,   &  ! 15
                           vmax0,  &  ! 16
                           effcon, &  ! 17
                           gradm,  &  ! 18
                           binter, &  ! 19
                           atheta, &  ! 20
                           btheta, &  ! 21
                           trda,   &  ! 22
                           trdm,   &  ! 23
                           trop,   &  ! 24
                           respcp, &  ! 25
                           slti,   &  ! 26
                           shti,   &  ! 27
                           hltii,  &  ! 28
                           hhti,   &  ! 29
                           soref1, &  ! 30
                           soref2, &  ! 31
                           bee,    &  ! 32
                           phsat,  &  ! 33
                           satco,  &  ! 34 
                           poros,  &  ! 35
                           slope,  &  ! 36
                           wopt,   &  ! 37
                           wsat,   &  ! 38
                           zm )       ! 39

! !USES:
  use ESMF
  use LDT_coreMod,     only : LDT_rc, LDT_domain
  use LDT_logMod,      only : LDT_logunit, LDT_getNextUnitNumber, &
          LDT_releaseUnitNumber, LDT_endrun
  use SiB2_parmsMod

  implicit none

! !ARGUMENTS: 
  integer,   intent(in) :: n
  integer,   intent(in) :: nvegtypes
  real,   intent(inout) :: z2(LDT_rc%lnc(n),LDT_rc%lnr(n),nvegtypes)
  real,   intent(inout) :: z1(LDT_rc%lnc(n),LDT_rc%lnr(n),nvegtypes)
  real,   intent(inout) :: vcover(LDT_rc%lnc(n),LDT_rc%lnr(n),nvegtypes)
  real,   intent(inout) :: chil(LDT_rc%lnc(n),LDT_rc%lnr(n),nvegtypes)
  real,   intent(inout) :: sodep(LDT_rc%lnc(n),LDT_rc%lnr(n),nvegtypes)
  real,   intent(inout) :: rootd(LDT_rc%lnc(n),LDT_rc%lnr(n),nvegtypes)
  real,   intent(inout) :: phc(LDT_rc%lnc(n),LDT_rc%lnr(n),nvegtypes)
  real,   intent(inout) :: tran1(LDT_rc%lnc(n),LDT_rc%lnr(n),nvegtypes)
  real,   intent(inout) :: tran2(LDT_rc%lnc(n),LDT_rc%lnr(n),nvegtypes)
  real,   intent(inout) :: tran3(LDT_rc%lnc(n),LDT_rc%lnr(n),nvegtypes)
  real,   intent(inout) :: tran4(LDT_rc%lnc(n),LDT_rc%lnr(n),nvegtypes)
  real,   intent(inout) :: ref1(LDT_rc%lnc(n),LDT_rc%lnr(n),nvegtypes)
  real,   intent(inout) :: ref2(LDT_rc%lnc(n),LDT_rc%lnr(n),nvegtypes)
  real,   intent(inout) :: ref3(LDT_rc%lnc(n),LDT_rc%lnr(n),nvegtypes)
  real,   intent(inout) :: ref4(LDT_rc%lnc(n),LDT_rc%lnr(n),nvegtypes)
  real,   intent(inout) :: vmax0(LDT_rc%lnc(n),LDT_rc%lnr(n),nvegtypes)
  real,   intent(inout) :: effcon(LDT_rc%lnc(n),LDT_rc%lnr(n),nvegtypes)
  real,   intent(inout) :: gradm(LDT_rc%lnc(n),LDT_rc%lnr(n),nvegtypes)
  real,   intent(inout) :: binter(LDT_rc%lnc(n),LDT_rc%lnr(n),nvegtypes)
  real,   intent(inout) :: atheta(LDT_rc%lnc(n),LDT_rc%lnr(n),nvegtypes)
  real,   intent(inout) :: btheta(LDT_rc%lnc(n),LDT_rc%lnr(n),nvegtypes)
  real,   intent(inout) :: trda(LDT_rc%lnc(n),LDT_rc%lnr(n),nvegtypes)
  real,   intent(inout) :: trdm(LDT_rc%lnc(n),LDT_rc%lnr(n),nvegtypes)
  real,   intent(inout) :: trop(LDT_rc%lnc(n),LDT_rc%lnr(n),nvegtypes)
  real,   intent(inout) :: respcp(LDT_rc%lnc(n),LDT_rc%lnr(n),nvegtypes)
  real,   intent(inout) :: slti(LDT_rc%lnc(n),LDT_rc%lnr(n),nvegtypes)
  real,   intent(inout) :: shti(LDT_rc%lnc(n),LDT_rc%lnr(n),nvegtypes)
  real,   intent(inout) :: hltii(LDT_rc%lnc(n),LDT_rc%lnr(n),nvegtypes)
  real,   intent(inout) :: hhti(LDT_rc%lnc(n),LDT_rc%lnr(n),nvegtypes)
  real,   intent(inout) :: soref1(LDT_rc%lnc(n),LDT_rc%lnr(n),nvegtypes)
  real,   intent(inout) :: soref2(LDT_rc%lnc(n),LDT_rc%lnr(n),nvegtypes)
  real,   intent(inout) :: bee(LDT_rc%lnc(n),LDT_rc%lnr(n),nvegtypes)
  real,   intent(inout) :: phsat(LDT_rc%lnc(n),LDT_rc%lnr(n),nvegtypes)
  real,   intent(inout) :: satco(LDT_rc%lnc(n),LDT_rc%lnr(n),nvegtypes)
  real,   intent(inout) :: poros(LDT_rc%lnc(n),LDT_rc%lnr(n),nvegtypes)
  real,   intent(inout) :: slope(LDT_rc%lnc(n),LDT_rc%lnr(n),nvegtypes)
  real,   intent(inout) :: wopt(LDT_rc%lnc(n),LDT_rc%lnr(n),nvegtypes)
  real,   intent(inout) :: wsat(LDT_rc%lnc(n),LDT_rc%lnr(n),nvegtypes)
  real,   intent(inout) :: zm(LDT_rc%lnc(n),LDT_rc%lnr(n),nvegtypes)
!
! !DESCRIPTION:
!  This subroutine retrieves the greenness fraction climatology for the 
!  specified month and returns the values in the latlon projection
!  
!  The arguments are:
!  \begin{description}
!  \item[n]
!   index of the nest
!  \item[array]
!   output field with the retrieved greenness fraction
!  \end{description}
!
!EOP      
  integer        :: c,r,v
  integer        :: nc, nr
  integer        :: ftn
  logical        :: file_exists
  character(140) :: vegfile
  character(2)   :: numvegchar
  integer        :: ip,nperiod,nfiles,i,j,nf
  real           :: vegtyp(LDT_rc%lnr(n),LDT_rc%lnc(n),nvegtypes)
! _______________________________________________________

!  parameter (nr=490,nc=1160,nperiod = 46,nfiles=11)
   nr = LDT_rc%lnr(n)
   nc = LDT_rc%lnc(n)

! _____________________________________________________

! Check for SiB2 file directory:
!  inquire(file=trim(LDT_rc%sib2parmsdir(n)), exist=file_exists)
!  if(.not.file_exists) then 
!     write(LDT_logunit,*) "SiB2 Directory: ",trim(LDT_rc%sib2parmsdir(n))," not found."
!     write(LDT_logunit,*) "Program stopping ..."
!     call LDT_endrun
!  endif

! Initialize inputs:
!  input : time independent variables from mapper.
  z2 = LDT_rc%udef
  z1 = LDT_rc%udef
  vcover = LDT_rc%udef
  chil = LDT_rc%udef
  sodep = LDT_rc%udef
  rootd = LDT_rc%udef
  phc = LDT_rc%udef
  tran1 = LDT_rc%udef
  tran2 = LDT_rc%udef
  tran3 = LDT_rc%udef
  tran4 = LDT_rc%udef
  ref1 = LDT_rc%udef
  ref2 = LDT_rc%udef
  ref3 = LDT_rc%udef
  ref4 = LDT_rc%udef
  vmax0 = LDT_rc%udef
  effcon = LDT_rc%udef
  gradm = LDT_rc%udef
  binter = LDT_rc%udef
  atheta = LDT_rc%udef
  btheta = LDT_rc%udef
  trda = LDT_rc%udef
  trdm = LDT_rc%udef
  trop = LDT_rc%udef
  respcp = LDT_rc%udef
  slti = LDT_rc%udef
  shti = LDT_rc%udef
  hltii = LDT_rc%udef
  hhti = LDT_rc%udef
  soref1 = LDT_rc%udef
  soref2 = LDT_rc%udef
  bee = LDT_rc%udef
  phsat = LDT_rc%udef
  satco = LDT_rc%udef
  poros = LDT_rc%udef
  slope = LDT_rc%udef
  wopt = LDT_rc%udef
  wsat = LDT_rc%udef
  zm = LDT_rc%udef


! Loop over the ISA/SiB2 veg classes (1 file per class):
  do v=1,nvegtypes

     write( numvegchar, '(i2.2)' ) v 
     vegfile = trim(SiB2_struc(n)%sib2parmsdir)//"/mapper_"//numvegchar//"_inv.grd"


     inquire(file=trim(vegfile), exist=file_exists)
     if( file_exists ) then 

       write(LDT_logunit,*) "Opening SiB2 file: ",trim(vegfile)
       ftn = LDT_getNextUnitNumber()
       open(ftn, file=trim(vegfile), status='unknown', &
          form='unformatted',access='direct',recl=4*nr*nc)
  
       read(ftn,rec=1)((vegtyp  (j,i,v),j=1,nc),i=1,nr)     ! 0
       read(ftn,rec=2)((z2      (j,i,v),j=1,nc),i=1,nr)     ! 1
       read(ftn,rec=3)((z1      (j,i,v),j=1,nc),i=1,nr)     ! 2
       read(ftn,rec=4)((vcover  (j,i,v),j=1,nc),i=1,nr)     ! 3
       read(ftn,rec=5)((chil    (j,i,v),j=1,nc),i=1,nr)     ! 4
       read(ftn,rec=6)((sodep   (j,i,v),j=1,nc),i=1,nr)     ! 5
       read(ftn,rec=7)((rootd   (j,i,v),j=1,nc),i=1,nr)     ! 6
       read(ftn,rec=8)((phc     (j,i,v),j=1,nc),i=1,nr)     ! 7
       read(ftn,rec=9)((tran1   (j,i,v),j=1,nc),i=1,nr)     ! 8
       read(ftn,rec=10)((tran2  (j,i,v),j=1,nc),i=1,nr)     ! 9
       read(ftn,rec=11)((tran3  (j,i,v),j=1,nc),i=1,nr)     ! 10
       read(ftn,rec=12)((tran4  (j,i,v),j=1,nc),i=1,nr)     ! 11
       read(ftn,rec=13)((ref1   (j,i,v),j=1,nc),i=1,nr)     ! 12
       read(ftn,rec=14)((ref2   (j,i,v),j=1,nc),i=1,nr)     ! 13
       read(ftn,rec=15)((ref3   (j,i,v),j=1,nc),i=1,nr)     ! 14
       read(ftn,rec=16)((ref4   (j,i,v),j=1,nc),i=1,nr)     ! 15
       read(ftn,rec=17)((vmax0  (j,i,v),j=1,nc),i=1,nr)    ! 16
       read(ftn,rec=18)((effcon (j,i,v),j=1,nc),i=1,nr)    ! 17
       read(ftn,rec=19)((gradm  (j,i,v),j=1,nc),i=1,nr)    ! 18
       read(ftn,rec=20)((binter (j,i,v),j=1,nc),i=1,nr)    ! 19
       read(ftn,rec=21)((atheta (j,i,v),j=1,nc),i=1,nr)    ! 20
       read(ftn,rec=22)((btheta (j,i,v),j=1,nc),i=1,nr)    ! 21
       read(ftn,rec=23)((trda   (j,i,v),j=1,nc),i=1,nr)    ! 22
       read(ftn,rec=24)((trdm   (j,i,v),j=1,nc),i=1,nr)    ! 23
       read(ftn,rec=25)((trop   (j,i,v),j=1,nc),i=1,nr)    ! 24
       read(ftn,rec=26)((respcp (j,i,v),j=1,nc),i=1,nr)    ! 25
       read(ftn,rec=27)((slti   (j,i,v),j=1,nc),i=1,nr)    ! 26
       read(ftn,rec=28)((shti   (j,i,v),j=1,nc),i=1,nr)    ! 27
       read(ftn,rec=29)((hltii  (j,i,v),j=1,nc),i=1,nr)    ! 28
       read(ftn,rec=30)((hhti   (j,i,v),j=1,nc),i=1,nr)    ! 29
       read(ftn,rec=31)((soref1 (j,i,v),j=1,nc),i=1,nr)    ! 30
       read(ftn,rec=32)((soref2 (j,i,v),j=1,nc),i=1,nr)    ! 31
       read(ftn,rec=33)((bee    (j,i,v),j=1,nc),i=1,nr)    ! 32
       read(ftn,rec=34)((phsat  (j,i,v),j=1,nc),i=1,nr)    ! 33
       read(ftn,rec=35)((satco  (j,i,v),j=1,nc),i=1,nr)    ! 34
       read(ftn,rec=36)((poros  (j,i,v),j=1,nc),i=1,nr)    ! 35
       read(ftn,rec=37)((slope  (j,i,v),j=1,nc),i=1,nr)    ! 36
       read(ftn,rec=38)((wopt   (j,i,v),j=1,nc),i=1,nr)    ! 37
       read(ftn,rec=39)((wsat   (j,i,v),j=1,nc),i=1,nr)    ! 38
       read(ftn,rec=40)((zm     (j,i,v),j=1,nc),i=1,nr)    ! 39

       call LDT_releaseUnitNumber(ftn)

  ! If file is missing:
    else 
       write(LDT_logunit,*) "SiB2 file: ",trim(vegfile)," not found."
       write(LDT_logunit,*) " File not read in for veg type :: ",v 
!       write(LDT_logunit,*) "Program stopping ..."
!       call LDT_endrun
     endif

  enddo  ! End veg type loop

end subroutine read_SiB2_parms

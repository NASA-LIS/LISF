!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #######################
      MODULE  MODN_IO_OFFLINE
!     #######################
!
!!****  *MODN_IO_OFFLINE* define the variables and namelist for SURFEX
!                         offline programs (pgd, prep, offline)
!!
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!      V. Masson   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    04/2008
!!      P. Lemoigne 04/2013 Add XDELTA_OROG to fix the maximum difference allowed between
!!                          forcing and surface file orographies if LSET_FORC_ZS=.F
!!      M. Dumont 12/2016 spectral calculation for Crocus CSPECSNOW
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
IMPLICIT NONE
!
!*    Types of files
!     --------------
!
 CHARACTER(LEN=6) :: CSURF_FILETYPE       = 'ASCII ' ! type of SURFEX surface files
!                                                   ! 'NETDF '
!                                                   ! 'FA    '
!                                                   ! 'ASCII '
!                                                   ! 'LFI   '
 CHARACTER(LEN=6) :: CTIMESERIES_FILETYPE = 'NONE  ' ! type of the files contining the
!                                                   ! output diagnostic time series
!                                                   ! 'NETCDF ', 'TEXTE '
 CHARACTER(LEN=6) :: CFORCING_FILETYPE    = 'NETCDF' ! type of atmospheric FORCING files
!                                                   ! 'NETDF', 'BINARY', or 'ASCII '
!
!*    Names of files
!     --------------
!
 CHARACTER(LEN=28):: CPGDFILE  ='PGD'          ! name of the PGD file
 CHARACTER(LEN=28):: CPREPFILE ='PREP'         ! name of the INITIAL file
 CHARACTER(LEN=28):: CSURFFILE ='SURFOUT'      ! name of the final output CSURFEX file
 CHARACTER(LEN=28):: CNAMELIST ='OPTIONS.nam'  ! name of namelist file
!
!
!*    General flags defining output options
!     -------------------------------------
!
LOGICAL          :: LPRINT      = .FALSE.  ! write some information on screen 
LOGICAL          :: LRESTART    = .FALSE.  ! write restart file
LOGICAL          :: LRESTART_2M = .FALSE.  ! write restart file
LOGICAL          :: LINQUIRE    = .FALSE.  ! inquiry mode
!      
LOGICAL          :: LWRITE_COORD = .FALSE. ! write lat/lon of the target grid
LOGICAL          :: LWRITE_TOPO  = .FALSE. ! write topography of the target grid
!
LOGICAL          :: LOUT_TIMENAME = .FALSE.! change the name of output file at the end of a day
                                           ! (ex: 19860502_00h00 -> 19860501_24h00)
!
LOGICAL          :: LDIAG_FA_NOCOMPACT = .FALSE. ! fa compaction for diagnostic files
!
 LOGICAL           :: LALLOW_ADD_DIM   = .FALSE. ! allow multi-dimensional output 
                                                 ! if IO scheme can deal with- e.g. XIOS
!
!*    Variables
!     ---------
!
INTEGER          :: NSCAL = 0                 ! Number of scalar species
INTEGER          :: NHALO = 0
!
!*    Time steps
!     ----------
!
REAL             :: XTSTEP_SURF   = 300.   ! time step of the surface 
REAL             :: XTSTEP_OUTPUT = 1800.  ! time step of the output time-series
INTEGER          :: NB_READ_FORC  = 0      ! subdivisions of the reading of forcings
!
!*    Allow the simulation to start from a different time step than the first record of a netcdf file
!     ----------
LOGICAL          :: LDELAYEDSTART_NC = .FALSE.
INTEGER,DIMENSION(4) :: NDATESTOP=(/0,0,0,0/) ! Year month day time (sec) to stop the simulation before the end of the netcdf forcing file
!
!*    General flag for coherence between forcing file orography and surface file orography
!     ----------
!
LOGICAL          :: LSET_FORC_ZS =.FALSE.  ! .T. : the orography of the
!                                          !  forcing file is
!                                          !  automatically set to the same
!                                          !  value as in the surface file
!                                          ! .F. : the orography of the
!                                          !  forcing file is kept as it is
REAL             :: XDELTA_OROG   = 200. ! maximum difference allowed between
!                                          ! forcing and surface file
!                                          ! orographies if LSET_FORC_ZS=.F.
!
!*    General flag for coherence between forcing Qair and calculated Qsat(Tair)
!     ----------
!
LOGICAL          :: LLIMIT_QAIR = .FALSE. ! .T. : Qair always <= Qsat(Tair)
                                          ! .F. : No limitation
!
!*    General flag for using land use scheme
!     ----------
!
LOGICAL          :: LLAND_USE = .FALSE.
!
!*    General flag for using simple coherence between solar zenithal angle and radiation
!     ----------
!
LOGICAL          :: LADAPT_SW = .FALSE.
LOGICAL          :: LINTERP_SW = .FALSE.
!
!*    General flag to modify direct solar radiation due to slopes and shadows.
!     ----------
!
LOGICAL          :: LSHADOWS_SLOPE = .FALSE.
LOGICAL          :: LSHADOWS_OTHER = .FALSE.
!
LOGICAL          :: LWR_VEGTYPE = .FALSE.
!
! * For offline driver with openMP
INTEGER         :: NPROMA                 ! Size of openMP packets
INTEGER         :: NI,NJ                  ! Domain size
!
REAL            :: XIO_FRAC = 1.          ! fraction of ISIZE deduced to I/O
!
CHARACTER(LEN=4) :: YALG_MPI = "LIN "     ! type of distribution algorithm for MPI
!
! * autorize spectral caculation for snow 
LOGICAL     :: CSPECSNOW=.FALSE.

! * forcing of impurity deposit coefficients
LOGICAL     :: LFORCIMP=.FALSE.
!
INTEGER     :: NIMPUROF=0  !Number of impurity types in the run, similar to nimpur variable but available in Offline
!
! * autorize forcing of total aerosol optical depth and ozone column
LOGICAL     :: LFORCATMOTARTES=.FALSE.
!
!-------------------------------------------------------------------------------
!
!*       1.    NAMELISTS
!              ---------
!
NAMELIST/NAM_IO_OFFLINE/CSURF_FILETYPE, CTIMESERIES_FILETYPE, CFORCING_FILETYPE, &
                        CPGDFILE, CPREPFILE, CSURFFILE, LRESTART_2M,             &
                        LPRINT, LRESTART, LINQUIRE, NSCAL, NHALO,                &
                        XTSTEP_SURF, XTSTEP_OUTPUT, LDIAG_FA_NOCOMPACT,          &
                        LSET_FORC_ZS, LWRITE_COORD, LWRITE_TOPO,                 &
                        LOUT_TIMENAME, LLIMIT_QAIR,                              &
                        LSHADOWS_SLOPE,LSHADOWS_OTHER, LWR_VEGTYPE,              &
                        NB_READ_FORC, LLAND_USE, NPROMA, NI, NJ, XIO_FRAC,       &
                        YALG_MPI, XDELTA_OROG, LADAPT_SW, LINTERP_SW,            &
                        LALLOW_ADD_DIM, LDELAYEDSTART_NC, NDATESTOP, CSPECSNOW,  &
			                  LFORCIMP,NIMPUROF, LFORCATMOTARTES              
                        
!
!-------------------------------------------------------------------------------
END MODULE MODN_IO_OFFLINE

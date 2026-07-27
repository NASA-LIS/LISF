!=====================================================================
! Program: mergefiles
! Author: Fadji Maina
! Modified by: Kristen Whitney
! Modification Date: January 16, 2025
! Description:
!   This Fortran program processes hourly LIS output files to
!   generate daily NetCDF files that include all hours of the day.
!   The program allows processing for specific date ranges, defined
!   via command-line arguments, to enable efficient batch processing.
!   The daily outputs are aligned to the same hourly structure and
!   attributes as MERRA-2 data, ensuring compatibility for analysis.
!=====================================================================
!
! Required Command-Line Arguments:
!   1. <start_year>   - Start year of the processing range (e.g., 2000)
!   2. <start_month>  - Start month of the processing range (e.g., 1 for January)
!   3. <start_day>    - Start day of the processing range (e.g., 1)
!   4. <end_year>     - End year of the processing range (e.g., 2000)
!   5. <end_month>    - End month of the processing range (e.g., 1 for January)
!   6. <end_day>      - End day of the processing range (e.g., 31)
!
! Command line arguments to compile this script:
!   module purge
!   module use --append $HOME/privatemodules
!   module load lisf_7_intel_19_1_0_166
!   ifort -c -g -check all -traceback -names lowercase -convert big_endian -assume byterecl -I /discover/nobackup/projects/lis/libs/netcdf/4.5.0_intel-19.1.0.166_sles12/include merge_files_for_date_range.f90
!   ifort -o merge_files_for_date_range merge_files_for_date_range.o -L /discover/nobackup/projects/lis/libs/netcdf/4.5.0_intel-19.1.0.166_sles12/lib/ -lnetcdf -lnetcdff
!
!
! Usage:
!   ./merge_files_for_date_range <start_year> <start_month> <start_day> <end_year> <end_month> <end_day>
!
! Example:
!   ./mergefiles 2000 1 1 2000 1 30
!   This processes files for the date range January 1–30, 2000.
!
! Steps Performed:
!   1. Parses command-line arguments to retrieve the date range.
!   2. Validates the number and format of the input arguments.
!   3. Loops through years, months, and days to process hourly LIS files
!      within the specified date range.
!   4. Reads data for Tair, Qair, Psurf, and LWdown from hourly files.
!   5. Constructs daily NetCDF files with 24-hour data arrays.
!   6. Ensures time variables and attributes match MERRA-2 specifications.
!      - Time is centered at 00:30, 01:30, ..., 23:30 GMT.
!      - Attributes include "minutes since YYYY-MM-DD 00:30:00".
!   7. Outputs daily NetCDF files named with the MERRA-2 convention.
!=====================================================================

Program mergefiles
  use netcdf
  implicit none

  ! Date range variables (to be set via command-line arguments)
  integer :: start_year, start_month, start_day,end_year, end_month, end_day
  integer :: start_date, end_date, current_date
  integer :: month_start, month_end
  integer :: day_start, day_end
  character(len=10) :: start_year_arg, start_month_arg, start_day_arg
  character(len=10) :: end_year_arg, end_month_arg, end_day_arg

  ! Declare other variables
  integer i,j,n,iyear,imonth,iday,ihour,id,jd,iaa,iaa1,ihoura,nday(12),dd
  integer, parameter:: NX=11700,NY=6500
  Character (LEN=200) fmonth,fday,fhour,filename1a,filename1,filenameo
  Character (LEN=200) varname_in_tair,varname_out_tair,longname_out_tair,units_tair
  Character (LEN=200) varname_in_qair,varname_out_qair,longname_out_qair,units_qair
  Character (LEN=200) varname_in_psurf,varname_out_psurf,longname_out_psurf,units_psurf
  Character (LEN=200) varname_in_lwdown,varname_out_lwdown,longname_out_lwdown,units_lwdown
  Character (LEN=200) indir, outdir  
  Character(LEN=50) :: year_str
  character(len=6) :: year_month
  Integer, parameter :: ndims=2,ndims1=3,ndt=24
  integer :: dimids(ndims),numLons,numLats,dimids1(ndims1)
  real*4, allocatable, dimension (:,:) :: tairlis,qairlis,psurflis,lwdownlis  ! Updated to 4-byte precision (REAL*4) 
  real*4, allocatable, dimension (:,:,:) :: tairh,qairh,psurfh,lwdownh  ! Updated to 4-byte precision (REAL*4)  
  real*4, allocatable, dimension (:) :: lat,lon,t_dt  ! Updated to 4-byte precision (REAL*4)
  integer :: ncid, tair_id,qair_id,psurf_id,lwdown_id,lat_dimid,lon_dimid,lat_id,lon_id
  integer :: x, y, ndimsp, nvarsp, nattsp, unlimdimidp,dt_dimid,dt_id
  integer, allocatable, dimension (:) :: start,count
  LOGICAL :: file_exists, force, all_hours_exist

  ! Retrieve command-line arguments
  call get_command_argument(1, start_year_arg)
  call get_command_argument(2, start_month_arg)
  call get_command_argument(3, start_day_arg)
  call get_command_argument(4, end_year_arg)
  call get_command_argument(5, end_month_arg)
  call get_command_argument(6, end_day_arg)
  call get_command_argument(7, indir)
  call get_command_argument(8, outdir)

  ! Validate the number of arguments
  if (command_argument_count() < 8) then
     print *, "Usage: ./mergefiles <start_year> <start_month> <start_day> <end_year> <end_month> <end_day> <indir> <outdir>"
     stop "Insufficient arguments provided."
  end if

  ! Validate argument format
  ! (Assumes an is_numeric function is defined in the contains section)
  if (.NOT. is_numeric(start_year_arg)) then
     print *, "Error: Start year must be numeric. Provided: ", start_year_arg
     stop "Invalid start year."
  end if
  if (.NOT. is_numeric(start_month_arg)) then
     print *, "Error: Start month must be numeric. Provided: ", start_month_arg
     stop "Invalid start month."
  end if
  if (.NOT. is_numeric(start_day_arg)) then
     print *, "Error: Start day must be numeric. Provided: ", start_day_arg
     stop "Invalid start day."
  end if
  if (.NOT. is_numeric(end_year_arg)) then
     print *, "Error: End year must be numeric. Provided: ", end_year_arg
     stop "Invalid end year."
  end if
  if (.NOT. is_numeric(end_month_arg)) then
     print *, "Error: End month must be numeric. Provided: ", end_month_arg
     stop "Invalid end month."
  end if
  if (.NOT. is_numeric(end_day_arg)) then
     print *, "Error: End day must be numeric. Provided: ", end_day_arg
     stop "Invalid end day."
  end if
  if (trim(indir) == '') then
     print *, "Error: Input directory (indir) not provided."
     stop
  end if
  if (trim(outdir) == '') then
     print *, "Error: Output directory (indir) not provided."
     stop
  end if

  ! Convert string arguments to integers
  read(start_year_arg, *) start_year
  read(start_month_arg, *) start_month
  read(start_day_arg, *) start_day
  read(end_year_arg, *) end_year
  read(end_month_arg, *) end_month
  read(end_day_arg, *) end_day

  ! Add trailing slashes to input/output directories need be
  if (indir(len_trim(indir):len_trim(indir)) /= '/') indir = trim(indir)//'/'
  if (outdir(len_trim(outdir):len_trim(outdir)) /= '/') outdir = trim(outdir)//'/'

  ! Calculate start_date and end_date in YYYYMMDD format
  start_date = start_year * 10000 + start_month * 100 + start_day
  end_date = end_year * 10000 + end_month * 100 + end_day

  ! Allocate arrays
  allocate (start(ndims),count(ndims),tairlis(NX,NY),t_dt(24),qairlis(NX,NY),psurflis(NX,NY),lwdownlis(NX,NY))
  allocate (tairh(NX,NY,24))
  allocate (qairh(NX,NY,24))
  allocate (psurfh(NX,NY,24))
  allocate (lwdownh(NX,NY,24))
  allocate (lat(NY),lon(NX))

  ! Initialize arrays to default values
  tairh=0.0
  qairh=0.0
  psurfh=0.0
  lwdownh=0.0

  ! Initialize time array
  do i=1,24
     t_dt(i) = (i-1) * 60  ! Time in minutes since 00:30 (centered)
  end do

  ! Set the strings
  force=.FALSE. ! FALSE skips existing output files, TRUE re-creates/overwrites existing output files
  varname_in_tair='Tair_f_tavg'
  varname_out_tair='Tair'
  longname_out_tair='2-m air temperature'
  units_tair='K'
  varname_in_qair='Qair_f_tavg'
  varname_out_qair='Qair'
  longname_out_qair='specific humidity'
  units_qair='kg kg-1'
  varname_in_psurf='Psurf_f_tavg'
  varname_out_psurf='Psurf'
  longname_out_psurf='surface pressure'
  units_psurf='Pa'
  varname_in_lwdown='LWdown_f_tavg'
  varname_out_lwdown='LWdown'
  longname_out_lwdown='surface downward longwave radiation'
  units_lwdown='W m-2'

  ! allocate dimensions
  nday(1)=31	;	nday(2)=29	;	nday(3)=31   ;  nday(4)=30
  nday(5)=31        ;       nday(6)=30      ;       nday(7)=31   ;  nday(8)=31 
  nday(9)=30        ;       nday(10)=31     ;       nday(11)=30  ;  nday(12)=31
  do i=1,nx
     do j=1,ny
        lat(j)=7.005+(j-1)*0.01
        lon(i)=-168.995+(i-1)*0.01
        !lat1(i,j)=7.005+(j-1)*0.01
        !lon1(i,j)=-168.995+(i-1)*0.01
     end do
  end do
  dd=0

  ! Begin loops through dates/times
  DO iyear = start_year, end_year
     ! Handle leap years
     IF (MOD(iyear, 4) == 0 .AND. (MOD(iyear, 100) /= 0 .OR. MOD(iyear, 400) == 0)) THEN
        nday(2) = 29  ! February has 29 days in a leap year
     ELSE
        nday(2) = 28  ! February has 28 days in a non-leap year
     END IF
     WRITE(year_str, "(I4)") iyear

     ! Determine month bounds
     IF (iyear == start_year) THEN 
        month_start = start_month
     ELSE
        month_start = 1
     END IF
     IF (iyear == end_year) THEN
        month_end = end_month
     ELSE
        month_end = 12
     END IF

     ! Loop through month
     DO imonth = month_start, month_end
        WRITE(fmonth, "(I2.2)") imonth
        WRITE(year_month, "(I4,I2.2)") iyear, imonth

        ! Determine day bounds
        IF (iyear == start_year .AND. imonth == start_month) THEN
           day_start = start_day
        ELSE
           day_start = 1
        END IF
        IF (iyear == end_year .AND. imonth == end_month) THEN
           day_end = MIN(end_day, nday(imonth))  ! Ensure day_end does not exceed
        ELSE
           day_end = nday(imonth)
        END IF

        ! Loop through the days
        DO iday = day_start, day_end
           ! Write day string
           WRITE(fday, "(I2.2)") iday
           PRINT *, "Processing date: Year=", iyear, "Month=", imonth, "Day=", iday

           ! Construct output filename
           filename1 = 'LIS_HIST_' // TRIM(year_str) // TRIM(fmonth) // TRIM(fday) // '0030.d01.nc'
           filenameo = TRIM(outdir) // filename1

           ! Check if output file exists
           PRINT *, "Checking existance of output File: ", filenameo
           INQUIRE(FILE=TRIM(filenameo), EXIST=file_exists)
           IF (file_exists .AND. .NOT. force) THEN
              PRINT *, "File exists and force is FALSE. Skipping creation of file: ", filenameo
              CYCLE  ! Skip this day
           END IF

           ! Loop to open hourly files and populate output array
           all_hours_exist = .TRUE.
           DO ihoura = 1, 24
              ihour = ihoura - 1  ! Adjust hour for 0-based indexing

              ! Write hour string
              WRITE(fhour, "(I2.2)") ihour

              ! Construct input hourly filename
              filename1 = 'LIS_HIST_' // TRIM(year_str) // TRIM(fmonth) // TRIM(fday) // TRIM(fhour) // '00.d01.nc'
              filename1a = TRIM(indir) // TRIM(year_month) // "/" // filename1

              ! Check if file exists and if not, skip the day
              INQUIRE(FILE=TRIM(filename1a), EXIST=file_exists)
              IF (.NOT. file_exists) THEN
                 PRINT *, "Skipping day, b/c hourly file does not exist: ", TRIM(filename1a)
                 all_hours_exist = .FALSE.
                 EXIT  ! Exit the hour loop and go to the next day
              END IF
              IF (.NOT. all_hours_exist) CYCLE ! Skip the rest of the day loop

              call check(NF90_OPEN(filename1a, NF90_NOWRITE, ncid) )
              call check( nf90_inquire (ncid, ndimsp, nvarsp, nattsp, unlimdimidp))
              call check( nf90_inq_varid(ncid, varname_in_tair, tair_id) )
              call check( nf90_inquire_variable(ncid, tair_id, dimids = dimIDs))
              call check( nf90_inq_varid(ncid, varname_in_qair, qair_id) )
              call check( nf90_inquire_variable(ncid, qair_id, dimids = dimIDs))
              call check( nf90_inq_varid(ncid, varname_in_psurf, psurf_id) )
              call check( nf90_inquire_variable(ncid, psurf_id, dimids = dimIDs))
              call check( nf90_inq_varid(ncid, varname_in_lwdown, lwdown_id) )
              call check( nf90_inquire_variable(ncid, lwdown_id, dimids = dimIDs))
              call check (nf90_inquire_dimension(ncid, dimIDs(1), len = numLons))
              call check (nf90_inquire_dimension(ncid, dimIDs(2), len = numLats))
              dimids = (/numLons, numLats /)
              call check( nf90_inq_varid(ncid, varname_in_tair, tair_id))
              call check( nf90_inq_varid(ncid, varname_in_qair, qair_id))
              call check( nf90_inq_varid(ncid, varname_in_psurf, psurf_id))
              call check( nf90_inq_varid(ncid, varname_in_lwdown, lwdown_id))
              count = 1
              start =1
              call check( nf90_get_var(ncid, tair_id, tairlis))
              call check( nf90_get_var(ncid, qair_id, qairlis))
              call check( nf90_get_var(ncid, psurf_id, psurflis))
              call check( nf90_get_var(ncid, lwdown_id, lwdownlis))
              call check( nf90_close(ncid))
              do i=1,nx
                 do j=1,ny
                    tairh(i,j,ihoura)=tairlis(i,j)
                    qairh(i,j,ihoura)=qairlis(i,j)
                    psurfh(i,j,ihoura)=psurflis(i,j)
                    lwdownh(i,j,ihoura)=lwdownlis(i,j)
                 end do
              end do
           END DO !HOUR

           IF (.NOT. all_hours_exist) CYCLE ! Skip the rest of the day loop
           PRINT *, "Creating file: ", filenameo

           ! Initialize the file
           call check(NF90_create(filenameo,IOR(NF90_NETCDF4, NF90_CLOBBER),ncid))

           ! Define lat and lon dimensions
           call check(nf90_def_dim(ncid, 'latitude', NY, lat_dimid))  ! Latitude dimension (1D)
           call check(nf90_def_dim(ncid, 'longitude', NX, lon_dimid))  ! Longitude dimension (1D)
           call check(nf90_def_var(ncid, 'latitude', NF90_FLOAT, (/lat_dimid/), lat_id))  ! Updated to 1D latitude
           call check(nf90_def_var(ncid, 'longitude', NF90_FLOAT, (/lon_dimid/), lon_id))  ! Updated to 1D longitude
           dimids=(/lon_dimid,lat_dimid/)
           CALL check(nf90_put_att(ncid,lat_id,"units","degree_north"))
           call check(nf90_put_att(ncid,lat_id,"_FillValue", -9999.00))
           CALL check(nf90_put_att(ncid,lat_id,"long_name","latitude"))
           CALL check(nf90_put_att(ncid,lon_id,"units","degree_east"))
           call check(nf90_put_att(ncid,lon_id,"_FillValue", -9999.00))
           CALL check(nf90_put_att(ncid,lon_id,"long_name","longitude"))
           call check(nf90_put_var(ncid, lat_id, lat))
           call check(nf90_put_var(ncid, lon_id, lon))

           ! Define the time dimension and variable
           call check(nf90_def_dim(ncid, 'time', ndt, dt_dimid))  ! Define time dimension
           call check(nf90_def_var(ncid, 'time', NF90_INT, (/dt_dimid/), dt_id))  ! Time variable

           ! Add time variable attributes
           call check(nf90_put_att(ncid, dt_id, "long_name", "time"))
           call check(nf90_put_att(ncid, dt_id, "units", "minutes since "//TRIM(year_str)//"-"//TRIM(fmonth(14:15))//"-"//TRIM(fday(16:17))//" 00:30:00"))
           call check(nf90_put_att(ncid, dt_id, "time_increment", 10000))
           call check(nf90_put_att(ncid, dt_id, "begin_date", (iyear) * 10000 + imonth * 100 + iday))  ! Correct YYYYMMDD
           call check(nf90_put_att(ncid, dt_id, "begin_time", 3000))  ! 00:30 in HHMM

           ! Write time variable data
           call check(nf90_put_var(ncid, dt_id, t_dt))

           ! Define the variables upfront
           dimids1=(/lon_dimid,lat_dimid,dt_dimid/)
           call check(nf90_def_var(ncid,varname_out_tair,NF90_FLOAT,dimids1,tair_id,shuffle=.TRUE.,deflate_level=1))  ! Ensured single precision (NF90_FLOAT)
           call check(nf90_def_var(ncid,varname_out_qair,NF90_FLOAT,dimids1,qair_id,shuffle=.TRUE.,deflate_level=1))  ! Ensured single precision (NF90_FLOAT)
           call check(nf90_def_var(ncid,varname_out_psurf,NF90_FLOAT,dimids1,psurf_id,shuffle=.TRUE.,deflate_level=1))  ! Ensured single precision (NF90_FLOAT)
           call check(nf90_def_var(ncid,varname_out_lwdown,NF90_FLOAT,dimids1,lwdown_id,shuffle=.TRUE.,deflate_level=1))  ! Ensured single precision (NF90_FLOAT)

           ! Add attributes for variables
           CALL check(nf90_put_att(ncid,tair_id,"units",units_tair))
           call check(nf90_put_att(ncid,tair_id,"_FillValue", -9999.00))
           CALL check(nf90_put_att(ncid,tair_id,"long_name",longname_out_tair))
           CALL check(nf90_put_att(ncid,qair_id,"units",units_qair))
           call check(nf90_put_att(ncid,qair_id,"_FillValue", -9999.00))
           CALL check(nf90_put_att(ncid,qair_id,"long_name",longname_out_qair))
           CALL check(nf90_put_att(ncid,psurf_id,"units",units_psurf))
           call check(nf90_put_att(ncid,psurf_id,"_FillValue", -9999.00))
           CALL check(nf90_put_att(ncid,psurf_id,"long_name",longname_out_psurf))
           CALL check(nf90_put_att(ncid,lwdown_id,"units",units_lwdown))
           call check(nf90_put_att(ncid,lwdown_id,"_FillValue", -9999.00))
           CALL check(nf90_put_att(ncid,lwdown_id,"long_name",longname_out_lwdown))

           ! Exit define mode
           CALL check(nf90_enddef(ncid))

           ! Write data to each variable
           call check(nf90_put_var(ncid,tair_id,tairh))
           call check(nf90_put_var(ncid,qair_id,qairh))
           call check(nf90_put_var(ncid,psurf_id,psurfh))
           call check(nf90_put_var(ncid,lwdown_id,lwdownh))

           ! Close the file
           CALL check(nf90_close(ncid))  
        END DO	!DAY
     END DO	!MONTH
  END DO	!YEAR

14 FORMAT(2000(1x,e14.7))

contains

  subroutine check(status)
    integer, intent ( in) :: status

    if(status /= nf90_noerr) then
       print *, trim(nf90_strerror(status))
       stop "Stopped"
    end if
  end subroutine check

  logical function is_numeric(str)
    character(len=*), intent(in) :: str
    integer :: num
    logical :: read_success

    read_success = .FALSE.
    read(str, '(I)', IOSTAT=num) num
    is_numeric = (num == 0)  ! If IOSTAT = 0, the string is numeric
  end function is_numeric

end Program mergefiles

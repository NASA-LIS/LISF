module ac_project_input

use ac_kinds, only: dp, &
                    int8, &
                    int32
use ac_utils, only: assert
use iso_fortran_env, only: iostat_end
implicit none


type ProjectInput_type
    !! Container for project file input data
    real(dp) :: VersionNr
        !! AquaCrop version number (common for all runs)
    character(len=:), allocatable :: Description
        !! Project description (common for all runs)
    integer(int8)  :: Simulation_YearSeason
        !! Year number of cultivation (1 = seeding/planting year)
    integer(int32) :: Simulation_DayNr1
        !! First day of simulation period
    integer(int32) :: Simulation_DayNrN
        !! Last day of simulation period
    integer(int32) :: Crop_Day1
        !! First day of cropping period
    integer(int32) :: Crop_DayN
        !! Last day of cropping period
    character(len=:), allocatable :: Climate_Info
        !! Climate info
    character(len=:), allocatable :: Climate_Filename
        !! Climate file name
    character(len=:), allocatable :: Climate_Directory
        !! Climate file directory
    character(len=:), allocatable :: Temperature_Info
        !! Temperature info
    character(len=:), allocatable :: Temperature_Filename
        !! Temperature file name
    character(len=:), allocatable :: Temperature_Directory
        !! Temperature file directory
    character(len=:), allocatable :: ETo_Info
        !! ETo info
    character(len=:), allocatable :: ETo_Filename
        !! ETo file name
    character(len=:), allocatable :: ETo_Directory
        !! ETo file directory
    character(len=:), allocatable :: Rain_Info
        !! Rain info
    character(len=:), allocatable :: Rain_Filename
        !! Rain file name
    character(len=:), allocatable :: Rain_Directory
        !! Rain file directory
    character(len=:), allocatable :: CO2_Info
        !! CO2 info
    character(len=:), allocatable :: CO2_Filename
        !! CO2 file name
    character(len=:), allocatable :: CO2_Directory
        !! CO2 file directory
    character(len=:), allocatable :: Calendar_Info
        !! Calendar info
    character(len=:), allocatable :: Calendar_Filename
        !! Calendar file name
    character(len=:), allocatable :: Calendar_Directory
        !! Calendar file directory
    character(len=:), allocatable :: Crop_Info
        !! Crop info
    character(len=:), allocatable :: Crop_Filename
        !! Crop file name
    character(len=:), allocatable :: Crop_Directory
        !! Crop file directory
    character(len=:), allocatable :: Irrigation_Info
        !! Irrigation info
    character(len=:), allocatable :: Irrigation_Filename
        !! Irrigation file name
    character(len=:), allocatable :: Irrigation_Directory
        !! Irrigation file directory
    character(len=:), allocatable :: Management_Info
        !! Management info
    character(len=:), allocatable :: Management_Filename
        !! Management file name
    character(len=:), allocatable :: Management_Directory
        !! Management file directory
    character(len=:), allocatable :: GroundWater_Info
        !! Groundwater info
    character(len=:), allocatable :: GroundWater_Filename
        !! Groundwater file name
    character(len=:), allocatable :: GroundWater_Directory
        !! Groundwater file directory
    character(len=:), allocatable :: Soil_Info
        !! Soil info
    character(len=:), allocatable :: Soil_Filename
        !! Soil file name
    character(len=:), allocatable :: Soil_Directory
        !! Soil file directory
    character(len=:), allocatable :: SWCIni_Info
        !! SWCIni info
    character(len=:), allocatable :: SWCIni_Filename
        !! SWCIni file name
    character(len=:), allocatable :: SWCIni_Directory
        !! SWCIni file directory
    character(len=:), allocatable :: OffSeason_Info
        !! OffSeason info
    character(len=:), allocatable :: OffSeason_Filename
        !! OffSeason file name
    character(len=:), allocatable :: OffSeason_Directory
        !! OffSeason file directory
    character(len=:), allocatable :: Observations_Info
        !! Observations info
    character(len=:), allocatable :: Observations_Filename
        !! Observations file name
    character(len=:), allocatable :: Observations_Directory
        !! Observations file directory
    contains
    procedure :: read_project_file
end type ProjectInput_type


interface get_project_input
    module procedure get_project_input_dp
    module procedure get_project_input_int8
    module procedure get_project_input_int32
    module procedure get_project_input_string
end interface get_project_input


interface set_project_input
    module procedure set_project_input_dp
    module procedure set_project_input_int8
    module procedure set_project_input_int32
    module procedure set_project_input_string
end interface set_project_input


type(ProjectInput_type), dimension(:), allocatable :: ProjectInput
    !! Project file input data for every run in a project


contains


subroutine allocate_project_input(NrRuns)
    !! Simply (re)allocates the ProjectInput module variable,
    integer(int32), intent(in) :: NrRuns
        !! total number of runs

    if (allocated(ProjectInput)) then
        deallocate(ProjectInput)
    end if

    allocate(ProjectInput(NrRuns))
end subroutine allocate_project_input


subroutine initialize_project_input(filename, NrRuns)
    !! Initializes the ProjectInput module variable,
    !! if it has not yet been allocated.
    character(len=*), intent(in) :: filename
        !! PRM or PRO file name
    integer(int32), intent(in), optional :: NrRuns
        !! total number of runs (if known beforehand)

    integer :: i, NrRuns_local

    if (present(NrRuns)) then
        NrRuns_local = NrRuns
    else
        call ReadNumberSimulationRuns(filename, NrRuns_local)
    end if
    call allocate_project_input(NrRuns_local)

    do i = 1, NrRuns_local
        call ProjectInput(i)%read_project_file(filename, i)
    end do
end subroutine initialize_project_input


subroutine ReadNumberSimulationRuns(TempFileNameFull, NrRuns)
    !! Reads the project file to get the total number of runs.
    character(len=*), intent(in) :: TempFileNameFull
        !! PRM or PRO file name
    integer(int32), intent(out) :: NrRuns
        !! total number of runs

    integer :: fhandle
    integer(int32) :: NrFileLines, rc, i

    NrRuns = 1

    open(newunit=fhandle, file=trim(TempFileNameFull), status='old', &
         action='read', iostat=rc)
    read(fhandle, *, iostat=rc)  ! Description
    read(fhandle, *, iostat=rc)  ! AquaCrop version Nr

    do i = 1, 5
        read(fhandle, *, iostat=rc) ! Type year and Simulation and Cropping period Run 1
    end do

    NrFileLines = 42 ! Clim(15),Calendar(3),Crop(3),Irri(3),Field(3),Soil(3),Gwt(3),Inni(3),Off(3),FieldData(3)
    do i = 1, NrFileLines
        read(fhandle, *, iostat=rc) ! Files Run 1
    end do

    read_loop: do
        i = 0
        do while (i < (NrFileLines+5))
            read(fhandle, *, iostat=rc)
            if (rc == iostat_end) exit read_loop
            i = i + 1
        end do

        if (i == (NrFileLines+5)) then
            NrRuns = NrRuns + 1
        end if
    end do read_loop
    close(fhandle)
end subroutine ReadNumberSimulationRuns


function GetNumberSimulationRuns() result(NrRuns)
    !! Returns the total number of runs.
    integer :: NrRuns

    NrRuns = size(ProjectInput)
end function GetNumberSimulationRuns


subroutine read_project_file(self, filename, NrRun)
    !! Reads in the project file contents that apply to the given run index.
    class(ProjectInput_type), intent(out) :: self
    character(len=*), intent(in) :: filename
        !! PRM or PRO file name
    integer, intent(in) :: NrRun
        !! Run index (should be 1 in the case of a PRO file)

    integer :: fhandle, i, rc, Runi
    character(len=1024) :: buffer

    open(newunit=fhandle, file=trim(filename), status='old', action='read', &
         iostat=rc)
    read(fhandle, '(a)', iostat=rc) buffer
    self%Description = trim(buffer)
    read(fhandle, *, iostat=rc) self%VersionNr ! AquaCrop version Nr

    if (NrRun > 1) then
        ! Skip sections belonging to previous runs
        do Runi = 1, (NrRun - 1)
            do i = 1, 47
                read(fhandle, *, iostat=rc) ! 5 + 42 lines with files
            end do
        end do
    end if

    ! 0. Year of cultivation and Simulation and Cropping period
    read(fhandle, *, iostat=rc) self%Simulation_YearSeason
    read(fhandle, *, iostat=rc) self%Simulation_DayNr1
    read(fhandle, *, iostat=rc) self%Simulation_DayNrN
    read(fhandle, *, iostat=rc) self%Crop_Day1
    read(fhandle, *, iostat=rc) self%Crop_DayN

    ! 1. Climate
    read(fhandle, '(a)', iostat=rc) buffer
    self%Climate_Info = trim(buffer)
    read(fhandle, *, iostat=rc) buffer
    self%Climate_Filename = trim(buffer)
    read(fhandle, *, iostat=rc) buffer
    self%Climate_Directory = trim(buffer)

    ! 1.1 Temperature
    read(fhandle, '(a)', iostat=rc) buffer
    self%Temperature_Info = trim(buffer)
    read(fhandle, *, iostat=rc) buffer
    self%Temperature_Filename = trim(buffer)
    read(fhandle, *, iostat=rc) buffer
    self%Temperature_Directory = trim(buffer)

    ! 1.2 ETo
    read(fhandle, '(a)', iostat=rc) buffer
    self%ETo_Info = trim(buffer)
    read(fhandle, *, iostat=rc) buffer
    self%ETo_Filename = trim(buffer)
    read(fhandle, *, iostat=rc) buffer
    self%ETo_Directory = trim(buffer)

    ! 1.3 Rain
    read(fhandle, '(a)', iostat=rc) buffer
    self%Rain_Info = trim(buffer)
    read(fhandle, *, iostat=rc) buffer
    self%Rain_Filename = trim(buffer)
    read(fhandle, *, iostat=rc) buffer
    self%Rain_Directory = trim(buffer)

    ! 1.4 CO2
    read(fhandle, '(a)', iostat=rc) buffer
    self%CO2_Info = trim(buffer)
    read(fhandle, *, iostat=rc) buffer
    self%CO2_Filename = trim(buffer)
    read(fhandle, *, iostat=rc) buffer
    self%CO2_Directory = trim(buffer)

    ! 2. Calendar
    read(fhandle, '(a)', iostat=rc) buffer
    self%Calendar_Info = trim(buffer)
    read(fhandle, *, iostat=rc) buffer
    self%Calendar_Filename = trim(buffer)
    read(fhandle, *, iostat=rc) buffer
    self%Calendar_Directory = trim(buffer)

    ! 3. Crop
    read(fhandle, '(a)', iostat=rc) buffer
    self%Crop_Info = trim(buffer)
    read(fhandle, *, iostat=rc) buffer
    self%Crop_Filename = trim(buffer)
    read(fhandle, *, iostat=rc) buffer
    self%Crop_Directory = trim(buffer)

    ! 4. Irrigation
    read(fhandle, '(a)', iostat=rc) buffer
    self%Irrigation_Info = trim(buffer)
    read(fhandle, *, iostat=rc) buffer
    self%Irrigation_Filename = trim(buffer)
    read(fhandle, *, iostat=rc) buffer
    self%Irrigation_Directory = trim(buffer)

    ! 5. Field Management
    read(fhandle, '(a)', iostat=rc) buffer
    self%Management_Info = trim(buffer)
    read(fhandle, *, iostat=rc) buffer
    self%Management_Filename = trim(buffer)
    read(fhandle, *, iostat=rc) buffer
    self%Management_Directory = trim(buffer)

    ! 6. Soil Profile
    read(fhandle, '(a)', iostat=rc) buffer
    self%Soil_Info = trim(buffer)
    read(fhandle, *, iostat=rc) buffer
    self%Soil_Filename = trim(buffer)
    read(fhandle, *, iostat=rc) buffer
    self%Soil_Directory = trim(buffer)

    ! 7. GroundWater
    read(fhandle, '(a)', iostat=rc) buffer
    self%GroundWater_Info = trim(buffer)
    read(fhandle, *, iostat=rc) buffer
    self%GroundWater_Filename = trim(buffer)
    read(fhandle, *, iostat=rc) buffer
    self%GroundWater_Directory = trim(buffer)

    ! 8. Initial conditions
    read(fhandle, '(a)', iostat=rc) buffer
    self%SWCIni_Info = trim(buffer)
    read(fhandle, *, iostat=rc) buffer
    self%SWCIni_Filename = trim(buffer)
    read(fhandle, *, iostat=rc) buffer
    self%SWCIni_Directory = trim(buffer)

    ! 9. Off-season conditions
    read(fhandle, '(a)', iostat=rc) buffer
    self%OffSeason_Info = trim(buffer)
    read(fhandle, *, iostat=rc) buffer
    self%OffSeason_Filename = trim(buffer)
    read(fhandle, *, iostat=rc) buffer
    self%OffSeason_Directory = trim(buffer)

    ! 10. Field data
    read(fhandle, '(a)', iostat=rc) buffer
    self%Observations_Info = trim(buffer)
    read(fhandle, *, iostat=rc) buffer
    self%Observations_Filename = trim(buffer)
    read(fhandle, *, iostat=rc) buffer
    self%Observations_Directory = trim(buffer)

    close(fhandle)
end subroutine read_project_file


function get_project_input_dp(index, key, mold) result(value)
    !! Returns the chosen double precision attribute of the
    !! ProjectInput module variable.
    integer, intent(in) :: index
    character(len=*), intent(in) :: key
    real(dp), intent(in) :: mold
    real(dp) :: value

    if (key == 'VersionNr') then
        value = ProjectInput(index)%VersionNr
    else
        call assert(.false., 'Unknown dp key: ' // key)
    end if
end function get_project_input_dp


subroutine set_project_input_dp(index, key, value)
    !! Sets the chosen double precision attribute of the
    !! ProjectInput module variable.
    integer, intent(in) :: index
    character(len=*), intent(in) :: key
    real(dp), intent(in) :: value

    if (key == 'VersionNr') then
        ProjectInput(index)%VersionNr = value
    else
        call assert(.false., 'Unknown dp key: ' // key)
    end if
end subroutine set_project_input_dp


function get_project_input_int8(index, key, mold) result(value)
    !! Returns the chosen int8 attribute of the
    !! ProjectInput module variable.
    integer, intent(in) :: index
    character(len=*), intent(in) :: key
    integer(int8), intent(in) :: mold
    integer(int8) :: value


    if (key == 'Simulation_YearSeason') then
        value = ProjectInput(index)%Simulation_YearSeason
    else
        call assert(.false., 'Unknown int8 key: ' // key)
    end if
end function get_project_input_int8


subroutine set_project_input_int8(index, key, value)
    !! Sets the chosen int8 attribute of the
    !! ProjectInput module variable.
    integer, intent(in) :: index
    character(len=*), intent(in) :: key
    integer(int8), intent(in) :: value

    if (key == 'Simulation_YearSeason') then
        ProjectInput(index)%Simulation_YearSeason = value
    else
        call assert(.false., 'Unknown int8 key: ' // key)
    end if
end subroutine set_project_input_int8


function get_project_input_int32(index, key, mold) result(value)
    !! Returns the chosen int32 attribute of the
    !! ProjectInput module variable.
    integer, intent(in) :: index
    character(len=*), intent(in) :: key
    integer(int32), intent(in) :: mold
    integer(int32) :: value

    if (key == 'Simulation_DayNr1') then
        value = ProjectInput(index)%Simulation_DayNr1
    elseif (key == 'Simulation_DayNrN') then
        value = ProjectInput(index)%Simulation_DayNrN
    elseif (key == 'Crop_Day1') then
        value = ProjectInput(index)%Crop_Day1
    elseif (key == 'Crop_DayN') then
        value = ProjectInput(index)%Crop_DayN
    else
        call assert(.false., 'Unknown int32 key: ' // key)
    end if
end function get_project_input_int32


subroutine set_project_input_int32(index, key, value)
    !! Sets the chosen int32 attribute of the
    !! ProjectInput module variable.
    integer, intent(in) :: index
    character(len=*), intent(in) :: key
    integer(int32), intent(in) :: value

    if (key == 'Simulation_DayNr1') then
        ProjectInput(index)%Simulation_DayNr1 = value
    elseif (key == 'Simulation_DayNrN') then
        ProjectInput(index)%Simulation_DayNrN = value
    elseif (key == 'Crop_Day1') then
        ProjectInput(index)%Crop_Day1 = value
    elseif (key == 'Crop_DayN') then
        ProjectInput(index)%Crop_DayN = value
    else
        call assert(.false., 'Unknown int32 key: ' // key)
    end if
end subroutine set_project_input_int32


function get_project_input_string(index, key) result(value)
    !! Returns the chosen string attribute of the
    !! ProjectInput module variable.
    integer, intent(in) :: index
    character(len=*), intent(in) :: key
    character(len=:), allocatable :: value

    if (key == 'Description') then
        value = ProjectInput(index)%Description
    elseif (key == 'Climate_Info') then
        value = ProjectInput(index)%Climate_Info
    elseif (key == 'Climate_Filename') then
        value = ProjectInput(index)%Climate_Filename
    elseif (key == 'Climate_Directory') then
        value = ProjectInput(index)%Climate_Directory
    elseif (key == 'Temperature_Info') then
        value = ProjectInput(index)%Temperature_Info
    elseif (key == 'Temperature_Filename') then
        value = ProjectInput(index)%Temperature_Filename
    elseif (key == 'Temperature_Directory') then
        value = ProjectInput(index)%Temperature_Directory
    elseif (key == 'ETo_Info') then
        value = ProjectInput(index)%ETo_Info
    elseif (key == 'ETo_Filename') then
        value = ProjectInput(index)%ETo_Filename
    elseif (key == 'ETo_Directory') then
        value = ProjectInput(index)%ETo_Directory
    elseif (key == 'Rain_Info') then
        value = ProjectInput(index)%Rain_Info
    elseif (key == 'Rain_Filename') then
        value = ProjectInput(index)%Rain_Filename
    elseif (key == 'Rain_Directory') then
        value = ProjectInput(index)%Rain_Directory
    elseif (key == 'CO2_Info') then
        value = ProjectInput(index)%CO2_Info
    elseif (key == 'CO2_Filename') then
        value = ProjectInput(index)%CO2_Filename
    elseif (key == 'CO2_Directory') then
        value = ProjectInput(index)%CO2_Directory
    elseif (key == 'Calendar_Info') then
        value = ProjectInput(index)%Calendar_Info
    elseif (key == 'Calendar_Filename') then
        value = ProjectInput(index)%Calendar_Filename
    elseif (key == 'Calendar_Directory') then
        value = ProjectInput(index)%Calendar_Directory
    elseif (key == 'Crop_Info') then
        value = ProjectInput(index)%Crop_Info
    elseif (key == 'Crop_Filename') then
        value = ProjectInput(index)%Crop_Filename
    elseif (key == 'Crop_Directory') then
        value = ProjectInput(index)%Crop_Directory
    elseif (key == 'Irrigation_Info') then
        value = ProjectInput(index)%Irrigation_Info
    elseif (key == 'Irrigation_Filename') then
        value = ProjectInput(index)%Irrigation_Filename
    elseif (key == 'Irrigation_Directory') then
        value = ProjectInput(index)%Irrigation_Directory
    elseif (key == 'Management_Info') then
        value = ProjectInput(index)%Management_Info
    elseif (key == 'Management_Filename') then
        value = ProjectInput(index)%Management_Filename
    elseif (key == 'Management_Directory') then
        value = ProjectInput(index)%Management_Directory
    elseif (key == 'GroundWater_Info') then
        value = ProjectInput(index)%GroundWater_Info
    elseif (key == 'GroundWater_Filename') then
        value = ProjectInput(index)%GroundWater_Filename
    elseif (key == 'GroundWater_Directory') then
        value = ProjectInput(index)%GroundWater_Directory
    elseif (key == 'Soil_Info') then
        value = ProjectInput(index)%Soil_Info
    elseif (key == 'Soil_Filename') then
        value = ProjectInput(index)%Soil_Filename
    elseif (key == 'Soil_Directory') then
        value = ProjectInput(index)%Soil_Directory
    elseif (key == 'SWCIni_Info') then
        value = ProjectInput(index)%SWCIni_Info
    elseif (key == 'SWCIni_Filename') then
        value = ProjectInput(index)%SWCIni_Filename
    elseif (key == 'SWCIni_Directory') then
        value = ProjectInput(index)%SWCIni_Directory
    elseif (key == 'OffSeason_Info') then
        value = ProjectInput(index)%OffSeason_Info
    elseif (key == 'OffSeason_Filename') then
        value = ProjectInput(index)%OffSeason_Filename
    elseif (key == 'OffSeason_Directory') then
        value = ProjectInput(index)%OffSeason_Directory
    elseif (key == 'Observations_Info') then
        value = ProjectInput(index)%Observations_Info
    elseif (key == 'Observations_Filename') then
        value = ProjectInput(index)%Observations_Filename
    elseif (key == 'Observations_Directory') then
        value = ProjectInput(index)%Observations_Directory
    else
        call assert(.false., 'Unknown string key: ' // key)
    end if
end function get_project_input_string


subroutine set_project_input_string(index, key, value)
    !! Sets the chosen string attribute of the
    !! ProjectInput module variable.
    integer, intent(in) :: index
    character(len=*), intent(in) :: key
    character(len=*), intent(in) :: value

    if (key == 'Description') then
        ProjectInput(index)%Description = value
    elseif (key == 'Climate_Info') then
        ProjectInput(index)%Climate_Info = value
    elseif (key == 'Climate_Filename') then
        ProjectInput(index)%Climate_Filename = value
    elseif (key == 'Climate_Directory') then
        ProjectInput(index)%Climate_Directory = value
    elseif (key == 'Temperature_Info') then
        ProjectInput(index)%Temperature_Info = value
    elseif (key == 'Temperature_Filename') then
        ProjectInput(index)%Temperature_Filename = value
    elseif (key == 'Temperature_Directory') then
        ProjectInput(index)%Temperature_Directory = value
    elseif (key == 'ETo_Info') then
        ProjectInput(index)%ETo_Info = value
    elseif (key == 'ETo_Filename') then
        ProjectInput(index)%ETo_Filename = value
    elseif (key == 'ETo_Directory') then
        ProjectInput(index)%ETo_Directory = value
    elseif (key == 'Rain_Info') then
        ProjectInput(index)%Rain_Info = value
    elseif (key == 'Rain_Filename') then
        ProjectInput(index)%Rain_Filename = value
    elseif (key == 'Rain_Directory') then
        ProjectInput(index)%Rain_Directory = value
    elseif (key == 'CO2_Info') then
        ProjectInput(index)%CO2_Info = value
    elseif (key == 'CO2_Filename') then
        ProjectInput(index)%CO2_Filename = value
    elseif (key == 'CO2_Directory') then
        ProjectInput(index)%CO2_Directory = value
    elseif (key == 'Calendar_Info') then
        ProjectInput(index)%Calendar_Info = value
    elseif (key == 'Calendar_Filename') then
        ProjectInput(index)%Calendar_Filename = value
    elseif (key == 'Calendar_Directory') then
        ProjectInput(index)%Calendar_Directory = value
    elseif (key == 'Crop_Info') then
        ProjectInput(index)%Crop_Info = value
    elseif (key == 'Crop_Filename') then
        ProjectInput(index)%Crop_Filename = value
    elseif (key == 'Crop_Directory') then
        ProjectInput(index)%Crop_Directory = value
    elseif (key == 'Irrigation_Info') then
        ProjectInput(index)%Irrigation_Info = value
    elseif (key == 'Irrigation_Filename') then
        ProjectInput(index)%Irrigation_Filename = value
    elseif (key == 'Irrigation_Directory') then
        ProjectInput(index)%Irrigation_Directory = value
    elseif (key == 'Management_Info') then
        ProjectInput(index)%Management_Info = value
    elseif (key == 'Management_Filename') then
        ProjectInput(index)%Management_Filename = value
    elseif (key == 'Management_Directory') then
        ProjectInput(index)%Management_Directory = value
    elseif (key == 'GroundWater_Info') then
        ProjectInput(index)%GroundWater_Info = value
    elseif (key == 'GroundWater_Filename') then
        ProjectInput(index)%GroundWater_Filename = value
    elseif (key == 'GroundWater_Directory') then
        ProjectInput(index)%GroundWater_Directory = value
    elseif (key == 'Soil_Info') then
        ProjectInput(index)%Soil_Info = value
    elseif (key == 'Soil_Filename') then
        ProjectInput(index)%Soil_Filename = value
    elseif (key == 'Soil_Directory') then
        ProjectInput(index)%Soil_Directory = value
    elseif (key == 'SWCIni_Info') then
        ProjectInput(index)%SWCIni_Info = value
    elseif (key == 'SWCIni_Filename') then
        ProjectInput(index)%SWCIni_Filename = value
    elseif (key == 'SWCIni_Directory') then
        ProjectInput(index)%SWCIni_Directory = value
    elseif (key == 'OffSeason_Info') then
        ProjectInput(index)%OffSeason_Info = value
    elseif (key == 'OffSeason_Filename') then
        ProjectInput(index)%OffSeason_Filename = value
    elseif (key == 'OffSeason_Directory') then
        ProjectInput(index)%OffSeason_Directory = value
    elseif (key == 'Observations_Info') then
        ProjectInput(index)%Observations_Info = value
    elseif (key == 'Observations_Filename') then
        ProjectInput(index)%Observations_Filename = value
    elseif (key == 'Observations_Directory') then
        ProjectInput(index)%Observations_Directory = value
    else
        call assert(.false., 'Unknown string key: ' // key)
    end if
end subroutine set_project_input_string

end module ac_project_input

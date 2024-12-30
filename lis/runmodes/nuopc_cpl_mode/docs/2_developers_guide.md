# Developer's Guide
This document is intended for LIS NUOPC Cap developers. The User's Guide can
be found [here](1_users_guide.md).

## NUOPC Model Phases

### Register
During the register phase subroutines are added to the ESMF registry table.
The NUOPC driver then calls each subroutine at defined points during
initialization, run, and finalize.
| Phase     | Cap Subroutines   | Description                                
|-----------|-------------------|--------------------------------------------
| Register  | ::setservices     | Register NUOPC component subroutines       

### Initialize
During the initialization phase subroutines are called to configure the
model settings, advertise fields, realize fields, initalize field data,
and initialize the model clock.
| Phase     | Cap Subroutines   | Description
|-----------|-------------------|--------------------------------------------
| Initalize | ::initializep0    | Set the IPD and configure model
| Initalize | ::initializep1    | Initialize model and advertise fields
| Initalize | ::initializep3    | Realize fields
| Initalize | ::datainitialize  | Initialize field data
| Initalize | ::setclock        | Initialize model clock

### Run
During the run phase subroutines are called to check the import field
timestamps and advance the model.
| Phase     | Cap Subroutines   | Description
|-----------|-------------------|--------------------------------------------
| Run       | ::checkimport     | Check timestamp on import data
| Run       | ::modeladvance    | Advance model by a timestep

### Finalize
Durint the finalize phase subroutines are called to destroy objects and
release memory.
| Phase     | Cap Subroutines   | Description
|-----------|-------------------|--------------------------------------------
| Finalize  | ::modelfinalize   | Release memory

## Model Fields

The following tables list the import and export fields.

### Import Fields

Import fields arelisted in the import_list parameter.

| Standard Name  | Units  | Model Variable  | Description                                | Notes
| ---------------|--------|-----------------|--------------------------------------------|--------------------------------------
| dummy_field_1  | Pa     | forcing_1       | field description for first import field   | |
| dummy_field_2  | kg     | forcing_2       | field description for second import field  | |
| dummy_field_3  | W m-2  | forcing_3       | field description for third import field   | field notes

###Export Fields

Export fields are listed in the export_list parameter.

| Standard Name  | Units   | Model Variable  | Description                               | Notes
| ---------------|---------|-----------------|-------------------------------------------|---------------------------
| dummy_field_1  | m       | output_1        | field description for first export field  | field notes
| dummy_field_2  | kg      | output_2        | field description for second export field | |
| dummy_field_3  | m s-1   | output_3        | field description for third export field  | field notes

## Memory Management

Model configuration is stored in a custom internal state data type. A
pointer to the custom internal state data type is stored in the component.

The cap allocates new memory for each field so that 2-D coordinate points
can be translated into the LIS tiled field points.


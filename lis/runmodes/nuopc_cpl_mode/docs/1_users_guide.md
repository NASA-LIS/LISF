# User's Guide
This document is intended for LIS NUOPC Cap users. The Developer's Guide can
be found [here](2_developers_guide.md).

## Building and Installing
Environment Variables
- ESMFMKFILE

NUOPC Makefile Targets
- nuopc
- nuopcinstall
- nuopcclean

The build system in <a href="Makefile">./Makefile</a> wraps the LIS build system
and adds the nuopc, nuopcinstall, and nuopcclean targets. Before building
make sure to configure the internal model.

To build and install into the current directory run:
   $ make nuopc

To install into an alternative directory run:
   $ make nuopcinstall DESTDIR=\<INSTALL\_DIR\> INSTDIR=\<SUBDIR\>

To build with debugging information run:
   $ make nuopc DEBUG=on

## Model Configuration
Model attributes are used to configure the model. The cap configuration is
stored in the interal state structure see \ref type_internalstatestruct

| Attribute            | Default         | Description
|----------------------|-----------------|-------------------------------------------------------------------------------------
| Verbosity            | 0               | String, converted into an integer. Bit 16: LIS cap information will be logged.
| Diagnostic           | 0               | String, converted into an integer. Bit 16: LIS cap diagnostics will be written.
| config\_file         | lis.config      | Set the LIS configuration file.
| nest\_to\_nest       | false           | Turn on nest to nest coupling. Each nest will be identified with an integer.
| coupled\_ensemble    | false           | Couple ensemble members using 3D fields with ungridded dimension.
| import\_dependency   | false           | Data initialization will loop until all import field dependencies are satisfied.
| output\_directory    | [CNAME]\_OUTPUT | Configure the LIS Cap output directory.

## Input and Output
Cap diagnostic output is written to the ESMF PET Logs. Cap diagnostic
output can be increased or decreased by setting the Verbosity attribute.

NUOPC state restart write files are written depending on the
RestartInterval attribute. If set to 0 then NUOPC state restart write files
will never be written.

LIS diagnostics output is written to the LIS logs configured in the LIS
configuration file.

LIS output files are written to the output directory configured in the LIS
configuration file.  LIS output includes LIS history files and LIS restart
files.


# FAO AquaCrop

<img src="https://encrypted-tbn1.gstatic.com/images?q=tbn:ANd9GcSEmcLG0zbWXIaLClg09v77ZccbFH_zuDVRbH-eBLxAdmcZ4nZ7" align="right" width="200px">

AquaCrop v7.0 and higher versions are released as open-source Fortran code, 
developed at KU Leuven and the Food and Agriculture Organization (FAO) 
of the United Nations (FAO and KU Leuven copyright),
and based on the original AquaCrop v6.0 (FAO copyright). 
Compared to AquaCrop v6.0, the AquaCrop v7.0 and higher versions
feature bug fixes, performance improvements and internal restructuring,
a translation from Pascal to Fortran, and a range of new and/or updated 
scientific features.

The following applications are publicly distributed:
* AquaCrop version-controlled open **source code** (this GitHub page, and zip-file under [Releases](
  https://github.com/KUL-RSDA/AquaCrop/releases))
* AquaCrop standard Windows **graphical user interface** (zip-file under [Releases](
  https://github.com/KUL-RSDA/AquaCrop/releases))
* AquaCrop **standalone executable** (zip-file under [Releases](
  https://github.com/KUL-RSDA/AquaCrop/releases)) for
  * Windows
  * Linux
  * MacOS

From v7.0 onwards, it is also possible to use AquaCrop
as a crop model within [NASA’s Land Information System (LIS)](
https://github.com/NASA-LIS/LISF). More information can be found in
the LIS section below.

## Documentation

Online documentation and contact information are available at the [FAO website](https://www.fao.org/aquacrop/en/). The AquaCrop core team is small and answers will be found fastest in the release notes, training handbooks and youtube videos provided by FAO.

Please also visit our [Discussions](https://github.com/KUL-RSDA/AquaCrop/discussions) forum for FAQ, or to contribute.

## Running the executable

Download the ZIP file with the AquaCrop executable (v7.0 or higher) for
Windows, Linux or MacOS from the release page. 
Consult the reference manual (FAO website) for details about the AquaCrop stand-alone program.

Optionally, it can be verified if the executable produces the expected output on the user's system, by running a self-contained testcase for which reference output is provided (zip-file under [Releases](
  https://github.com/KUL-RSDA/AquaCrop/releases)). Please carefully read the README file to run the testcase.

## Building the executable

Either unzip the ZIP file with the source code from the release page, 
or if you wish to contribute to the code, then fork the repository and locally clone your fork.

The source code is under the `src` directory. Building the Aquacrop executable on a Linux system requires:

* GNU Make (>= v3.82)
* a GNU or Intel Fortran compiler (GNU Fortran >= v6.4.0 and ifort >= v18.0.1).
  MinGW can be used to (cross)compile for Windows. In that case `make` needs
  to be called with an additional `CPPFLAGS=-D_WINDOWS` option.

```bash
cd AquaCrop/src
make
```

The main `make` targets are `bin` (producing an `aquacrop` executable),
`lib` (producing a `libaquacrop.so` library). The default target is
`all`, which combines the `bin` and `lib` targets.

## Optional contributing and support

Please follow good practices. New features, enhancements or suggestions will only be considered and reviewed once a year by the core AquaCrop developers.

We encourage scientific (only) exchanges via our [Discussions](https://github.com/KUL-RSDA/AquaCrop/discussions) forum. Only if the wealth of [documentation](https://www.fao.org/aquacrop/en/) or the FAO Contact (aquacrop@fao.org) did not provide sufficient help, or if you have a good suggestion, then start a new "Discussion" with the information you already gathered from the FAO Contact or documentation. Please do not open an "Issue" to ask your question and do not offer a "Pull Request" without any prior "Discussion" with the AquaCrop core team.


## LIS integration

The distribution of AquaCrop v7.0 and higher versions within NASA's LIS is currently being reviewed
and will be accessible after approval of an upcoming pull request to NASA via the [LIS GitHub page](https://github.com/NASA-LIS/LISF).

## Citation

A wide range of publications is available to refer to AquaCrop in the GUI or standalone version.
The users can refer to any publication of their choice, when using these AquaCrop assets. 

In case the open source code of AquaCrop (from v7.0 onwards) is used, please refer to one of the following papers, covering precursory work in preparation for the open source release. Publications using AquaCrop v7.0 and higher versions in NASA's LIS will be added for reference as soon as they become available.
* de Roos, S., De Lannoy, G.J.M., Raes, D. (2021). Performance analysis of regional AquaCrop (v6.1) biomass and surface soil moisture simulations using satellite and in situ observations. Geoscientific Model Development, 14(12), 7309-7328, 10.5194/gmd-14-7309-2021.
* Busschaert, L., de Roos, S., Thiery, W., Raes, D., De Lannoy, G.J.M. (2022). Net irrigation requirement under different climate scenarios using AquaCrop over Europe. Hydrology and Earth System Sciences, 26, 3731–3752, 10.5194/hess-26-3731-2022.


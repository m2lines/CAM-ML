# CAM-ML: A Fork of The Community Atmosphere Model Implementing a Machine Learning Convection Parameterisation

This Fork of the CAM model implements a [new convection parameterisation, YOG](https://github.com/m2lines/convection-parameterization-in-CAM).
The parameterisation is a machine learning implementation using a neural net trained
on high-resolution cloud-resolving simulations in the SAM model as described in:\

- Yuval, J., O’Gorman, P.A.
  _Stable machine-learning parameterization of subgrid processes for climate modeling at
  a range of resolutions_
  Nat Commun 11, 3295 (2020). DOI: [10.1038/s41467-020-17142-3](https://doi.org/10.1038/s41467-020-17142-3)
- Yuval, J., O'Gorman, P.A., Hill, C.N.
  _Use of Neural Networks for Stable, Accurate and Physically Consistent Parameterization of Subgrid Atmospheric Processes With Good Performance at Reduced Precision_
  Geophysical Research Letters, 48, e2020GL091363 (2021). DOI: [10.1029/2020GL091363](https://doi.org/10.1029/2020GL091363)

The work is contained in a `CAM-ML` branch which is based off the `cam_cesm2_1_rel_60` tag.

Note that the `doc` folder and `README_EXTERNALS` are not updated, but just artifacts from the fork that we are keeping to be able to merge back to CESM later.
Information specific to running this fork of CAM within CESM is included in this README.

## Using this model in a CESM Run

This section describes how to use CAM-ML as the atmospheric component in a normal CESM run. In its current build pipeline, the code assumes that it is run on NCAR's supercomputer Derecho.

### Obtaining CESM

Clone a copy of CESM from git and checkout the `cesm2.1.5` tag on which this work is based:
```
git clone https://github.com/escomp/cesm.git my_cesm_sandbox_2_1
cd my_cesm_sandbox_2_1/
git checkout cesm2.1.5
```

### Setting CAM version in CESM

To use this model in a CESM run you need to modify the `Externals.cfg` file in the
main CESM directory to replace the CAM entry with:
```toml
[cam]
branch = CAM-ML
protocol = git
repo_url = https://github.com/m2lines/CAM-ML
local_path = components/cam
externals = Externals_CAM.cfg
required = True
```
This will pull the `CAM-ML` branch of this repo in as the CAM component.

You can now run, from within the CESM root directory,
```
./manage_externals/checkout_externals
```
to fetch the external components.

**Note**   
If you want to change the externals, or have made a mistake in this step, you have to delete the newly created `components` folder (or the respective subfolders therein) in the base directory before you rerun `checkout_externals`.)

### Creating and running a case

Details on creating a case can be found
[here](https://ncar.github.io/CAM/doc/build/html/CAM6.0_users_guide/building-and-running-cam.html) on the NCAR website.
For this work we are using the gate III testcase which can be set up by running:
```
./create_newcase --case <path_to_testcase_directory> --compset FSCAM --res T42_T42 --user-mods-dir ../../components/cam/cime_config/usermods_dirs/scam_gateIII --project NCGD0054
```
from `<cesm_root>/cime/scripts/`.

The `<testcase_directory>` should be a separate directory outside of the code directory, to avoid cluttering up the local repository.

Once this has been done then edit `user_nl_cam` for the case as detailed below.
This is a CAM namelist generated from the default for the case.
Add the following lines:

1. `deep_scheme = 'off'`\
   `yog_scheme = 'on'` \
    If running a comparison to the ZM scheme also add `run_deep_comp = 'on'`.
2. `yog_nn_weights = '<PATH/TO/WEIGHTS.nc>'`\
    The path to the NN weights. There are some weights in `src/physics/cam/' of this respository (CAM-ML) which can be used, or they can be generated from the [standalone model](https://github.com/m2lines/convection-parameterization-in-CAM).
3. `SAM_sounding = '<PATH/TO/SAM/SOUNDING.nc>'`\
    The path to the SAM sounding for the NN.\
    This file is generated using the `sounding_to_netcdf.py` script in the resources of the [standalone NN code](https://github.com/m2lines/convection-parameterization-in-CAM), and it can be just copied over to a suitable place.

All the paths have to be absolute, as they will be used when the code is run on the compute nodes.

We can then run `./case.setup` and `./case.build`to setup and build the test case, and `./case.submit` to submit the job to the scheduler.

**Note:**  
By default, CESM will place outputs and logs/restart files on Derecho in `/glade/derecho/scratch/<user>/archive/<case>/`.
Outputs of failed runs can be found in`/glade/derecho/scratch/<user>/<case>/`.

To place all output with logs in `archive/case` switch 'short term archiving' on by editing `env_run.xml` in the case directory to change `DOUT_S` from `FALSE` to `TRUE` -- you can do this quickly by running `﻿./xmlchange DOUT_S=FALSE`.

### Expected Output

The run will create a `bld` and a `run` directory in the output case folder on `/glade/derecho/scratch/<user>/. The `bld` directory contains build files, logs and executables. The `run` directory contains the run logs, namelists and output NetCDF files.

A successful model run will have the line \
`******* END OF MODEL RUN *******` \
at the bottom of the run log file for the atmosphere run, i.e.,  `atm.log.<stuff>.gz`.

## CAM Documentation

CAM Documentation - https://ncar.github.io/CAM/doc/build/html/index.html

CAM6 namelist settings - http://www.cesm.ucar.edu/models/cesm2/settings/current/cam_nml.html

Please see the [wiki](https://github.com/ESCOMP/CAM/wiki) for information.

## Contributing

Contributions to the repository are welcome from members of [M2Lines](https://m2lines.github.io/) and [ICCS](https://iccs.cam.ac.uk/).

Open tickets can be viewed at [**Issues**](https://github.com/m2lines/CAM-ML/issues).

To contribute, find a relevant issue or open a new one and assign yourself to work on it. Then create a branch in which to add your contribution and open a pull request. Once ready assign a reviewer and request a code review. Merging should only be performed once a reviewer has approved the changes.

Interested contributors from outside M2Lines are invited to comment on issues to propose solutions.

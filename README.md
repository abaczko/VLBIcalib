# VLBIcalib
This is a set of functions to ease the VLBI calibration by using ParselTongue.

By Anne-Kathrin Baczko (baczko@mpifr-bonn.mpg.de).

If you like to use these scripts for any publication, please let me know.

It has only been tested on a small amount of datasets, so you may have to change them according to your needs.

For the imported modules to be found the repository has to been added to the local PYTHONPATH variable.

## Requirenments
- astropy, numpy
- ParselTongue: for information about installation please see https://github.com/kernsuite-debian/parseltongue 
- AIPS: in case you do not want to use AIPSLite (which comes with ParselTongue) but a local installation. Please see http://www.aips.nrao.edu/index.shtml for installation instructions

## Submodules and functions
Below is a short description of the individual files:

### observation_parameters.py
Ideally this should be the only file that has to be changed for each observation. It contains paths to calibration files, other data, settings for the different calibration steps (as manual phase cal parameters timer,refant...) and variables that should be set either to 'True' or 'False'. These define what steps should be run. From loading the correlated fits files, writing Flag tables, doing manual phasecal, bandpass or applying SPLIT and FITTP.

### PT_calibration_pipe.py: 
This is the main script to be run in a terminal as
   ```
   $ ParselTongue PT_calibration_pipe.py
   ```
Depending on the variables given in observation_parameters.py it runs different calibration steps. There is the option to use AIPSLite, in case no local aips installation is available or you do not want to use it for any reason. In that case folders will be created to hold the AIPS files and tasks in the same directory the script had been invoken.
### logging_config.json
A simple Java script that must reside in the directory where ParselTongue is called. It has some settings for the logger that is implemented in all the scripts. Everything that will be done durign calibration with theses scripts is logged to a file PT_inputs_xxx.log. For each day a new logging script will be made.

### aips_tasks.py
A collection of functions to make it easier to start aips tasks, with less input parameters to give to the functions. Also some things are done automatically, as downloading the usno.finals file for EOP. In addition there are some functions, that make the life easier as returning the header of a uvdatafile to a python dictionary, writing the header to a text file, get all keyword parameters for a given disk and catalog number or having the same functionality as used in aips (e.g. extd).

### functions.py
Additional functions that combine some python functionality with AIPSTasks or define a set of things that are normally done in a row. E.g. do_global_fring translates antennas names given in strings to numbers as input to antennas or refant, do a global fring, apply snsmo to the sntable if desired and applies CLCAL. Another example is get_scans() which returns an astropy table of all scans of the experiment as written in LISTR(opty='scan').

### helper_functions.py
A collection of functions needed in aips_tasks, e.g. downloading a file or catching output from the command line.

### extract_info_sn.py
Returns information saved in an SN table produced by fring. This function is needed by derive_solint() in fucntions.py

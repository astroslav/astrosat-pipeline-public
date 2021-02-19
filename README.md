# README #

### What is this repository for? ###

* This repository is for post-processing of level 1 AstroSAT data from the LAXPC instrument to level 2 and analyzing the data.

### LAXPC data analysis using LAXPC software available [here](http://astrosat-ssc.iucaa.in/?q=laxpcData)

The steps to obtain the following four "start-up" files:

- eventfiles
- filterfiles
- level2.event.fits
- usergti.fits

can be found on that website, at the bottom as a link called "**LAXPC_walkthrough**". These files will be necessary for the next steps. I installed Format A. 



## "Quick" Guide to post-processing using this repository

1. After you pull the repository, make a directory called _start_files_ and copy the four files mentioned above into it.
2. Create a different working directory within which you will do your analysis.
3. Initialize HEASoft in your terminal.
4. Run `laxpc_make_backlc.py` in your working directory. This script calls on many others to generate spectra, back spectra, and back lightcurve. The spectra are necessary for the LAXPC software to make the back lightcurve.
5. Run `plot_simple_lc.py` to generate a simple lightcurve of your data in the working directory.
... more steps to continue...



## Detailed guide

### Making spectra, back spectra, and back lightcurves

`laxpc_make_backlc.py` is the python script . It calls upon several other scripts in the directory to make spectra, background spectra and background lightcurves based on your inputs. I will describe this script step-by-step below, including the functions called in the order that `laxpc_make_backlc.py` calls them

1. `symlink_startfiles.py` - This script symbolically copies the four start-up files to your working directory so they can be used in the analysis. **You should change this file to point to the path where your four start-up files are.**

2. The script asks for the usergti file and the level 2 file. The defaults assumed are `usergti.fits` and `level2.event.fits`, respectively.

3. `get_PCUs_dt_layer.py` - This script contains a function called **get_PCUs_dt_layer** which asks for the PCUs used, time resolution you want, and the layers of the instrument. The defaults used for PCUs are the LAXPC 10 and 20 (LAXPC 30 is known to have hardware issues), the default time resolution is 100 ms, and the default layers is all of them (with the other option being '1', where only the first layer is used).

4. `split_gti.py` - This script contains a function called **split_gti** which takes a time resolution and a usergti.fits file as inputs and writes usergti_XXXX.fits files, where the XXXX is the number of the file, as outputs based on the time resolution requested. The function **laxpc_make_backlightcurve** defined by the LAXPC software has an arbitrary maximum limit of one million points, so this circumvents that by creating many fits files that all have less points than that. 

5. `context_manager.py` - This script contains a class called **cd**, which changes the current working directory. Used in the case where multiple usergti_XXXX.fits files are created and each script has to be ran in each of those working directories.

6. `call_laxpc_software.py` - This script contains three functions, each calling the LAXPC software functions.

    * **call_spectra** - Function for calling **laxpc_make_spectra**, where the inputs are a usergti.fits file, layer choice, and a level2event.fits file
	
	* **call_backspectra** - Function for calling **laxpc_make_backspectra**, where the inputs are a usergti.fits file and layer choice.
	
	* **call_backlc** - Function for calling **laxpc_make_backlightcurve**, where the inputs are PCUs, dt, usergti.fits file, layer choice, and "backlc_choice". The "backlc_choice" will ask you how you want to generate the back lightcurve file(s). There are two possible ways right now: (i) one file with an energy range of 3-80 keV or (ii) 77 files, each with an energy range of 1 keV (e.g. 3-4 keV, 4-5 keV, and so on). 

7. `get_ebounds_from_rmf.py` - This script creates the ebounds file necessary for **laxpc_make_backlightcurve** from the rmf files. This script is particularly involved, so I will describe it in steps.

    * It (the script) searches for any files named "lx10cshp*" in the directory (rmf files), extracts the "EBOUNDS" data from those fits files, and writes them to numpy arrays for each PCU.
	
	* It asks you whether to include the extremes of the energy (the bins including the edges of the energy range, for example 3.0 keV and 80.0 keV for the 3.0-80.0 keV range).
	
	* It creates three arrays (for the three PCUs) that contain the lower and upper energy bounds and combines them into one full vertically stacked array that contains all the energy bounds.
	
	* The script contains a function called **check_arr_counts** which checks to make sure that all counts in all energy bounds defined are counted only once. If any ebounds counts a single event more than once or misses a single event, then this checkpoint would stop the rest of the script.

8. `recombine_backlc.py` - This script contains a function called **recombine_backlc** which creates a combined background lightcurve text in the case that multiple usergti_XXXX.fits were generated. The input is, by default, the background lightcurve filename "Back_lightcurve.txt", which is the default name created by **laxpc_make_backlightcurve**.



Others:

* `find_nearest_idx.py` - Used for finding the index in the array closest to a given value. Used in `get_ebounds_from_rmf.py`

* `create_gti_from_ascii.sh` - Bash script for reading in the original user gti fits data file and creating a new ascii "lcusergti.fits" file. HEASoft is initialized within this script, so either change the directory in the script so that it automatically initializes it, or remove the lines and manually initialize HEASoft in the terminal before running the script.

* `backlc1keV_save_combined.py` - Python script used for reading a series of files named "Back_lightcurve_%.1f_%.1fkeV.lc", which are generated when you choose to generate background lightcurves at 1keV intervals during the "call_laxpc_software" step of "laxpc_make_backlc.py". The inputs are the desired energy ranges for the combined background lightcurve. For example, if you have 77 background lightcurves at 1keV intervals (e.g. 3-4 keV, 4-5 keV, and so on), you can run this script to generate a combined background lightcurve for a user-defined range. If your inputs are, for example, "3.0" for the lower energy bound and "10.0" for the upper energy bound, then you would generate a combined background lightcurve for the energy range of 3.0-10.0 keV.



### Making lightcurves

#### Simple lightcurve

`plot_simple_lc.py` generates a simple lightcurve from the usergti fits and level2data file. You can choose the time resolution of the lightcurve, with the default set to 100 ms. The default mask for creating the lightcurve is for the energy range of 3.0-80.0 keV, using LAXPC 10 and 20, and all five layers of each PCU. The script also calls on the following scripts:

* `get_gtis.py` - A script that reads a gti fits file and returns two matrices: (i) "gtis", which is the full matrix of START and STOP times for all good time intervals, and (ii) "gtis_min", which is only the START of the first GTI and the STOP of the last GTI, i.e. the full time spanned by all GTI's.

* `get_level2data.py` - A script that reads an event fits file and returns three things: (i) a numpy matrix of the level2data, (ii) the number of counts, and (iii) the Modified Julian Date reference. The level2data matrix contains five columns - time, channel, layer, laxpc_no, and energy, respectively.

* `EnLxLa_mask.py` - A script that defines a mask on the data based on the energy, laxpc PCU number, and laxpc PCU layers. The default mask for the energy, laxpc PCU number, and laxpc PCU layers are 3.0-80.0 keV, laxpc 10 and 20, and all five layers.

##### The simple lightcurve script also uses the Stingray Python package. More information about it can be found [here](https://stingraysoftware.github.io/stingray/index.html).

#### Simple background lightcurve

`plot_simple_bkg_lc.py` generates a simple background lightcurve from the "Back_lightcurve.txt" background lightcurve file generated from earlier. It assumes that the time and the count_rate columns are the fourth and second, respectively. Very simple.



### Working with Stingray

`stingray_import.py` - A script used for importing the data from the usergti fits and level2data file, as well as making and returning a lightcurve file.

`stingray_binned_lc.py` - A script which calls the lightcurve data using stingray and bins the lightcurve to a given number of counts per bin. It plots the full lightcurve and the binned lightcurve for comparison.

**Any file with the prefix of `test_` is an incomplete file or currently in-progress.**

`test_stingray_powerspectrum.py` - A script which tests the functionality of the Powerspectrum and AveragedPowerspectrum classes from Stingray.

`test_stingray_bkg_subtract_and_crosscorrelation.py` - A script which explores a background subtracted lightcurve and the CrossCorrelation class from Stingray. It contains several different useful functions. Here is a list of the functions:

* **create_gti_from_time_arr**: You can use this to create your own GTI. You need to define the time array you want as the input and the time resolution of the GTI.

* **self_interpolate**: This function interpolates in between a given array. So if you have a time array that has a time resolution of 0.1s, you can make it smaller, e.g. 0.01s. 

* **add_poisson_to_lc**: This function adds poisson noise to the Stingray Lightcurve object passed as an input.

* **interp_time_series**: You can use this function to interpolate the time and countrates of a given input Lightcurve Stingray object. It returns an interpolated Lightcurve object. 

* **subtract_background**: Inputs are two Lightcurve objects: the lightcurve and its background lightcurve, and the output is a single Lightcurve object: a background-subtracted lightcurve.

* **In_Po_Sub**: This function combines the interpolation steps, adding poisson noise, and subtracting the background into one function.


### Working with PyXspec

##### The documentation on PyXspec's functions can be found [here](https://heasarc.gsfc.nasa.gov/docs/xanadu/xspec/python/html/index.html)

`test_pyxspec.py` - This is the main script which contains everything I've done with PyXspec. A lot of the inner workings are explained in comments blocks within the file.
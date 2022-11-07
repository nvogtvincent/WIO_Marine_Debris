# WIO_Marine_Debris
Repository for scripts used to reproduce figures in "Sources of marine debris for Seychelles and other remote islands in the western Indian Ocean"

## Directory overview
```
ANALYSIS/         Scripts to reproduce figures
FIGURES/          Figures are saved here
FISHERIES/        Fisheries input scripts and postprocessed data 
GRID_DATA/        Contains analysis grid (with source/sink sites)
MATRICES/         Postprocessed terrestrial and marine matrices
PLASTIC_DATA/     List of input and sink countries used in analysis 
REFERENCE/        Reference datasets for plotting/analyses
SHIPPING/         Reference shipping density file used for Figure 5
SIM/              Scripts to rerun particle-tracking analyses (relies on ocean/wind/wave data as described in manuscript)
TRAJ/             Postprocessing scripts to generate matrices from model output  
```

## User guide

Firstly, set up a conda environment with the `environment.yml` file in the root directory. Alternatively, the packages needed for these analyses are numpy, matplotlib, xarray, xskillscore, pandas, cmasher, scipy, scikit-image, gdal, cartopy, resterio, geopandas, and tqdm.  

### I want to generate results using a new combination of beaching and sinking rate, and windage
1. Retrieve the trajectory files for the physical scenario of interest from the BODC repository, and put them in in the relevant subdirectory in `TRAJ`. For instance, for the terrestrial CS3 scenario, you should put all 1056 terrestrial CS3 netcdf files in `TRAJ/TER/CS3/`; for the marine C0 scenario, you should put all 960 marine C0 netcdf files in `TRAJ/MAR/C0`, etc.
2. Run the relevant script to generate a processed matrix (note that this may take several hours):
    1. **Terrestrial**: Run `land_process.py [us_d] [ub_d] [scenario]` where `[us_d]` is the sinking timescale in days, `[us_b]` is the beaching timescale in days, and `[scenario]` is the scenario code (e.g. C0, CS0, etc.)
    2. **Marine climatology**: Run `marine_process_monclim.py [us_d] [ub_d] [scenario]` where `[us_d]` is the sinking timescale in days, `[us_b]` is the beaching timescale in days, and `[scenario]` is the scenario code (e.g. C0, CS0, etc.). By default, this script generates results for debris beaching at the Aldabra Group. If you want to generate results for debris beaching at any other site combination in Seychelles, modify line 43 to the relevant list of site numbers (see the `sink_id_psi` variable in `GRID_DATA/griddata_land.nc`). For instance, to generate results for Cosmoledo and Astove, this line should read `'sites': np.array([3,4])`
    3. **Marine source/sink time**: Run `marine_process.py [us_d] [ub_d] [scenario] [source/sink]` where `[us_d]` is the sinking timescale in days, `[us_b]` is the beaching timescale in days, `[scenario]` is the scenario code (e.g. C0, CS0, etc.), and `[source/sink]` specifies whether you want the matrix to give fluxes as a function of when the debris entered the ocean (source), or when the debris beached (sink). As above, this script generates results for the Aldabra Group by default. If you want to analyse another island/s, change line 46 (as described above).

### I want to reproduce a figure, or produce a figure using a new set of parameters
1. Put the prerequisite processed matrix (either as supplied, or using a matrix you have generated from one of the processing scripts) in `MATRICES/`
2. Run the relevant script in `ANALYSIS/`

**Script-specific notes:**
* `fig_4.py`: Run `fig4.py [us_d] [us_b] [scenario] [name]` where `[us_d]` is the sinking timescale in days, `[us_b]` is the beaching timescale in days, `[scenario]` is the scenario code (e.g. C0, CS0, etc.), and `[name]` is an optional debris class name used for plotting. For example, to reproduce Figure 4(c) in the main text, run `fig_4.py 360 30 CS3 'Class C'`. Please note that, if you want to reset the colour scheme (which you may wish to do if you use very different parameters to those used in our study), you can do this by setting the value of the `write_cmap` key on line 30 to `True`.
* `fig_5.py`: Run `fig5.py [us_d] [us_b] [scenario] [name]` where `[us_d]` is the sinking timescale in days, `[us_b]` is the beaching timescale in days, `[scenario]` is the scenario code (e.g. C0, CS0, etc.), and `[name]` is an optional debris class name used for plotting. For example, to reproduce Figure 5 in the main text, run `fig_5.py 360 30 CS3 'Class C'`. This script plots results for the Aldabra Group by default. To change this, edit line 45.
* `fig_6.py`: Run `fig6.py [us_d] [us_b] [scenario] [name]` where `[us_d]` is the sinking timescale in days, `[us_b]` is the beaching timescale in days, `[scenario]` is the scenario code (e.g. C0, CS0, etc.), and `[name]` is an optional debris class name used for plotting. For example, to reproduce Figure 6 in the main text, run `fig_6.py 360 30 CS3 'Class C'`. This script plots results for the Aldabra Group versus Seychelles Plateau by default. To change this, change the list of names in the lists on lines 22-24.
* `fig_7.py`: Run `fig7.py [us_d] [us_b] [scenario] [source/sink] [name]` where `[us_d]` is the sinking timescale in days, `[us_b]` is the beaching timescale in days, `[scenario]` is the scenario code (e.g. C0, CS0, etc.), [source/sink] specifies whether you want the matrix to give fluxes as a function of when the debris entered the ocean (source), or when the debris beached (sink), and `[name]` is an optional debris class name used for plotting. For example, to reproduce Figure 7 in the main text, run `fig_7.py 360 30 CS3 sink 'Class C'`. This script plots results for the Aldabra Group by default. To change this, edit line 38.
* `fig_8.py`: Please note that this script relies on proprietary dFAD data, which we do not have the permissions to make publicly available. Please contact the author (Noam Vogt-Vincent) if you would like to discuss this script. 
* `fig_9a.py`: Run `fig_9a.py [us_d] [us_b] [scenario] [source/sink] [mode] [delay] [name]` where `[us_d]` is the sinking timescale in days, `[us_b]` is the beaching timescale in days, `[scenario]` is the scenario code (e.g. C0, CS0, etc.), `[source/sink]` specifies whether you want the matrix to give fluxes as a function of when the debris entered the ocean (source), or when the debris beached (sink), `[mode]` is the ocean mode to carry out correlations with, `[delay]` is a delay to include for plotting (in months), and `[name]` is an optional debris class name used for plotting. For example, to reproduce Figure 9(a) in the main text, run `fig_9a.py 360 30 CS3 sink DMI 0 'Class C'`. This script plots results for the Aldabra Group by default. To change this, edit line 38.
* `fig_9b.py`: Just run the script. You can change the ocean mode to regress against in line 40.

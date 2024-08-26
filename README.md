# Global-extreme-events-tool
A flexible toolbox to study such concurrent temperature extremes, with adjustable parameters that different users can tailor to their specific needs

Usage of the tool.

Setting up the folders

The analysis tool expects to find the data organized in folders. This is defined in the file “file_details.py”. The root path where the data will be stored is specified in the global variable “root_path” that must be changed accordingly. Inside the root path, a folder “Heatwaves” must be present with the daily statistics together with an additional folder (in the example generically named “ERA5_mask”) with the spatial mask to distinguish between sea and land downloaded from ERA5 for instance. All the paths to be used by the tool are specified in the “file_details.py” script.

In particular:
1) filePath_absolute: set the global path where the daily statistics should be stored and found.
2) filePath_trend_climatology: set where the climatology and the year trends are stored
3) filePath_mask: set the path where the spatial mask is stored
4) filePath_mask_ERA5: function that provides paths of the mask, the resolutions (original and desired) as well as the folder where the mask is located
5) start_calendar_day: this function returns the number of seconds between the dataset time origin and a reference day (say the central date) of the used time period and can be computed as
(datetime.datetime(year,month,day)-datetime.datetime(year0,month0,day0)).days


Download the dataset

The dataset to be used in the analysis can be downloaded from a generic repository and stored as daily statistics with any spatial resolution. The file “download_ERA5_eliminate.py” downloads ERA5 single level data with hourly resolution for each month, computes the daily statistics and saves then as a pickle file for later analyses. Next, the script removes the downloaded netCDF file to save disk space (although this last step is not essential). Having available a dataset that was already downloaded, it is possible to modify the script to just create the pickle file and store the daily statistics. 

Managing the functions

The main functions are located in the "methods.py" and "clustering.py" scripts, related to the data organization (climatology and detrending) in the former and the clustering routines in the latter. An example about the usage of some key functions is reported in "main.py" while "plot_stuff.py" makes a plot of the cold spells between 2022 and 2023 (the pickle file obtained for the example is also among the files to give an idea). 

import numpy as np
from methods import save_detrended_and_anomalies , save_mask
from clustering import t2m_extreme
import pickle

# actual resolution of the data (0.5 deg x 0.5 deg)
latitudeb  = -np.linspace(-90,90,361)                       # original latitude [-90,90]

#############################################################################

# save_detrended_and_anomalies(2022,2023)

#############################################################################

# save_mask()

#############################################################################

# # basic usage to retrieve the dataset time series
# from methods import ssn_data_load
# import matplotlib.pyplot as plt
# month_list=[1,2]
# yr_list=[2023]
# latind=np.where((latitudeb>=30) & (latitudeb<=70))[0]
# data,time,lat,lon=ssn_data_load( month_list, yr_list, latind,  data_type='absolute', data_det=False, mask=True)
# plt.plot(data[:,20,200])
# plt.xlabel('Time [days]')
# plt.ylabel('Temperature [K]')
# plt.title('Example of time series at a given location: lat '+str(lat[20])+' lon '+str(lon[200]))
# plt.figure()
# plt.contourf(lon,lat,data[30,:,:],30)
# plt.xlabel('Longitude')
# plt.ylabel('Latitude')
# plt.colorbar()
# plt.show()

#############################################################################
# 
# # this is the main script to run the clustering

yr_list=np.arange(2022,2024)    # list of the years to analyse
data_type='anomaly'             # type of data to analyse
data_det=False                  # if True, the data are detrended
mask=True                       # if True, the land-sea mask is applied
duration=4                      # minimum duration of the event
extent_flag=True                # if True, the extent of the event is used to filter the events
cluster_ex=1000                 
dist_type='centroid'            # type of distance to use in the clustering
extent=2e5                      

# example
month_list=[1,2]; percentile=5; Name='ColdSpell_NH_5'       # example for cold spells in the NH at the 5th percentile
latind=np.where((latitudeb>=30) & (latitudeb<=70))[0]       # latitude indices of the NH band of interest

# extreme events identification
time, time_events, Areas, average_variable_all, LATLON = t2m_extreme(month_list, yr_list, latind, percentile, data_type, data_det, mask, duration, extent_flag, cluster_ex, dist_type, extent)

# saving the results
with open(Name+'.pkl','wb') as f:
    pickle.dump({'LATLON':LATLON,'time_original':time,'time_events':time_events,'Area': Areas,'Tmean':average_variable_all},f)
 

import numpy as np
import os

root_path='/home/antonio/Documents/Messori/final_code'

####################################################################

def filePath_absolute(month,year):
    # file and path where the daily ERA5 data are stored in .pkl format
    datapath = os.path.join(root_path,'Heatwaves')
    return os.path.join(datapath,'ERA5_t2m_'+str(month).zfill(2)+'_'+str(year)+'.pkl')

####################################################################

def filePath_trend_climatology():
    # file and path where the trend and climatology are stored
    datapath = os.path.join(root_path,'Heatwaves')
    return os.path.join(datapath,'ERA5_t2m_trend_climatology.pkl')

####################################################################

def filePath_mask():
    # file and path where the land-sea mask is stored
    datapath=os.path.join(root_path,'ERA5_mask','surface')
    return os.path.join(datapath,'lsm_0.5x0.5.pkl')

####################################################################

def filePath_mask_ERA5():
    # file and path where the land-sea mask is stored
    datapath=os.path.join(root_path,'ERA5_mask','surface','lsm_1279l4_0.1x0.1.grb_v4_unpack.nc') # original mask file
    original_resolution=0.1                                                                 # resolution of the original mask
    desired_resolution=0.5                                                                  # resolution of the desired mask
    root_path_loc=os.path.join(root_path,'ERA5_mask','surface')                                  # path where the mask is stored
    return datapath , original_resolution , desired_resolution , root_path_loc

####################################################################

def start_calendar_day():
    # number of seconds at 1/1/1980 (central day of the analysis) from 1/1/1970 (origin of the ERA5 dataset)
    # obtained from (dt.datetime(1980,1,1)-dt.datetime(1970,1,1)).total_seconds()
    return 315532800    
    
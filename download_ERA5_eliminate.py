import cdsapi
from netCDF4 import Dataset
import numpy as np
import pickle
from file_details import *
import os
c = cdsapi.Client()

for year in range(2023,2024):
    for month in range(1,13):
        fileName = 't2m_ERA5_1hr_05x05_'+str(year)+str(month).zfill(2)+'.nc'
        c.retrieve(    'reanalysis-era5-single-levels',
        {
        'product_type': 'reanalysis',
        'format': 'netcdf',
        'variable': '2m_temperature',
        'year': str(year),
        'month': str(month).zfill(2),
        'day': [
            '01', '02', '03',
            '04', '05', '06',
            '07', '08', '09',
            '10', '11', '12',
            '13', '14', '15',
            '16', '17', '18',
            '19', '20', '21',
            '22', '23', '24',
            '25', '26', '27',
            '28', '29', '30',
            '31',
        ],
        'time': [
            '00:00', '01:00', '02:00',
            '03:00', '04:00', '05:00',
            '06:00', '07:00', '08:00',
            '09:00', '10:00', '11:00',
            '12:00', '13:00', '14:00',
            '15:00', '16:00', '17:00',
            '18:00', '19:00', '20:00',
            '21:00', '22:00', '23:00',
        ],
        },  fileName )

        # additional part that decimates the data, makes daily averages and saves them in a .pkl file
        fh = Dataset(fileName, mode='r')
        decimation_lat,decimation_lon=2,2   # 0.5,0.5 degrees
        latval = fh.variables['latitude'][::decimation_lat]
        lonval = fh.variables['longitude'][::decimation_lon]
        time_orig = fh.variables['time'][:]
        Ndays=len(time_orig)//24
        t_list,var_data=[],[]
        for i in range(Ndays):
            t_list.append(time_orig[i*24])
            var_data.append( np.mean(fh.variables['t2m'][i*24:(i+1)*24,::decimation_lat,::decimation_lon],axis=0) ) # Daily variable (lat,lon)

        fh.close()
        var_data=np.stack(var_data,axis=0) # (time,lat,lon)
        with open(filePath_absolute(month,year),'wb') as f:
            pickle.dump({'time':np.array(t_list),'longitude':lonval,'latitude':latval,'data':var_data},f, pickle.HIGHEST_PROTOCOL)
        os.remove(fileName)
        
        
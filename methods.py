import numpy as np
import pickle , os
from file_details import *

#############################################################################

def haversine(P1, P2):
    # Haversine formula (P1 is a single point, P2 is an array of points [lat,lon])
    lat1 = np.deg2rad(P1[0])
    lat2 = np.deg2rad(P2[:,0])
    dlat = np.deg2rad(P2[:,0] - P1[0])
    dlon = np.deg2rad(P2[:,1] - P1[1])
    a = np.sin(dlat/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2)**2
    return 2 * 6371 * np.arcsin(np.sqrt(a))

#############################################################################

def areaquad(P1,P2): # Checked 20240404
    # spherical area between two points (P1 and P2 are corner points)
    Az  = np.sin(np.deg2rad(P2[0])) - np.sin(np.deg2rad(P1[0]))
    Az *= np.deg2rad(P2[1]-P1[1])
    Az *= 6371**2 
    return np.abs(Az)

#############################################################################

def Area_rectangle_unit(latitude,longitude): # Checked 20240404
    Area=np.array([areaquad([latitude[i],longitude[0]],[latitude[i+1],longitude[1]]) for i in range(len(latitude)-1)])
    return np.concatenate( (np.zeros(1) , (Area[1:]+Area[:-1])/2 , np.zeros(1)) )

#############################################################################

def moving_average_periodic_FFT(x, w, axis=0):
    """
    Routine designed to smooth periodic data with an odd moving window with size w
    """
    x_hat=np.fft.fft(np.moveaxis(x,axis,0),axis=0)
    
    filter_hat=np.fft.fft(np.ones(w)/w,n=x.shape[axis])   

    x_hat[:]=(x_hat.T * filter_hat).T

    x_hat=np.fft.ifft(x_hat,axis=0).real
    
    return np.moveaxis(np.roll(x_hat,-(w-1)//2,axis=0),0,axis)

#############################################################################

def save_detrended_and_anomalies(startyr, endyr): # Checked 20240404
    """
    This routine assumes that the data have been already saved in daily format
    The trends are computed from start to end of the analysed period
    The climatology is computed from start to end and it includes the 29th of February
    The outputs are stored in a climatology file that is retrieved when needed
    """

    d0 = 29219    # number of days at 1/1/1980
    days = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31] 
    sum_days = np.concatenate((np.zeros(1,dtype=np.int32),np.cumsum(days)),axis=0)  
        
    for year in range(startyr,endyr+1):
        for month in range(1,13):
            print(month,year)

            if os.path.isfile(filePath_absolute(month,year)):
                with open(filePath_absolute(month,year),'rb') as f:
                    F=pickle.load(f)
                Nd=F['data'].shape[0]   #actual number of days

                if month==1 and year==startyr:   # initialisation

                    A,B,X=np.zeros_like(F['data'][0]),np.zeros_like(F['data'][0]),np.arange(0)

                    climatology=[np.zeros((day,*A.shape)) for day in days]
                    climatology_time=[np.zeros(day) for day in days] # days climatology to accelerate the algorithm
                    NN=[0 for _ in days]
                    NNleap=0                    # counter for the 29th of February
                    
                    latitude=F['latitude']
                    longitude=F['longitude']

                # detrending algorithm (this is to avoid to load the entire dataset but rather ony a month at time)               
                s=F['time']//24-d0  # number of days from the 1/1/1980
                A+=np.sum((F['data'].T*s).T,axis=0)
                B+=np.sum(F['data'],axis=0)
                X=np.concatenate((X,s))

                #estimating climatology
                climatology[month-1][:Nd]+=F['data']
                climatology_time[month-1][:Nd]+=F['time']//24
                NN[month-1]+=1
                if Nd==29:     
                    NNleap+=1
            
    X1=np.sum(X)
    X2=np.sum(X**2)
    N=len(X)
    den=N*X2-X1**2
    m=(N*A-X1*B)/den    # slope of trend (K/days)
    q=(X2*B-X1*A)/den   # intercept of trend (K)

    q-=m*d0             # back after the normalisation
    Ymean=B/N           # mean of data
    q-=Ymean            # intercept of detrended data (once we subtract the trend, the mean remains)

    for month in [ 1,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12]:
        climatology[month-1]/=NN[month-1]
        climatology_time[month-1]/=NN[month-1]

    month=2
    climatology[month-1][:-1]/=NN[month-1]
    climatology_time[month-1][:-1]/=NN[month-1]
    if NNleap>0:
        climatology[month-1][-1]/=NNleap
        climatology_time[month-1][-1]/=NNleap
    else: # no 29th of February in the analysed period
        climatology[1][-1]=(climatology[1][-2]+climatology[2][0])/2
        climatology_time[1][-1]=(climatology_time[1][-2]+climatology_time[2][0])/2
    
    # 15 days moving average (with periodic correction)
    climatology=moving_average_periodic_FFT(np.concatenate(climatology,axis=0),15,axis=0)
    climatology_time=moving_average_periodic_FFT(np.concatenate(climatology_time),15,axis=0)

    with open(filePath_trend_climatology(),'wb') as f:
        pickle.dump({'slope':m,'intercept':q,'data_mean':Ymean,
                     'latitude':latitude,'longitude':longitude,
                     'climatology':[climatology[sum_days[i]:sum_days[i+1]] for i in range(12)],
                     'climatology_time':[climatology_time[sum_days[i]:sum_days[i+1]] for i in range(12)] },f, pickle.HIGHEST_PROTOCOL)

#############################################################################

def save_mask(): # Checked 20240404
    from netCDF4 import Dataset
    fileName , original_resolution , desired_resolution , root_path = filePath_mask_ERA5()
    delta=int(desired_resolution/original_resolution) # original resolution 0.1 degrees
    fh = Dataset(fileName, mode='r')
    LSM = np.squeeze(fh.variables['lsm'][:,::delta,::delta])
    fh.close()
    LSM[LSM>0]=1
    LSM[LSM==0]=np.nan
    with open(os.path.join(root_path,'lsm_'+str(desired_resolution)+'x'+str(desired_resolution)+'.pkl'),'wb') as f:
        pickle.dump(LSM,f, pickle.HIGHEST_PROTOCOL)
    
#############################################################################

def ssn_data_load( month_list, yr_list, latind,  data_type='absolute', data_det=False, mask=True): # Checked 20231209
    """
    Main function which defines areas and characteristics of extreme 
    temperatures. Takes as input monthly files of daily 2m temperature and
    temperature anomalies on a regular lat-lon grid.

    To think about, need to include in output info on lat and lon of each
    region to be able to couple e.g. duration to specific regions, as well as
    temp anomalies/abs values

    Input
    month_list: season. Can be provided as a row vector of month using values 1 to 12.
    yr_list: list of years in data.
    latind: latitude indices of selected domain.
    data_type: 'absolute' (absolute values) or 'anomaly' (anomalies).
    data_det: True or False. If True, applies a linear detrend to the raw temperature data.
    mask: True or False. If True, applies a mask where the points over sea are replaced by np.nan

    Output
    var_data: required temperature data, in a single timexlonxlat matrix
    time: time indices (hours from 1/1/1900)
    latitude: Latitude values restricted to latind [deg]
    longitude: Longitude values [deg]
    
    """
    
    if data_type=='anomaly' or data_det:
        with open(filePath_trend_climatology(),'rb') as f:
            G=pickle.load(f)

    if mask:
        with open(filePath_mask(),'rb') as f:
            MASK=pickle.load(f)

    #############################################################

    def get_data(month,year):
        # Script to load one month and apply detrending/anomaly
        
        with open(filePath_absolute(month,year),'rb') as f:
            F=pickle.load(f)

        Nd=F['data'].shape[0]   # number of days in the month
        s=F['time']//24         # number of days from 1/1/1900
        OUTPUT=F['data']        # regular data
        
        if data_type=='anomaly':
            # Regular anomaly
            OUTPUT-=G['climatology'][month-1][:Nd]
            if data_det:
                # detrended anomaly
                OUTPUT-=(np.repeat(G['slope'][np.newaxis,:],Nd,axis=0).T*(s-G['climatology_time'][month-1][:Nd])).T 
        elif data_det:
            # detrended absolute
            OUTPUT-=(np.repeat(G['slope'][np.newaxis,:],Nd,axis=0).T*s).T + G['intercept']

        if mask:
            OUTPUT*=MASK

        return OUTPUT[:,latind],F['time'],F['latitude'][latind],F['longitude']
    
    #############################################################

    var_data=[]
    time_data=[]
    for year in np.sort(yr_list):
        for month in np.sort(month_list):
            if os.path.isfile(filePath_absolute(month,year)):
                U,tempo,lat,lon=get_data(month,year)
                var_data.append(U)
                time_data.append(tempo)
    
    return np.concatenate(var_data,axis=0),np.concatenate(time_data),lat,lon

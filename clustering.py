import numpy as np
from methods import haversine, Area_rectangle_unit , ssn_data_load
from skimage.measure import label, regionprops, find_contours
import multiprocessing as mpr
from concurrent.futures import ProcessPoolExecutor
from functools import partial

#############################################################################

def min_dur_tex(latitude,longitude,tex_gridpoints,duration): # Checked 20231209
    # function to eliminate heatwaves/cold spells shorter than duration days
    Ntime = np.max(tex_gridpoints[0])+1
    Nlat  = len(latitude)
    Nlon  = len(longitude)

    if duration>Ntime:
        return tex_gridpoints
    
    filter=np.ones(duration)
    total=[]
    t=np.zeros((Ntime,Nlon))
    for i in range(Nlat):
        print('Duration identification for latitude index: ',i)
        t[:]=0
        s=np.where(tex_gridpoints[1]==i)[0]
        if len(s)>0:
            t[tex_gridpoints[0][s],tex_gridpoints[2][s]]=1            
            for j in range(Nlon):            
                z=np.where(np.convolve(t[:,j], filter, mode='valid')==duration)[0]
                if len(z)>0:
                    indices = np.unique(np.concatenate([z + q for q in range(duration)]))
                    for index in indices:
                        total.append([index,i,j])
    total=np.vstack(total).T

    return (total[0],total[1],total[2])

#############################################################################

def interpolate_points(P,latitude,longitude): # Checked 20231209
    if len(P)>0:
        # interpolation to get latitude/longitude
        lax, lox=np.floor(P[:,0]).astype(int), np.floor(P[:,1]).astype(int)
        lax[lax==len(latitude)-1]=len(latitude)-2
        lox[lox==len(longitude)-1]=len(longitude)-2
        P[:,0]=latitude[lax] + (latitude[lax+1]-latitude[lax]) * (P[:,0]-lax)
        P[:,1]=longitude[lox] + (longitude[lox+1]-longitude[lox]) * (P[:,1]-lox)
        
#############################################################################

def cluster_analysis_centroid_edges(latitude , longitude , extent_flag , cluster_type , cluster_distance , \
                                    extent , Area , gridpoint): # Checked 20231209
    """
    Function that perform the cluster analysis of the heatwave gridboxes for a single time step (parallelised)
    """

    print('time index: ',gridpoint[0]) 

    Nlat,Nlong=len(latitude),len(longitude)

    A=np.full((Nlat,Nlong),False)
    A[gridpoint[1],gridpoint[2]]=True

    #################################### cluster analysis ####################################

    label_A=label(A)        
    regions_A= regionprops(label_A)
    
    if extent_flag:        # remove small clusters
        for r in regions_A:
            if r.area==1:
                label_A[label_A == r.label] = 0    
        regions_A= regionprops(label_A)

    NA=len(regions_A)
    
    if cluster_type=='centroid':
        
        centroids_A=np.array([r.centroid for r in regions_A],dtype=float)
        interpolate_points(centroids_A,latitude,longitude)
        
        distance_matrix=np.zeros((NA,NA))+1e6
        for ii in range(NA-1):
            distance_matrix[ii,ii+1:]=haversine(centroids_A[ii],centroids_A[ii+1:])

        rows, cols = np.unravel_index(np.argsort(distance_matrix, axis=None), distance_matrix.shape)
        sorted_values = distance_matrix[rows, cols]
        rows=rows[sorted_values<cluster_distance]
        cols=cols[sorted_values<cluster_distance]

    else:

        perimeter_A=[]
        for r in regions_A:
            label_B=np.zeros_like(label_A)
            label_B[label_A == r.label] = 1
            perimeter_A.append(np.vstack(find_contours(label_B,0.5,fully_connected='high'),dtype=float))
            interpolate_points(perimeter_A[-1],latitude,longitude)
            
        rows,cols=[],[]
        for ii in range(NA-1):
            for kk in range(len(perimeter_A[ii])): # search throughout the perimeter points
                for jj in range(ii+1,NA): # search throughout the other regions
                    dist=haversine(perimeter_A[ii][kk],perimeter_A[jj])
                    if np.min(dist)<cluster_distance:
                        rows.append(ii)
                        cols.append(jj)
                        break
        rows,cols=np.array(rows),np.array(cols)

    ####################################################################################################
    
    if len(rows)>0:
        for ii in range(len(rows)):
            i1,i2=rows[ii],cols[ii]
            if i1>i2:
                i1,i2=i2.copy(),i1.copy()
            label_A[label_A == regions_A[i2].label] = regions_A[i1].label
            rows[rows==i2]=i1
            cols[cols==i2]=i1

    regions_A= regionprops(label_A)

    # Area filtering
    B=(label_A.T*Area).T
    for r in regions_A:
        if np.sum(B[label_A == r.label])/r.label<extent:
            label_A[label_A == r.label] = 0
    
    # Cleaning time
    Q=np.unique(label_A)[1:]
    for i in range(len(Q)):
        label_A[label_A == Q[i]] = i+1        
    regions_A=regionprops(label_A) 

    B[:]=(label_A.T*Area).T  
    Areas=np.array([np.sum(B[label_A == r.label])/r.label for r in regions_A])
    LATLON=[np.where(label_A == r.label) for r in regions_A]

    return [gridpoint[0],Areas,LATLON]

#############################################################################

def cluster_tex(time, latitude, longitude, tex_gridpoints, extent_flag=True, cluster_distance=250, \
                dist_type='centroid', extent=10000): # Checked 20231209
    """
    Function which implements clustering of heatwave gridboxes (it just organises the parallelisation)

    Input:
    time:       time values.
    latitude:   latitude values.
    longitude:  longitude values.
    tex_gridpoints: position indeces of temperature extremes satisfying the minimum duration criterion
    extent_flag: True or False. If True, ignores heatwaves detected at single pixels.
    cluster_distance: numerical distance in km. Clusters heatwave gridboxes within "distance" km 
                of each other into a single heatwave region even if they are not connected. Note that clustering 
                occurs before the minimum extent threshold is applied.
    dist_type:  type of distance used for cluster_ex or sep. Can be "centroid" or "edges". 
                If "centroid", then distance is computed between heatwave centroids. 
                If "edges", then distance is computed between closest points of heatwaves.
    extent: minimum extent of heatwave region in km^2. Heatwave regions smaller than this are discarded.
    

    Output:
    time_cluster: time position of the heatwave gridboxes satisfying the minimum duration criterion
    Areas:      List of arrays (N x 1) of heatwave region areas. Each array corresponds to a time step.
    LATLON:     List of arrays (N x 2) of heatwave region gridbox positions (indices). Each array corresponds 
                to a time step and an heat wave.
    """
    ll=np.unique(tex_gridpoints[0])     # list of temporal indices 
    
    Area_unit=Area_rectangle_unit(latitude,longitude) # area of each gridbox in km^2
    
    #############################################################################

    gridpoints=[]
    # conversion of indices for the parallelisation
    for i in range(len(ll)):
        s=np.where(tex_gridpoints[0]==ll[i])[0]
        gridpoints.append([ll[i],tex_gridpoints[1][s],tex_gridpoints[2][s]])

    # parallelisation
    with ProcessPoolExecutor(max_workers=min([mpr.cpu_count(),30])) as executor:
        f = partial(cluster_analysis_centroid_edges, latitude, longitude, extent_flag, dist_type, cluster_distance, extent, Area_unit)
        results = executor.map(f, gridpoints)
            
    time_cluster_index , Areas , LATLON = [],[],[]
    for result in results:
        time_cluster_index.append(result[0])
        Areas.append(result[1])
        LATLON.append(result[2])

    return time[time_cluster_index] , Areas , LATLON

#############################################################################

def t2m_extreme(month_list, yr_list, latind, percentile=5, data_type='anomaly', data_det=False, \
                mask=True, duration=4, extent_flag=True, cluster_ex=1000, \
                dist_type='centroid', extent=2e5): # Checked 20231209
    """
    Main function which defines areas and characteristics of extreme 
    temperatures. Takes as input montly files of daily 2m temperature and
    temperature anomalies on a regular lat-lon grid.

    Input
    month_list: Can be provided as a list of month using values 1 to 12
    yr_list:    Can be provided as a list of years 
    latind:     Indices of latitudes to be considered
    percentile: integer between 0 and 100, fixing percentile used to define
                heatwaves/cold spells. Assumes that if percentile is larger than
                the median, then one is looking at heatwaves and that if the
                percentile is smaller than the median then one is looking at cold
                spells.
    data_type:  'absolute' or 'anomaly'
    data_det:   True or False. If True, applies a linear detrend to the raw
                temperature data.
    mask:       True or False. If True, applies a mask of NaN above water
    duration:   minimum duration (consecutive percentile exceedances) of the
                heatwaves/cold spells
    extent_flag: True or False. If True, ignores heatwaves detected at single 
                pixels. This is useful to speed up the calculation. 
    cluster_ex: numerical distance in km . It clusters heatwave gridboxes 
                within "distance" km of each other into a single heatwave region 
                even if they are not connected. Note that clustering occurs 
                before the minimum extent threshold is applied.
    dist_type:  type of distance used for cluster_ex. Can be "centroid" 
                or "edges". If "centroid", then distance is computed between 
                heatwave centroids. If "edges", then distance is computed 
                between closest points of heatwaves.'centroid'
    extent:     minimum areal extent of the heatwave/cold spell region, in km2

    Output:
    time_cluster:   Time position
    Areas:          List (along time) of areas lists (heatwave/cold spell regions) for each time
    average_variable_all: Area averaged variable (heatwave/cold spell regions) for each time (list of lists)
    LATLON:         List (along time) of LATLON points lists (heatwave/cold spell regions) for each time
    """

    var_data,time,latitude,longitude=ssn_data_load(month_list, yr_list, latind, data_type, data_det, mask)
    
    pct_data = np.percentile(var_data, percentile, axis=0)

    if percentile>50:       # Heatwave
        tex_gridpoints = np.where(var_data > pct_data)    
    else:                   # Cold Spell
        tex_gridpoints = np.where(var_data < pct_data)

    tex_gridpoints = min_dur_tex(latitude, longitude, tex_gridpoints, duration)    # duration filter

    time_cluster, Areas, LATLON = cluster_tex(time, latitude, longitude, tex_gridpoints, extent_flag, cluster_ex, dist_type, extent)
    
    # attention here since LATLON are the indeces of the gridpoints and not the latitude/longitude

    Area_unit = Area_rectangle_unit(latitude,longitude)

    average_variable_all=[]
    for ii in range(len(time_cluster)):
        average_variable=[]
        if len(Areas[ii])>0:
            Q=np.squeeze(var_data[time==time_cluster[ii]])
            for jj in range(len(Areas[ii])):
                S=LATLON[ii][jj]
                SS=Area_unit[S[0]]
                average_variable.append(np.sum(np.array(Q[S])*SS)/ np.sum(SS))
                LATLON[ii][jj]=np.array([S[0],S[1]],dtype=float).T # conversion in latitude/longitude
                interpolate_points(LATLON[ii][jj],latitude,longitude)
            
        average_variable_all.append(average_variable)

    return time, time_cluster, Areas, average_variable_all, LATLON

#############################################################################

# if __name__=='__main__':
#     import matplotlib.pyplot as plt
#     import time

#     #############################################################################

#     # # test for the duration sorting
#     # A=[0,0,0,1,1,1,1,1,0,0,0,1,1,0,1,0,0,0,1,1,1]
#     # duration=4
#     # Nx,Ny=1,2
#     # lst=np.array(A)
#     # lst=np.repeat(lst[...,np.newaxis],Nx,axis=-1)
#     # lst=np.repeat(lst[...,np.newaxis],Ny,axis=-1)
#     # Q=np.ones((Nx,Ny))*0.1
#     # tex_gridpoints=np.where(lst>Q)
#     # print(tex_gridpoints)
#     # p=min_dur_tex(tex_gridpoints,duration)
#     # print(p)
    
#     #############################################################################

#     # var_data,time,latitude,longitude=ssn_data_load([9], [2020], np.arange(361), 'anomaly', False, False)

#     # X=var_data[6]
#     # plt.subplot(131)
#     # plt.contourf(longitude,latitude,X,30);plt.colorbar()
#     # plt.contour(longitude,latitude,X,levels=[-4],colors='r')
#     # Y=np.zeros_like(X)
#     # Y[X<-4]=1
#     # plt.subplot(132)
#     # plt.contourf(longitude,latitude,Y,30);plt.colorbar()
#     # plt.subplot(133)
#     # plt.contour(longitude,latitude,X,levels=[-4],colors='r')
#     # Y=np.repeat(Y[np.newaxis],1,axis=0)
#     # tex_gridpoints=np.where(Y==1)

#     # gridpoint=[]
#     # # conversion of indices for the parallelisation
#     # s=np.where(tex_gridpoints[0]==0)[0]
#     # gridpoint.append([0,tex_gridpoints[1][s],tex_gridpoints[2][s]])

#     # Area=Area_rectangle_unit(latitude,longitude)
#     # HH=cluster_analysis_centroid_edges(latitude , longitude , True , 'centroid' , 2000 , 10000 , Area , gridpoint[0])

#     # istante=HH[0]
#     # Areas=HH[1]
#     # LATLON=HH[2]
    
#     # for i in range(len(LATLON)):
#     #     LATLON[i]=np.array([LATLON[i][0],LATLON[i][1]]).T # conversion in latitude/longitude
#     #     interpolate_points(LATLON[i],latitude,longitude)
#     #     LL=LATLON[i]
#     #     plt.plot(LL[:,1],LL[:,0],'o')
#     #     plt.pause(0.1)
#     #     plt.title(i)

#     # plt.show()
    
#     #############################################################################   
            
#     tt=time.time()
#     X=t2m_extreme([1,2,12], np.arange(1980,1990), np.arange(361),  10, 'anomaly', True, True, 3, True, 400,'centroid',10000)
#     print(time.time()-tt)

#     # pos=15
#     # print(X[0][pos],'\n')
#     # print(X[1][pos],'\n')
#     # print(X[2][pos],'\n')
#     # for i in range(len(X[3][pos])):
#     #     LL=X[3][pos][i]
#     #     plt.plot(LL[:,1],LL[:,0],'o')
#     #     plt.pause(1)
#     # plt.show()
    
#     #############################################################################
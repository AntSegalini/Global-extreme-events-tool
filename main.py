import numpy as np
import datetime as dt
from methods import Area_rectangle_unit , save_detrended_and_anomalies , save_mask
from clustering import t2m_extreme
import pickle
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from scipy.stats import linregress
from scipy.ndimage import gaussian_filter

#############################################################################

# actual resolution of the data (0.5 deg x 0.5 deg)
latitudeb  = -np.linspace(-90,90,361)                       # original latitude [-90,90]
longitude  = np.arange(720)/2                               # original longitude [0,360]   
Area_unitb = Area_rectangle_unit(latitudeb,longitude)       # area of the grid cells
longitude  = np.fft.fftshift(np.mod(longitude+180,360)-180) # shift the longitude to [-180,180]

#############################################################################

# useful functions

def index_latlon_NH(lat,lon):
    return (((70-lat)*2).astype(int),(lon*2).astype(int))

def index_latlon_SH(lat,lon):
    return (((-30-lat)*2).astype(int),(lon*2).astype(int))

def plot_map(Nx,Ny):
    fig,axs=plt.subplots(Nx,Ny,figsize=(6*Ny,3*Nx), subplot_kw={'projection': ccrs.PlateCarree()})
    axs=axs.flatten()
    for ax in axs:
        ax.set_extent([-180, 180, -90, 90], crs=ccrs.PlateCarree())
        ax.coastlines(resolution='110m', color='black', linewidth=1, linestyle='solid', zorder=10)

        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linestyle='--', linewidth=0.5, color='gray', alpha=0.7)
        gl.top_labels, gl.bottom_labels = False,False
        gl.right_labels , gl.left_labels = False,False

        # Set ticks
        ax.set_xticks(np.arange(-180, 181, 60), crs=ccrs.PlateCarree())
        ax.set_yticks(range(-90, 91, 30), crs=ccrs.PlateCarree())
        ax.set_xticklabels(['180°W','120°W','60°W','0°','60°E','120°E','180°E'], fontsize=12)
    return fig,axs

#############################################################################

save_detrended_and_anomalies(1940,2023)

#############################################################################

save_mask()

#############################################################################

# this is the main script to run the clustering

yr_list=np.arange(1940,2024)    # list of the years to analyse
data_type='anomaly'             # type of data to analyse
data_det=False                  # if True, the data are detrended
mask=True                       # if True, the land-sea mask is applied
duration=4                      # minimum duration of the event
extent_flag=True                # if True, the extent of the event is used to filter the events
cluster_ex=1000                 
dist_type='centroid'            # type of distance to use in the clustering
extent=2e5                      

for ii in range(8):
    # Change in percentile 5/95
    # Northern hemisphere
    if ii==0:
        month_list=[1,2,12]; percentile=5; Name='ColdSpell_NH_5'
        latind=np.where((latitudeb>=30) & (latitudeb<=70))[0]
    elif ii==1:
        month_list=[6,7,8]; percentile=95; Name='HeatWave_NH_5'
        latind=np.where((latitudeb>=30) & (latitudeb<=70))[0]

    # Southern hemisphere
    elif ii==2:
        month_list=[6,7,8]; percentile=5; Name='ColdSpell_SH_5'
        latind=np.where((latitudeb>=-70) & (latitudeb<=-30))[0]
    elif ii==3:
        month_list=[1,2,12]; percentile=95; Name='HeatWave_SH_5'
        latind=np.where((latitudeb>=-70) & (latitudeb<=-30))[0]   

    # Change in percentile 10/90
    # Northern hemisphere
    if ii==4:
        month_list=[1,2,12]; percentile=10; Name='ColdSpell_NH_10'
        latind=np.where((latitudeb>=30) & (latitudeb<=70))[0]
    elif ii==5:
        month_list=[6,7,8]; percentile=90; Name='HeatWave_NH_10'
        latind=np.where((latitudeb>=30) & (latitudeb<=70))[0]

    # Southern hemisphere
    elif ii==6:
        month_list=[6,7,8]; percentile=10; Name='ColdSpell_SH_10'
        latind=np.where((latitudeb>=-70) & (latitudeb<=-30))[0]
    elif ii==7:
        month_list=[1,2,12]; percentile=90; Name='HeatWave_SH_10'
        latind=np.where((latitudeb>=-70) & (latitudeb<=-30))[0]   

    time, time_events, Areas, average_variable_all, LATLON = t2m_extreme(month_list, yr_list, latind, percentile, data_type, data_det, mask, duration, extent_flag, cluster_ex, dist_type, extent)
    
    with open(Name+'.pkl','wb') as f:
        pickle.dump({'LATLON':LATLON,'time_original':time,'time_events':time_events,'Area': Areas,'Tmean':average_variable_all},f)

  
#############################################################################

for ii in range(4):
    # Northern hemisphere
    if ii==0:
        month_list=[1,2,12]; percentile=5; Name='ColdSpell_NH_5'
        latind=np.where((latitudeb>=30) & (latitudeb<=70))[0]
    elif ii==1:
        month_list=[6,7,8]; percentile=95; Name='HeatWave_NH_5'
        latind=np.where((latitudeb>=30) & (latitudeb<=70))[0]
    # Southern hemisphere
    elif ii==2:
        month_list=[6,7,8]; percentile=5; Name='ColdSpell_SH_5'
        latind=np.where((latitudeb>=-70) & (latitudeb<=-30))[0]
    elif ii==3:
        month_list=[1,2,12]; percentile=95; Name='HeatWave_SH_5'
        latind=np.where((latitudeb>=-70) & (latitudeb<=-30))[0] 
    var_data,time,latitude,longitude=ssn_data_load(month_list, yr_list, latind, data_type, data_det, mask)
    print('ciao')
    with open(Name+'_variance.pkl','wb') as f:
        pickle.dump({'VAR':np.var(var_data,axis=0),'latitude':latitude,'longitude':longitude},f)
    var_data=None
exit()
#############################################################################

# this part creates Figures 2, S1, S3

# for co_occurrence,percentile,Name_fig in [True,5,'Figure_2.png'],[False,5,'Figure_S3.png']:#,[True,10,'Figure_S1.png']:

#     fig,axs=plot_map(3,2)
#     colore={'HeatWave_NH':'Reds','ColdSpell_NH':'Blues','HeatWave_SH':'Reds','ColdSpell_SH':'Blues'}

#     jj=0
#     for Name in ['HeatWave_NH','ColdSpell_NH','HeatWave_SH','ColdSpell_SH']:
        
#         with open(Name+'_'+str(percentile)+'.pkl','rb') as f:
#             F=pickle.load(f)
            
#         if 'NH' in Name:    latind=np.where((latitudeb>=30) & (latitudeb<=70))[0]
#         else:               latind=np.where((latitudeb>=-70) & (latitudeb<=-30))[0]
#         latitude=latitudeb[latind]

#         A=np.zeros((len(latitude),len(longitude)))
#         Nevents=np.zeros_like(A)
#         T=np.zeros_like(A)
#         S0=[]

#         for Nj in range(len(F['time_events'])):
#             NN=len(F['LATLON'][Nj])

#             S=[]
#             if NN>0: #there are heatwaves
#                 j=0
#                 for PP in F['LATLON'][Nj]:
#                     if 'NH' in Name:       S.append(index_latlon_NH(PP[:,0],PP[:,1]))
#                     else:                  S.append(index_latlon_SH(PP[:,0],PP[:,1]))
#                     if NN>co_occurrence: A[S[-1]]+=1; T[S[-1]]+=F['Tmean'][Nj][j]   # concurrent events    
#                     j+=1      
                
#                 S=[tuple(row) for row in np.hstack(S).T] # S is a list of tuples (lat,lon) with the indices of the events

#                 S0=np.array(list(set(S)-set(S0)))
#                 if len(S0)>0: Nevents[S0[:,0],S0[:,1]]+=1 
#             S0=S 

#         # relative frequency of events
#         A[A==0]=np.nan
#         PO=axs[jj].contourf(longitude,latitude,np.fft.fftshift(A/(0.01*percentile*len(F['time_original'])),axes=1),np.linspace(0,0.7,8),transform=ccrs.PlateCarree(),cmap=plt.get_cmap(colore[Name]))
#         if 'NH' in Name: plt.colorbar(PO, ax=axs[jj], orientation='vertical', fraction=0.05, pad=0.05, aspect=30, shrink=0.8)
#         if 'HeatWave_NH' in Name: axs[jj].set_title('(a) Relative frequency heatwaves',fontsize=16)
#         if 'ColdSpell_NH' in Name: axs[jj].set_title('(b) Relative frequency cold spells',fontsize=16)

#         A[Nevents==0]=np.nan
#         Nevents[A==np.nan]=np.nan   
#         # duration
#         PO=axs[jj+2].contourf(longitude,latitude,np.fft.fftshift(A/Nevents,axes=1),np.arange(2,8),transform=ccrs.PlateCarree(),cmap=plt.get_cmap(colore[Name]),extend='max')
#         if 'NH' in Name: plt.colorbar(PO, ax=axs[jj+2], orientation='vertical', fraction=0.05, pad=0.05, aspect=30, shrink=0.8)
#         if 'HeatWave_NH' in Name: axs[jj+2].set_title('(c) Duration heatwaves [days]',fontsize=16)
#         if 'ColdSpell_NH' in Name: axs[jj+2].set_title('(d) Duration cold spells [days]',fontsize=16)
        
#         # temperature anomaly
#         if 'HeatWave' in Name: cc=np.arange(5,16);qq=''
#         if 'ColdSpell' in Name: cc=np.arange(-20,-4);qq='_r'
#         PO=axs[jj+4].contourf(longitude,latitude,np.fft.fftshift(T/A,axes=1),cc,transform=ccrs.PlateCarree(),cmap=plt.get_cmap(colore[Name]+qq))
#         if 'NH' in Name: plt.colorbar(PO, ax=axs[jj+4], orientation='vertical', fraction=0.05, pad=0.05, aspect=30, shrink=0.8)
#         if 'HeatWave_NH' in Name: axs[jj+4].set_title('(e) Temperature anomaly heatwaves [K]',fontsize=16)
#         if 'ColdSpell_NH' in Name: axs[jj+4].set_title('(f) Temperature anomaly cold spells [K]',fontsize=16)
            
#         jj=np.mod(jj+1,2)
        
#     # plt.subplot_tool()
#     plt.subplots_adjust(left=0.033,bottom=0.052,right=0.952,top=0.955,wspace=0.348, hspace=0.207)
#     # plt.tight_layout(rect=[0, 0, 1, 0.95])
#     plt.savefig(Name_fig,dpi=300)
#     plt.close()

#############################################################################               

# this part creates Figures 3, S4, S5, S6, S7, S8

# for co_occurrence,Name_fig,EM in [True,'Figure_3.png','NH_5'],[True,'Figure_S4.png','SH_5'],[False,'Figure_S7.png','NH_5'],[False,'Figure_S8.png','SH_5']:#,[True,'Figure_S5.png','NH_10'],[True,'Figure_S6.png','SH_10']:
#     plt.figure(1,figsize=(12,8))
#     plt.figure(2,figsize=(12,8))
#     jj=0
#     for Name in ['HeatWave_'+EM,'ColdSpell_'+EM]:
        
#         jj=1+np.mod(jj,2)

#         if 'NH' in Name:    latind=np.where((latitudeb>=30) & (latitudeb<=70))[0]
#         else:               latind=np.where((latitudeb>=-70) & (latitudeb<=-30))[0]
        
#         latitude=latitudeb[latind]
#         Area_unit=Area_unitb[latind]
        
#         A=np.zeros((len(latitude),len(longitude)))
#         Ab=np.zeros_like(A)
#         Nevents=np.zeros_like(A)

#         with open(Name+'.pkl','rb') as f:
#             F=pickle.load(f)
#         year=[(dt.datetime(1900,1,1)+dt.timedelta(days=int(i/24))).year for i in F['time_events']]

#         NUMBER , DURATION , ANOMALY , AREA=[],[],[],[]
#         NUMBER_all , DURATION_all , ANOMALY_all , AREA_all=[],[],[],[]
        
#         for yy in np.unique(year):
#             z=np.where(np.array(year)==yy)[0]

#             NTOT_part , ANOMALY_part , AREA_part , NUMBER_part = [],[],[],0
#             NTOT_partb , ANOMALY_partb , AREA_partb , NUMBER_partb = [],[],[],0

#             A[:]=0 
#             Ab[:]=0     
#             Nevents[:]=0
#             S0=[]
#             for i in z: # all the events in one year
#                 NN=len(F['Area'][i]) # number of heatwaves active at that day
#                 if NN>0:
#                     NUMBER_partb+=1       # number of days with heatwaves
#                     ANOMALY_partb.append(np.sum(F['Tmean'][i]*F['Area'][i])) # sum of the anomalies weighted by the area
#                     AREA_partb.append(np.sum(F['Area'][i]))
#                 if NN>co_occurrence:
#                     # NTOT_part.append(NN) # number of detected heatwaves
#                     NUMBER_part+=1       # number of days with heatwaves
#                     ANOMALY_part.append(np.sum(F['Tmean'][i]*F['Area'][i])) # sum of the anomalies weighted by the area
#                     AREA_part.append(np.sum(F['Area'][i]))
                
#                 #duration analysis
#                 S=[]
#                 if NN>0: # there are heatwaves
#                     j=0
#                     for PP in F['LATLON'][i]:
#                         if 'NH' in Name:       S.append(index_latlon_NH(PP[:,0],PP[:,1]))
#                         else:                  S.append(index_latlon_SH(PP[:,0],PP[:,1]))
#                         if NN>co_occurrence: A[S[-1]]+=1    # number of concurrent events
#                         if NN>0: Ab[S[-1]]+=1    # number of events
#                         j+=1
                    
#                     S=[tuple(row) for row in np.hstack(S).T] # S is a list of tuples (lat,lon) with the indices of the events

#                     S0=np.array(list(set(S)-set(S0)))
#                     if len(S0)>0: Nevents[S0[:,0],S0[:,1]]+=1 # Number of events 
#                 S0=S 

#             NUMBER.append(NUMBER_part)
#             if np.sum(AREA_part)>0:
#                 ANOMALY.append(np.sum(ANOMALY_part)/np.sum(AREA_part))
#             else:
#                 ANOMALY.append(np.nan)
#                 # print('No events: ',Name, yy,NUMBER_part)
#             AREA.append(np.sum(AREA_part))

#             NUMBER_all.append(NUMBER_partb)
#             if np.sum(AREA_partb)>0:
#                 ANOMALY_all.append(np.sum(ANOMALY_partb)/np.sum(AREA_partb))
#             else:
#                 ANOMALY_all.append(np.nan)
#                 # print('No events: ',Name, yy,NUMBER_part)
#             AREA_all.append(np.sum(AREA_partb))

#             # relative frequency of events
#             Nevents[Nevents==0]=np.nan
#             A[np.isnan(Nevents)]=np.nan 
#             A[np.isfinite(Nevents)]/=Nevents[np.isfinite(Nevents)]
#             DURATION.append(np.nansum(A.T*Area_unit))
#             A[~np.isnan(A)]=1
#             DURATION[-1]/=np.nansum(np.repeat(Area_unit[:,np.newaxis],len(longitude),axis=1)*A)            
            
#             # relative frequency of events
#             Ab[np.isnan(Nevents)]=np.nan 
#             Ab[np.isfinite(Nevents)]/=Nevents[np.isfinite(Nevents)]
#             DURATION_all.append(np.nansum(Ab.T*Area_unit))
#             Ab[~np.isnan(Ab)]=1
#             DURATION_all[-1]/=np.nansum(np.repeat(Area_unit[:,np.newaxis],len(longitude),axis=1)*Ab)            
            
#             if yy==2023 and 'Cold' in Name: 
#                 NUMBER[-1]=np.nan
#                 DURATION[-1]=np.nan
#                 ANOMALY[-1]=np.nan
#                 AREA[-1]=np.nan 

#                 NUMBER_all[-1]=np.nan
#                 DURATION_all[-1]=np.nan
#                 ANOMALY_all[-1]=np.nan
#                 AREA_all[-1]=np.nan
                
#         year=np.unique(year)

#         def plot_error_bar(year,K,print_stuff=False):
#             fx=np.array(year)[np.isfinite(K)]
#             fy=np.array(K)[np.isfinite(K)]
#             c=np.polyfit(fx,fy,1)
#             if print_stuff:    print(co_occurrence,Name_fig,EM,c, np.polyval(c,2010)/np.polyval(c,1950))
#             sy=np.sqrt(np.sum((fy-np.polyval(c,fx))**2)/(len(fx)-2))
#             sxx=np.sum(fx**2)-np.sum(fx)**2/len(fx)
#             plt.plot(year,np.polyval(c,year),'-r')
#             plt.plot(year,np.polyval(c,year)+1.96*sy*np.sqrt(1/len(fx)+((year-np.mean(fx))**2)/sxx),':r')
#             plt.plot(year,np.polyval(c,year)-1.96*sy*np.sqrt(1/len(fx)+((year-np.mean(fx))**2)/sxx),':r')
#             result=linregress(fx,fy)
#             if result.pvalue<0.05: return True,c[0]
#             else: return False,c[0]
        
#         plt.figure(1)
#         plt.subplot(4,2,jj)
#         K=NUMBER
#         plt.plot(year,K,'ob')
#         dd,slope=plot_error_bar(year,K,True);    hh=''
#         if dd: hh=' *'
#         cx,cy=plt.gca().get_xlim(),plt.gca().get_ylim()
#         plt.text(cx[0]+(0.05+(jj-1)/1.8)*(cx[1]-cx[0]),cy[0]+0.8*(cy[1]-cy[0]),'slope='+str(np.round(slope,3))+' [1/year]',fontsize=10)
#         if jj==1: plt.title('(a) Number of '+EM[:2]+' heatwaves'+hh,fontsize=12)
#         else: plt.title('(b) Number of '+EM[:2]+' cold spells'+hh,fontsize=12)
#         # plt.ylim(0,95)    

#         plt.subplot(4,2,jj+2)
#         K=[]
#         for i in DURATION:
#             if np.isfinite(i):
#                 if i>0.1:                    K.append(i)
#                 else:                    K.append(np.nan)
#             else:                K.append(np.nan)
        
#         plt.plot(year,K,'ob')
#         dd,slope=plot_error_bar(year,K);    hh=''
#         if dd: hh=' *'
#         cx,cy=plt.gca().get_xlim(),plt.gca().get_ylim()
#         plt.text(cx[0]+(0.05+(jj-1)/1.8)*(cx[1]-cx[0]),cy[0]+0.8*(cy[1]-cy[0]),'slope='+str(np.round(slope,3))+' [days/year]',fontsize=10)
#         if jj==1: plt.title('(c) Duration of '+EM[:2]+' heatwaves [days]'+hh,fontsize=12)
#         else: plt.title('(d) Duration of '+EM[:2]+' cold spells [days]'+hh,fontsize=12)

#         plt.subplot(4,2,jj+4)
#         K=ANOMALY
#         plt.plot(year,K,'ob')
#         dd,slope=plot_error_bar(year,K);    hh=''
#         if dd: hh=' *'
#         cx,cy=plt.gca().get_xlim(),plt.gca().get_ylim()
#         plt.text(cx[0]+(0.05+(jj-1)/1.8)*(cx[1]-cx[0]),cy[0]+0.8*(cy[1]-cy[0]),'slope='+str(np.round(slope,3))+' [K/year]',fontsize=10)
#         if jj==1: plt.title('(e) Temp anom of '+EM[:2]+' heatwaves [K]'+hh,fontsize=12)
#         else: plt.title('(f) Temp anom of '+EM[:2]+' cold spells [K]'+hh,fontsize=12)

#         plt.subplot(4,2,jj+6)
#         K=[]
#         for i in AREA:
#             if np.isfinite(i):
#                 if i>1e4:                    K.append(i)
#                 else:                    K.append(np.nan)
#             else:                K.append(np.nan)

#         plt.plot(year,K,'ob')
#         dd,slope=plot_error_bar(year,K);    hh=''
#         if dd: hh=' *'
#         cx,cy=plt.gca().get_xlim(),plt.gca().get_ylim()
#         plt.text(cx[0]+(0.05+(jj-1)/2)*(cx[1]-cx[0]),cy[0]+0.8*(cy[1]-cy[0]),'slope='+str(np.round(slope))+' [km$^2$/year]',fontsize=10)
#         if jj==1: plt.title('(g) Area of '+EM[:2]+' heatwaves [km$^2$]'+hh,fontsize=12)
#         else: plt.title('(h) Area of '+EM[:2]+' cold spells [km$^2$]'+hh,fontsize=12)
#         # plt.ylim(0,3.5e8)

#         # figure with ratio
#         plt.figure(2)
#         plt.subplot(4,2,jj)
#         K=np.array(NUMBER)/np.array(NUMBER_all)
#         plt.plot(year,K,'ob')
#         dd,slope=plot_error_bar(year,K);    hh=''
#         if dd: hh=' *'
#         if jj==1: plt.title('(a) Number of '+EM[:2]+' heatwaves'+hh,fontsize=12)
#         else: plt.title('(b) Number of '+EM[:2]+' cold spells'+hh,fontsize=12)
#         # plt.ylim(0,95)    

#         plt.subplot(4,2,jj+2)
#         K=[]
#         for i in np.array(DURATION)/np.array(DURATION_all):
#             if np.isfinite(i):
#                 if i>0.001:                    K.append(i)
#                 else:                    K.append(np.nan)
#             else:                K.append(np.nan)
        
#         plt.plot(year,K,'ob')
#         dd,slope=plot_error_bar(year,K);    hh=''
#         if dd: hh=' *'
#         if jj==1: plt.title('(c) Duration of '+EM[:2]+' heatwaves'+hh,fontsize=12)
#         else: plt.title('(d) Duration of '+EM[:2]+' cold spells'+hh,fontsize=12)

#         plt.subplot(4,2,jj+4)
#         K=np.array(ANOMALY)/np.array(ANOMALY_all)
#         plt.plot(year,K,'ob')
#         dd,slope=plot_error_bar(year,K);    hh=''
#         if dd: hh=' *'
#         if jj==1: plt.title('(e) Temp anom of '+EM[:2]+' heatwaves'+hh,fontsize=12)
#         else: plt.title('(f) Temp anom of '+EM[:2]+' cold spells'+hh,fontsize=12)

#         plt.subplot(4,2,jj+6)
#         K=[]
#         for i in np.array(AREA)/np.array(AREA_all):
#             if np.isfinite(i):
#                 if i>0:                    K.append(i)
#                 else:                    K.append(np.nan)
#             else:                K.append(np.nan)

#         plt.plot(year,K,'ob')
#         dd,slope=plot_error_bar(year,K);    hh=''
#         if dd: hh=' *'
#         if jj==1: plt.title('(g) Area of '+EM[:2]+' heatwaves [km$^2$]'+hh,fontsize=12)
#         else: plt.title('(h) Area of '+EM[:2]+' cold spells [km$^2$]'+hh,fontsize=12)
#     plt.figure(1)
#     plt.tight_layout()
#     plt.savefig(Name_fig,dpi=300)
#     plt.figure(2)
#     plt.tight_layout()
#     plt.savefig('ratio'+Name_fig,dpi=300)
#     plt.close()    

#############################################################################               

for co_occurrence,percentile,Name_fig in [True,5,'Figure_4.png'],[False,5,'Figure_S10.png']:#,[True,10,'Figure_S9.png']:
    fig,axs=plot_map(3,2)
    jj=0
    for Name in ['HeatWave_NH','ColdSpell_NH','HeatWave_SH','ColdSpell_SH']:
        Ne,inte,dura=[],[],[]
        
        with open(Name+'_'+str(percentile)+'.pkl','rb') as f:
            F=pickle.load(f)

        year=[(dt.datetime(1900,1,1)+dt.timedelta(days=int(i/24))).year for i in F['time_events']]
            
        if 'NH' in Name:    latind=np.where((latitudeb>=30) & (latitudeb<=70))[0]
        else:               latind=np.where((latitudeb>=-70) & (latitudeb<=-30))[0]
        latitude=latitudeb[latind]

        for yy in np.unique(year):
            A=np.zeros((len(latitude),len(longitude)))
            Nevents=np.zeros_like(A)
            T=np.zeros_like(A)
            S0=[]
            z=np.where(np.array(year)==yy)[0]
            
            for Nj in z:
                NN=len(F['LATLON'][Nj])
                S=[]
                if NN>0: 
                    j=0
                    for PP in F['LATLON'][Nj]:
                        if 'NH' in Name:          S.append(index_latlon_NH(PP[:,0],PP[:,1]))
                        else:                     S.append(index_latlon_SH(PP[:,0],PP[:,1]))
                        if NN>co_occurrence: A[S[-1]]+=1; T[S[-1]]+=F['Tmean'][Nj][j]   # concurrent events  
                        j+=1      
                    S=np.hstack(S).T
                    S=[tuple(row) for row in S] # S is a list of tuples (lat,lon) with the indices of the events

                    S0=np.array(list(set(S)-set(S0)))
                    if len(S0)>0: Nevents[S0[:,0],S0[:,1]]+=1 # Number of events                     
                S0=S 
                
            # relative frequency of events
            A[A==0]=np.nan
            Ne.append(A.copy())
            A[A==0]=np.nan
            A[Nevents==0]=np.nan
            Nevents[A==np.nan]=np.nan   
            dura.append(A/Nevents)
            inte.append(T/A)
        Ne=np.stack(Ne)
        dura=np.stack(dura)
        inte=np.stack(inte)

        A[:]=0
        B=np.zeros_like(A)
        Nevents[:]=0
        Ax=np.full(A.shape,False)
        Bx=np.full(A.shape,False)
        Neventsx=np.full(A.shape,False)

        for i in range(A.shape[0]):
            for j in range(A.shape[1]):
                s=Ne[:,i,j]
                d=dura[:,i,j]
                ii=inte[:,i,j]
                if len(s[np.isfinite(s)])>1:
                    c=linregress((np.arange(1940,2024))[np.isfinite(s)],s[np.isfinite(s)])
                    Nevents[i,j]=c.slope*10
                    if c.pvalue<0.05: Neventsx[i,j]=True
                else:
                    Nevents[i,j]=np.nan
                
                if len(d[np.isfinite(d)])>1:
                    c=linregress((np.arange(1940,2024))[np.isfinite(d)],d[np.isfinite(d)])
                    A[i,j]=c.slope*10
                    if c.pvalue<0.05: Ax[i,j]=True
                else:
                    A[i,j]=np.nan

                if len(ii[np.isfinite(ii)])>1:
                    c=linregress((np.arange(1940,2024))[np.isfinite(ii)],ii[np.isfinite(ii)])
                    B[i,j]=c.slope*10
                    if c.pvalue<0.05: Bx[i,j]=True
                else:
                    B[i,j]=np.nan
        # relative frequency of events
        # plt.rcParams['hatch.color'] = 'green'
        Nevents[~np.isnan(Nevents)]=gaussian_filter(Nevents[~np.isnan(Nevents)],sigma=1)
        PO=axs[jj].contourf(longitude,latitude,np.fft.fftshift(Nevents,axes=1),np.linspace(-1,1,21),transform=ccrs.PlateCarree(),extend='both',cmap=plt.get_cmap('coolwarm'))
        if 'NH' in Name: plt.colorbar(PO, ax=axs[jj], orientation='vertical', fraction=0.05, pad=0.05, aspect=30, shrink=0.8)
        if 'HeatWave_NH' in Name: axs[jj].set_title('(a) Number of heatwaves [events/decade]',fontsize=16)
        if 'ColdSpell_NH' in Name: axs[jj].set_title('(b) Number of cold spells [events/decade]',fontsize=16)
        axs[jj].contourf(longitude,latitude,np.fft.fftshift(Neventsx,axes=1), levels=[0.1, 1.5], colors='none', hatches=['/////////', '/////////'],transform=ccrs.PlateCarree())
        
        # duration
        A[~np.isnan(A)]=gaussian_filter(A[~np.isnan(A)],sigma=1)
        PO=axs[jj+2].contourf(longitude,latitude,np.fft.fftshift(A,axes=1),np.linspace(-1,1,21),transform=ccrs.PlateCarree(),extend='both',cmap=plt.get_cmap('coolwarm'))
        if 'NH' in Name: plt.colorbar(PO, ax=axs[jj+2], orientation='vertical', fraction=0.05, pad=0.05, aspect=30, shrink=0.8)
        if 'HeatWave_NH' in Name: axs[jj+2].set_title('(c) Heatwaves duration [days/decade]',fontsize=16)
        if 'ColdSpell_NH' in Name: axs[jj+2].set_title('(d) Cold spells duration [days/decade]',fontsize=16)
        axs[jj+2].contourf(longitude,latitude,np.fft.fftshift(Ax,axes=1), levels=[0.1, 1.5], colors='none', hatches=['/////////', '/////////'],transform=ccrs.PlateCarree())
        
        # temperature anomaly
        B[~np.isnan(B)]=gaussian_filter(B[~np.isnan(B)],sigma=1)
        PO=axs[jj+4].contourf(longitude,latitude,np.fft.fftshift(B,axes=1),np.linspace(-1,1,21),transform=ccrs.PlateCarree(),extend='both',cmap=plt.get_cmap('coolwarm'))
        if 'NH' in Name: plt.colorbar(PO, ax=axs[jj+4], orientation='vertical', fraction=0.05, pad=0.05, aspect=30, shrink=0.8)
        if 'HeatWave_NH' in Name: axs[jj+4].set_title('(e) Heatwave temperature [K/decade]',fontsize=16)
        if 'ColdSpell_NH' in Name: axs[jj+4].set_title('(f) Cold spell temperature [K/decade]',fontsize=16)
        axs[jj+4].contourf(longitude,latitude,np.fft.fftshift(Bx,axes=1), levels=[0.1, 1.5], colors='none', hatches=['/////////', '/////////'],transform=ccrs.PlateCarree())
        
        jj=np.mod(jj+1,2)

    # plt.subplot_tool()
    plt.subplots_adjust(left=0.033,bottom=0.052,right=0.952,top=0.955,wspace=0.348, hspace=0.207)
    # plt.tight_layout(rect=[0, 0, 1, 0.95])

    plt.savefig(Name_fig,dpi=300)
    
    plt.close()

#############################################################################

fig,axs=plot_map(1,2)
with open('HeatWave_NH_5_variance.pkl','rb') as f:
    F=pickle.load(f)
h=F['VAR']
PO=axs[0].contourf(longitude,F['latitude'],np.fft.fftshift(h,axes=1),np.arange(21),transform=ccrs.PlateCarree(),extend='both',cmap=plt.get_cmap('Reds'))
plt.colorbar(PO, ax=axs[0], orientation='vertical', fraction=0.05, pad=0.05, aspect=30, shrink=0.8)
axs[0].set_title('(a) Temp summertime variance [K$^2$]',fontsize=16)
with open('HeatWave_SH_5_variance.pkl','rb') as f:
    F=pickle.load(f)
h=F['VAR']
axs[0].contourf(longitude,F['latitude'],np.fft.fftshift(h,axes=1),np.arange(21),transform=ccrs.PlateCarree(),extend='both',cmap=plt.get_cmap('Reds'))

with open('ColdSpell_NH_5_variance.pkl','rb') as f:
    F=pickle.load(f)
h=F['VAR']
PO=axs[1].contourf(longitude,F['latitude'],np.fft.fftshift(h,axes=1),np.arange(91),transform=ccrs.PlateCarree(),extend='both',cmap=plt.get_cmap('Blues'))
plt.colorbar(PO, ax=axs[1], orientation='vertical', fraction=0.05, pad=0.05, aspect=30, shrink=0.8)
axs[1].set_title('(b) Temp wntertime variance [K$^2$]',fontsize=16)
with open('ColdSpell_SH_5_variance.pkl','rb') as f:
    F=pickle.load(f)
h=F['VAR']
axs[1].contourf(longitude,F['latitude'],np.fft.fftshift(h,axes=1),np.arange(0,91,21),transform=ccrs.PlateCarree(),extend='both',cmap=plt.get_cmap('Blues'))
# plt.subplot_tool()
plt.subplots_adjust(left=0.033,bottom=0.052,right=0.952,top=0.955)
# plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.savefig('Figure_S2.png',dpi=300)
plt.close()         
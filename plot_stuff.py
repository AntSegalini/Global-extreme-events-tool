import numpy as np
import datetime as dt
import pickle
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from scipy.stats import linregress
from scipy.ndimage import gaussian_filter
from methods import Area_rectangle_unit

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

# this part creates Figures 2
with open('ColdSpell_NH_5.pkl','rb') as f:
    F=pickle.load(f)
percentile=5
fig,axs=plot_map(3,1) # plotting only the cold spells
    
latind=np.where((latitudeb>=30) & (latitudeb<=70))[0]
latitude=latitudeb[latind]

A=np.zeros((len(latitude),len(longitude)))
Nevents=np.zeros_like(A)
T=np.zeros_like(A)
S0=[]

for Nj in range(len(F['time_events'])):
    NN=len(F['LATLON'][Nj])

    S=[]
    if NN>0: #there are heatwaves
        j=0
        for PP in F['LATLON'][Nj]:
            S.append(index_latlon_NH(PP[:,0],PP[:,1]))
            A[S[-1]]+=1; T[S[-1]]+=F['Tmean'][Nj][j]   # concurrent events    
            j+=1      
        
        S=[tuple(row) for row in np.hstack(S).T] # S is a list of tuples (lat,lon) with the indices of the events

        S0=np.array(list(set(S)-set(S0)))
        if len(S0)>0: Nevents[S0[:,0],S0[:,1]]+=1 
    S0=S 

# relative frequency of events
A[A==0]=np.nan
PO=axs[0].contourf(longitude,latitude,np.fft.fftshift(A/(0.01*percentile*len(F['time_original'])),axes=1),np.linspace(0,0.7,8),transform=ccrs.PlateCarree(),cmap=plt.get_cmap('Blues'),extend='both')
plt.colorbar(PO, ax=axs[0], orientation='vertical', fraction=0.05, pad=0.05, aspect=30, shrink=0.8)
axs[0].set_title('(a) Relative frequency cold spells',fontsize=16)

A[Nevents==0]=np.nan
Nevents[A==np.nan]=np.nan   
# duration
PO=axs[1].contourf(longitude,latitude,np.fft.fftshift(A/Nevents,axes=1),np.linspace(2,7,6),transform=ccrs.PlateCarree(),cmap=plt.get_cmap('Blues'),extend='both')
plt.colorbar(PO, ax=axs[1], orientation='vertical', fraction=0.05, pad=0.05, aspect=30, shrink=0.8)
axs[1].set_title('(b) Duration cold spells [days]',fontsize=16)

# temperature anomaly
PO=axs[2].contourf(longitude,latitude,np.fft.fftshift(T/A,axes=1),np.linspace(-20,-5,16),transform=ccrs.PlateCarree(),cmap=plt.get_cmap('Blues'))
plt.colorbar(PO, ax=axs[2], orientation='vertical', fraction=0.05, pad=0.05, aspect=30, shrink=0.8)
axs[2].set_title('(f) Temperature anomaly cold spells [K]',fontsize=16)
  
plt.tight_layout()
plt.show()
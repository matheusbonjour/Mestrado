import matplotlib.pyplot as plt
import numpy as np
import math as mat
import scipy.io
import xarray as xr 
import sys, glob
import matplotlib
import cartopy, cartopy.crs as ccrs 

from celluloid import Camera


lista_ccmp = glob.glob("CCMP*.nc")
lista_oisst = glob.glob("oisst*.nc")


#ds_vento = ds_vento.mean(dim="time")
avgfile_wind = xr.open_mfdataset(lista_ccmp,  concat_dim='time',combine='nested')
avgfile_sst = xr.open_mfdataset(lista_oisst,  concat_dim='time',combine='nested')
#   avgfile_sst = avgfile_sst['sst'].sel(lat=slice(-45,0),lon=slice(300,330)).squeeze()
vavg_u_wind = np.array( [:])
vavg_v_wind = np.array(avgfile_wind.resample(time='1D').mean('time').variables['vwnd'][:])
vavg_sst=np.array(avgfile_sst.resample(time='1D').mean('time').variables['sst'][:])


# camera = Camera(fig)
for dado_ccmp, dado_oisst in zip(lista_ccmp,lista_oisst):
    ds_vento = xr.open_dataset('C:/Users/matheus/Desktop/labsin/mono/{}'.format(dado_ccmp))
    ds_tsm = xr.open_dataset('C:/Users/matheus/Desktop/labsin/mono/{}'.format(dado_oisst))
    #ds_vento = xr.open_dataset('C:/Users/matheus/Desktop/labsin/mono/CCMP_Wind_Analysis_20180103_V02.0_L3.0_RSS.nc',decode_times=False)
    data2 = dado_ccmp[19:27]
    #plt.rcParams["figure.autolayout"] = True

    ds2_tsm = ds_tsm['sst'].sel(zlev = 0.0, lat=slice(-45, 0), lon=slice(300, 330)).squeeze() 
    lat = ds_tsm['lat'].sel(lat=slice(-45,0))
    lon = ds_tsm['lon'].sel(lon=slice(300,330))
    extent = [-60, -45, -30, 0]
        
    clim_sst = []
    clim_sst.append(ds2_tsm)

    lat = ds_vento['latitude'].sel(latitude=slice(-45,0))
    lon = ds_vento['longitude'].sel(longitude=slice(300,330))
    ds_vento = ds_vento.mean(dim="time")
    u_vento = ds_vento['uwnd'].sel(latitude=slice(-45,0), longitude=slice(300,330)).squeeze()
    v_vento = ds_vento['vwnd'].sel(latitude=slice(-45,0), longitude=slice(300,330)).squeeze()
    extent = [-60, -45, -30, 0]
    #sys.exit()
    ws = np.sqrt(u_vento**2 + v_vento**2) # Calcular vel do vento 
    sys.exit()
    
    clim_u_vento = []
    clim_u_vento.append(u_vento)
    clim_v_vento = []
    clim_v_vento.append(v_vento)
    clim_ws =[]
    clim_ws.append(ws)
            
  


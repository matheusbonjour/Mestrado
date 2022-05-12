import matplotlib.pyplot as plt
import numpy as np
import math as mat
import scipy.io
import xarray as xr 
import cartopy, cartopy.crs as ccrs 
import matplotlib

import sys, glob

lista_ccmp = glob.glob("CCMP*.nc")
lista_oisst = glob.glob("oisst*.nc")


for dado_ccmp, dado_oisst in zip(lista_ccmp,lista_oisst):
    ds_vento = xr.open_dataset('C:/Users/matheus/Desktop/labsin/mono/{}'.format(dado_ccmp),decode_times=False)
    ds_tsm = xr.open_dataset('C:/Users/matheus/Desktop/labsin/mono/{}'.format(dado_oisst),decode_times=False)
    #ds_vento = xr.open_dataset('C:/Users/matheus/Desktop/labsin/mono/CCMP_Wind_Analysis_20180103_V02.0_L3.0_RSS.nc',decode_times=False)
    data2 = dado_ccmp[19:27]
    
    ds2_tsm = ds_tsm['sst'].sel(zlev = 0.0, lat=slice(-45, 0), lon=slice(300, 330)).squeeze() 
    lat = ds_tsm['lat'].sel(lat=slice(-45,0))
    lon = ds_tsm['lon'].sel(lon=slice(300,330))
    extent = [-60, -45, -30, 0]
    
    # GERAR FIG
    plt.figure(figsize=(2,2))
    
    ax1 = plt.subplot(2,2,1,projection=ccrs.PlateCarree())
    
    #img_extent = [extent[0], extent[2], extent[1], extent[3]] # [min. lon, max. lon, min. lat, max. lat]
    ax1.coastlines(resolution='50m', color='black', linewidth=0.8)
    ax1.add_feature(cartopy.feature.BORDERS, edgecolor='black', linewidth=0.5)
    gl = ax1.gridlines(crs=ccrs.PlateCarree(), color='gray', alpha=1.0, linestyle='--', linewidth=0.25, xlocs=np.arange(-180, 180, 5), ylocs=np.arange(-90, 90, 5), draw_labels=True)
    gl.top_labels = False
    gl.right_labels = False 
        # Define de contour interval
    data_min = ds2_tsm.min()
    data_max = ds2_tsm.max()
    interval = 1
    levels = np.arange(data_min,data_max,interval)

    # Plot the contours
    img1 = ax1.contourf(lon, lat, ds2_tsm, cmap='jet', levels=levels, extend='both')    
    img2 = ax1.contour(lon, lat, ds2_tsm, colors='white', linewidths=0.3, levels=levels)
    ax1.clabel(img2, inline=2, inline_spacing=0, fontsize='10',fmt = '%1.0f', colors= 'blue')
    plt.colorbar(img1, label='Temperatura da Superfície do Mar(°C)', orientation='vertical', pad=0.01, fraction=0.05)
    ax1.set_title('OISST' + data2 , fontweight='bold', fontsize=10, loc='left')
    #plt.savefig('tsm_oisst.png')
    
    


    lat = ds_vento['latitude'].sel(latitude=slice(-45,0))
    lon = ds_vento['longitude'].sel(longitude=slice(300,330))
    ds_vento = ds_vento.mean(dim="time")
    u_vento = ds_vento['uwnd'].sel(latitude=slice(-45,0), longitude=slice(300,330)).squeeze()
    v_vento = ds_vento['vwnd'].sel(latitude=slice(-45,0), longitude=slice(300,330)).squeeze()
    extent = [-60, -45, -30, 0]
    #sys.exit()
    ws = np.sqrt(u_vento**2 + v_vento**2) # Calcular vel do vento 
    
    ax2 = plt.subplot(2,2,2,projection=ccrs.PlateCarree())
    
    
    #img_extent = [extent[0], extent[2], extent[1], extent[3]]
    ax2.set_extent([extent[0], extent[2], extent[1], extent[3]], ccrs.PlateCarree())
    ax2.coastlines(resolution='10m', color='black', linewidth=0.8)
    ax2.add_feature(cartopy.feature.BORDERS, edgecolor='black', linewidth=0.5)
    g2 = ax2.gridlines(crs=ccrs.PlateCarree(), color='gray', alpha=1.0, linestyle='--', linewidth=0.25, xlocs=np.arange(-180, 180, 5), ylocs=np.arange(-90, 90, 5), draw_labels=True)
    g2.top_labels = False
    g2.right_labels = False
    data_min = 0
    data_max = 20
    interval = 2
    levels = np.arange(data_min,data_max,interval) # Define intervalo do contorno 
    # Paleta de cores customizada # 
    colors = ["#e7f2f4", "#ceeaee", "#b6e2e8", "#abdcff", "#a4d685", "#9cd04e", 
              "#abcf2a", "#c9d21b", "#e8d50c", "#ffd100", "#ffba00", "#ffa200"]
    cmap = matplotlib.colors.ListedColormap(colors)
    cmap.set_over('#ff8c00')
    cmap.set_under('#fffafa')
    # Plotando contornos 
    img4 = ax2.contourf(lon, lat, ws, cmap=cmap, levels=levels, extend='both')    
    img5 = ax2.contour(lon, lat, ws, colors='white', linewidths=0.3, levels=levels)
    ax2.clabel(img5, inline=1, inline_spacing=0, fontsize='10',fmt = '%1.0f', colors= 'black')
    # Plotando vetores de vento 
    img6 = ax2.quiver(lon[::3], lat[::3], u_vento[::3,::3], v_vento[::3,::3])
    qk = ax2.quiverkey(img6, 0.85, 0.89, 15, '15 m/s', labelpos='E', coordinates='figure')
    plt.colorbar(img4, label='Isotacas (m/s)', orientation='vertical', pad=0.05, fraction=0.05)
    ax2.set_title('CCMP: Direção e intensidade do vento' , fontweight='bold', fontsize=10, loc='left')
    ax2.set_title('Data: ' + data2 , fontsize=10, loc='right')
    plt.show()
    sys.exit()






#fig, axs = plt.subplots(1,2)
#axs[0,0].contourf(ds2_tsm)
#axs[0,0].set_title('TSM',fontsize=16)
#axs[0,0].colorbar()
#axs[0,0].axis('equal')
#axs[0,1].quiver(u_vento,v_vento)
# axs[0,1].set_title('TSM',fontsize=16)
# axs[0,1].colorbar()
# axs[0,1].axis('equal')




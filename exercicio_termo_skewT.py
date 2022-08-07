# Importando m´odulos necess´arios
import matplotlib.pyplot as plt
import matplotlib.axis as maxis
from matplotlib.ticker import (MultipleLocator, NullFormatter, ScalarFormatter)
import pandas as pd
import numpy as np
import metpy as mt
from metpy.plots import SkewT
from metpy.units import pandas_dataframe_to_unit_arrays, units
#from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from metpy.plots import Hodograph, SkewT
import metpy.calc as mpcalc
from metpy.units import units
import math
class SkewXAxis(maxis.XAxis):
def _get_tick(self, major):
return SkewXTick(self.axes, None, major=major)

# Definindo constantes
g=9.81
cp=1004
Rd=287
Lv=2.501*10**6
eps=0.622
gamad=g/cp

# Lendo dados
gamad=g/cp
df=pd.read_csv(’/matheus/sondagem.csv’,sep=’;’)
p=df[’PRES’].values
t=df[’TEMP’].values
td=df[’DWPT’].values
rv=df[’MIXR’].values
tetak=df[’THTA’].values
tetaek=df[’THTE’].values
z=df["HGHT"].values
UR=df[’RELH’].values
tk=df[’TEMP’].values+273.15
tdk=df[’DWPT’].values+273.15
i500 = df[df[’PRES’] == 500.0].index[0]
i850 = df[df[’PRES’] == 850.0].index[0]

e=[]
es=[]
tvk=[]
tv=[]
q=[]
rs=[]
rvc=[]
rvg=[]
rsg=[]
ro=[]
teta=[]
tetae=[]
tetakc=[]
tetaekc=[]
tetac=[]
tetaec=[]

# Calculando todos os niveis:
for i in range(len(p)):
    e.append(6.11*math.exp(17.67*td[i]/(td[i]+243.5)))
    es.append(6.11*math.exp(17.67*t[i]/(t[i]+243.5)))
    q.append(rv[i]/1000/(1+rv[i]/1000))
    tvk.append(tk[i]*(1+0.61*q[i]))
    tv.append(tvk[i]-273.15)
    rs.append(eps*(es[i]/(p[i]-es[i])))
    rvc.append(eps*(e[i]/(p[i]-e[i])))
    rvg.append(rvc[i]*1000)
    rsg.append(rs[i]*1000)
    ro.append(p[i]*100/Rd/(tvk[i]))
    # ro.append((p[i]-e[i])*100/Rd/(tvk[i]))
    teta.append(tetak[i]-273.15)
    tetae.append(tetaek[i]-273.15)
    tetakc.append((tk[i])*(1000/p[i])**(Rd/cp))
    tetaekc.append(tetakc[i]*math.exp(Lv*rs[i]/cp/(tk[i])))
    tetac.append(tetakc[i]-273.15)
    tetaec.append(tetaekc[i]-273.15)
    

# plotando o skew T
fig = plt.figure(figsize=(9, 11))
skew = SkewT(fig, rotation=45)
skew.plot(p, t, ’red’)
skew.plot(p, td, ’blue’)
# definindo os limites e nomes dos eixos
skew.ax.set_xlim(-50, 70)
skew.ax.set_ylim(1050, 50)
plt.xlabel(’Temperatura (°C)’, fontsize=14)
plt.ylabel(’Press˜ao (hPa)’, fontsize=14)
plt.title(’Perfil vertical da temperatura (°C)\nSBMT C. Marte 05/03/22 12Z’,
fontsize=15, ha=’center’)
plt.legend([’T’, ’Td’])
# adicionando as linhas de adiabaticas e raz˜ao de mistura
skew.plot_dry_adiabats(t0=np.arange(233, 533, 10) * units.K,
alpha=0.25, color=’orangered’)
skew.plot_moist_adiabats(t0=np.arange(233, 400, 5) * units.K,
alpha=0.25, color=’tab:green’)
skew.plot_mixing_lines(pressure=np.arange(1000, 20, -5) * units.hPa,
linestyle=’dotted’, color=’tab:blue’)
#skew.shade_cin(p, t, tp, td)
#skew.shade_cape(p, t, tp, **kwargs)
r=rvc[0]/rs[0]
#Temperatura no NCL
tncl=1/(1/(t[0]+273.15-55)-(math.log(r)/2840))+55
tnclk=tncl+273.15
kd=0.286
kp=kd*(1-0.26*q[0])
pncl=p[0]*(tncl/(t[0]+273.15))**(1/kp) #Tsonis pag 137
tetancl=tncl*(1000/pncl)**(Rd/cp)
zncl=z[0]+Rd/g*(tvk[0]+tncl)/2*math.log(p[0]/pncl)

gamasp=[]
gamadc=[]
gamaspkm=[]
a=Lv/Rd
c=eps*(Lv**2)
d=cp*Rd
for i in range(len(z)):
    b=rs[i]/tvk[i]
    e=rs[i]/(tvk[i]**2)
    gamasp.append(gamad*((1+a*b)/(1+c/d*e)))
    gamaspkm.append(gamasp[i]*1000)
    gamadc.append(9.8*z[i]/1000)
    
    
tpk=[tk[0]]
tp=[t[0]]
tp=[]
for i in range(len(z)-1):
    #tpk.append(tpk[i]-gamasp[i]*(z[i+1]-z[i]))
    if z[i]<zncl:
    tpk.append(tpk[i]-gamad*(z[i+1]-z[i]))
    #tpk.append(tpk[i]-gamasp[i]*(z[i+1]-z[i]))
    else:
    tpk.append(tpk[i]-gamasp[i]*(z[i+1]-z[i]))
    
    
for i in range(len(tpk)):
    tp.append(round(tpk[i]-273.15,5))
    

# plotando o skew T
fig = plt.figure(figsize=(9, 11))
skew = SkewT(fig, rotation=45)
skew.plot(p, t, ’red’)
skew.plot(p, tp, ’blue’)
skew.plot(p, td, ’black’)
# definindo os limites e nomes dos eixos
skew.ax.set_xlim(-50, 70)
skew.ax.set_ylim(1050, 50)
plt.xlabel(’Temperatura (°C)’, fontsize=14)
plt.ylabel(’Press˜ao (hPa)’, fontsize=14)
plt.title(’Perfil vertical da temperatura (°C)’,
fontsize=15, ha=’center’)
plt.legend([’T’, ’$T_p$’,’$T_D$’])
# adicionando as linhas de adiabaticas e raz˜ao de mistura
skew.plot_dry_adiabats(t0=np.arange(233, 533, 10) * units.K,
alpha=0.25, color=’orangered’)
skew.plot_moist_adiabats(t0=np.arange(233, 400, 5) * units.K,
alpha=0.25, color=’tab:green’)
skew.plot_mixing_lines(pressure=np.arange(1000, 20, -5) * units.hPa,
linestyle=’dotted’, color=’tab:blue’)
skew.shade_cin(np.array(p[0:45]), np.array(t[0:45]), np.array(tp[0:45]))
skew.shade_cape(p, t, tp)
r=rvc[i850]/rs[i850]
#Temperatura no NCL
tncl=1/(1/(t[i850]+273.15-55)-(math.log(r)/2840))+55
tnclk=tncl+273.15
kd=0.286
kp=kd*(1-0.26*q[i850])
pncl=p[i850]*(tncl/(t[i850]+273.15))**(1/kp) #Tsonis pag 137
tetancl=tncl*(1000/pncl)**(Rd/cp)
zncl=z[i850]+Rd/g*(tvk[i850]+tncl)/2*math.log(p[i850]/pncl)
pncl,tetancl,zncl
gamasp=[]
gamadc=[]
gamaspkm=[]
a=Lv/Rd
c=eps*(Lv**2)
d=cp*Rd
for i in range(len(z)):
    b=rs[i]/tvk[i]
    e=rs[i]/(tvk[i]**2)
    gamasp.append(gamad*((1+a*b)/(1+c/d*e)))
    gamaspkm.append(gamasp[i]*1000)
    gamadc.append(9.8*z[i]/1000)
    
   
tpk=[tk[i850]]
tp=[t[i850]]
tp=[]
for i in range(i850, len(z)-1):
#tpk.append(tpk[i]-gamasp[i]*(z[i+1]-z[i]))
    if z[i]<zncl:
        pk.append(tpk[i-i850]-gamad*(z[i+1]-z[i]))
#tpk.append(tpk[i]-gamasp[i]*(z[i+1]-z[i]))
    else:
        tpk.append(tpk[i-i850]-gamasp[i]*(z[i+1]-z[i]))
        
for i in range(len(tpk)):
    tp.append(round(tpk[i]-273.15,5))
ps=p
ts=t
zs=z
tl = np.array(tp)
tl = tl[i500]
t5=t[i500]
# Showalter
S=t5-tl


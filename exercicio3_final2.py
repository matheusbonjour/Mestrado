import numpy as np
import math
import matplotlib.pyplot as plt 
from celluloid import Camera
import matplotlib.colors as colors
from matplotlib import gridspec
import pylab as pl
from scipy.sparse.linalg import spsolve
from scipy.sparse import spdiags
import matplotlib.colors as colors
import time
from matplotlib import cm
import matplotlib.pyplot as plt 
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LightSource
import cmocean.cm as cmo
from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation
import metpy.calc
from metpy.units import units
import sys 
import cmocean.cm as cm





 # Salvar tempo inicial    
result_u = []
result_v = []
result_h = []
# Para conferir a conservação da massa e energia 
Vol = [] # volume
vort = [] # volume
Ep = []  # energia potencial
Ek = []  # energia cinética
conser = []
# dado de um ponto da grade para ser salvo ao longo da integração 
histgrid1 = [] 
histgrid2 = []  
histgrid3 = []   
histgrid4 = []    



# -------------------------------
# Parametros da Grade 
nx = 100 #  numero de pontos da grade
dx = 40000 #   espaçamento da grade (m)
dy = dx
x_axis = np.arange(0,(nx*dx/1000),dx/1000) # para plot
dt = 40 # tamanho do passo de tempo (s)

delta = 100e3 #  grid spacing (m) 
xmax,xmin = 4e6,-4e6
ymax,ymin = 4e6,-4e6
nx = int((xmax+delta-xmin)/delta) # (numero de pontos da grade em x)
ny = int((ymax+delta-ymin)/delta) # (numero de pontos da grade em y)
# Grade para u
xu = np.arange(-40, 42,1)*delta
yu = np.arange(40, -41,-1)*delta
lonsu,latsu = np.meshgrid(xu,yu) 
# Grade para csi 
xz = np.arange(-40, 42.,1.)*delta
yz = (np.arange(41, -41,-1)-.5)*delta
lonsz,latsz = np.meshgrid(xz,yz) 
# Grade para v   
xv = (np.arange(-40, 41.,1.)+.5)*delta
yv = (np.arange(41, -41,-1)-.5)*delta
lonsv,latsv = np.meshgrid(xv,yv)

lons_plotnl = (lonsz[1:,1:] + lonsz[:-1,:-1])/2
lats_plotnl = (latsz[1:,1:] + latsz[:-1,:-1])/2

#-------------------------
# Constants
H=250
g = 9.8 # gravidade(m/s-2)
mu = 0.2
dt = (mu*delta/(np.sqrt(g*H)))/2
tmax= 1501
crls=1
forcante = 2
# ----------------------------
if crls == 1:
    # Calculando o parâmetro de Coriolis utilizando aproximação de plano Beta 
    beta = (2*7.29*1e-5*np.cos(0*np.pi/180)/6400000)
    
if crls == 0:
    # Desliga Coriolis 
    beta = 0
# Coriolis force exerced by each wind component
# at v
# yv = np.arange(ymax+(delta/2), ymin-1.5*delta,-delta)
# xv = np.arange(xmin, xmax+1*delta,delta)
# lonsv,latsv = np.meshgrid(xv,yv) 
fv = latsv*beta
# # at u
# yu = np.arange(ymax, ymin-1*delta,-delta)
# xu = np.arange(xmin, xmax+2*delta,delta)
# lonsu,latsu = np.meshgrid(xu,yu) 
fu = latsu*beta
fz = latsz*beta
a =0.25
cx = 0
cy = -2e6 # -2e6
nrx = 4
nry = 4
dx = delta 
xgauss = np.linspace(xmin, xmax,nx)
ygauss = np.linspace(ymax, ymin,ny) #centrado
#ygauss = np.linspace(7e6, ymin,ny)    
########################################################

longa,latga = np.meshgrid(xgauss,ygauss) 

gauss_space = 1* a * np.exp(- (((longa-cx)**2)/(nrx*dx)**2) - (((latga-cy)**2)/(nry*dx)**2))
Nt = tmax 
u = np.zeros([Nt,ny,nx+1])
v = np.zeros([Nt,ny+1,nx])
h = np.zeros([Nt,ny,nx])

frecplot = 50
alpha = 0.02
t = 0
salvafig = 0 
propriedades = 1
salvaframes = 0
for t in range(tmax):

    if forcante == 2:
        
        # A pertubação é a gaussiana no espaço * seu avanço temporal             
        decay = (H/2)*(alpha**3)*(t**2)*math.exp(-alpha*t)
        # gauss_time = ((500/2)*alpha**3)*(t/3600)*np.exp(-alpha*t/3600)
        frc = (gauss_space * decay)
        
    if t == 0:        
        # Grades iniciais 
        un = np.zeros((ny,nx+1))  # (81,82)
        vn =  np.zeros((ny+1,nx)) # (82,81)
        hn =  np.zeros((ny,nx))   # (81,81)

       
        
        # Atualiza matrizes do tempo 
        unm1,vnm1,hnm1 = un*np.nan,vn*np.nan,hn*np.nan
        unp1,vnp1,hnp1 = un*np.nan,vn*np.nan,hn*np.nan
        u[t,:,:] = un; v[t,:,:] = vn; h[t,:,:]=hn
        
        # Condição Inicial (Euler Avançado no tempo) 
    elif t == 1:
        
        # --------------------------------------------
        # Efeito de rotação para os campos de velocidade
        # Coriolis v 
        corv = 0.25  * (fv[:-1,:-1] * (vn[:-1,1:] + vn[:-1,:-1]) + fv[1:,:-1] * (vn[1:,:-1] + vn[1:,1:]))
        # Coriolis u 
        coru = 0.25 * (fu[1:,:-1] * (un[1:,:-1] + un[1:,1:]) + fu[:-1,:-1] * (un[:-1,1:] + un[:-1,:-1])) 
        
        # --------------------------------------------
        # Usando Euler Avançado no tempo para t = 1:        
        # Prognóstico do campo de u 
        unp1[:,1:-1] = un[:,1:-1] + dt * ((-g * ((hn[:,1:] - hn[:,:-1]) / delta)) + corv) 
        # Prognóstico do campo de v
        vnp1[1:-1] = vn[1:-1] + dt * ((-g * ((hn[:-1] - hn[1:]) / delta)) - coru) 
        # Prognóstico do campo de h 
        hnp1[:] = hn[:] - (dt * H / delta) * (un[:,1:] - un[:,:-1] + (vn[:-1] - vn[1:]))
        hnp1[:,:] =  hnp1[:,:] + frc[:,:] # Adicionando Forçante! 

             
        # Restante da integração (t > 1) (Leapfrog) 
    else:  
        # --------------------------------------------
        # Efeito de rotação para os campos de velocidade
        # Coriolis v            
        corv = 0.25  * (fv[:-1,:-1] * (vn[:-1,1:] + vn[:-1,:-1]) + fv[1:,:-1] * (vn[1:,:-1] + vn[1:,1:]))
        # Coriolis u 
        coru = 0.25 * (fu[1:,:-1] * (un[1:,:-1] + un[1:,1:]) + fu[:-1,:-1] * (un[:-1,1:] + un[:-1,:-1]))
        
        # --------------------------------------------
        # Usando Leap-Frog t > 1:        
        # Prognóstico do campo de u 
        unp1[:,1:-1] =  unm1[:,1:-1] + 2*dt * ((-g * ((hn[:,1:] - hn[:,:-1]) / delta)) + corv) 
        # Prognóstico do campo de u 
        vnp1[1:-1] = vnm1[1:-1] + 2*dt * ((-g * ((hn[:-1] - hn[1:]) / delta)) - coru) 
        # Prognóstico do campo de h 
        hnp1 = hnm1 - (2*dt*H/delta) * (un[:,1:] - un[:,:-1] + (vn[:-1] - vn[1:]))       
        hnp1[:,:] =  hnp1[:,:] + frc[:,:] # Adicionando Forçante! 

        # Condições de Contorno Radiacionais 
    if t != 0:
    
        # --------------------------------------------
        # Contorno de u: Radiacional à esquerda (Oeste) 
        for i in range(len(unp1[:,0])):
            c = np.sqrt(g*(hn[i,1]+H))
            # if c > 0 or np.isnan(c):
                # c = np.sqrt(g*H)
            # elif c > dx/dt:
                # c = dx/dt
            unp1[i,0] = un[i,0] + dt * ((c*(un[i,1]-un[i,0])/delta) + (0.5 * (fv[i,0]*vn[i,0] + fv[i+1,0]*vn[i+1,0])))
        # --------------------------------------------
        # Contorno de u: Radiacional à direita (Leste)

        # unp1[:,-1] = 0
        for i in range(len(unp1[:,-1])):
            c = np.sqrt(g*(hn[i,-2]+H))
            # if c < 0 or np.isnan(c):
                # c = np.sqrt(g*H)
            # elif c > dx/dt:
                # c = dx/dt
            unp1[i,-1] = un[i,-1] + dt * ((-c*(un[i,-1]-un[i,-2])/delta) + (0.5 * (fv[i,-1] * vn[i,-1] + fv[i+1,-1] * vn[i+1,-1])))
        # --------------------------------------------
        # Contorno de v: Radiacional parte superior (Norte)
        for i in range(len(vnp1[0])):
            c = np.sqrt(g*(hn[1,i]+H))
            # if c < 0 or np.isnan(c):
                # c = np.sqrt(g*H)
            # elif c > dx/dt:
                # c = dx/dt                        
            vnp1[0,i] = vn[0,i] + dt * ((-c*(vn[0,i]-vn[1,i])/delta))
        # --------------------------------------------
        # Contorno de v: Radiacional parte inferior (Sul)
        for i in range(len(vnp1[-1])):
            c = np.sqrt(g*(hn[-2,i]+H))
            # if c > 0 or np.isnan(c):
                # c = np.sqrt(g*H)
            # elif c > dx/dt:
                # c = dx/dt                        
            vnp1[-1,i] = vn[-1,i] + dt * ((c*(vn[-2,i]-vn[-1,i])/delta) - (0.5 * (fu[-1,i] * un[-1,i] + fu[-1,i+1] * un[-1,i+1])))  
            

        
    # Atualizando as matrizes no tempo         
        unm1,vnm1,hnm1 = un,vn,hn         
        un,vn,hn = unp1,vnp1,hnp1 
        unp1,vnp1,hnp1 = un*np.nan,vn*np.nan,hn*np.nan
        u[t,:,:] = un; v[t,:,:] = vn; h[t,:,:]=hn

    #salvar resultados desses tempos nas matrizes
    
    # Calculo de vorticidade relativa e absoluta 
    vort_rel = ((-vn[:-1]+vn[1:])/delta) - ((-un[:,:-1]+un[:,1:])/delta) 
    vort_abs = (vort_rel + beta)/(hn+H)
    tmp = (vort_abs)*(delta**2)/2
    # Salva na vorticidade 
    vort.append((tmp*(hn+H)).sum())
    
    # Volume total        
    Vol.append(((hn+H)*delta**2).sum())
    
    # Energia Potencial 
    Ep.append(((hn+H)**2*delta**2).sum()*g/2)


    # Energia Cinética 
    EKu = (un[:,:-1]+un[:,1:])/2
    EKv = (vn[:-1]+vn[1:])/2
    Ek.append((((EKu**2 + EKv**2)*delta*delta).sum())*H/2)
    ç=10e10
    
    # ----------------------------------
    # Salva os pontos centrais das grades 
    histgrid1.append(hn[40,40])
    histgrid2.append(un[40,41])
    histgrid3.append(vn[40,40])
    histgrid4.append(vn[41,40])
    # --------------------------
    # Salvar os resultados de acordo com a frequência de plot 
    if t % frecplot == 0:        
        result_u.append(un)
        result_v.append(vn)
        result_h.append(hn)
        
    # Plota e salva o momento inicial e os momentos em que t%x == 0:
    # | x = frequencia de plot 
    if salvafig == 1:
        if t == 1 or t%40 ==0:  
            # -----------------------------------
            # Interpolação de un e vn para plot 
            uplot = (un[:,1:]+un[:,0:-1])/2
            vplot = (vn[1:,:]+vn[0:-1,:])/2
            # -----------------------------------
            # Inicia o plot
            plt.clf()
            # Define subplots, tamanho e resolução da figura 
            fig, ax1 =plt.subplots(figsize=(9, 8),dpi=300)
            # ---------------------------------------------------------------
            # Visualização de H 
            # Intervalo de valores de h 
            hmax = 3
            hmin = -hmax 
            # ----------------------------------------------------------------
            # Plota contorno preenchido e isolinhas de elevação da superfície (h) 
            cf=plt.contourf(lons_plot*1e-3,lats_plot*1e-3,hn,np.linspace(hmin,hmax,51),cmap=cm.balance,extend='both')
            plt.contour(lons_plot*1e-3,lats_plot*1e-3,hn,np.linspace(hmin,hmax,30),linewidths=1,colors='k',extend='both')
            # ===TESTE===TESTE===TESTE===TESTE===TESTE                 
            # cf=plt.contourf(lons_plot*1e-3,lats_plot*1e-3,hn-H,np.linspace(-0.06,0.06,51),cmap=cm.balance,extend='both')
            # plt.contour(lons_plot*1e-3,lats_plot*1e-3,hn-H,np.linspace(-0.06,0.06,30),linewidths=1,colors='k',extend='both')
            # ===TESTE===TESTE===TESTE===TESTE===TESTE 
            # ------------------------------------------------------------------
            # Visualização do campo de velocidade 
            # Definir parâmetros para os vetores (intervalo, largura e escala)
            step=3
            width=0.002
            scale=2.5
            # Plota o campo de velocidade 
            qv = ax1.quiver(lons_plot[::step,::step]*1e-3,lats_plot[::step,::step]*1e-3,uplot[::step,::step],vplot[::step,::step],scale=scale,width=width)
            # Insere vetor de escala no local desejado 
            ax1.quiverkey(qv, 0.85, 0.87,0.1,r'$0.1 \frac{m}{s}$', labelpos='E', coordinates='figure')
            # -----------------------------------------------------------------
            # Configurações de colorbar, labels e legendas
            cbar=plt.colorbar(cf)
            plt.xlabel('x [km]',weight='bold',fontsize=12)
            plt.ylabel('y [km]',weight='bold',fontsize=12)
            cbar.set_label('h [m]',weight='bold',fontsize=12)
            plt.savefig('a4e6'+str(10000+t)+'.png',dpi=300,bbox_inches='tight')
            # -----------------------------------------------------------------
        
         
        
# hmax = np.abs([h.min(),h.max()]).max()
# hmin= -hmax


Vol =  np.array(Vol) #- Vol[0]
Enstrophy =  np.array(vort)
Enstrophy = Enstrophy*Enstrophy/ç

Ep= np.array(Ep) #- Ep[0]
Ek= np.array(Ek)

conser= Ek+Ep

tx = np.arange(0,(len(histgrid1))*dt/60/60/24,dt/60/60/24)

 

fig,axs = plt.subplots(5,figsize=(10, 12), constrained_layout=True)

axs[0].plot(tx[1:],Vol,linewidth=1,color='k')
# axs[0].set_ylim(1.6400*1e16,1.6600*1e16)
axs[0].set_title('Volume',color='k', fontsize = 22)

# axs[0].plot(tx, zorooo,linewidth=5,color='k')
axs[0].tick_params(axis='both', which='major', labelsize=16)    
axs[1].plot(tx[1:],Ep,linewidth=1,color='k')
# axs[1].set_ylim(2.000*1e19,2.0600*1e19)
axs[1].set_title('Energia Potencial',color='k', fontsize = 22)
axs[1].tick_params(axis='both', which='major', labelsize=16)
axs[4].set_xlabel('Tempo (dias)', fontsize=18)       # plot energy and mass
axs[1].grid() 
axs[0].grid() 

axs[2].plot(tx[1:],Ek,linewidth=1,color='k')
# axs[2].set_ylim(0,2*1e16)
axs[2].set_title('Energia Cinética',color='k', fontsize = 22)      
axs[2].tick_params(axis='both', which='major', labelsize=16)    
axs[2].grid() 

axs[3].plot(tx[1:],conser,linewidth=1,color='k')
# axs[3].set_ylim(2.0000*1e19,2.0600*1e19)
axs[3].set_title('Energia Total',color='k', fontsize = 22)      
axs[3].tick_params(axis='both', which='major', labelsize=16)  
axs[3].grid() 

axs[4].plot(tx[1:],Enstrophy,linewidth=1,color='k')

axs[4].set_title('Enstrofia',color='k', fontsize = 22)      
axs[4].tick_params(axis='both', which='major', labelsize=16)  
axs[4].grid() 

zorooo = Enstrophy
zorooo = np.zeros(np.shape(Enstrophy))
zorooo= zorooo + Vol[0]
axs[0].fill_between(tx[1:], Vol, Vol[0], where=Vol >= zorooo, facecolor='0.7', interpolate=True, alpha=0.4)
axs[0].fill_between(tx[1:], Vol, Vol[0], where=Vol <= zorooo, facecolor='0.2', interpolate=True, alpha=0.4)
zorooo = np.zeros(np.shape(Enstrophy))
zorooo= zorooo + Ep[0]
axs[1].fill_between(tx[1:], Ep, Ep[0], where=Ep >= zorooo, facecolor='0.7', interpolate=True, alpha=0.4)
axs[1].fill_between(tx[1:], Ep, Ep[0], where=Ep <= zorooo, facecolor='0.2', interpolate=True, alpha=0.4)
zorooo = np.zeros(np.shape(Enstrophy))
zorooo= zorooo + Ek[0]
axs[2].fill_between(tx[1:], Ek, Ek[0], where=Ek >= zorooo, facecolor='0.7', interpolate=True, alpha=0.4)
axs[2].fill_between(tx[1:], Ek, Ek[0], where=Ek <= zorooo, facecolor='0.2', interpolate=True, alpha=0.4)
zorooo = np.zeros(np.shape(Enstrophy))
zorooo= zorooo + conser[0]
axs[3].fill_between(tx[1:], conser, conser[0], where=conser >= zorooo, facecolor='0.7', interpolate=True, alpha=0.4)
axs[3].fill_between(tx[1:], conser, conser[0], where=conser <= zorooo, facecolor='0.2', interpolate=True, alpha=0.4)
zorooo = np.zeros(np.shape(Enstrophy))
zorooo= zorooo + Enstrophy[0]
axs[4].fill_between(tx[1:], Enstrophy, Enstrophy[0], where=Enstrophy >= zorooo, facecolor='0.7', interpolate=True, alpha=0.4)
axs[4].fill_between(tx[1:], Enstrophy, Enstrophy[0], where=Enstrophy <= zorooo, facecolor='0.2', interpolate=True, alpha=0.4)

pl.savefig('Linear_Deslocado1.png')     
 

      
      


if salvaframes == 1:

    # ============ t1 ============= #
    t = 20
    uplot = (u[t,:,1:]+u[t,:,0:-1])/2
    vplot = (v[t,1:,:]+v[t,0:-1,:])/2

    # -----------------------------------
    # Inicia o plot
    plt.clf()
    # Define subplots, tamanho e resolução da figura 
    fig2, ax2 =plt.subplots(4,2,figsize=(8, 12),sharex=True,sharey=True,dpi=450)
    hmax = 3
    hmin = -hmax 
    titlestr='tempo = '+str(int(t*dt/3600)) + ' horas ($\mu$='+str(mu)+')'  
    ax2[0,0].set_title(titlestr, fontsize = 12, fontweight = 'bold')
    # ---------------------------------------------------------------
    # Visualização de H 
    # Intervalo de valores de h 

    # ----------------------------------------------------------------
    # Plota contorno preenchido e isolinhas de elevação da superfície (h) 
    cf=ax2[0,0].contourf(lons_plotnl*1e-3,lats_plotnl*1e-3,h[t],np.linspace(hmin,hmax,51),cmap=cm.balance,extend='both')
    ax2[0,0].contour(lons_plotnl*1e-3,lats_plotnl*1e-3,h[t],np.linspace(hmin,hmax,30),linewidths=0.7,colors='k',extend='both')
    
    # ------------------------------------------------------------------
    # Visualização do campo de velocidade 
    # Definir parâmetros para os vetores (intervalo, largura e escala)
    step=3
    width=0.002
    scale=2.5
    # Plota o campo de velocidade 
    qv = ax2[0,0].quiver(lons_plotnl[::step,::step]*1e-3,lats_plotnl[::step,::step]*1e-3,uplot[::step,::step],vplot[::step,::step],scale=scale,width=width)
    # Insere vetor de escala no local desejado 
    ax2[0,0].quiverkey(qv, 0.85, 0.87,0.1,r'$0.1 \frac{m}{s}$', labelpos='E', coordinates='figure')
    # -----------------------------------------------------------------
    # Configurações de colorbar, labels e legendas
    # cbar=plt.colorbar(cf)
    # ax2[0,0].set_xlabel('x [km]',weight='bold',fontsize=12)
    ax2[0,0].set_ylabel('y [km]',weight='bold',fontsize=12)
    # cbar.set_label('h [m]',weight='bold',fontsize=12)
    
    
    # ========= t2 ============================================================ #
    
    t = 120
    uplot = (u[t,:,1:]+u[t,:,0:-1])/2
    vplot = (v[t,1:,:]+v[t,0:-1,:])/2
    
    titlestr='tempo = '+str(int(t*dt/3600)) + ' horas ($\mu$='+str(mu)+')'  
    ax2[0,1].set_title(titlestr, fontsize = 12, fontweight = 'bold')
        
    # ----------------------------------------------------------------
    # Plota contorno preenchido e isolinhas de elevação da superfície (h) 
    cf=ax2[0,1].contourf(lons_plotnl*1e-3,lats_plotnl*1e-3,h[t],np.linspace(hmin,hmax,51),cmap=cm.balance,extend='both')
    ax2[0,1].contour(lons_plotnl*1e-3,lats_plotnl*1e-3,h[t],np.linspace(hmin,hmax,30),linewidths=1,colors='k',extend='both')
    
    # ------------------------------------------------------------------
    # Plota o campo de velocidade 
    qv = ax2[0,1].quiver(lons_plotnl[::step,::step]*1e-3,lats_plotnl[::step,::step]*1e-3,uplot[::step,::step],vplot[::step,::step],scale=scale,width=width)
    # Insere vetor de escala no local desejado 
    ax2[0,1].quiverkey(qv, 0.85, 0.87,0.1,r'$0.1 \frac{m}{s}$', labelpos='E', coordinates='figure')
    # -----------------------------------------------------------------
    # Configurações de colorbar, labels e legendas
    # cbar=plt.colorbar(cf)
    # ax2[0,1].set_xlabel('x [km]',weight='bold',fontsize=12)
    # ax2[0,1].set_ylabel('y [km]',weight='bold',fontsize=12)
    # cbar.set_label('h [m]',weight='bold',fontsize=12)
    # ========= t3 ============ # 
    
    t = 240
    uplot = (u[t,:,1:]+u[t,:,0:-1])/2
    vplot = (v[t,1:,:]+v[t,0:-1,:])/2
    
    titlestr='tempo = '+str(int(t*dt/3600)) + ' horas ($\mu$='+str(mu)+')'  
    ax2[1,0].set_title(titlestr, fontsize = 12, fontweight = 'bold')
        
    # ----------------------------------------------------------------
    # Plota contorno preenchido e isolinhas de elevação da superfície (h) 
    cf=ax2[1,0].contourf(lons_plotnl*1e-3,lats_plotnl*1e-3,h[t],np.linspace(hmin,hmax,51),cmap=cm.balance,extend='both')
    ax2[1,0].contour(lons_plotnl*1e-3,lats_plotnl*1e-3,h[t],np.linspace(hmin,hmax,30),linewidths=1,colors='k',extend='both')
    
    # ------------------------------------------------------------------
    # Visualização do campo de velocidade 
    # Plota o campo de velocidade 
    qv = ax2[1,0].quiver(lons_plotnl[::step,::step]*1e-3,lats_plotnl[::step,::step]*1e-3,uplot[::step,::step],vplot[::step,::step],scale=scale,width=width)
    # Insere vetor de escala no local desejado 
    ax2[1,0].quiverkey(qv, 0.85, 0.87,0.1,r'$0.1 \frac{m}{s}$', labelpos='E', coordinates='figure')
    # -----------------------------------------------------------------
    # Configurações de colorbar, labels e legendas
    # cbar=plt.colorbar(cf)
    # ax2[1,0].set_xlabel('x [km]',weight='bold',fontsize=12)
    ax2[1,0].set_ylabel('y [km]',weight='bold',fontsize=12)
    # cbar.set_label('h [m]',weight='bold',fontsize=12)
    
        # ========= t4 ============ # 
    
    t = 360
    uplot = (u[t,:,1:]+u[t,:,0:-1])/2
    vplot = (v[t,1:,:]+v[t,0:-1,:])/2
    
    titlestr='tempo = '+str(int(t*dt/3600)) + ' horas ($\mu$='+str(mu)+')'  
    ax2[1,1].set_title(titlestr, fontsize = 12, fontweight = 'bold')
        
    # ----------------------------------------------------------------
    # Plota contorno preenchido e isolinhas de elevação da superfície (h) 
    cf=ax2[1,1].contourf(lons_plotnl*1e-3,lats_plotnl*1e-3,h[t],np.linspace(hmin,hmax,51),cmap=cm.balance,extend='both')
    ax2[1,1].contour(lons_plotnl*1e-3,lats_plotnl*1e-3,h[t],np.linspace(hmin,hmax,30),linewidths=1,colors='k',extend='both')
    
    # ------------------------------------------------------------------
    # Visualização do campo de velocidade 
    # Plota o campo de velocidade 
    qv = ax2[1,1].quiver(lons_plotnl[::step,::step]*1e-3,lats_plotnl[::step,::step]*1e-3,uplot[::step,::step],vplot[::step,::step],scale=scale,width=width)
    # Insere vetor de escala no local desejado 
    ax2[1,1].quiverkey(qv, 0.85, 0.87,0.1,r'$0.1 \frac{m}{s}$', labelpos='E', coordinates='figure')
    # -----------------------------------------------------------------
    # Configurações de colorbar, labels e legendas
    # cbar=plt.colorbar(cf)
    # ax2[1,1].set_xlabel('x [km]',weight='bold',fontsize=12)
    # ax2[1,1].set_ylabel('y [km]',weight='bold',fontsize=12)
    # cbar.set_label('h [m]',weight='bold',fontsize=12)
    
    # ========= t5 ============ # 
    
    t = 480
    uplot = (u[t,:,1:]+u[t,:,0:-1])/2
    vplot = (v[t,1:,:]+v[t,0:-1,:])/2
    
    titlestr='tempo = '+str(int(t*dt/3600)) + ' horas ($\mu$='+str(mu)+')' 
    ax2[2,0].set_title(titlestr, fontsize = 12, fontweight = 'bold')
        
    # ----------------------------------------------------------------
    # Plota contorno preenchido e isolinhas de elevação da superfície (h) 
    cf=ax2[2,0].contourf(lons_plotnl*1e-3,lats_plotnl*1e-3,h[t],np.linspace(hmin,hmax,51),cmap=cm.balance,extend='both')
    ax2[2,0].contour(lons_plotnl*1e-3,lats_plotnl*1e-3,h[t],np.linspace(hmin,hmax,30),linewidths=1,colors='k',extend='both')
    
    # ------------------------------------------------------------------
    # Visualização do campo de velocidade 
    # Plota o campo de velocidade 
    qv = ax2[2,0].quiver(lons_plotnl[::step,::step]*1e-3,lats_plotnl[::step,::step]*1e-3,uplot[::step,::step],vplot[::step,::step],scale=scale,width=width)
    # Insere vetor de escala no local desejado 
    ax2[2,0].quiverkey(qv, 0.85, 0.87,0.1,r'$0.1 \frac{m}{s}$', labelpos='E', coordinates='figure')
    # -----------------------------------------------------------------
    # Configurações de colorbar, labels e legendas
    # cbar=plt.colorbar(cf)
    # ax2[2,0].set_xlabel('x [km]',weight='bold',fontsize=12)
    ax2[2,0].set_ylabel('y [km]',weight='bold',fontsize=12)
    # cbar.set_label('h [m]',weight='bold',fontsize=12)
    
    # ========= t6 ============ # 
    
    t = 600
    uplot = (u[t,:,1:]+u[t,:,0:-1])/2
    vplot = (v[t,1:,:]+v[t,0:-1,:])/2
    
    titlestr='tempo = '+str(int(t*dt/3600)) + ' horas ($\mu$='+str(mu)+')'  
    ax2[2,1].set_title(titlestr, fontsize = 12, fontweight = 'bold')
        
    # ----------------------------------------------------------------
    # Plota contorno preenchido e isolinhas de elevação da superfície (h) 
    cf=ax2[2,1].contourf(lons_plotnl*1e-3,lats_plotnl*1e-3,h[t],np.linspace(hmin,hmax,51),cmap=cm.balance,extend='both')
    ax2[2,1].contour(lons_plotnl*1e-3,lats_plotnl*1e-3,h[t],np.linspace(hmin,hmax,30),linewidths=1,colors='k',extend='both')
    
    # ------------------------------------------------------------------
    # Visualização do campo de velocidade 
    # Plota o campo de velocidade 
    qv = ax2[2,1].quiver(lons_plotnl[::step,::step]*1e-3,lats_plotnl[::step,::step]*1e-3,uplot[::step,::step],vplot[::step,::step],scale=scale,width=width)
    # Insere vetor de escala no local desejado 
    ax2[2,1].quiverkey(qv, 0.85, 0.87,0.1,r'$0.1 \frac{m}{s}$', labelpos='E', coordinates='figure')
    # -----------------------------------------------------------------
    # Configurações de colorbar, labels e legendas
    # cbar=plt.colorbar(cf)
    # ax2[2,1].set_xlabel('x [km]',weight='bold',fontsize=12)
    #ax2[2,1].set_ylabel('y [km]',weight='bold',fontsize=12)
    # cbar.set_label('h [m]',weight='bold',fontsize=12)

    # ========= t5 ============ # 
    
    t = 720
    uplot = (u[t,:,1:]+u[t,:,0:-1])/2
    vplot = (v[t,1:,:]+v[t,0:-1,:])/2
    
    titlestr='tempo = '+str(int(t*dt/3600)) + ' horas ($\mu$='+str(mu)+')'  
    ax2[3,0].set_title(titlestr, fontsize = 12, fontweight = 'bold')
        
    # ----------------------------------------------------------------
    # Plota contorno preenchido e isolinhas de elevação da superfície (h) 
    cf=ax2[3,0].contourf(lons_plotnl*1e-3,lats_plotnl*1e-3,h[t],np.linspace(hmin,hmax,51),cmap=cm.balance,extend='both')
    ax2[3,0].contour(lons_plotnl*1e-3,lats_plotnl*1e-3,h[t],np.linspace(hmin,hmax,30),linewidths=1,colors='k',extend='both')
    
    # ------------------------------------------------------------------
    # Visualização do campo de velocidade 
    # Plota o campo de velocidade 
    qv = ax2[3,0].quiver(lons_plotnl[::step,::step]*1e-3,lats_plotnl[::step,::step]*1e-3,uplot[::step,::step],vplot[::step,::step],scale=scale,width=width)
    # Insere vetor de escala no local desejado 
    ax2[3,0].quiverkey(qv, 0.85, 0.87,0.1,r'$0.1 \frac{m}{s}$', labelpos='E', coordinates='figure')
    # -----------------------------------------------------------------
    # Configurações de colorbar, labels e legendas
    # cbar=plt.colorbar(cf)
    ax2[3,0].set_xlabel('x [km]',weight='bold',fontsize=12)
    ax2[3,0].set_ylabel('y [km]',weight='bold',fontsize=12)
    # cbar.set_label('h [m]',weight='bold',fontsize=12)

    # ========= t6 ============ # 
    
    t = 1200
    uplot = (u[t,:,1:]+u[t,:,0:-1])/2
    vplot = (v[t,1:,:]+v[t,0:-1,:])/2
    
    titlestr='tempo = '+str(int(t*dt/3600)) + ' horas ($\mu$='+str(mu)+')'  
    ax2[3,1].set_title(titlestr, fontsize = 12, fontweight = 'bold')
        
    # ----------------------------------------------------------------
    # Plota contorno preenchido e isolinhas de elevação da superfície (h) 
    cf=ax2[3,1].contourf(lons_plotnl*1e-3,lats_plotnl*1e-3,h[t],np.linspace(hmin,hmax,51),cmap=cm.balance,extend='both')
    ax2[3,1].contour(lons_plotnl*1e-3,lats_plotnl*1e-3,h[t],np.linspace(hmin,hmax,30),linewidths=1,colors='k',extend='both')
    
    # ------------------------------------------------------------------
    # Visualização do campo de velocidade 
    # Plota o campo de velocidade 
    qv = ax2[3,1].quiver(lons_plotnl[::step,::step]*1e-3,lats_plotnl[::step,::step]*1e-3,uplot[::step,::step],vplot[::step,::step],scale=scale,width=width)
    # Insere vetor de escala no local desejado 
    ax2[3,1].quiverkey(qv, 0.85, 0.87,0.1,r'$0.1 \frac{m}{s}$', labelpos='E', coordinates='figure')
    # -----------------------------------------------------------------
    # Configurações de colorbar, labels e legendas
    # cbar=plt.colorbar(cf)
    ax2[3,1].set_xlabel('x [km]',weight='bold',fontsize=12)
    #ax2[2,1].set_ylabel('y [km]',weight='bold',fontsize=12)
    # cbar.set_label('h [m]',weight='bold',fontsize=12)
    
    # ======= FINAL ====== #
    # ax2[1,0].set_xticks([])
    # ax2[1,1].set_xticks([])
    # ax2[0,0].set_xticks([])
    # ax2[0,1].set_xticks([])

    #label
    fig2.subplots_adjust(right=0.8)
    cbar_ax = fig2.add_axes([0.85, 0.15, 0.05, 0.7])
    fig2.colorbar(cf,cax=cbar_ax) 
        
    pl.savefig('Frames_Linear_Centrado.png',bbox_inches='tight',dpi=300)   



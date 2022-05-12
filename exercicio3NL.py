# Import modules
import numpy as np
import time
import matplotlib.pyplot as plt 
import pylab as pl
import os
import glob
import math
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.mplot3d import Axes3D 
from matplotlib.colors import LightSource
from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation
import metpy.calc
from metpy.units import units
import cmocean.cm as cm
import sys 

g = 9.8
delta = 100e3
xmax,xmin = 40e5,-40e5
ymax,ymin = 40e5,-40e5
nx = int((xmax+delta-xmin)/delta) # (número de pontos da grade em x)
ny = int((ymax+delta-ymin)/delta) # (número de pontos da grade em y)  
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
mu = 0.2
H = 250
dt = (mu*delta/(np.sqrt(g*H)))/2
tmax = 1501
Nt = tmax 
# Iniciando integração # 
u = np.zeros([Nt,ny,nx+1])
v = np.zeros([Nt,ny+1,nx])
h = np.zeros([Nt,ny,nx])

 
# Salvar tempo inicial    
result_u = []
result_v = []
result_h = []
result_csi = []
# Para conferir a conservação da massa e energia 
Vol = [] # volume
vort = [] # volume
Ep = []  # energia potencial
Ek = []  # energia cinética
ape = [] # enstrofia potencial absoluta 
# dado de um ponto da grade para ser salvo ao longo da integração 
histgrid1 = [] 
histgrid2 = []  
histgrid3 = []   
histgrid4 = []    

crls = 1
if crls == 1:
    # Calculando o parâmetro de Coriolis utilizando aproximação de plano Beta 
    beta = (2*7.29*1e-5*np.cos(0*np.pi/180)/6400000)
if crls == 0:
    # Desliga Coriolis
    beta = 0

dx = delta 
dy = delta 


fu = latsu*beta
fv = latsv*beta
fz = latsz*beta

xgauss = np.linspace(xmin, xmax,nx)
#ygauss = np.linspace(4e6, ymin,ny) # centrado 4e6
ygauss = np.linspace(ymax, ymin,ny)  

longa,latga = np.meshgrid(xgauss,ygauss)

a = 0.25 # Variar 
cx = 0 # -4e6 (pulso na borda oeste)
cy = 0  #  -4e6 (pulso na borda sul) 
nrx = 4
nry = 4
    
gauss =  1 * a * np.exp(- (((longa-cx)**2)/(nrx*dx)**2) -  (((latga-cy)**2)/(nry*dx)**2))
    
alpha = 0.02 # Variar 

# Grades iniciais
unp1 =  np.zeros((ny,nx+1))
vnp1 =  np.zeros((ny+1,nx))
hnp1 =  np.zeros((ny,nx)) + H
U_np1 = np.zeros((ny,nx+1))
V_np1 = np.zeros((ny+1,nx)) 
Bnp1 = np.zeros((ny,nx))        
csinp1 = np.zeros((ny+1,nx+1))*np.nan

# Atualiza matrizes do tempo 
unm1,vnm1,hnm1 = unp1*np.nan,vnp1*np.nan,hnp1*np.nan              
un,vn,hn = unp1,vnp1,hnp1 
unp1,vnp1,hnp1 = un*np.nan,vn*np.nan,hn*np.nan 
t=0                                                                                  
U_nm1,V_nm1 = unp1*np.nan,vnp1*np.nan                               
U_n,V_n = U_np1,V_np1 
U_np1,V_np1 = un*np.nan,vn*np.nan 

n2 = 0 
n2dt = 5 
u[n2,:,:] = un; v[n2,:,:] = vn; h[n2,:,:]=hn
#=== ANIMAR OU SALVAR FIGURA ===# 
animacao = 0  # gera animação após rodar o código
salvafig = 0  # salva frames dos campos gerados pela integração
propriedades = 1 # calcula energia, volume, enstrofia
salvaframes = 0 # salva frames para o artigo
for t in range(1,tmax):
    # dt = 40
    decay = (H/2)*(alpha**3)*(t**2)*math.exp(
    -alpha*t)
    
    frc = gauss*decay
    
    # Calculo das matrizes de diagnóstico
    U_np1[:,1:-1] = un[:,1:-1] * .5 * (hn[:,:-1] + hn[:,1:])
    U_np1[:,0] = un[:,0] * hn[:,0]
    U_np1[:,-1] = un[:,-1] * hn[:,-1]   
    
    V_np1[1:-1] = vn[1:-1] * .5 * (hn[1:] + hn[:-1])
    V_np1[0] = vn[0] * hn[0]
    V_np1[-1] = vn[-1] * hn[-1]

    u_mean = .5 * (un[:,1:]**2 + un[:,:-1]**2)
    v_mean = .5 * (vn[:-1]**2 + vn[1:]**2)
    Bnp1 = g*(hn) + .5*(u_mean+v_mean)
    
    delta_v = (vn[1:-1,1:]-vn[1:-1,:-1])/dx
    delta_u = (un[:-1,1:-1] - un[1:,1:-1])/dy
    h_mean = .25*(hn[1:,:-1] + hn[1:,1:] + hn[:-1,:-1] + hn[:-1,1:])
    csinp1[1:-1,1:-1] = (fz[1:-1,1:-1] + delta_v - delta_u)/h_mean
    
    #NORTE E SUL
    csinp1[0,1:-1] = (fz[0,1:-1] + (vn[0,1:]-vn[0,0:-1])/dx)/((hn[0,:-1]+hn[0,1:])/2) 
    csinp1[-1,1:-1] = (fz[-1,1:-1] + (vn[-1,1:]-vn[-1,0:-1])/dx)/((hn[-1,0:-1]+hn[-1,1:])/2) 
    #LESTE E OESTE
    csinp1[1:-1,0] = (fz[1:-1,0] - (un[1:,0]-un[0:-1,0])/dy)/((hn[0:-1,0]+hn[1:,0])/2) 
    csinp1[1:-1,-1] = (fz[1:-1,-1] - (un[1:,-1]-un[0:-1,-1])/dy)/((hn[0:-1,-1]+hn[1:,-1])/2)    
    
    csinp1[np.isnan(csinp1)]=0
         

    # Primeiro passo de tempo:
    # discretização usando Euler avançado no tempo e centrado no espaço: 
    if t == 1:
                                    
        # calcula o campo de u
        V_mean = .5 * ((csinp1[:-1,1:-1] * .5 * (V_np1[:-1,:-1] + V_np1[:-1,1:])) + (csinp1[1:,1:-1] * .5 * (V_np1[1:,:-1] + V_np1[1:,1:])))
        delta_Bu = (Bnp1[:,1:] - Bnp1[:,:-1])/delta
        unp1[:,1:-1] = un[:,1:-1] + dt * (V_mean-delta_Bu)           
        
        # calcula o campo de v
        U_mean = .5 * ((csinp1[1:-1,1:] * .5 * (U_np1[1:,1:] + U_np1[:-1,1:])) +  (csinp1[1:-1,:-1] * .5 * (U_np1[1:,:-1] + U_np1[:-1,:-1])))
        delta_Bv = (Bnp1[:-1] - Bnp1[1:])/delta            
        vnp1[1:-1] = vn[1:-1] + dt * (-U_mean-delta_Bv)
        
        # calcula o campo de h
        delta_U =  (U_np1[:,1:] - U_np1[:,:-1])/delta
        delta_V = (V_np1[:-1,:] - V_np1[1:,:])/delta            
        hnp1 = hn + dt *  (-delta_U-delta_V)
    
    # Outros passos de tempo usando o esquema Leapfrog 
    # Restante da integração (t > 1) 
    else:           
    
        # Prognóstico do campo de u
        V_mean = .5 * ((csinp1[:-1,1:-1] * .5 * (V_np1[:-1,:-1] + V_np1[:-1,1:])) + (csinp1[1:,1:-1] * .5 * (V_np1[1:,:-1] + V_np1[1:,1:])))
        delta_Bu = (Bnp1[:,1:] - Bnp1[:,:-1])/delta
        unp1[:,1:-1] =  unm1[:,1:-1] + 2*dt * (V_mean-delta_Bu)  
        # Prognóstico do campo de v
        U_mean = .5 * ((csinp1[1:-1,1:] * .5 * (U_np1[1:,1:] + U_np1[:-1,1:])) + (csinp1[1:-1,:-1] * .5 * (U_np1[1:,:-1] + U_np1[:-1,:-1])))
        delta_Bv = (Bnp1[:-1] - Bnp1[1:])/delta
        vnp1[1:-1] = vnm1[1:-1] + 2*dt * (-U_mean-delta_Bv)     
        # Prognóstico do campo de h
        delta_U =  (U_np1[:,1:] - U_np1[:,:-1])/delta
        delta_V = (V_np1[:-1,:] - V_np1[1:,:])/delta
        hnp1 = hnm1 + 2*dt *  (-delta_U-delta_V)

    # ---------------------------------------------    
    # Adicionando a forçante 
    hnp1[1:-1,1:-1] =  hnp1[1:-1,1:-1] + (frc[1:-1,1:-1])     
    # ---------------------------------------------    

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
        for i in range(len(unp1[:,-1])):
            c = np.sqrt(g*(hn[i,-2]+H))
            # if c < 0 or np.isnan(c):
                # c = np.sqrt(g*H)
                # unp1[i,-1] = 0

            # elif c >= dx/dt:
                # c = dx/dt
                # unp1[i,-1] = un[i,-1] + dt * ((-c*(un[i,-1]-un[i,-2])/delta) + (0.5 * (fv[i,-1] * vn[i,-1] + fv[i+1,-1] * vn[i+1,-1])))
            
            unp1[i,-1] = un[i,-1] + dt * ((-c*(un[i,-1]-un[i,-2])/delta) + (0.5 * (fv[i,-1] * vn[i,-1] + fv[i+1,-1] * vn[i+1,-1])))
        
        # --------------------------------------------
        # Contorno de v: Radiacional parte superior (Norte)
        for i in range(len(vnp1[0])):
            c = np.sqrt(g*(hn[1,i]+H))
            # if c < 0 or np.isnan(c):
                # c = np.sqrt(g*H)
            # elif c > dx/dt:
                # c = dx/dt                        
            vnp1[0,i] = vn[0,i] + dt * ((-c*(vn[0,i]-vn[1,i])/delta) - (0.5 * (fu[0,i] * un[0,i] + fu[0,i+1] * un[0,i+1])))
        
        # --------------------------------------------
        # Contorno de v: Radiacional parte inferior (Sul)
        for i in range(len(vnp1[-1])):
            c = np.sqrt(g*(hn[-2,i]+H))
            # if c > 0 or np.isnan(c):
                # c = np.sqrt(g*H)
            # elif c > dx/dt:
                # c = dx/dt                        
            vnp1[-1,i] = vn[-1,i] + dt * ((c*(vn[-2,i]-vn[-1,i])/delta) - (0.5 * (fu[-1,i] * un[-1,i] + fu[-1,i+1] * un[-1,i+1])))  
            
            ç=1000000
            
    # ---------------------------------------------                 
    # Atualiza matrizes         
    unm1,vnm1,hnm1 = un,vn,hn           
    un,vn,hn = unp1,vnp1,hnp1
    csin = csinp1 
    unp1,vnp1,hnp1 = un*np.nan,vn*np.nan,hn*np.nan
    # 
    U_nm1,V_nm1 = U_n,V_n                   
    U_n,V_n,csin = U_np1,V_np1,csinp1        
    U_np1,V_np1,csinp1 = U_n*np.nan,V_n*np.nan,csin*np.nan
    u[t,:,:] = un; v[t,:,:] = vn; h[t,:,:]=hn
    if propriedades == 1:
        # ------------------------------------------------
        # Calculo das propriedades conservativas 
        # Volume total 
        Vol.append(((hn)*delta**2).sum())
        
        # Energia Potencial
        Ep.append(((hn)**2*delta**2).sum()*g/2)
        
        # Energia Cinética
        EKu = (un[:,:-1]+un[:,1:])/2
        EKv = (vn[:-1]+vn[1:])/2
        h_meanxy = hn.mean()
        Ek.append((((EKu**2 + EKv**2)*h_meanxy/2*(delta**2)).sum()))

        # Enstrofia # CONFERIR COM PROF
        h_meanx = (hn[:,:-1] + hn[:,1:])/2
        h_meanxy = (h_meanx[:-1] + h_meanx[1:])/2 
        tmp = csin*np.nan
        tmp[1:-1,1:-1] = (csin[1:-1,1:-1]**2)*h_meanxy
        tmp[0,:-1] = (csin[0,:-1]**2)*hn[0]
        tmp[-1,1:] = (csin[-1,1:]**2)*hn[-1]
        tmp[1:,0] = (csin[1:,0]**2)*hn[:,0]
        tmp[:-1,-1] = (csin[:-1,-1]**2)*hn[:,-1]
        ape.append((tmp*(delta**2)).sum()/2*ç)


        # ---------------------------------------------            
        # Salva os pontos centrais das grades
        histgrid1.append(hn[40,40])
        histgrid2.append(U_n[40,41])
        histgrid3.append(V_n[41,40])
        histgrid4.append(csin[41,41])
    
    
    # ================ ANIMAÇÃO ====================== # 
    # n2dt é o intervalo entre os frames de animação 
    if animacao == 1:
        if t % n2dt == 0:
            n2 = n2 + 1
            u[n2,:,:] = un; v[n2,:,:] = vn; h[n2,:,:]=hn
            
    # =============== SALVAR CAMPOS =================  # 
    # Plota e salva o momento inicial e os momentos em que t%x == 0:
    # | x = frequencia de plot 
    if salvafig == 1: 
        xxt = 40
        if t == 1 or t%xxt ==0:  
            # -----------------------------------
            # Interpolação de un e vn para plot 
            uplot = (un[:,1:]+un[:,0:-1])/2
            vplot = (vn[1:,:]+vn[0:-1,:])/2
            # -----------------------------------
            # Inicia o plot
            plt.clf()
            # Define subplots, tamanho e resolução da figura 
            fig, ax1 =plt.subplots(figsize=(9, 8),dpi=300)
            titlestr='tempo = '+str(int(t*dt/3600)) + 'horas ($\mu$='+str(mu)+')'  
            ax1.set_title(titlestr, fontsize = 12, fontweight = 'bold')
            # ---------------------------------------------------------------
            # Visualização de H 
            # Intervalo de valores de h 
            hmax = 3
            hmin = -hmax 
            # ----------------------------------------------------------------
            # Plota contorno preenchido e isolinhas de elevação da superfície (h) 
            cf=plt.contourf(lons_plot*1e-3,lats_plot*1e-3,hn-H,np.linspace(hmin,hmax,51),cmap=cm.balance,extend='both')
            plt.contour(lons_plot*1e-3,lats_plot*1e-3,hn-H,np.linspace(hmin,hmax,30),linewidths=1,colors='k',extend='both')
            
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
            plt.savefig('nl_sul_'+str(10000+t)+'.png',dpi=300,bbox_inches='tight')
        # -----------------------------------------------------------------
    
if propriedades == 1:
    Ep= np.array(Ep) #- Ep[0]
    Ek= np.array(Ek)
    conser= Ek+Ep
    
    
    
    # ============= Visualização Propriedades Conservativas ============ # 
    fig,axs = plt.subplots(5,figsize=(10, 12), constrained_layout=True)
    tx = np.arange(0,len(histgrid1)*dt/60/60/24,dt/60/60/24)
    axs[0].plot(tx,Vol,linewidth=1,color='k')
    #axs[0].set_ylim(1.6400*1e16,1.6600*1e16)
    axs[0].set_title('Volume',color='k', fontsize = 22)
    # axs[0].plot(tx, zorooo,linewidth=5,color='k')
    axs[0].tick_params(axis='both', which='major', labelsize=16)    
    axs[1].plot(tx,Ep,linewidth=1,color='k')
    #axs[1].set_ylim(2.000*1e19,2.0600*1e19)
    
    axs[1].set_title('Energia Potencial',color='k', fontsize = 22)
    axs[1].tick_params(axis='both', which='major', labelsize=16)
    axs[4].set_xlabel('Tempo (dias)', fontsize=18)       
    axs[1].grid() 
    axs[0].grid() 
    
    axs[2].plot(tx,Ek,linewidth=1,color='k')
    #axs[2].set_ylim(0,2*1e16)

    axs[2].set_title('Energia Cinética',color='k', fontsize = 22)      
    axs[2].tick_params(axis='both', which='major', labelsize=16)    
    axs[2].grid() 

    axs[3].plot(tx,conser,linewidth=1,color='k')
    #axs[3].set_ylim(2.0000*1e19,2.0600*1e19)
    axs[3].set_title('Energia Total',color='k', fontsize = 22)      
    axs[3].tick_params(axis='both', which='major', labelsize=16)  
    axs[3].grid() 

    axs[4].plot(tx,ape,linewidth=1,color='k')
    #axs[4].set_ylim(3.8500*1e8,3.9040*1e8)

    axs[4].set_title('Enstrofia',color='k', fontsize = 22)      
    axs[4].tick_params(axis='both', which='major', labelsize=16)  
    axs[4].grid() 
    

    
    zorooo = ape
    zorooo = np.zeros(np.shape(ape))
    zorooo= zorooo + Vol[0]
    axs[0].fill_between(tx, Vol, Vol[0], where=Vol >= zorooo, facecolor='0.7', interpolate=True, alpha=0.4)
    axs[0].fill_between(tx, Vol, Vol[0], where=Vol <= zorooo, facecolor='0.2', interpolate=True, alpha=0.4)
    zorooo = np.zeros(np.shape(ape))
    zorooo= zorooo + Ep[0]
    axs[1].fill_between(tx, Ep, Ep[0], where=Ep >= zorooo, facecolor='0.7', interpolate=True, alpha=0.4)
    axs[1].fill_between(tx, Ep, Ep[0], where=Ep <= zorooo, facecolor='0.2', interpolate=True, alpha=0.4)
    zorooo = np.zeros(np.shape(ape))
    zorooo= zorooo + Ek[0]
    axs[2].fill_between(tx, Ek, Ek[0], where=Ek >= zorooo, facecolor='0.7', interpolate=True, alpha=0.4)
    axs[2].fill_between(tx, Ek, Ek[0], where=Ek <= zorooo, facecolor='0.2', interpolate=True, alpha=0.4)
    zorooo = np.zeros(np.shape(ape))
    zorooo= zorooo + conser[0]
    axs[3].fill_between(tx, conser, conser[0], where=conser >= zorooo, facecolor='0.7', interpolate=True, alpha=0.4)
    axs[3].fill_between(tx, conser, conser[0], where=conser <= zorooo, facecolor='0.2', interpolate=True, alpha=0.4)
    zorooo = np.zeros(np.shape(ape))
    zorooo= zorooo + ape[0]
    axs[4].fill_between(tx, ape, ape[0], where=ape >= zorooo, facecolor='0.7', interpolate=True, alpha=0.4)
    axs[4].fill_between(tx, ape, ape[0], where=ape <= zorooo, facecolor='0.2', interpolate=True, alpha=0.4)
    
    pl.savefig('N_Linear_Centrado1.png')     



if animacao == 1:   
    #Determinar valores maximos e minimos para plotar:
    hmax = 15
    hmin= -hmax


    i = 0 #frame zero da animacao
    uplot = np.zeros([Nt,ny,nx])
    vplot = np.zeros([Nt,ny,nx])
    titlestr = 'tempo = '+str(int(i*dt/3600))+' horas ($\mu$='+str(mu)+')'

    # Define subplots, tamanho e resolução da figura 
    fig, ax1 =plt.subplots(figsize=(9, 8))
    hplot = ax1.contourf(lons_plot*1e-3,lats_plot*1e-3,h[i]-H,np.linspace(hmin,hmax,51),cmap=cm.balance,extend='both')
    # Interpolação para plot de u e v 
    uplot[i,:,:] = (u[i,:,1:]+u[i,:,0:-1])/2
    vplot[i,:,:] = (v[i,1:,:]+v[i,0:-1,:])/2
    # Parâmetros do quiver 
    step=3
    width=0.002
    scale=2.5
    qv = ax1.quiver(lons_plot[::step,::step]*1e-3,lats_plot[::step,::step]*1e-3,uplot[i,::step,::step],vplot[i,::step,::step],scale=scale,width=width)
    ax1.quiverkey(qv, 0.95, 0.97,0.1,r'$0.1 \frac{m}{s}$', labelpos='E', coordinates='figure')
         
    # Configurações 
    cbar=fig.colorbar(hplot)
    ax1.set_xlabel('x [km]',weight='bold',fontsize=12)
    ax1.set_ylabel('y [km]',weight='bold',fontsize=12)
    cbar.set_label('h [m]',weight='bold',fontsize=12)
    ax1.set_title(titlestr, fontsize = 12, fontweight = 'bold')
    # Função para animar ao longo do tempo 
    def animate(i):
        titlestr='tempo = '+str(int(i*dt*n2dt/3600)) + 'horas ($\mu$='+str(mu)+')'        
        ax1.collections = []
        cf=plt.contourf(lons_plot*1e-3,lats_plot*1e-3,h[i]-H,np.linspace(hmin,hmax,51),cmap=cm.balance,extend='both')
        plt.contour(lons_plot*1e-3,lats_plot*1e-3,h[i]-H,np.linspace(hmin,hmax,30),linewidths=1,colors='k',extend='both')
        uplot[i,:,:] = (u[i,:,1:]+u[i,:,0:-1])/2
        vplot[i,:,:] = (v[i,1:,:]+v[i,0:-1,:])/2
        qv = ax1.quiver(lons_plot[::step,::step]*1e-3,lats_plot[::step,::step]*1e-3,uplot[i,::step,::step],vplot[i,::step,::step],scale=scale,width=width)
        ax1.quiverkey(qv, 0.95, 0.97,0.1,r'$0.1 \frac{m}{s}$', labelpos='E', coordinates='figure')
        ax1.set_title(titlestr, fontsize = 12, fontweight = 'bold')
    
    anim = FuncAnimation(
        fig, animate, interval = 100, frames = Nt
    )

    plt.tight_layout()
    plt.show()




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
    cf=ax2[0,0].contourf(lons_plotnl*1e-3,lats_plotnl*1e-3,h[t]-H,np.linspace(hmin,hmax,51),cmap=cm.balance,extend='both')
    ax2[0,0].contour(lons_plotnl*1e-3,lats_plotnl*1e-3,h[t]-H,np.linspace(hmin,hmax,30),linewidths=0.7,colors='k',extend='both')
    
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
    cf=ax2[0,1].contourf(lons_plotnl*1e-3,lats_plotnl*1e-3,h[t]-H,np.linspace(hmin,hmax,51),cmap=cm.balance,extend='both')
    ax2[0,1].contour(lons_plotnl*1e-3,lats_plotnl*1e-3,h[t]-H,np.linspace(hmin,hmax,30),linewidths=1,colors='k',extend='both')
    
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
    cf=ax2[1,0].contourf(lons_plotnl*1e-3,lats_plotnl*1e-3,h[t]-H,np.linspace(hmin,hmax,51),cmap=cm.balance,extend='both')
    ax2[1,0].contour(lons_plotnl*1e-3,lats_plotnl*1e-3,h[t]-H,np.linspace(hmin,hmax,30),linewidths=1,colors='k',extend='both')
    
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
    cf=ax2[1,1].contourf(lons_plotnl*1e-3,lats_plotnl*1e-3,h[t]-H,np.linspace(hmin,hmax,51),cmap=cm.balance,extend='both')
    ax2[1,1].contour(lons_plotnl*1e-3,lats_plotnl*1e-3,h[t]-H,np.linspace(hmin,hmax,30),linewidths=1,colors='k',extend='both')
    
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
    cf=ax2[2,0].contourf(lons_plotnl*1e-3,lats_plotnl*1e-3,h[t]-H,np.linspace(hmin,hmax,51),cmap=cm.balance,extend='both')
    ax2[2,0].contour(lons_plotnl*1e-3,lats_plotnl*1e-3,h[t]-H,np.linspace(hmin,hmax,30),linewidths=1,colors='k',extend='both')
    
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
    cf=ax2[2,1].contourf(lons_plotnl*1e-3,lats_plotnl*1e-3,h[t]-H,np.linspace(hmin,hmax,51),cmap=cm.balance,extend='both')
    ax2[2,1].contour(lons_plotnl*1e-3,lats_plotnl*1e-3,h[t]-H,np.linspace(hmin,hmax,30),linewidths=1,colors='k',extend='both')
    
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
    cf=ax2[3,0].contourf(lons_plotnl*1e-3,lats_plotnl*1e-3,h[t]-H,np.linspace(hmin,hmax,51),cmap=cm.balance,extend='both')
    ax2[3,0].contour(lons_plotnl*1e-3,lats_plotnl*1e-3,h[t]-H,np.linspace(hmin,hmax,30),linewidths=1,colors='k',extend='both')
    
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
    cf=ax2[3,1].contourf(lons_plotnl*1e-3,lats_plotnl*1e-3,h[t]-H,np.linspace(hmin,hmax,51),cmap=cm.balance,extend='both')
    ax2[3,1].contour(lons_plotnl*1e-3,lats_plotnl*1e-3,h[t]-H,np.linspace(hmin,hmax,30),linewidths=1,colors='k',extend='both')
    
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
        
    pl.savefig('Frames_NL_Deslocado.png',bbox_inches='tight',dpi=300)    

# calculate the rotation sum for the seismograms received from seissol output
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from lib_seissol_data import data_processing
from _utils import phases,filter,gauss_window
import re
import os, sys

path='/home/djamel/PHD_projects/force_on_hill/results_seismogram/'

seissol_arr=['model_real_2_f_peak']

d_a=19000   #distance first receiver
d_r=500

dt=0.01
n_circ=28   #circular receiver
n_rad=63    # radial receiver
ts=700      # timestep





####

data_seissol = {}
for nu1 in seissol_arr:
     data_seissol[nu1] =np.load( path+nu1+'.npz')



############## choose parameters for plotting

ampl=np.arange(25, 88)*75000
ampl[0:1]=ampl[0:1]*np.arange(2, 3)[::-1]*1.2


comp='R'                # R or T component

if comp=='T':
    # Love modes ######
    phase_y=[7500, 49500, 49500]
    phase1_x=[3, 23.0, 32.5]
    phase2_x=[12.0, 33.5,34.0]
    hase4_x = [12.0, 33.5, 34.0]
    phase3_x=[-1.0, 11.5, 23.0]

    fm=str(int((phase_y[1]-phase_y[0]) / (phase1_x[1]-phase1_x[0])))            # fundamental Love wave
    m1=str(int((phase_y[1]-phase_y[0]) / (phase2_x[1]-phase2_x[0])))             # 1. Love mode
    m2=str(int((phase_y[1]-phase_y[0]) / (phase3_x[1]-phase3_x[0])))             # 1. Love mode



if comp=='R':
    # Love modes ######
    phase_y=[7500, 51000, 51000]
    phasefm_x=[0, 21.0, 27.5]

    phase1_x=[0.0, 17.5,30.0]
    phase2_x = [0.0, 14.8, 30.0]
    phase3_x = [0.0, 12.8, 30.0]
    # phase4_x = [0.0, 16.5, 30.0]
    #
    #
    # phase3_x=[0.0, 12.5, 21.0]

    fm=str(int((phase_y[1]-phase_y[0]) / (phasefm_x[1]-phasefm_x[0])))            # fundamental Love wave
    m1=str(int((phase_y[1]-phase_y[0]) / (phase1_x[1]-phase1_x[0])))             # 1. Love mode
    m2 = str(int((phase_y[1] - phase_y[0]) / (phase2_x[1] - phase1_x[0])))  # 1. Love mode
    m3 = str(int((phase_y[1] - phase_y[0]) / (phase3_x[1] - phase1_x[0])))  # 1. Love mode
    # m3 = str(int((phase_y[1] - phase_y[0]) / (phase4_x[1] - phase2_x[0])))  # 1. Love mode
    # m2=str(int((phase_y[1]-phase_y[0]) / (phase3_x[1]-phase3_x[0])))             # 1. Love mode



########################
global col
col = ['darkred', 'b', 'g', 'r', 'm', 'c']

def legend_plot(data_seissol):

    global  seissol_arr
    red={}
    strii=[]
    for ii in range(len(seissol_arr)):
        red[ii] = mlines.Line2D([], [], color=col[ii], markersize=10, linewidth=2,label='$F_t/F_n$= '+str(data_seissol[seissol_arr[ii]]['ttn']))
        strii.append(red[ii])

    plt.axis([0, 27, 15000, 54000])
    #plt.legend(handles=strii[0:len(seissol_arr)], loc='lower right')
    #plt.title(comp+'-component' + ', $\phi$= ' + str(int(phii[kk])), fontsize=17)
    plt.ylabel('distance [m]')
    plt.xlabel('t [s]')
    plt.rc('xtick', labelsize=10)



####################### filter!!!
dt=0.05
f_n=1/(dt*2)

f = 0.01
fe = 5.2
Wn = f / f_n
Wn = [f / f_n, fe / f_n]


#data_filt=filter(data.data_R, 4, Wn)
####################


phii=np.arange(0,360,360/float(data_seissol[seissol_arr[0]]['n_circ']),dtype=float)
for kk in range(5,6):
    fig=plt.figure(figsize=(24, 16))
    #plt.subplot(121)
    ax=fig.add_subplot(1, 2, 1)
    #ax.patch.set_alpha(0.2)


    hh = 0
    for ii in range(0,n_rad):
        if ii % 4== 0:
            ll = 0
            for ari in seissol_arr:
                plt.plot(data_seissol[ari]['t'],d_a + ii * d_r + ampl[ii] * gauss_window(filter(data_seissol[ari][comp][:,kk,ii],6,Wn),data_seissol[ari]['t'],ii*5,-275) , col[ll], linewidth=2)
                #plt.plot(data_seissol[ari]['t'],d_a + ii * d_r + ampl[ii] * data_seissol[ari][comp][:, kk+10, ii], 'g',linewidth=2)
                #plt.plot(data_seissol[ari]['t'],d_a + ii * d_r + ampl[ii] * data_seissol[ari][comp][:, kk, ii], 'b',linewidth=2)
                ll=ll+1

        # plt.plot(phasefm_x, phase_y, color='wheat', linewidth=3, alpha=1.0)
        # plt.plot(phase1_x, phase_y, color='wheat', linewidth=3, alpha=1.0)
        # plt.plot(phase2_x, phase_y, color='blue', linewidth=3, alpha=1.0)
        # plt.plot(phase3_x, phase_y, color='green', linewidth=3, alpha=1.0)
        # plt.plot(phase1_x, phase_y,'chocolate', linewidth=3, alpha=1.0)

        hh=hh+1

    plt.fill_between(phase3_x, phase_y,color='wheat', alpha=1.0)
    #plt.fill_between(phase2_x, phase_y,color='wheat', alpha=0.8)
    plt.fill_between(phase1_x, phase_y, color='ivory', alpha=0.8)
    plt.fill_between(phasefm_x, phase_y, color='tan', alpha=0.8)


    #plt.text(21.0, 51500, '0. mode', fontsize=18,bbox={'facecolor': 'chocolate', 'alpha': 1})
    plt.text(17.8, 51500, '1./2./3. mode ', fontsize=18, bbox={'facecolor': 'ivory', 'alpha': 1})
    #plt.text(19.0, 51500, '2. mode ' + str(m2), fontsize=18, bbox={'facecolor': 'chocolate', 'alpha': 1})
    #plt.text(13.0, 51500, 'higher modes ', fontsize=18, bbox={'facecolor': 'wheat', 'alpha': 1})
    #plt.text(21.0, 51500, '0. mode ', fontsize=18, bbox={'facecolor': 'tan', 'alpha': 1})

    #plt.text(12.5, 51500, 'overtones', fontsize=18, bbox={'facecolor': 'wheat', 'alpha': 1})


    #plt.text(25.5, 55500, comp+'-component' + ', $\phi$= ' + str(int(phii[kk])), fontsize=23,weight='semibold')

    legend_plot(data_seissol)



    # stri='wiggle, models real 2, f peak, '+comp+'-component'+str(360/26*kk)+'.png'
    # fig.savefig(stri,format='png')      # save figure
    #
    #
    #
    # data_str='/home/djamel/PHD_projects/force_on_hill/results_seismogram/model_real_2_f_peak.npy'
    # data = np.load(data_str)
    # data=data_processing(data,0.05,700,28,63)
    # data.radial_transversal()
    # swi='R'      # choose which phase 'R' or 'T'
    # v_r=[kk,kk+1]      # range arround phi
    #
    #
    # data.disperison(100,4000,1500,500,swi,v_r,data_str)


    data_str='/home/djamel/PHD_projects/force_on_hill/results_seismogram/model_real_2_f_peak.npy'
    data = np.load(data_str)

    dt=0.05
    f_n=1/(dt*2)

    f=0.1
    fe=2
    Wn=f/f_n
    Wn = [f / f_n, fe/f_n]



    data=data_processing(data,0.05,700,28,63)
    data.radial_transversal()
    data.filter_but(6,Wn)


    swi='R'      # choose which phase 'R' or 'T'
    v_r=[kk,kk+1]      # range arround phi

    data.disperison(100,4000,1500,500,swi,v_r,data_str)



#
    # stri='wiggle, models real 2, f peak, higher modes, '+comp+'-component'+str(360/26*kk)+'.png'
    # fig.savefig(stri,format='png')      # save figure



plt.show()
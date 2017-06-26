# class to change data from 2D to 3D

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import inspect
import os.path
import re



from _utils import max_2D_array, _theo_disp

class data_processing(object):

    '''
    this class reads data from npy array and rearrange it to process it; several methods to analyis the data
    dt -> delta t
    t_n-> number of time samples
    n_circ-> number of receivers in circular dimensions
    n_r -> numbers of receivers in radial direction
    '''

    def __init__(self,data,dt,t_n,n_circ,n_r):
        self.data=data
        self.dt = dt
        self.t_n = t_n
        self.n_circ = n_circ
        self.n_r = n_r

        # brauche ich fuer jeden vierten output
        self.data_t = np.zeros((self.t_n))
        self.data_X = np.zeros((self.t_n, self.n_r, self.n_circ))
        self.data_Y = np.zeros((self.t_n, self.n_r, self.n_circ))
        self.data_Z = np.zeros((self.t_n, self.n_r, self.n_circ))
        self.rot_X = np.zeros((self.t_n, self.n_r, self.n_circ))
        self.rot_Y = np.zeros((self.t_n, self.n_r, self.n_circ))
        self.rot_Z = np.zeros((self.t_n, self.n_r, self.n_circ))

        self.data_t[:] = self.data[:, 0, 0]
        ll = 0
        for ii in range(0, self.n_circ):
            for jj in range(0, self.n_r):

                self.data_X[:, jj, ii] = self.data[:, 1, ll]
                self.data_Y[:, jj, ii] = self.data[:, 2, ll]
                self.data_Z[:, jj, ii] = self.data[:, 3, ll]
                self.rot_X[:,jj,ii]=self.data[:, 4, ll]
                self.rot_Y[:,jj,ii]=self.data[:, 5, ll]
                self.rot_Z[:,jj,ii]=self.data[:, 6, ll]
                ll = ll + 1


    # diese funktion berechnet aus X,Y,Z die radiale und transversal kompnente aus (Z bleicht gleich)

    #die parameter sind die gleichen wie oben

    def radial_transversal(self):

        phii = np.linspace(0.0, 2.0 * np.pi * (self.n_circ - 1) / self.n_circ, self.n_circ)

        self.data_R = np.zeros((self.t_n, self.n_r, self.n_circ))  # data for RADIAL rotation
        self.data_T = np.zeros((self.t_n, self.n_r, self.n_circ))  # data for TRANSVERSAL rotation

        for ii in range(0, self.n_circ):
            for jj in range(0, self.n_r):
                self.data_R[:, jj, ii] = np.cos(phii[ii]) * self.data_X[:, jj, ii] + np.sin(phii[ii]) * self.data_Y[:, jj,ii]       # radial direction
                self.data_T[:, jj, ii] = -np.sin(phii[ii]) * self.data_X[:, jj, ii] + np.cos(phii[ii]) * self.data_Y[:, jj,ii]      # transversal direction









    def disperison(self,incr,v_a,v_e,dist,swi,vec_phi,data_str):
        '''swi: switch -> 0 to calculate Raylegih disp; 1 to calculate Love disp
        PARAMETERS:
            incr: icrement of velocity points
            v_a: velocity anfang
            v_e: velocity ende
            dist: receiver distance: only possible for the same intervall between the receivers
            swi: switch between T and R
            vec_phi: calculate disp curves in phi direction
            data: receiver for the calculation of the disp curves
            data_max: max value for the normalisation of tranversal component -> R compoment is allways 1; therefore you have for the default value
            a data_max=1; then you switch from 1 to the array that incluedes the max values for all R components
        '''

        v = np.linspace(v_a, v_e, incr)  # velocity vector
        freq = np.fft.fftfreq(len(self.data_R[:, 0, 0]), self.dt)
        w = 2 * np.pi * freq
        x = np.linspace(0, (self.n_r - 1) * dist, self.n_r)
        phii = np.arange(0, 360, 360 / float(len(self.data_R[0,0,:])), dtype=float)

        maxx = np.zeros((self.n_circ, 1))  # array fuer maximale werte kannst du nachher brauchen um max value vs. phi zu plotten

        ####### strings for saving
        part_str = re.split(r"/", data_str)

        part_str_2 = re.split(r"\.n", part_str[-1])


        str_max = 'max_' + part_str_2[0]


        if (swi == 'R'):
            str_max = str_max + '_R.npy'
        elif (swi == 'T'):
            str_max = str_max + '_T.npy'


        # string for max values saving
        if (swi=='T'):
            data_max = np.load('max_' + part_str_2[0]+ '_R.npy')
        elif(swi=='R'):
            data_max = 1

        # string for plot saving
        if(swi=='R'):
            stri_save = 'dispersion_'+part_str_2[0]+'_R'
        elif(swi=='T'):
            stri_save = 'dispersion_'+part_str_2[0]+'_T'

        # create folder
        if(swi=='R'):
            stri_fold = part_str_2[0]+'_R'
        elif(swi=='T'):
            stri_fold = part_str_2[0]+'_T'
        if not os.path.exists(stri_fold):
            os.makedirs(stri_fold)


        #################################


        for iii in range(vec_phi[0], vec_phi[1]):
            DATA = np.zeros((self.t_n, self.n_r), dtype=complex)

            if swi=='R':
                for ii in range(0, self.n_r):
                    DATA[:, ii] = np.fft.fft(self.data_R[:, ii, iii])

            if swi=='T':
                for ii in range(0, self.n_r):
                    DATA[:, ii] = np.fft.fft(self.data_T[:, ii, iii])

            DISP = np.zeros((len(self.data_R[:, 1, 0]), len(v)), dtype=complex)
            for ii in range(0, len(freq)):
                for jj in range(0, len(v)):
                    for kk in range(0, 85):
                        DISP[ii, jj] = DISP[ii, jj] + np.exp(1j * w[ii] * x[kk] / v[jj]) * DATA[ii, kk]

            maxx[iii, 0] = max_2D_array(DISP)


            for ii in range(0,700):
                DISP[ii,:]=DISP[ii,:]/max(abs(DISP[ii,:]))

            if swi=='T':
                DISP = DISP * maxx[iii] / data_max[iii]

            fig = plt.figure(figsize=(9, 10))
            _theo_disp(swi)
            plt.imshow(abs(DISP[0:180, :]).T, aspect='auto', extent=(freq[0], freq[180], 1500, 4000))
            plt.colorbar()
            plt.xlabel('frequency [Hz]', fontsize=17)
            plt.ylabel('velocity [m/s]', fontsize=17)

            if swi=='R':
                plt.title('R-component, $F_t/F_n$=' + str(np.round(0.583, 2)) + ', $\phi$= ' + str(int(phii[iii])),fontsize=17)
            elif swi=='T':
                plt.title('T-component, $F_t/F_n$=' + str(np.round(0.583, 2)) + ', $\phi$= ' + str(int(phii[iii])),fontsize=17)

            str_save=stri_save+ ', phi= ' + str(int(phii[iii]))+'.png'
            print(str_save)
            fig.savefig(stri_fold+'/'+str_save,format='png')      # save figure

            plt.show()


        np.save(str_max, maxx)

data_str='/home/djamel/PHD_projects/force_on_hill/results_seismogram/model_3_f_0.2.npy'
data = np.load(data_str)
data=data_processing(data,0.05,700,28,85)
data.radial_transversal()
swi='R'      # choose which phase 'R' or 'T'
v_r=[0,1]      # range arround phi

data.disperison(100,4000,1500,500,swi,v_r,data_str)


print('hallo wiggle')
print('qiabudaliz')
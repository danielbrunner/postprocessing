# mit diesem script will ich die first arrival berchnen
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines


def phases(d_r,dr_1,num_r):


    # dr_1=7600   #distance first receiver
    # # d_r=500    # distance receivers
    # # num_r=85    # number of receivers

    vp_1=3398.0   # velocity erste layer
    vp_2=4016       # velocity zweite layer
    vp_3=5009       # velocity  dritte layer
    vp_4 = 5818  # velocity  dritte layer

    vs_1=2274   # velocity erste layer
    vs_2=2723       # velocity zweite layer
    vs_3=3393       # velocity  dritte layer
    vs_4 =3964  # velocity  dritte layer


    h1=1000          # thickness of the first and second layer
    h2=1000             # thickness of the third layer
    h3=2000  # thickness of the third layer

    #########
    ts=0# sollte eigentlich 0 sein ---> habe es so eingestellt das first arrival uebereinstimmt
    ########

    dist_r=np.zeros((num_r,1))      # distance between receivers
    for ii in range(0,num_r):
        dist_r[ii,0]=dr_1+d_r*ii


    t_dir_p=np.zeros((num_r,1))
    t_dir_s=np.zeros((num_r,1))

    for ii in range(0,num_r):
        t_dir_p[ii,0]=ts+dist_r[ii,0]/vp_1
        t_dir_s[ii,0]=ts+dist_r[ii,0]/vs_1


    # reflected phase
    t_refl_p=np.zeros((num_r,1))
    t_refl_s=np.zeros((num_r,1))

    for ii in range(0,num_r):
        t_refl_p[ii,0]=ts+2.0/vp_1*np.sqrt((dist_r[ii,0]/2.0)**2.0+h1**2.0)
        t_refl_s[ii,0]=ts+2.0/vs_1*np.sqrt((dist_r[ii,0]/2.0)**2.0+h1**2.0)



    # refrected wave at the first layer


    t_refr_p_1=np.zeros((num_r,1))
    t_refr_s_1=np.zeros((num_r,1))

    for ii in range(0,num_r):
        t_refr_p_1[ii,0]=ts+2*h1*np.sqrt((1/vp_1)**2-(1/vp_2)**2)+dist_r[ii]/vp_2
        t_refr_s_1[ii,0]=ts+2*h1*np.sqrt((1/vs_1)**2-(1/vs_2)**2)+dist_r[ii]/vs_2




    # refrected wave at the second layer

    t_refr_p_2=np.zeros((num_r,1))
    t_refr_s_2=np.zeros((num_r,1))

    for ii in range(0,num_r):
        t_refr_p_2[ii,0]=ts+2*h1*np.sqrt((1/vp_1)**2-(1/vp_2)**2)+2*h2*np.sqrt((1/vp_2)**2-(1/vp_3)**2)+dist_r[ii]/vp_3
        t_refr_s_2[ii,0]=ts+2*h1*np.sqrt((1/vs_1)**2-(1/vs_2)**2)+2*h2*np.sqrt((1/vs_2)**2-(1/vs_3)**2)+dist_r[ii]/vs_3






    # refrected wave at the third layer

    t_refr_p_3=np.zeros((num_r,1))
    t_refr_s_3=np.zeros((num_r,1))

    for ii in range(0,num_r):
        t_refr_p_3[ii,0]=ts+2*h1*np.sqrt((1/vp_1)**2-(1/vp_2)**2)+2*h2*np.sqrt((1/vp_2)**2-(1/vp_3)**2)+2*h3*np.sqrt((1/vp_3)**2-(1/vp_4)**2)+dist_r[ii]/vp_4
        t_refr_s_3[ii,0]=ts+2*h1*np.sqrt((1/vs_1)**2-(1/vs_2)**2)+2*h2*np.sqrt((1/vs_2)**2-(1/vs_3)**2)+2*h3*np.sqrt((1/vs_3)**2-(1/vs_4)**2)+dist_r[ii]/vs_4



    return dist_r,t_dir_p,t_dir_s,t_refl_p,t_refl_s,t_refr_p_1,t_refr_s_1, t_refr_p_2,t_refr_s_2,t_refr_p_3,t_refr_s_3



def max_2D_array(A):
    '''
    find maximum in a 2D array
    '''
    len_x=len(A[0,:])
    ma = np.zeros((1, len_x))
    for kk in range(0, len_x):
        ma[0, kk] = max(A[:,kk])
    maa=max(ma[0,:])
    return maa



def _theo_disp(swi):
    'diese funktion ist nur fuer dispersion methode zu verwenden'

    green_line = mlines.Line2D([], [], color='chocolate', markersize=10, linewidth=2, label='Fundamental Mode')
    red_line = mlines.Line2D([], [], color='wheat', markersize=10, linewidth=2, label='overtones')
    black_line = mlines.Line2D([], [], color='c', markersize=10, linewidth=2, label='2. overtone')
    cyan_line = mlines.Line2D([], [], color='c', markersize=10, linewidth=2, label='3. overtone')
    magenta_line = mlines.Line2D([], [], color='c', markersize=10, linewidth=2, label='4. overtone')
    yellow_line = mlines.Line2D([], [], color='c', markersize=10, linewidth=2, label='5. overtone')

    prov=[green_line, red_line, black_line, cyan_line, magenta_line,yellow_line]

    file = open("/home/djamel/GEOPSY/gpdc/four_layer_model_5_mode.disp", "r")
    line = file.readlines()

    if swi=='R':
        swi=0
    elif swi=='T':
        swi=1



    b=[]
    ll=0
    for ii in range(0,len(line)):

        if line[ii][0]=='#':
            if line[ii+1][0] != '#':
                b.append(ll)
            elif line[ii+1][0] == '#':
                if line[ii+2][0] == '#':
                    continue
                elif line[ii+2][0] != '#':
                    b.append(ll)
        elif line[ii][0]!='#':
            ll=ll+1
        last=ii
    b.append(last)

    b=[x-1 for x in b[1:]]
    b[0]=0
    b=np.reshape(b,(2,len(b)/2))






    data = np.loadtxt('/home/djamel/GEOPSY/gpdc/four_layer_model_5_mode.disp', skiprows=8)




    for ii in range(0,len(b[0,:])-1):
        if ii==0:
            plt.plot(data[b[swi,ii]+1:b[swi,ii+1]-2, 0], 1 / data[b[swi,ii]+1:b[swi,ii+1]-2, 1], linewidth=4,color='chocolate', alpha=0.9)
        if ii != 0:
            plt.plot(data[b[swi, ii] + 1:b[swi, ii + 1] - 2, 0], 1 / data[b[swi, ii] + 1:b[swi, ii + 1] - 2, 1],linewidth=4, color='wheat', alpha=0.9)
        ###############!!!!!!! -2 bei array fuer 20Hz getan weil sonst kurve udberlappt vorher war -2 nicht da!!!!!

    plt.legend(handles=prov[0:2], loc=3, fontsize=16)


def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx


def _mode_max(DISP,freq,v,swi):
    '''hier wird die maximum values der einzelnen moden gefunden und die freq und v dieses ortes zurueckgegegen'''
    file = open("/home/djamel/GEOPSY/gpdc/four_layer_model_5_mode.disp", "r")
    line = file.readlines()
    #swi='T'
    if swi == 'R':
        swi = 0
    elif swi == 'T':
        swi = 1

    b = []
    ll = 0
    for ii in range(0, len(line)):

        if line[ii][0] == '#':

            if line[ii + 1][0] != '#':

                b.append(ll)

            elif line[ii + 1][0] == '#':

                if line[ii + 2][0] == '#':

                    continue
                elif line[ii + 2][0] != '#':
                    b.append(ll)
        elif line[ii][0] != '#':
            ll = ll + 1
        last = ii
    b.append(last)

    b = [x - 1 for x in b[1:]]
    b[0] = 0
    b = np.reshape(b, (2, len(b) / 2))

    data = np.loadtxt('/home/djamel/GEOPSY/gpdc/four_layer_model_5_mode.disp', skiprows=8)

    fr=[]
    vv=[]
    ma=0.0
    maa=[]
    int_h = 1
    int_l=3

    incr_a=[17,15,10,4]
    incr_e=[70,50,41,33]
    ll=0

    for kk in range(0,4):
        fr=[]
        vv=[]
        for jj in range(incr_a[ll],incr_e[ll],1):
            for ii in range(b[swi,kk]+1,b[swi,kk+1]-2):
                if ii==b[swi,kk]+jj:
                    ww=[data[ii, 0], 1 / data[ii, 1]]
                    fr.append(find_nearest(freq,ww[0]))
                    vv.append(find_nearest(v, ww[1]))


            for ii in range(0, len(vv)):
                prov=max_2D_array(abs(DISP[fr[ii]-int_h:fr[ii]+int_h,vv[ii]-int_l:vv[ii]+int_l]))
                if prov>ma:
                    ma=prov
        maa=np.append(maa,ma)
        ma = 0.0

        for ii in range(0,len(vv)):
            plt.imshow(abs(DISP[fr[ii]-int_h:fr[ii]+int_h,vv[ii]-int_l:vv[ii]+int_l]).T,
                       aspect='auto', extent=(freq[fr[ii]-int_l], freq[fr[ii]+int_l], v[vv[ii]+int_h], v[vv[ii]-int_h]))



        plt.plot(data[b[swi, kk] + 1:b[swi, kk + 1] - 2, 0], 1 / data[b[swi, kk] + 1:b[swi, kk + 1] - 2, 1], linewidth=3)


        ll=ll+1


    #plt.show()

    return maa



    # plt.plot(ww[:,0])
    # plt.show()
        ###############!!!!!!! -2 bei array fuer 20Hz getan weil sonst kurve udberlappt vorher war -2 nicht da!!!!!

#_mode_max()

    # def theo_dispersion_curve_GEOPSY():
#     # f=open('/home/djamel/PHD_projects/force_on_hill/theoretische_dispersion_curve/model_1.disp','r')
#     # for ii in range(0,8):
#     #    trash=f.readline()
#
#
#     # data=np.zeros((80,2))
#     data = np.loadtxt('/home/djamel/PHD_projects/force_on_hill/theoretische_dispersion_curve/model_1.disp', skiprows=8)
#
#     # Rayleigh wave 1 mode
#     plt.plot(data[1:99, 0], 1 / data[1:99, 1], 'r')
#
#     # Love wave 1 mode
#     plt.plot(data[101:180, 0], 1 / data[101:180, 1], 'b')
#     plt.axis([0, 5, 2000, 4000])
#     plt.show()
# test=np.load('./data python/data max polar/model 1; max value in dependence of phi; transversal component; frequency=1 Hz.npy')
# test2=np.load('./data python/data max polar/model 1; max value in dependence of phi; radial component; frequency=1 Hz.npy')
# test2=np.take(test2, range(0, len(test2)+1), mode='wrap')
# test=np.take(test, range(0, len(test)+1), mode='wrap')
# testt1=test/test2
#
#
#
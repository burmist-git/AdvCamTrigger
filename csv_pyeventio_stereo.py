from eventio import SimTelFile

import numpy as np
import pickle as pkl
import pandas as pd
import sys
import time

def test(datafilein = "/scratch/snx3000/lburmist/simtel_data/proton/data/corsika_run1.simtel.gz"):
    #
    print("test")
    print("datafilein : ",datafilein)
    #
    sf = SimTelFile(datafilein)
    #
    nevloop = 0
    #
    for ev in sf:
        print("------------------------------")
        print("event_id     ",ev['event_id'])
        print("keys         ",ev['photoelectrons'].keys())    
        print("energy       ",ev['mc_shower']['energy'])
        print("azimuth      ",ev['mc_shower']['azimuth'])
        print("altitude     ",ev['mc_shower']['altitude'])
        print("h_first_int  ",ev['mc_shower']['h_first_int'])
        print("xmax         ",ev['mc_shower']['xmax'])
        print("hmax         ",ev['mc_shower']['hmax'])
        print("emax         ",ev['mc_shower']['emax'])
        print("cmax         ",ev['mc_shower']['cmax'])
        print("xcore        ",ev['mc_event']['xcore'])
        print("ycore        ",ev['mc_event']['ycore'])
        ev_time=[0,0,0,0]
        LSTID=ev['telescope_events'].keys()
        for i in LSTID :
            ev_time[int(i-1)] = ev['telescope_events'][i]['header']['readout_time']
        #
        print("ev_time_LST1 ",ev_time[0])
        print("ev_time_LST2 ",ev_time[1])
        print("ev_time_LST3 ",ev_time[2])
        print("ev_time_LST4 ",ev_time[3])
        #
        nphotons=[0,0,0,0]
        n_pe=[0,0,0,0]
        n_pixels=[0,0,0,0]
        #
        telNUM=ev['photoelectrons'].keys()
        for i in telNUM :
            nphotons[i]=len(ev['photons'][i])
            n_pe[i]=ev['photoelectrons'][i]['n_pe']
            n_pixels[i]=(ev['photoelectrons'][i]['n_pixels']-np.sum(ev['photoelectrons'][i]['photoelectrons']==0))
        #
        print("nphotons_LST1 ",nphotons[0])
        print("nphotons_LST2 ",nphotons[1])
        print("nphotons_LST3 ",nphotons[2])
        print("nphotons_LST4 ",nphotons[3])
        #
        print("n_pe_LST1 ",n_pe[0])
        print("n_pe_LST2 ",n_pe[1])
        print("n_pe_LST3 ",n_pe[2])
        print("n_pe_LST4 ",n_pe[3])
        #
        print("n_pixels_LST1 ",n_pixels[0])
        print("n_pixels_LST2 ",n_pixels[1])
        print("n_pixels_LST3 ",n_pixels[2])
        print("n_pixels_LST4 ",n_pixels[3])
        #
        if nevloop > 0 :
            break
        nevloop = nevloop + 1
        
    sf.close()

def loop_header_pe(datafilein = "/scratch/snx3000/lburmist/simtel_data/proton/data/corsika_run1.simtel.gz", max_ev = 100, headerfilename = 'header.csv',
                   pefilename_LST1 = 'pe_info_LST1.csv',
                   pefilename_LST2 = 'pe_info_LST2.csv',
                   pefilename_LST3 = 'pe_info_LST3.csv',
                   pefilename_LST4 = 'pe_info_LST4.csv'):
    #
    print("loop_header_pe")
    #
    sf = SimTelFile(datafilein)
    it_cout = 0
    it_cout_LST1 = 0
    it_cout_LST2 = 0
    it_cout_LST3 = 0
    it_cout_LST4 = 0    
    file_counter = 0;
    #
    tot_list=[]
    tot_arr_LST1=np.array([])
    tot_arr_LST2=np.array([])
    tot_arr_LST3=np.array([])
    tot_arr_LST4=np.array([])
    #
    tic = time.time()
    toc = time.time()
    #
    for ev in sf:
        #
        if (it_cout%1000==0) :
            toc = time.time()
            print('{:10d} {:10d} {:10.2f} s'.format(it_cout, ev['event_id'], toc - tic))
            tic = time.time()                    
        #
        ev_time=[0,0,0,0]
        nphotons=[0,0,0,0]
        n_pe=[0,0,0,0]
        n_pixels=[0,0,0,0]
        #
        telNUM=ev['photoelectrons'].keys()
        LSTID=ev['telescope_events'].keys()
        for i in LSTID :
            ev_time[(i-1)] = float(ev['telescope_events'][i]['header']['readout_time'])
        #
        for i in telNUM :
            nphotons[i]=int(len(ev['photons'][i]))
            n_pe[i]=int(ev['photoelectrons'][i]['n_pe'])
            n_pixels[i]=int(ev['photoelectrons'][i]['n_pixels']-np.sum(ev['photoelectrons'][i]['photoelectrons']==0))
        #
        #
            
        if (it_cout == 0):
            tot_list = []
        #
        #
        tot_list.append([ev['event_id'],
                         ev['mc_shower']['energy'],
                         ev['mc_shower']['azimuth'],
                         ev['mc_shower']['altitude'],
                         ev['mc_shower']['h_first_int'],
                         ev['mc_shower']['xmax'],
                         ev['mc_shower']['hmax'],
                         ev['mc_shower']['emax'],
                         ev['mc_shower']['cmax'],
                         ev['mc_event']['xcore'],
                         ev['mc_event']['ycore'],
                         ev_time[0],
                         ev_time[1],
                         ev_time[2],
                         ev_time[3],
                         nphotons[0],
                         nphotons[1],
                         nphotons[2],
                         nphotons[3],
                         n_pe[0],
                         n_pe[1],
                         n_pe[2],
                         n_pe[3],
                         n_pixels[0],
                         n_pixels[1],
                         n_pixels[2],
                         n_pixels[3]])
        #
        #
        #
        for i in range(len(n_pe)) :
            n_pe_per_tel = n_pe[i]
            if (n_pe_per_tel > 0) :
                npe = ev['photoelectrons'][i]['n_pe']
                new_ev = np.concatenate((np.ones((npe,1))*ev['event_id'],
                                         np.reshape(np.array(ev['photoelectrons'][i]['pixel_id']),(npe,1)),
                                         np.reshape(np.array(ev['photoelectrons'][i]['time']),(npe,1))), axis=1)
                #
                if (i == 0) :
                    if(it_cout_LST1 == 0):
                        tot_arr_LST1 = new_ev
                    else :
                        tot_arr_LST1 = np.concatenate((tot_arr_LST1,new_ev), axis=0)
                    #
                    it_cout_LST1 = it_cout_LST1 + 1
                #
                if (i == 1) :
                    if(it_cout_LST2 == 0):
                        tot_arr_LST2 = new_ev
                    else :
                        tot_arr_LST2 = np.concatenate((tot_arr_LST2,new_ev), axis=0)
                    #
                    it_cout_LST2 = it_cout_LST2 + 1
                #
                if (i == 2) :
                    if(it_cout_LST3 == 0):
                        tot_arr_LST3 = new_ev
                    else :
                        tot_arr_LST3 = np.concatenate((tot_arr_LST3,new_ev), axis=0)
                    #
                    it_cout_LST3 = it_cout_LST3 + 1
                #
                if (i == 3) :
                    if(it_cout_LST4 == 0):
                        tot_arr_LST4 = new_ev
                    else :
                        tot_arr_LST4 = np.concatenate((tot_arr_LST4,new_ev), axis=0)
                    #
                    it_cout_LST4 = it_cout_LST4 + 1
        #
        #
        it_cout = it_cout + 1

        #if(it_cout == 10000):
        #    pkl.dump(np.array(tot_list), open(str(headrefilename + "_" +str(file_counter)), "wb"), protocol=pkl.HIGHEST_PROTOCOL)    
        #    pkl.dump(tot_arr, open(str(pefilename + "_" +str(file_counter)), "wb"), protocol=pkl.HIGHEST_PROTOCOL)
        #    file_counter = file_counter + 1
        #
        #    it_cout = 0
        
        if (it_cout>=max_ev) :
            break
    #
    #
    #pkl.dump(np.array(tot_list), open(str(headrefilename + "_" +str(file_counter)), "wb"), protocol=pkl.HIGHEST_PROTOCOL)    
    #pkl.dump(tot_arr, open(str(pefilename + "_" +str(file_counter)), "wb"), protocol=pkl.HIGHEST_PROTOCOL)
    #
    #
    df_pe_LST1 = pd.DataFrame({'event_id': tot_arr_LST1[:,0], 
                               'pixel_id': tot_arr_LST1[:,1],
                               'time': tot_arr_LST1[:,2]})
    df_pe_LST1.to_csv(pefilename_LST1, sep=' ',header=False)
    #
    df_pe_LST2 = pd.DataFrame({'event_id': tot_arr_LST2[:,0], 
                               'pixel_id': tot_arr_LST2[:,1],
                               'time': tot_arr_LST2[:,2]})
    df_pe_LST2.to_csv(pefilename_LST2, sep=' ',header=False)
    #
    df_pe_LST3 = pd.DataFrame({'event_id': tot_arr_LST3[:,0], 
                               'pixel_id': tot_arr_LST3[:,1],
                               'time': tot_arr_LST3[:,2]})
    df_pe_LST3.to_csv(pefilename_LST3, sep=' ',header=False)
    #
    df_pe_LST4 = pd.DataFrame({'event_id': tot_arr_LST4[:,0], 
                               'pixel_id': tot_arr_LST4[:,1],
                               'time': tot_arr_LST4[:,2]})
    df_pe_LST4.to_csv(pefilename_LST4, sep=' ',header=False)
    #
    header=np.array(tot_list)
    df = pd.DataFrame({'event_id': header[:,0], 
                       'energy': header[:,1],
                       'azimuth': header[:,2],
                       'altitude': header[:,3],
                       'h_first_int': header[:,4],
                       'xmax': header[:,5],
                       'hmax': header[:,6],
                       'emax': header[:,7],
                       'cmax': header[:,8],
                       'xcore': header[:,9],
                       'ycore': header[:,10],                       
                       'ev_time_LST1': header[:,11],
                       'ev_time_LST2': header[:,12],
                       'ev_time_LST3': header[:,13],
                       'ev_time_LST4': header[:,14],
                       'nphotons_LST1': header[:,15],
                       'nphotons_LST2': header[:,16],
                       'nphotons_LST3': header[:,17],
                       'nphotons_LST4': header[:,18],
                       'n_pe_LST1': header[:,19],
                       'n_pe_LST2': header[:,20],
                       'n_pe_LST3': header[:,21],
                       'n_pe_LST4': header[:,22],
                       'n_pixels_LST1': header[:,23],
                       'n_pixels_LST2': header[:,24],
                       'n_pixels_LST3': header[:,25],
                       'n_pixels_LST4': header[:,26]})
    df.to_csv(headerfilename,sep=' ',header=False)
    #
    sf.close()
    
if __name__ == "__main__":

    #test()
    #loop_header_pe()
    
    if (len(sys.argv)==7):
        #
        #datafilein = "../scratch/data_nagaia/data/mono-lst-sipm-pmma-3ns-v1_triggerless/gamma_on_nsb_1x/output/corsika_run999.simtel.gz"
        #headerout = "corsika_run999.header.csv"
        #pe_info_out = "corsika_run999.pe_info.csv"
        #
        datafilein = str(sys.argv[1])
        headerout = str(sys.argv[2])
        pe_info_out_LST1 = str(sys.argv[3])
        pe_info_out_LST2 = str(sys.argv[4])
        pe_info_out_LST3 = str(sys.argv[5])
        pe_info_out_LST4 = str(sys.argv[6])
        #
        print("datafilein       = ", datafilein)
        print("headerout        = ", headerout)
        print("pe_info_out_LST1 = ", pe_info_out_LST1)
        print("pe_info_out_LST2 = ", pe_info_out_LST2)
        print("pe_info_out_LST3 = ", pe_info_out_LST3)
        print("pe_info_out_LST4 = ", pe_info_out_LST4)
        #
        tic = time.time()
        #test(datafilein)
        loop_header_pe(datafilein, 1000001, headerout, pe_info_out_LST1, pe_info_out_LST2, pe_info_out_LST3, pe_info_out_LST4)
        toc = time.time()
        print('{:.2f} s'.format(toc - tic))

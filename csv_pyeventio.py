from eventio import SimTelFile

import numpy as np
import pickle as pkl
import pandas as pd
import sys
import time

def loop_header_pe(datafilein = "./simtel_data/gamma/data/corsika_run307.simtel.gz", max_ev = 100, headerfilename = 'header.csv', pefilename = 'pe_info.csv'):
    #
    print("loop_header_pe")
    #
    sf = SimTelFile(datafilein)
    it_cout = 0
    file_counter = 0;
    tot_list=[]
    tot_arr=np.array([])
    #
    tic = time.time()
    toc = time.time()
    #
    for ev in sf:
        if (it_cout%1000==0) :
            toc = time.time()
            print('{:10d} {:10d} {:10.2f} s'.format(it_cout, ev['event_id'], toc - tic))
            tic = time.time()

        npe = ev['photoelectrons'][0]['n_pe']
        new_ev = np.concatenate((np.ones((npe,1))*ev['event_id'],
                                 np.reshape(np.array(ev['photoelectrons'][0]['pixel_id']),(npe,1)),
                                 np.reshape(np.array(ev['photoelectrons'][0]['time']),(npe,1))), axis=1)
        
        if(it_cout == 0):
            tot_arr = new_ev
            tot_list = []
        else :
            tot_arr = np.concatenate((tot_arr,new_ev), axis=0)                

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
                         ev['telescope_events'][1]['header']['readout_time'],
                         len(ev['photons'][0]),
                         ev['photoelectrons'][0]['n_pe'],
                        (ev['photoelectrons'][0]['n_pixels']-np.sum(ev['photoelectrons'][0]['photoelectrons']==0))])


        it_cout = it_cout + 1

        if(it_cout == 10000):
            pkl.dump(np.array(tot_list), open(str(headrefilename + "_" +str(file_counter)), "wb"), protocol=pkl.HIGHEST_PROTOCOL)    
            pkl.dump(tot_arr, open(str(pefilename + "_" +str(file_counter)), "wb"), protocol=pkl.HIGHEST_PROTOCOL)
            file_counter = file_counter + 1
            #
            it_cout = 0

        
        if (it_cout>=max_ev) :
            break
    #
    #
    #pkl.dump(np.array(tot_list), open(str(headrefilename + "_" +str(file_counter)), "wb"), protocol=pkl.HIGHEST_PROTOCOL)    
    #pkl.dump(tot_arr, open(str(pefilename + "_" +str(file_counter)), "wb"), protocol=pkl.HIGHEST_PROTOCOL)
    #
    #
    df_pe = pd.DataFrame({'event_id': tot_arr[:,0], 
                          'pixel_id': tot_arr[:,1],
                          'time': tot_arr[:,2]})
    df_pe.to_csv(pefilename, sep=' ',header=False)
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
                       'ev_time': header[:,11],
                       'nphotons': header[:,12],
                       'n_pe': header[:,13],
                       'n_pixels': header[:,14]})
    df.to_csv(headerfilename,sep=' ',header=False)
    #
    sf.close()
    
if __name__ == "__main__":

    if (len(sys.argv)==4):
        #
        #datafilein = "../scratch/data_nagaia/data/mono-lst-sipm-pmma-3ns-v1_triggerless/gamma_on_nsb_1x/output/corsika_run999.simtel.gz"
        #headerout = "corsika_run999.header.csv"
        #pe_info_out = "corsika_run999.pe_info.csv"
        #
        datafilein = str(sys.argv[1])
        headerout = str(sys.argv[2])
        pe_info_out = str(sys.argv[3])
        #
        print("datafilein  = ", datafilein)
        print("headerout   = ", headerout)
        print("pe_info_out = ", pe_info_out)
        #
        tic = time.time()
        #test(datafilein)
        loop_header_pe(datafilein, 1000001, headerout, pe_info_out)
        toc = time.time()
        print('{:.2f} s'.format(toc - tic))

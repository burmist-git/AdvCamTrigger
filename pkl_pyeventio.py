from eventio import SimTelFile
import numpy as np
import pickle as pkl

import time

def loop_wf(datafilein = "../simtel_data/gamma_diffuse/data/corsika_run307.simtel.gz", max_ev = 100, outfilename = 'wf_info'):
    #
    sf = SimTelFile(datafilein)
    it_cout = 0
    tot_list=[]
    #
    tic = time.time()
    toc = time.time()
    #
    file_counter=-1
    #
    for ev in sf :
        if (it_cout%1000==0) :
            toc = time.time()
            print('{:10d} {:10d} {:10.2f} s'.format(it_cout, ev['event_id'], toc - tic))
            tic = time.time()

            if ( file_counter > -1 ) :
                outfilename_i=str(outfilename+'_{0:05}'.format(file_counter)+'.pkl')
                pkl.dump(np.stack(tot_list, axis=2), open(outfilename_i, "wb"), protocol=pkl.HIGHEST_PROTOCOL)

            tot_list=[]
            file_counter = file_counter+1
            
        if (it_cout>max_ev) :
            break
        it_cout = it_cout + 1
        #
        tot_list.append(ev['telescope_events'][1]['adc_samples'][0])

    #print(tot_arr)

    if ( len(tot_list) > 0 ) :
        outfilename_i=str(outfilename+'_{0:05}'.format(file_counter)+'.pkl')
        pkl.dump(np.stack(tot_list, axis=2), open(outfilename_i, "wb"), protocol=pkl.HIGHEST_PROTOCOL)
    
    sf.close()
    

def loop_pe(datafilein = "../simtel_data/gamma_diffuse/data/corsika_run307.simtel.gz", max_ev = 100, headrefilename = 'pe_info.pkl'):
    #
    sf = SimTelFile(datafilein)
    it_cout = 0
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
        if (it_cout>max_ev) :
            break
        it_cout = it_cout + 1

        #
        #print("-->         "    , ev['photoelectrons'][0]['pixel_id'])
        #print("-->         "    , ev['photoelectrons'][0]['time'])

        npe = ev['photoelectrons'][0]['n_pe']

        #np.concatenate( (np.ones((npe,1))*ev['event_id'], np.array(ev['photoelectrons'][0]['pixel_id'])), axis=1)

        #print(type(np.ones((npe,1))*ev['event_id']))
        #print(type(np.array(ev['photoelectrons'][0]['pixel_id'])))
        #print((np.ones((npe,1))*ev['event_id']).shape)
        #print((np.reshape(np.array(ev['photoelectrons'][0]['pixel_id']),(npe,1))).shape)
        #print(np.concatenate((np.ones((npe,1))*ev['event_id'],
        #                      np.reshape(np.array(ev['photoelectrons'][0]['pixel_id']),(npe,1)),
        #                      np.reshape(np.array(ev['photoelectrons'][0]['time']),(npe,1))), axis=1))
        #
        #
        new_ev = np.concatenate((np.ones((npe,1))*ev['event_id'],
                                 np.reshape(np.array(ev['photoelectrons'][0]['pixel_id']),(npe,1)),
                                 np.reshape(np.array(ev['photoelectrons'][0]['time']),(npe,1))), axis=1)
        if(it_cout == 1):
            tot_arr = new_ev
        else :
            tot_arr = np.concatenate((tot_arr,new_ev), axis=0)                

    #print(tot_arr)
    
    pkl.dump(tot_arr, open(headrefilename, "wb"), protocol=pkl.HIGHEST_PROTOCOL)

    sf.close()

def loop_header(datafilein = "../simtel_data/gamma_diffuse/data/corsika_run307.simtel.gz", max_ev = 100, headrefilename = 'header.pkl'):
    #
    sf = SimTelFile(datafilein)
    it_cout = 0
    tot_list=[]
    #
    tic = time.time()
    toc = time.time()
    #
    for ev in sf:
        if (it_cout%1000==0) :
            toc = time.time()
            print('{:10d} {:10d} {:10.2f} s'.format(it_cout, ev['event_id'], toc - tic))
            tic = time.time()
        if (it_cout>max_ev) :
            break
        it_cout = it_cout + 1
        #
        #print("----------------------------------")
        #print("event_id         ", ev['event_id'])
        #print("energy           ", ev['mc_shower']['energy'])
        #print("xcore            ", ev['mc_event']['xcore'])
        #print("ycore            ", ev['mc_event']['ycore'])
        #print("ev_time          ", ev['telescope_events'][1]['header']['readout_time'])
        #print("nphotons         ", len(ev['photons'][0]))
        #print("n_pe             ", ev['photoelectrons'][0]['n_pe'])
        #print("n_pixels         ", (ev['photoelectrons'][0]['n_pixels']-np.sum(ev['photoelectrons'][0]['photoelectrons']==0)))
        #print("----------------------------------")
        tot_list.append([ev['event_id'],
                         ev['mc_shower']['energy'],
                         ev['mc_event']['xcore'],
                         ev['mc_event']['ycore'],
                         ev['telescope_events'][1]['header']['readout_time'],
                         len(ev['photons'][0]),
                         ev['photoelectrons'][0]['n_pe'],
                         (ev['photoelectrons'][0]['n_pixels']-np.sum(ev['photoelectrons'][0]['photoelectrons']==0))])
        
    pkl.dump(np.array(tot_list), open(headrefilename, "wb"), protocol=pkl.HIGHEST_PROTOCOL)    

    sf.close()
    
if __name__ == "__main__":
    #
    #datafilein = "../simtel_data/gamma/data/corsika_run307.simtel.gz"
    #headerout = "../compressed_data/gamma/corsika_run307.header.pkl"
    #pe_info_out = "../compressed_data/gamma/corsika_run307.pe_info.pkl"
    #wf_info_out = "../compressed_data/gamma/corsika_run307.wf_info"
    #
    datafilein = "../simtel_data/no_nsb_cut/gamma/data/corsika_run307.simtel.gz"
    headerout = "../compressed_data/no_nsb_cut/gamma/corsika_run307.header.pkl"
    pe_info_out = "../compressed_data/no_nsb_cut/gamma/corsika_run307.pe_info.pkl"
    wf_info_out = "../compressed_data/no_nsb_cut/gamma/corsika_run307.wf_info"
    #
    #datafilein = "../simtel_data/gamma_diffuse/data/corsika_run307.simtel.gz"
    #datafilein = "../simtel_data/electron/data/corsika_run307.simtel.gz"
    #datafilein = "../simtel_data/proton/data/corsika_run307.simtel.gz"
    tic = time.time()
    loop_header( datafilein, 100000, headerout)
    loop_pe( datafilein, 100000, pe_info_out)
    loop_wf( datafilein, 100000, wf_info_out)
    toc = time.time()
    print('{:.2f} s'.format(toc - tic))

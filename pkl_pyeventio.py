from eventio import SimTelFile

import numpy as np
import pickle as pkl
import sys
import time

def loop_wf(datafilein = "../simtel_data/gamma/data/corsika_run307.simtel.gz", max_ev = 100, outfilename = 'wf_info', csv_flag=True):
    #
    sf = SimTelFile(datafilein)
    tot_list=()
    #
    tic = time.time()
    toc = time.time()
    #
    it_cout = 0
    file_counter=0
    first_evne_in_file=True
    #
    for ev in sf :
        if (it_cout%1000==0) :
            toc = time.time()
            print('{:10d} {:10d} {:10.2f} s'.format(it_cout, ev['event_id'], toc - tic))
            tic = time.time()
            
            if ( file_counter > 0 ) :
                outfilename_i=str(outfilename+'_{0:05}'.format(int(file_counter-1)) + '.pkl')
                outfilename_csv_i=str(outfilename+'_{0:05}'.format(int(file_counter-1)) + '.csv')
                data_to_save=np.concatenate(tot_list)
                pkl.dump(data_to_save, open(outfilename_i, "wb"), protocol=pkl.HIGHEST_PROTOCOL)
                if csv_flag :
                    np.savetxt(outfilename_csv_i, data_to_save, delimiter=' ',fmt='%i')

                
            tot_list=(ev['telescope_events'][1]['adc_samples'][0],)

            file_counter = file_counter+1
            first_evne_in_file = True
            

        if first_evne_in_file == False :
            tot_list = tot_list + (ev['telescope_events'][1]['adc_samples'][0],)

        it_cout = it_cout + 1
        if (it_cout>=max_ev) :
            break

        first_evne_in_file = False 
            
        #
        #print(ev['telescope_events'][1]['adc_samples'][0].shape)
        #print(type(ev['telescope_events'][1]['adc_samples'][0]))
        #break
        
    #print(tot_arr)

    if ( len(tot_list) > 0 ) :
        outfilename_i=str(outfilename+'_{0:05}'.format(int(file_counter-1))+'.pkl')
        outfilename_csv_i=str(outfilename+'_{0:05}'.format(int(file_counter-1)) + '.csv')
        data_to_save=np.concatenate(tot_list)
        pkl.dump(data_to_save, open(outfilename_i, "wb"), protocol=pkl.HIGHEST_PROTOCOL)
        if csv_flag :
            np.savetxt(outfilename_csv_i, data_to_save, delimiter=' ',fmt='%i')
    
    sf.close()

def loop_wf_stack(datafilein = "../simtel_data/gamma/data/corsika_run307.simtel.gz", max_ev = 100, outfilename = 'wf_info'):
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
                outfilename_i=str(outfilename+'_{0:05}'.format(file_counter) + '.pkl')
                pkl.dump(np.stack(tot_list, axis=2), open(outfilename_i, "wb"), protocol=pkl.HIGHEST_PROTOCOL)

            tot_list=[]
            file_counter = file_counter+1
            
        it_cout = it_cout + 1
        if (it_cout>=max_ev) :
            break
        #
        #print(ev['telescope_events'][1]['adc_samples'][0].shape)
        #print(type(ev['telescope_events'][1]['adc_samples'][0]))
        tot_list.append(ev['telescope_events'][1]['adc_samples'][0])
        #break
        
    #print(tot_arr)

    if ( len(tot_list) > 0 ) :
        outfilename_i=str(outfilename+'_{0:05}'.format(file_counter)+'.pkl')
        pkl.dump(np.stack(tot_list, axis=2), open(outfilename_i, "wb"), protocol=pkl.HIGHEST_PROTOCOL)
    
    sf.close()
    

def loop_pe(datafilein = "../simtel_data/gamma/data/corsika_run307.simtel.gz", max_ev = 100, headrefilename = 'pe_info.pkl'):
    #
    print("loop_pe")
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
        if(it_cout == 0):
            tot_arr = new_ev
        else :
            tot_arr = np.concatenate((tot_arr,new_ev), axis=0)                

        it_cout = it_cout + 1
        if (it_cout>=max_ev) :
            break

            
    #print(tot_arr)
    
    pkl.dump(tot_arr, open(headrefilename, "wb"), protocol=pkl.HIGHEST_PROTOCOL)

    sf.close()

def loop_header(datafilein = "../simtel_data/gamma/data/corsika_run307.simtel.gz", max_ev = 100, headrefilename = 'header.pkl'):
    #
    print("loop_header")
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
        #'azimuth' 3.1415927410125732, 
        #'altitude' 1.2217304706573486, 
        #'h_first_int' 26537.1640625,
        #'xmax' 256.582275390625,
        #'hmax' 10661.2880859375,
        #'emax' 256.582275390625,
        #'cmax' 262.06817626953125,

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
        if (it_cout>=max_ev) :
            break
        
    pkl.dump(np.array(tot_list), open(headrefilename, "wb"), protocol=pkl.HIGHEST_PROTOCOL)    
    
    sf.close()


def loop_header_pe(datafilein = "../simtel_data/gamma/data/corsika_run307.simtel.gz", max_ev = 100, headrefilename = 'header.pkl', pefilename = 'pe_info.pkl'):
    #
    print("loop_header_pe")
    #
    sf = SimTelFile(datafilein)
    it_cout = 0
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

        npe = ev['photoelectrons'][0]['n_pe']
        new_ev = np.concatenate((np.ones((npe,1))*ev['event_id'],
                                 np.reshape(np.array(ev['photoelectrons'][0]['pixel_id']),(npe,1)),
                                 np.reshape(np.array(ev['photoelectrons'][0]['time']),(npe,1))), axis=1)
        
        if(it_cout == 0):
            tot_arr = new_ev
        else :
            tot_arr = np.concatenate((tot_arr,new_ev), axis=0)                
        
        it_cout = it_cout + 1
        if (it_cout>=max_ev) :
            break

    pkl.dump(np.array(tot_list), open(headrefilename, "wb"), protocol=pkl.HIGHEST_PROTOCOL)    
    pkl.dump(tot_arr, open(pefilename, "wb"), protocol=pkl.HIGHEST_PROTOCOL)
    
    sf.close()


def get_pixel_mapping(datafilein = "../simtel_data/gamma/data/corsika_run307.simtel.gz", outmap_csv = 'pixel_mapping.csv'):
    sf = SimTelFile(datafilein)
    #
    n_pixels=float(sf.telescope_descriptions[1]['camera_organization']['n_pixels'])
    n_drawers=float(sf.telescope_descriptions[1]['camera_organization']['n_drawers'])
    pixel_size=float(sf.telescope_descriptions[1]['camera_settings']['pixel_size'][0])
    #
    the_map=np.concatenate((sf.telescope_descriptions[1]['camera_settings']['pixel_x'].reshape(int(n_pixels),1),
                            sf.telescope_descriptions[1]['camera_settings']['pixel_y'].reshape(int(n_pixels),1),
                            sf.telescope_descriptions[1]['camera_organization']['drawer'].reshape(int(n_pixels),1)), axis=1)
    np.savetxt(outmap_csv, the_map, delimiter=' ',fmt='%f')
    #
    print('n_pixels   = ', int(n_pixels))
    print('n_drawers  = ', int(n_drawers))
    print('pixel_size = ', pixel_size)
    #
    # 0.024300
    # 0.023300
    #
    #
    sf.close()
    
if __name__ == "__main__":

    if (len(sys.argv)==4):
        filenameNpz = str(sys.argv[1])
        filenameRoot = str(sys.argv[2])
        #
        #datafilein = "../scratch/data_nagaia/data/mono-lst-sipm-pmma-3ns-v1_triggerless/gamma_on_nsb_1x/output/corsika_run999.simtel.gz"
        #headerout = "corsika_run999.header.pkl"
        #pe_info_out = "corsika_run999.pe_info.pkl"
        datafilein = str(sys.argv[1])
        headerout = str(sys.argv[2])
        pe_info_out = str(sys.argv[3])
        #
        print("datafilein  = ", datafilein)
        print("headerout   = ", headerout)
        print("pe_info_out = ", pe_info_out)
        #
        #datafilein = "../simtel_data/gamma/data/corsika_run307.simtel.gz"
        #headerout = "../compressed_data/gamma/corsika_run307.header.pkl"
        #pe_info_out = "../compressed_data/gamma/corsika_run307.pe_info.pkl"
        #wf_info_out = "../compressed_data/gamma/corsika_run307.wf_info"
        #
        #datafilein = "../simtel_data/no_nsb_cut/gamma/data/corsika_run307.simtel.gz"
        #headerout = "../compressed_data/no_nsb_cut/gamma/corsika_run307.header.pkl"
        #pe_info_out = "../compressed_data/no_nsb_cut/gamma/corsika_run307.pe_info.pkl"
        #wf_info_out = "../compressed_data/no_nsb_cut/gamma/corsika_run307.wf_info"
        #
        #datafilein = "../simtel_data/no_nsb_cut/proton/data/corsika_run307.simtel.gz"
        #headerout = "../compressed_data/no_nsb_cut/proton/corsika_run307.header.pkl"
        #pe_info_out = "../compressed_data/no_nsb_cut/proton/corsika_run307.pe_info.pkl"
        #wf_info_out = "../compressed_data/no_nsb_cut/proton/corsika_run307.wf_info"
        #
        #datafilein = "../simtel_data/gamma/data/corsika_run307.simtel.gz"
        #datafilein = "../simtel_data/electron/data/corsika_run307.simtel.gz"
        #datafilein = "../simtel_data/proton/data/corsika_run307.simtel.gz"
        #
        tic = time.time()
        loop_header_pe(datafilein, 1000001, headerout, pe_info_out)
        #loop_header( datafilein, 10000, headerout)
        #loop_pe( datafilein, 10000, pe_info_out)
        #loop_wf_stack( datafilein, 100000, wf_info_out)
        #loop_wf( datafilein, 500000, wf_info_out, True)
        #
        #loop_wf_test( datafilein, 100000, wf_info_out)
        #
        #get_pixel_mapping(datafilein, outmap_csv = 'pixel_mapping.csv')
        toc = time.time()
        print('{:.2f} s'.format(toc - tic))

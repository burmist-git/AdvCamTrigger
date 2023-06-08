import numpy as np
import pickle as pkl
import pandas as pd
import time
import subprocess

def conv_wf(fin):
    print(fin)
    wf_arr = pkl.load(open(fin, 'rb'))
    old_shape=wf_arr.shape
    df=pd.DataFrame(wf_arr.reshape((old_shape[0]*old_shape[2],old_shape[1])),
                    columns=[str(j) for j in range(old_shape[1])])
    csvfileout=str(fin[:-3]+'csv')
    print(csvfileout)
    df.to_csv(csvfileout,sep=' ',header=False)
    
def conv_pe(pklfilein = '../compressed_data/no_nsb_cut/gamma/corsika_run307.pe_info.pkl',
            csvfileout = '../compressed_data/no_nsb_cut/gamma/corsika_run307.pe_info.csv'):
    #
    pe_inf = pkl.load(open(pklfilein, 'rb'))
    df = pd.DataFrame({'event_id': pe_inf[:,0], 
                       'pixel_id': pe_inf[:,1],
                       'time': pe_inf[:,2]})
    df.to_csv(csvfileout,sep=' ',header=False)
    
def conv_header(pklfilein = '../compressed_data/no_nsb_cut/gamma/corsika_run307.header.pkl',
                csvfileout = '../compressed_data/no_nsb_cut/gamma/corsika_run307.header.csv'):
    #
    header = pkl.load(open(pklfilein, 'rb'))
    df = pd.DataFrame({'event_id': header[:,0], 
                       'energy': header[:,1],
                       'xcore': header[:,2],
                       'ycore': header[:,3],
                       'ev_time': header[:,4],
                       'nphotons': header[:,5],
                       'n_pe': header[:,6],
                       'n_pixels': header[:,7]})
    df.to_csv(csvfileout,sep=' ',header=False)
    
if __name__ == "__main__":
    #
    #header = pkl.load(open('../compressed_data/gamma/corsika_run307.header.pkl', 'rb'))
    #pe_inf = pkl.load(open('../compressed_data/gamma/corsika_run307.pe_info.pkl', 'rb'))
    #wf_inf = pkl.load(open('../compressed_data/gamma/corsika_run307.wf_info_00000.pkl', 'rb'))
    #pe_inf = pkl.load(open('../compressed_data/no_nsb_cut/gamma/corsika_run307.pe_info.pkl', 'rb'))
    #wf_inf = pkl.load(open('../compressed_data/no_nsb_cut/gamma/corsika_run307.wf_info_00000.pkl', 'rb'))
    #
    #fadc_sum_offset=15
    #fadc_amplitude=8.25
    #fadc_MHz=1024
    #fadc_sample_in_ns=1000.0/fadc_MHz
    #time_offset=fadc_sum_offset*fadc_sample_in_ns
    #fadc_bins=75
    #
    #print('fadc_sum_offset   = ', fadc_sum_offset)
    #print('fadc_amplitude    = ', fadc_amplitude)
    #print('fadc_MHz          = ', fadc_MHz)
    #print('fadc_sample_in_ns = ', fadc_sample_in_ns)
    #print('fadc_bins         = ', fadc_bins)
    #print('time_offset       = ', time_offset)
    #
    #wf_times=np.array(range(fadc_bins))*1000.0/1024.0    
    #
    tic = time.time()
    #conv_header(pklfilein = '../compressed_data/no_nsb_cut/gamma/corsika_run307.header.pkl',
    #            csvfileout = '../compressed_data/no_nsb_cut/gamma/corsika_run307.header.csv')
    conv_header(pklfilein = '../compressed_data/gamma/corsika_run307.header.pkl',
                csvfileout = '../compressed_data/gamma/corsika_run307.header.csv')
    #
    #conv_pe(pklfilein = '../compressed_data/no_nsb_cut/gamma/corsika_run307.pe_info.pkl',
    #        csvfileout = '../compressed_data/no_nsb_cut/gamma/corsika_run307.pe_info.csv')
    conv_pe(pklfilein = '../compressed_data/gamma/corsika_run307.pe_info.pkl',
            csvfileout = '../compressed_data/gamma/corsika_run307.pe_info.csv')
    #
    #wf_pkl_list=['../compressed_data/no_nsb_cut/gamma/corsika_run307.wf_info_00000.pkl',
    #             '../compressed_data/no_nsb_cut/gamma/corsika_run307.wf_info_00001.pkl',
    #             '../compressed_data/no_nsb_cut/gamma/corsika_run307.wf_info_00002.pkl']
    #wf_pkl_list=(subprocess.getoutput("ls ../compressed_data/no_nsb_cut/gamma/corsika_run307.wf_info_*.pkl")).split('\n')
    #[conv_wf(thefile) for thefile in wf_pkl_list]
    #
    toc = time.time()
    print('{:.2f} s'.format(toc - tic))

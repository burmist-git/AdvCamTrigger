#!/bin/bash

#'nightsky_background=all:0.001'
#inRootFile="../compressed_data/no_nsb_cut/gamma/corsika_run307.compressed.root"
#outHistF="./hist_corsika_run307_gamma_no_nsb_cut.root"
#inRootFile="../compressed_data/no_nsb_cut/proton/corsika_run307.compressed.root"
#outHistF="./hist_corsika_run307_proton_no_nsb_cut.root"
#'nightsky_background=all:0.386'
##
##inRootFile="../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/gamma_on_nsb_1x/root/0000/corsika_0000ID.root"
##outHistF="./hist_gamma_on_corsika_0000ID.root"
##
inRootFile="gamma_on_nsb_1x.list"
outHistF="./hist_gamma_on_tw.root"
#inRootFile="gamma_diffuse_nsb_1x.list"
#inRootFile="gamma_diffuse_nsb_1x.listshort"
#outHistF="./hist_gamma_diffuse_nsb_1x_tw.root"
#inRootFile="proton_nsb_1x.list"
#outHistF="./hist_proton_nsb_1x_tw.root"

#short
#'nightsky_background=all:0.001'
#inRootFile="../compressed_data/no_nsb_cut/gamma/corsika_run307.compressed.short.root"
#outHistF="./hist_short_corsika_run307_no_nsb_cut.root"
#'nightsky_background=all:0.386'
#inRootFile="../compressed_data/gamma/corsika_run307.compressed.root"
#outHistF="./hist_corsika_run307.root"

#make -f Makefileana clean; make -f Makefileana runana;
#make clean; make -f Makefileana clean ; make -j; make -f Makefileana -j

function printHelp {
    echo " --> ERROR in input arguments "
    echo " [0] -d     : single root file"
    echo " [0] -SiPM  : get SiPM distr."
    echo " [0] -ds    : single root file (short format)"
    echo " [0] -tw    : test waveforms"
    echo " [0] -th    : test cam sipmCameraHist"
    echo " [0] -th_flowerpixID     : flower pixID (out csv file with mapping out csv file with mapping)"
    echo " [0] -th_isolated_flower : test of camera hist: test_trigger_channel_mask_isolated_flower"
    echo " [0] -th_isolated_flower_plus_super : test of camera hist: test_trigger_channel_mask_isolated_flower_plus_super"
    echo " [0] -get_isolated_flower_mask : save trigger channel mask isolated flower"
    echo " [0] -sg    : gamma"
    echo " [0] -sgd   : gamma diffuse"
    echo " [0] -se    : electron"
    echo " [0] -sp    : proton"
    echo " [0] -fg    : fast gamma"
    echo " [0] -fgd   : fast gamma diffuse"
    echo " [0] -fe    : fast electron"
    echo " [0] -fp    : fast proton"
    echo " [0] -fall  : fast all"
    echo " [0] -sfall : super fast all"
    echo " [0] -pg    : PCA gamma"
    echo " [0] -pgd   : PCA gamma diffuse"
    echo " [0] -pp    : PCA proton"
    echo " [0] -pg_d  : Draw PCA gamma"
    echo " [0] -pgd_d : Draw PCA gamma diffuse"
    echo " [0] -pp_d  : Draw PCA proton"
    echo " [0] -pg_d_reco  : Draw reco. showers (PCA) gamma"
    echo " [0] -pgd_d_reco : Draw reco. showers (PCA) gamma diffuse"
    echo " [0] -pp_d_reco  : Draw reco. showers (PCA) proton"
    echo " [0] -trgp  : Trg proton"
    echo " [0] -trgg  : Trg gamma"
    echo " [0] -trggd : Trg gamma diffuse"
    echo " [0] -trgAg    : TrgA gamma"
    echo " [0] -trgAgd   : TrgA gamma diffuse"
    echo " [0] -trgAp    : TrgA proton"
    echo " [0] -trgA_NGB : TrgA NGB"
    echo " [0] -fscang        : fast scan gamma"
    echo " [0] -fscangd       : fast scan gamma diffuse"
    echo " [0] -fscane        : fast scan electron"
    echo " [0] -fscanp        : fast scan proton"
    echo " [0] -test_pe       : test single pe amplitude generator"
    echo " [0] -test_evstHist : test of evstHist"
    echo " [0] -rate_proton_binwise : calculate proton rate (binwise)"
    echo " [0] -rate_proton         : calculate proton rate"
    echo " [0] -rate_gamma          : calculate gamma rate"
    echo " [0] -EvPerEv             : event per event analysis"
    echo " [0] -evMAP               : event map"
    echo " [0] -c     : recompile"
    echo " [0] -h     : print help"
}

if [ $# -eq 0 ] 
then    
    printHelp
else
    if [ "$1" = "-d" ]; then
	inRootFile="../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/proton_nsb_1x/root/0000/corsika_0000ID.root"
	outHistF="./hist_corsika_0000ID.root"
	./runana 1 $inRootFile $outHistF
    elif [ "$1" = "-ds" ]; then
	./runana 3 $inRootFile $outHistF
    elif [ "$1" = "-trgp" ]; then
	inRootFile="../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/proton_nsb_1x/root/0000/corsika_0000ID.root"
	outHistF="./hist_trg_proton_nsb_1x_corsika_0000ID.root"
	./runana 11 $inRootFile $outHistF
    elif [ "$1" = "-trgg" ]; then
	inRootFile="../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/gamma_on_nsb_1x/root/0000/corsika_0000ID.root"
	outHistF="./hist_trg_gamma_on_nsb_1x_corsika_0000ID.root"
	./runana 11 $inRootFile $outHistF
    elif [ "$1" = "-trggd" ]; then
	inRootFile="../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/gamma_diffuse_nsb_1x/root/0000/corsika_0000ID.root"
	outHistF="./hist_trg_gamma_diffuse_nsb_1x_corsika_0000ID.root"
	./runana 11 $inRootFile $outHistF
    elif [ "$1" = "-SiPM" ]; then
	inRootFile="../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/gamma_diffuse_nsb_1x/root/0000/corsika_0000ID.root"
	outHistF="./hist_SiPM.root"
	./runana 100 $inRootFile $outHistF
    elif [ "$1" = "-trgAg" ]; then
	################################################
	#  E                                           #
	################################################
	#  0                    1              1.58489 #
	#  1              1.58489              2.51189 #
	#  2              2.51189              3.98107 #
	#  3              3.98107              6.30957 #
	#  4              6.30957                   10 #
	#  5                   10              15.8489 #
	#  6              15.8489              25.1189 #
	#  7              25.1189              39.8107 #
	#  8              39.8107              63.0957 #
	#  9              63.0957                  100 #
	# 10                  100              158.489 #
	################################################
	# theta                                        #
	################################################
	#  0                    0                    1 #
	#  1                    1                    2 #
	#  2                    2                    3 #
	#  3                    3                    4 #
	#  4                    4                    5 #
	#  5                    5                    6 #
	################################################
	#  r                                           #
	################################################
	#  0                    0                  150 #
	#  1                  150                  250 #
	#  2                  250                  350 #
	#  3                  350                  500 #
	#  4                  500                  650 #
	#  5                  650                  800 #
	#  6                  800                 1000 #
	#  7                 1000                 1300 #
	#  8                 1300                 1600 #
	#  9                 1600                 2000 #
	################################################
	inRootFile="../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/gamma_on_nsb_1x/root/0000/corsika_0000ID.root"
	outHistF="./hist_trg_gamma_on_nsb_1x_corsika_0000ID.root"
	binE=5
	binTheta=0
	binDist=0
	npe_min=40
	npe_max=200
	nEv_max=100
	rndseed=12312
	./runana 111 $inRootFile $outHistF $binE $binTheta $binDist $npe_min $npe_max $nEv_max $rndseed
    elif [ "$1" = "-trgAgd" ]; then
	inRootFile="../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/gamma_diffuse_nsb_1x/root/0000/corsika_0000ID.root"
	outHistF="./hist_trg_gamma_diffuse_nsb_1x_corsika_0000ID.root"
	binE=5
	binTheta=1
	binDist=2
	npe_min=40
	npe_max=200
	nEv_max=100	
	rndseed=12312
	./runana 111 $inRootFile $outHistF $binE $binTheta $binDist $npe_min $npe_max $nEv_max $rndseed
    elif [ "$1" = "-trgAp" ]; then
	inRootFile="../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/proton_nsb_1x/root/0000/corsika_0000ID.root"
	outHistF="./hist_trg_proton_nsb_1x_corsika_0000ID.root"
	binE=5
	binTheta=1
	binDist=1
	npe_min=40
	npe_max=200
	nEv_max=100
	rndseed=12312
	./runana 111 $inRootFile $outHistF $binE $binTheta $binDist $npe_min $npe_max $nEv_max $rndseed
    elif [ "$1" = "-trgA_NGB" ]; then
	inRootFile="../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/proton_nsb_1x/root/0000/corsika_0000ID.root"
	outHistF="./hist_trg_NGB.root"
	binE=5
	binTheta=5
	binDist=7
	npe_min=1
	npe_max=5
	nEv_max=100
	rndseed=12312
	./runana 111 $inRootFile $outHistF $binE $binTheta $binDist $npe_min $npe_max $nEv_max $rndseed
    elif [ "$1" = "-test_pe" ]; then
	inRootFile="../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/proton_nsb_1x/root/0000/corsika_0000ID.root"
	outHistF="./hist_test_single_pe_amplitude_generator.root"
	./runana 1111 $inRootFile $outHistF
    elif [ "$1" = "-tw" ]; then
	#./runana 4 $inRootFile $outHistF 35489 gamma # n_pe 80228  energy 3506.29 GeV
	#./runana 4 $inRootFile $outHistF 127417 gamma # n_pe 261729  energy 11278.8 GeV
	#./runana 4 $inRootFile $outHistF 144155 gamma # n_pe 251240  energy 11119.8 GeV
	#./runana 4 $inRootFile $outHistF 177092 gamma # n_pe 124584  energy 6220.12 GeV
	#./runana 4 $inRootFile $outHistF 218806 gamma # n_pe 55170  energy 3034.07 GeV
	#./runana 4 $inRootFile $outHistF 287758 gamma # n_pe 43400  energy 1676.92 GeV
	#./runana 4 $inRootFile $outHistF 58871 gamma # n_pe 24469  energy 1710.38 GeV
	#./runana 4 $inRootFile $outHistF 7677 gamma # n_pe 283606  energy 13327.7 GeV
	
	## azimuth  = 180 +/- 0.1
        ## altitude = 70  +/- 0.1
	#./runana 4 $inRootFile $outHistF 3192050  gamma # 4580
	#./runana 4 $inRootFile $outHistF 6792462  gamma # 3594
	#./runana 4 $inRootFile $outHistF 12118108 gamma # 5648 

	## azimuth  = 181 +/- 0.1
        ## altitude = 70  +/- 0.1
	#./runana 4 $inRootFile $outHistF 1724488  gamma # 5166
	#./runana 4 $inRootFile $outHistF 16219735 gamma # 3390 
	#./runana 4 $inRootFile $outHistF 27869493 gamma # 3522

	## azimuth  = 179 +/- 0.1
        ## altitude = 70  +/- 0.1
	#./runana 4 $inRootFile $outHistF 1979525 gamma  # 7011
	#./runana 4 $inRootFile $outHistF 5389323 gamma  # 5048
	#./runana 4 $inRootFile $outHistF 17040599 gamma # 4762

	## azimuth  = 177 +/- 0.1
        ## altitude = 70  +/- 0.1
	#./runana 4 $inRootFile $outHistF 2023749  gamma  # 3093
	#./runana 4 $inRootFile $outHistF 17269715 gamma  # 4986
	#./runana 4 $inRootFile $outHistF 35291254 gamma  # 4131

	## azimuth  = 170 +/- 0.1
        ## altitude = 70  +/- 0.1
	#./runana 4 $inRootFile $outHistF 524445 gamma    # 4104
	#./runana 4 $inRootFile $outHistF 20449567 gamma  # 4046
	#./runana 4 $inRootFile $outHistF 31214745 gamma  # 3717

	## azimuth  = 180 +/- 0.2
        ## altitude = 70  +/- 0.2
	#./runana 4 $inRootFile $outHistF 576379 gamma_diff        
	#./runana 4 $inRootFile $outHistF 783671 gamma_diff
	#./runana 4 $inRootFile $outHistF 1218953 gamma_diff	

	## azimuth  = 177 +/- 0.2
        ## altitude = 70  +/- 0.2
	#./runana 4 $inRootFile $outHistF 2301553 gamma_diff        
	#./runana 4 $inRootFile $outHistF 646949 gamma_diff
	#./runana 4 $inRootFile $outHistF 5712136 gamma_diff	

	## azimuth  = 174 +/- 0.2
        ## altitude = 70  +/- 0.2
	#./runana 4 $inRootFile $outHistF 474353 gamma_diff        
	#./runana 4 $inRootFile $outHistF 4529997 gamma_diff
	#./runana 4 $inRootFile $outHistF 7611489 gamma_diff	

	## azimuth  = 180 +/- 0.2
        ## altitude = 72  +/- 0.2
	#./runana 4 $inRootFile $outHistF 439530 gamma_diff        
	#./runana 4 $inRootFile $outHistF 2069890 gamma_diff
	#./runana 4 $inRootFile $outHistF 2481499 gamma_diff	

	## azimuth  = 180 +/- 0.2
        ## altitude = 68  +/- 0.2
	#./runana 4 $inRootFile $outHistF 1768962 gamma_diff        
	#./runana 4 $inRootFile $outHistF 2362795 gamma_diff
	#./runana 4 $inRootFile $outHistF 9353299 gamma_diff

	## azimuth  = 180 +/- 0.2
        ## altitude = 69  +/- 0.2
	#./runana 4 $inRootFile $outHistF 135723 gamma_diff        
	#./runana 4 $inRootFile $outHistF 845310 gamma_diff
	#./runana 4 $inRootFile $outHistF 2378773 gamma_diff
   
	#./runana 4 $inRootFile $outHistF 896 gamma    # 162 pe
	#./runana 4 $inRootFile $outHistF 929 gamma    # 220 pe
	#./runana 4 $inRootFile $outHistF 3140 gamma   # 180 pe
	#./runana 4 $inRootFile $outHistF 3139 gamma   # 210 pe
	#
	#./runana 4 $inRootFile $outHistF 39807 proton
	#./runana 4 $inRootFile $outHistF 66724 proton
	#./runana 4 $inRootFile $outHistF 71559 proton
	#./runana 4 $inRootFile $outHistF 124790 proton
	#./runana 4 $inRootFile $outHistF 143047 proton
	#./runana 4 $inRootFile $outHistF 145282 proton
	#./runana 4 $inRootFile $outHistF 154835 proton
	#./runana 4 $inRootFile $outHistF 338715 proton
	#./runana 4 $inRootFile $outHistF 380672 proton
	#./runana 4 $inRootFile $outHistF 499225 proton
	#./runana 4 $inRootFile $outHistF 628655 proton
	#./runana 4 $inRootFile $outHistF 663307 proton
	#./runana 4 $inRootFile $outHistF 714851 proton
	#./runana 4 $inRootFile $outHistF 736295 proton
	#./runana 4 $inRootFile $outHistF 762120 proton
	#./runana 4 $inRootFile $outHistF 809713 proton
	#./runana 4 $inRootFile $outHistF 890894 proton
	#
	#./runana 4 $inRootFile $outHistF 179699 proton
	#./runana 4 $inRootFile $outHistF 184478 proton
	#./runana 4 $inRootFile $outHistF 184579 proton
	#./runana 4 $inRootFile $outHistF 185077 proton
	#./runana 4 $inRootFile $outHistF 187373 proton
	#./runana 4 $inRootFile $outHistF 193092 proton
	#./runana 4 $inRootFile $outHistF 198048 proton
	#./runana 4 $inRootFile $outHistF 204414 proton
	#./runana 4 $inRootFile $outHistF 206725 proton
	#./runana 4 $inRootFile $outHistF 208503 proton
	#./runana 4 $inRootFile $outHistF 211130 proton
	#./runana 4 $inRootFile $outHistF 212463 proton
	#./runana 4 $inRootFile $outHistF 212589 proton
	#./runana 4 $inRootFile $outHistF 216551 proton
	#./runana 4 $inRootFile $outHistF 219463 proton
	#./runana 4 $inRootFile $outHistF 222736 proton
	#./runana 4 $inRootFile $outHistF 231021 proton
	#./runana 4 $inRootFile $outHistF 232611 proton
	#./runana 4 $inRootFile $outHistF 238618 proton
	#./runana 4 $inRootFile $outHistF 240388 proton
	#./runana 4 $inRootFile $outHistF 242426 proton
	#./runana 4 $inRootFile $outHistF 249587 proton
	#./runana 4 $inRootFile $outHistF 250125 proton
	#./runana 4 $inRootFile $outHistF 252482 proton
	#
	#./runana 4 $inRootFile $outHistF 1673 gamma   # 98 pe
	#./runana 4 $inRootFile $outHistF 1699 gamma   # 100 pe
	#./runana 4 $inRootFile $outHistF 1708 gamma   # 101 pe
	#./runana 4 $inRootFile $outHistF 375 gamma   # 43 pe
	#./runana 4 $inRootFile $outHistF 397 gamma   # 43 pe
	#./runana 4 $inRootFile $outHistF 3058 gamma  # 43 pe
	#./runana 4 $inRootFile $outHistF 2153 gamma  # 43 pe
	#./runana 4 $inRootFile $outHistF 3929 gamma  # 42 pe
        #./runana 4 $inRootFile $outHistF 3932 gamma  # 44 pe
	#./runana 4 $inRootFile $outHistF 4067 gamma   # 42 pe
	#./runana 4 $inRootFile $outHistF 4064 gamma   # 40 pe
	#./runana 4 $inRootFile $outHistF 1222 gamma   # 27 pe
	#./runana 4 $inRootFile $outHistF 1351 gamma    # 32 pe
	#./runana 4 $inRootFile $outHistF 1374 gamma    # 30 pe
	#./runana 4 $inRootFile $outHistF 2993 gamma    # 29 pe
	#./runana 4 $inRootFile $outHistF 3792 gamma    # 30 pe
	#./runana 4 $inRootFile $outHistF 3799 gamma    # 30 pe
	#./runana 4 $inRootFile $outHistF 4540 gamma    # 30 pe
	#./runana 4 $inRootFile $outHistF 4 gamma    # 20 pe
	#./runana 4 $inRootFile $outHistF 111 gamma # 20 pe
	#./runana 4 $inRootFile $outHistF 118 gamma # 20 pe
	#./runana 4 $inRootFile $outHistF 142 gamma # 20 pe
	#./runana 4 $inRootFile $outHistF 405 gamma # 20 pe
	#./runana 4 $inRootFile $outHistF 872 gamma # 20 pe
	#./runana 4 $inRootFile $outHistF 848 gamma # 20 pe
	#./runana 4 $inRootFile $outHistF 952 gamma # 20 pe
	#./runana 4 $inRootFile $outHistF 959 gamma # 20 pe
	#
	#./runana 4 $inRootFile $outHistF 214 gamma  # 55 pe
	#./runana 4 $inRootFile $outHistF 796 gamma  # 52 pe
	#./runana 4 $inRootFile $outHistF 811 gamma  # 51 pe
	#./runana 4 $inRootFile $outHistF 1113 gamma # 55 pe
	#./runana 4 $inRootFile $outHistF 1175 gamma # 52 pe
	#
	#./runana 4 $inRootFile $outHistF  588 gamma  #   68
	#./runana 4 $inRootFile $outHistF  590 gamma  #   88
	#./runana 4 $inRootFile $outHistF  824 gamma  #   70
	#./runana 4 $inRootFile $outHistF  862 gamma  #   87
	#./runana 4 $inRootFile $outHistF  1459 gamma  #   127
	#./runana 4 $inRootFile $outHistF  1486 gamma  #   133
	#./runana 4 $inRootFile $outHistF  1588 gamma  #   119
	#./runana 4 $inRootFile $outHistF  1693 gamma  #   72
	#./runana 4 $inRootFile $outHistF  1825 gamma  #   54
	#./runana 4 $inRootFile $outHistF  1918 gamma  #   81
	#./runana 4 $inRootFile $outHistF  2373 gamma  #   90
	#./runana 4 $inRootFile $outHistF  2405 gamma  #   73
	#./runana 4 $inRootFile $outHistF  2465 gamma  #   86
	#./runana 4 $inRootFile $outHistF  2619 gamma  #   58
	#./runana 4 $inRootFile $outHistF  3415 gamma  #   108
	#./runana 4 $inRootFile $outHistF  3509 gamma  #   90
	#./runana 4 $inRootFile $outHistF  3717 gamma  #   361
	#./runana 4 $inRootFile $outHistF  3718 gamma  #   187
	#./runana 4 $inRootFile $outHistF  4047 gamma  #   121
	#./runana 4 $inRootFile $outHistF  4531 gamma  #   265
	#./runana 4 $inRootFile $outHistF  4821 gamma  #   151
	#./runana 4 $inRootFile $outHistF  5061 gamma  #   184
	#
	#./runana 4 $inRootFile $outHistF 716 gamma    # n_pe 40  energy 11.4927 GeV
	#./runana 4 $inRootFile $outHistF 934 gamma    # n_pe 45  energy 13.5425 GeV
	#./runana 4 $inRootFile $outHistF 626 gamma    # n_pe 41  energy 15.969 GeV
	#
	#./runana 4 $inRootFile $outHistF 10206 gamma # n_pe 42  energy 14.6058 GeV
	#./runana 4 $inRootFile $outHistF 10391 gamma # n_pe 41  energy 13.9274 GeV
	#./runana 4 $inRootFile $outHistF 10580 gamma # n_pe 40  energy 18.939 GeV
	#./runana 4 $inRootFile $outHistF 11697 gamma # n_pe 42  energy 12.4828 GeV
	#./runana 4 $inRootFile $outHistF 11063 gamma # n_pe 40  energy 12.9968 GeV
	#./runana 4 $inRootFile $outHistF 12666 gamma # n_pe 48  energy 14.4241 GeV
	#./runana 4 $inRootFile $outHistF 13624 gamma # n_pe 44  energy 12.7906 GeV
	#./runana 4 $inRootFile $outHistF 13655 gamma # n_pe 41  energy 11.0314 GeV
	#
	#./runana 4 $inRootFile $outHistF 1026 gamma   # n_pe 46  energy 14.9026 GeV
	#./runana 4 $inRootFile $outHistF 1160 gamma   # n_pe 42  energy 15.0509 GeV
	#./runana 4 $inRootFile $outHistF 1444 gamma   # n_pe 45  energy 12.0146 GeV
	#./runana 4 $inRootFile $outHistF 1456 gamma   # n_pe 96  energy 11.5971 GeV
	#./runana 4 $inRootFile $outHistF 949715 gamma # n_pe 87  energy 12.5351 GeV
	#./runana 4 $inRootFile $outHistF 950077 gamma # n_pe 44  energy 11.9362 GeV

	#./runana 4 $inRootFile $outHistF 1430 gamma   # n_pe 50 energy  GeV
	#
	#./runana 4 $inRootFile $outHistF 6 gamma_diff
	#./runana 4 $inRootFile $outHistF 17 gamma_diff
	#./runana 4 $inRootFile $outHistF 19 gamma_diff
	#./runana 4 $inRootFile $outHistF 21 gamma_diff
	#./runana 4 $inRootFile $outHistF 29 gamma_diff
	#./runana 4 $inRootFile $outHistF 46 gamma_diff
	#./runana 4 $inRootFile $outHistF 48 gamma_diff
	#./runana 4 $inRootFile $outHistF 52 gamma_diff
	#./runana 4 $inRootFile $outHistF 53 gamma_diff
	#./runana 4 $inRootFile $outHistF 55 gamma_diff
	#./runana 4 $inRootFile $outHistF 58 gamma_diff
	#./runana 4 $inRootFile $outHistF 63 gamma_diff
	#./runana 4 $inRootFile $outHistF 66 gamma_diff
	#./runana 4 $inRootFile $outHistF 67 gamma_diff
	#./runana 4 $inRootFile $outHistF 68 gamma_diff
	#
	#./runana 4 $inRootFile $outHistF 201  gamma_diff
	#./runana 4 $inRootFile $outHistF 7617  gamma_diff
	#./runana 4 $inRootFile $outHistF 7830  gamma_diff
	#./runana 4 $inRootFile $outHistF 22248 gamma_diff
	#./runana 4 $inRootFile $outHistF 22461 gamma_diff
	#./runana 4 $inRootFile $outHistF 23475 gamma_diff
	#./runana 4 $inRootFile $outHistF 26666 gamma_diff
	#./runana 4 $inRootFile $outHistF 27676 gamma_diff
	#./runana 4 $inRootFile $outHistF 36010 gamma_diff
	#./runana 4 $inRootFile $outHistF 36176 gamma_diff
	#./runana 4 $inRootFile $outHistF 43005 gamma_diff
	#./runana 4 $inRootFile $outHistF 48574 gamma_diff
	#./runana 4 $inRootFile $outHistF 51743 gamma_diff
	#./runana 4 $inRootFile $outHistF 51744 gamma_diff
	#./runana 4 $inRootFile $outHistF 61300 gamma_diff
	#./runana 4 $inRootFile $outHistF 70822 gamma_diff
	#./runana 4 $inRootFile $outHistF 79945 gamma_diff
	#./runana 4 $inRootFile $outHistF 79963 gamma_diff
	#./runana 4 $inRootFile $outHistF 83105 gamma_diff
	#./runana 4 $inRootFile $outHistF 85595 gamma_diff
	#./runana 4 $inRootFile $outHistF 86650 gamma_diff
	#./runana 4 $inRootFile $outHistF 93711 gamma_diff
	#./runana 4 $inRootFile $outHistF 100131 gamma_diff

	#./runana 4 $inRootFile $outHistF 95611 gamma_diff   # 
	#./runana 4 $inRootFile $outHistF 4296421 gamma_diff #
	#./runana 4 $inRootFile $outHistF 728242 gamma_diff # 
	
	#./runana 4 $inRootFile $outHistF 4027 gamma   # n_pe 50 energy  GeV
	#./runana 4 $inRootFile $outHistF 4727 gamma   # n_pe 50 energy  GeV
	#./runana 4 $inRootFile $outHistF 6332 gamma   # n_pe 50 energy  GeV
	#./runana 4 $inRootFile $outHistF 7347 gamma   # n_pe 50 energy  GeV
	#./runana 4 $inRootFile $outHistF 7478 gamma   # n_pe 50 energy  GeV
	#./runana 4 $inRootFile $outHistF 8142 gamma   # n_pe 50 energy  GeV
	#./runana 4 $inRootFile $outHistF 9855 gamma   # n_pe 50 energy  GeV

	##
	#./runana 4 $inRootFile $outHistF 10392 gamma   # n_pe 50 energy  GeV

	#./runana 4 $inRootFile $outHistF 12409 gamma   # n_pe 50 energy  GeV

	##./runana 4 $inRootFile $outHistF 380672 gamma  # n_pe
	##
	
	#./runana 4 $inRootFile $outHistF 380672 proton # n_pe 652779  energy 77215.1 GeV


	#./runana 4 $inRootFile $outHistF 432957 proton # n_pe 12012  energy 12941.7 GeV azi. 179.938 alt. 69.0438
	#./runana 4 $inRootFile $outHistF 663468 proton # n_pe 19384  energy 4243.47 GeV azi. 179.988 alt. 69.0788
	#./runana 4 $inRootFile $outHistF 721973 proton # n_pe 28603  energy 19568.3 GeV azi. 180.181 alt. 68.8473
	#./runana 4 $inRootFile $outHistF 721975 proton # n_pe 28546  energy 19568.3 GeV azi. 180.181 alt. 68.8473
	#./runana 4 $inRootFile $outHistF 1842076 proton # n_pe 18990  energy 4119.99 GeV azi. 180.097 alt. 69.0628
	#./runana 4 $inRootFile $outHistF 2165305 proton # n_pe 51582  energy 2740.96 GeV azi. 180.045 alt. 69.0556
	#./runana 4 $inRootFile $outHistF 2592242 proton # n_pe 60448  energy 74998 GeV azi. 180.067 alt. 69.1491
	#./runana 4 $inRootFile $outHistF 2592340 proton # n_pe 40701  energy 5037.92 GeV azi. 179.902 alt. 69.0582
	#./runana 4 $inRootFile $outHistF 4851729 proton # n_pe 18199  energy 4219.2 GeV azi. 179.86 alt. 68.9217
	#./runana 4 $inRootFile $outHistF 6070811 proton # n_pe 26305  energy 11048.5 GeV azi. 180.182 alt. 68.9255
	#./runana 4 $inRootFile $outHistF 6591666 proton # n_pe 628094  energy 45336.2 GeV azi. 180.042 alt. 69.1757
	#./runana 4 $inRootFile $outHistF 7905149 proton # n_pe 54283  energy 4563.45 GeV azi. 179.881 alt. 68.8071
	#./runana 4 $inRootFile $outHistF 8799036 proton # n_pe 50236  energy 7297.92 GeV azi. 179.904 alt. 69.1129
	#./runana 4 $inRootFile $outHistF 9701121 proton # n_pe 13113  energy 3821.47 GeV azi. 179.902 alt. 69.1146

	#./runana 4 $inRootFile $outHistF 7014318 proton # n_pe 12438  energy 98388.8 GeV azi. 194.212 alt. 75.7073
	#./runana 4 $inRootFile $outHistF 7105497 proton # n_pe 12693  energy 95500.4 GeV azi. 164.026 alt. 75.1417

	#./runana 4 $inRootFile $outHistF 7705579 proton # n_pe 12485  energy 80426.4 GeV azi. 169.65 alt. 74.0754
	#./runana 4 $inRootFile $outHistF 16614505 proton # n_pe 12057  energy 76159.3 GeV azi. 164.719 alt. 71.2268
	#./runana 4 $inRootFile $outHistF 10023276 proton # n_pe 12255  energy 84172.5 GeV azi. 171.875 alt. 76.0566
	
	#./runana 4 $inRootFile $outHistF 10087523 proton # n_pe 12741  energy 85425 GeV azi. 174.655 alt. 62.3514
	#./runana 4 $inRootFile $outHistF 10336791 proton # n_pe 12561  energy 98093.2 GeV azi. 180.704 alt. 76.3203
	#./runana 4 $inRootFile $outHistF 12762628 proton # n_pe 12098  energy 83020.6 GeV azi. 190.514 alt. 68.0341
	#./runana 4 $inRootFile $outHistF 13362942 proton # n_pe 12211  energy 83828.8 GeV azi. 173.218 alt. 77.3532
	#./runana 4 $inRootFile $outHistF 15055958 proton # n_pe 12114  energy 94013.5 GeV azi. 153.116 alt. 73.9085
	#./runana 4 $inRootFile $outHistF 15539658 proton # n_pe 12065  energy 88918 GeV azi. 158.789 alt. 72.1095
	#./runana 4 $inRootFile $outHistF 15792857 proton # n_pe 12090  energy 82888.4 GeV azi. 179.127 alt. 76.8714
	#./runana 4 $inRootFile $outHistF 17010207 proton # n_pe 12775  energy 82958.2 GeV azi. 171.942 alt. 76.153
	#./runana 4 $inRootFile $outHistF 17251109 proton # n_pe 12111  energy 77332.2 GeV azi. 179.034 alt. 68.5323
	#./runana 4 $inRootFile $outHistF 19918811 proton # n_pe 12499  energy 76577.8 GeV azi. 198.538 alt. 66.4997
	
	#./runana 4 $inRootFile $outHistF 1955588 proton # n_pe 663  energy 237.659 GeV azi. 175.277 alt. 70.643
	#./runana 4 $inRootFile $outHistF 1956184 proton # n_pe 929  energy 499.373 GeV azi. 182.713 alt. 71.2581
	#./runana 4 $inRootFile $outHistF 1958417 proton # n_pe 632  energy 629.811 GeV azi. 174.318 alt. 70.4644
	#./runana 4 $inRootFile $outHistF 1961391 proton # n_pe 513  energy 370.986 GeV azi. 181.264 alt. 68.4045
	#./runana 4 $inRootFile $outHistF 1962023 proton # n_pe 527  energy 522.798 GeV azi. 183.21 alt. 70.7176
	#./runana 4 $inRootFile $outHistF 1963887 proton # n_pe 567  energy 53.1942 GeV azi. 180.077 alt. 68.9247


	#./runana 4 $inRootFile $outHistF 294147 gamma_diff # n_pe 744  energy 95.6272 GeV azi. 176.903 alt. 69.4121
	#./runana 4 $inRootFile $outHistF 295097 gamma_diff # n_pe 638  energy 104.954 GeV azi. 184.103 alt. 68.6591
	#./runana 4 $inRootFile $outHistF 297136 gamma_diff # n_pe 745  energy 631.211 GeV azi. 183.728 alt. 69.5613
	#./runana 4 $inRootFile $outHistF 297393 gamma_diff # n_pe 506  energy 104.338 GeV azi. 176.395 alt. 70.158
	#./runana 4 $inRootFile $outHistF 298833 gamma_diff # n_pe 884  energy 311.858 GeV azi. 175.135 alt. 70.8417
	#./runana 4 $inRootFile $outHistF 300017 gamma_diff # n_pe 788  energy 71.6998 GeV azi. 182.422 alt. 69.3299
	#./runana 4 $inRootFile $outHistF 300247 gamma_diff # n_pe 708  energy 120.475 GeV azi. 184.444 alt. 68.9793
	#./runana 4 $inRootFile $outHistF 301566 gamma_diff # n_pe 714  energy 237.02 GeV azi. 174.342 alt. 70.2848
	#./runana 4 $inRootFile $outHistF 302997 gamma_diff # n_pe 797  energy 184.711 GeV azi. 184.042 alt. 70.2937

	#./runana 4 $inRootFile $outHistF 309219 gamma # n_pe 188  energy 14.1998 GeV azi. 180 alt. 70
	#./runana 4 $inRootFile $outHistF 309248 gamma # n_pe 195  energy 38.9857 GeV azi. 180 alt. 70
	#./runana 4 $inRootFile $outHistF 309270 gamma # n_pe 142  energy 72.6081 GeV azi. 180 alt. 70
	#./runana 4 $inRootFile $outHistF 309372 gamma # n_pe 140  energy 121.083 GeV azi. 180 alt. 70
	#./runana 4 $inRootFile $outHistF 309375 gamma # n_pe 168  energy 121.083 GeV azi. 180 alt. 70
	#./runana 4 $inRootFile $outHistF 309407 gamma # n_pe 168  energy 23.2705 GeV azi. 180 alt. 70
	#./runana 4 $inRootFile $outHistF 309425 gamma # n_pe 197  energy 85.1306 GeV azi. 180 alt. 70
	#./runana 4 $inRootFile $outHistF 309443 gamma # n_pe 120  energy 29.1908 GeV azi. 180 alt. 70
	#./runana 4 $inRootFile $outHistF 309545 gamma # n_pe 138  energy 284.292 GeV azi. 180 alt. 70


	./runana 4 $inRootFile $outHistF 1253355 proton # n_pe 193  energy 607.424 GeV azi. 180.691 alt. 70.0067
	./runana 4 $inRootFile $outHistF 1253358 proton # n_pe 121  energy 607.424 GeV azi. 180.691 alt. 70.0067
	./runana 4 $inRootFile $outHistF 1259656 proton # n_pe 176  energy 323.335 GeV azi. 180.008 alt. 70.2926
	./runana 4 $inRootFile $outHistF 1276471 proton # n_pe 160  energy 367.73 GeV azi. 180.838 alt. 69.9816
	./runana 4 $inRootFile $outHistF 1276476 proton # n_pe 135  energy 367.73 GeV azi. 180.838 alt. 69.9816
	./runana 4 $inRootFile $outHistF 1276477 proton # n_pe 143  energy 367.73 GeV azi. 180.838 alt. 69.9816
	./runana 4 $inRootFile $outHistF 1289253 proton # n_pe 160  energy 1228.8 GeV azi. 180.513 alt. 70.079
	./runana 4 $inRootFile $outHistF 1304203 proton # n_pe 185  energy 123.244 GeV azi. 180.468 alt. 69.8788
	./runana 4 $inRootFile $outHistF 1311879 proton # n_pe 168  energy 175.343 GeV azi. 179.775 alt. 70.3453
	./runana 4 $inRootFile $outHistF 1313207 proton # n_pe 136  energy 36.6412 GeV azi. 180.166 alt. 69.8697

    elif [ "$1" = "-th" ]; then
	./runana 5
    elif [ "$1" = "-th_flowerpixID" ]; then
	./runana 55
    elif [ "$1" = "-th_isolated_flower" ]; then
	./runana 555
    elif [ "$1" = "-th_isolated_flower_plus_super" ]; then
	./runana 556
    elif [ "$1" = "-get_isolated_flower_mask" ]; then
	./runana 5550
    elif [ "$1" = "-sg" ]; then
	./runana 6 gamma_on_nsb_1x.list hist_gamma_on_nsb_1x.root
	#./runana 6 gamma_on_nsb_1x_short.list hist_gamma_on_nsb_1x.root
    elif [ "$1" = "-sgd" ]; then
	#./runana 6 gamma_diffuse_nsb_1x.list hist_gamma_diffuse_nsb_1x.root
	./runana 6 gamma_diffuse_nsb_1x.list hist_gamma_diffuse_nsb_1x.root
    elif [ "$1" = "-se" ]; then
	#./runana 6 electron_nsb_1x.list hist_electron_nsb_1x.root
	./runana 6 electron_nsb_1x.list hist_electron_nsb_1x.root
    elif [ "$1" = "-sp" ]; then
	#./runana 6 proton_nsb_1x.list hist_proton_nsb_1x.root
	./runana 6 proton_nsb_1x.list hist_proton_nsb_1x.root
    elif [ "$1" = "-pg" ]; then
	./runana 61 gamma_on_nsb_1x.list hist_gamma_on_nsb_1x_PCAp.root
    elif [ "$1" = "-pgd" ]; then
	./runana 61 gamma_diffuse_nsb_1x.list hist_gamma_diffuse_nsb_1x_PCAp.root
    elif [ "$1" = "-pp" ]; then
	./runana 61 proton_nsb_1x.list hist_proton_nsb_1x_PCAp.root
    elif [ "$1" = "-pg_d" ]; then
	./runana 62 gamma_on_nsb_1x.list hist_gamma_on_nsb_1x_PCAp_component.root
    elif [ "$1" = "-pgd_d" ]; then
	./runana 62 gamma_diffuse_nsb_1x.list hist_gamma_diffuse_nsb_1x_PCAp_component.root
    elif [ "$1" = "-pp_d" ]; then
	./runana 62 proton_nsb_1x.list hist_proton_nsb_1x_PCAp_component.root
    elif [ "$1" = "-pg_d_reco" ]; then
	nEvMax=$(more reco.cvs | wc -l)
	./runana 63 gamma_on_nsb_1x.list hist_gamma_on_nsb_1x_PCA_reco.root $nEvMax 'reco.cvs'
    elif [ "$1" = "-pgd_d_reco" ]; then
	nEvMax=$(more reco.cvs | wc -l)
	nEvMax=50
	./runana 63 gamma_diffuse_nsb_1x.list hist_gamma_diffuse_nsb_1x_PCA_reco.root $nEvMax 'reco.cvs'
    elif [ "$1" = "-pp_d_reco" ]; then
	nEvMax=$(more reco.cvs | wc -l)
	./runana 63 proton_nsb_1x.list hist_proton_nsb_1x_PCA_reco.root $nEvMax 'reco.cvs'
    elif [ "$1" = "-fg" ]; then
	#./runana 7 gamma_on_nsb_1x.list hist_fast_gamma_on_nsb_1x.root anaFast_g.conf | tee hist_fast_gamma_on_nsb_1x.log
	#./runana 7 gamma_on_nsb_1x.list hist_fast_gamma_on_nsb_1x_cut.root anaFast_g.conf | tee hist_fast_gamma_on_nsb_1x_cut.log
	#./runana 7 gamma_on_nsb_1x.list hist_fast_gamma_on_nsb_1x_cut_faint.root anaFast_g.conf | tee hist_fast_gamma_on_nsb_1x_cut_faint.log
	./runana 7 gamma_on_nsb_1x.list hist_fast_gamma_on_nsb_1x_cut_faint.root anaFast_g.conf | tee hist_fast_gamma_on_nsb_1x_cut.log
    elif [ "$1" = "-fgd" ]; then
	#./runana 7 gamma_diffuse_nsb_1x.list hist_fast_gamma_diffuse_nsb_1x.root anaFast_gd.conf | tee hist_fast_gamma_diffuse_nsb_1x.log
	#./runana 7 gamma_diffuse_nsb_1x.list hist_fast_gamma_diffuse_nsb_1x_cut.root anaFast_gd.conf | tee hist_fast_gamma_diffuse_nsb_1x_cut.log.ev
	./runana 7 gamma_diffuse_nsb_1x.list hist_fast_gamma_diffuse_nsb_1x_cut.root anaFast_gd.conf | tee hist_fast_gamma_diffuse_nsb_1x_cut.log
    elif [ "$1" = "-fe" ]; then
	#./runana 7 electron_nsb_1x.list hist_fast_electron_nsb_1x.root anaFast_e.conf | tee hist_fast_electron_nsb_1x.log
	./runana 7 electron_nsb_1x.list hist_fast_electron_nsb_1x.root anaFast_e.conf | tee hist_fast_electron_nsb_1x_cut.log
	#./runana 7 electron_nsb_1x.list hist_fast_electron_nsb_1x.root anaFast_e.conf | tee hist_fast_electron_nsb_1x_cut.log.ev
    elif [ "$1" = "-fp" ]; then
	./runana 7 proton_nsb_1x.list hist_fast_proton_nsb_1x.root anaFast_p.conf | tee hist_fast_proton_nsb_1x.log
	#./runana 7 proton_nsb_1x.list hist_fast_proton_nsb_1x.root anaFast_p.conf | tee hist_fast_proton_nsb_1x_cut.log
	#./runana 7 proton_nsb_1x.list hist_fast_proton_nsb_1x.root anaFast_p.conf | tee hist_fast_proton_nsb_1x_cut.log.ev
    elif [ "$1" = "-fall" ]; then	
	./runana 7 gamma_on_nsb_1x.list hist_fast_gamma_on_nsb_1x.root anaFast_g.conf | tee hist_fast_gamma_on_nsb_1x.log
	./runana 7 gamma_diffuse_nsb_1x.list hist_fast_gamma_diffuse_nsb_1x.root anaFast_gd.conf | tee hist_fast_gamma_diffuse_nsb_1x.log
	./runana 7 electron_nsb_1x.list hist_fast_electron_nsb_1x.root anaFast_e.conf | tee hist_fast_electron_nsb_1x.log
	./runana 7 proton_nsb_1x.list hist_fast_proton_nsb_1x.root anaFast_p.conf | tee hist_fast_proton_nsb_1x.log
    elif [ "$1" = "-sfall" ]; then
	time ./runana 77 gamma_on_nsb_1x.list hist_superfast_gamma_on_nsb_1x.root anaFast_g.conf | tee hist_superfast_gamma_on_nsb_1x.log
	#time ./runana 77 gamma_diffuse_nsb_1x.list hist_superfast_gamma_diffuse_nsb_1x.root anaFast_gd.conf | tee hist_superfast_gamma_diffuse_nsb_1x.log
	#time ./runana 77 electron_nsb_1x.list hist_superfast_electron_nsb_1x.root anaFast_e.conf | tee hist_superfast_electron_nsb_1x.log
	#time ./runana 77 proton_nsb_1x.list hist_superfast_proton_nsb_1x.root anaFast_p.conf | tee hist_superfast_proton_nsb_1x.log
    elif [ "$1" = "-fscang" ]; then
	time ./runana 8 gamma_on_nsb_1x.list hist_fast_scan_gamma_on_nsb_1x.root anaFast_g.conf ../cosmique_gamma_hadron_generator/gamma_on_axis_simtel.dat ../cosmique_gamma_hadron_generator/flux_gamma_crab.dat | tee hist_fast_scan_gamma_on_nsb_1x.log
    elif [ "$1" = "-fscangd" ]; then
	time ./runana 8 gamma_diffuse_nsb_1x.list hist_fast_scan_gamma_diffuse_nsb_1x.root anaFast_gd.conf ../cosmique_gamma_hadron_generator/gamma_diff_galactic_simtel.dat ../cosmique_gamma_hadron_generator/flux_gamma_diff_galactic.dat | tee hist_fast_scan_gamma_diffuse_nsb_1x.log
    elif [ "$1" = "-fscane" ]; then
	time ./runana 8 electron_nsb_1x.list hist_fast_scan_electron_nsb_1x.root anaFast_e.conf ../cosmique_gamma_hadron_generator/electron_simtel.dat ../cosmique_gamma_hadron_generator/flux_ele_pos_diff.dat | tee hist_fast_scan_electron_nsb_1x.log
    elif [ "$1" = "-fscanp" ]; then
	time ./runana 8 proton_nsb_1x.list hist_fast_scan_proton_nsb_1x.root anaFast_p.conf ../cosmique_gamma_hadron_generator/proton_diff_simtel.dat ../cosmique_gamma_hadron_generator/flux_diff_protons.dat | tee hist_fast_scan_proton_nsb_1x.log
    elif [ "$1" = "-test_evstHist" ]; then
	./runana 888
    elif [ "$1" = "-rate_proton_binwise" ]; then
	./runana 999 p ../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/proton_nsb_1x/trgA/ ./hist_rate_proton_diff.root
    elif [ "$1" = "-rate_proton" ]; then
	rsimulation=1500
	#./runana 9999 p ../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/proton_nsb_1x/trgA_20-10000pe/ ./hist_rate_proton_diff.root 1000 100 ../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/nsb_1x_268MHz/trgA/ $rsimulation
	#./runana 9999 p ../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/proton_nsb_1x/trgA_test/ ./hist_rate_proton_diff_test.root 15 15 ../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/nsb_1x_268MHz/trgA_test/ $rsimulation
	./runana 9999 p ../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/proton_nsb_1x/trgA_test_superwlower/ ./hist_rate_proton_diff_test_superwlower.root 15 15 ../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/nsb_1x_268MHz/trgA_test/ $rsimulation
    elif [ "$1" = "-rate_gamma" ]; then
	rsimulation=800
	./runana 9999 g ../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/gamma_on_nsb_1x/trgA_test_DBscan/ ./hist_rate_gamma_DBscan.root 10 15 ../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/nsb_1x_268MHz/trgA_test/ $rsimulation
	./runana 9999 g ../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/gamma_on_nsb_1x/trgA_test_superflower/ ./hist_rate_gamma_superflower.root 10 15 ../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/nsb_1x_268MHz/trgA_test/ $rsimulation
    elif [ "$1" = "-EvPerEv" ]; then
	inRootFiles="../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/gamma_on_nsb_1x/root/0000/corsika_0000ID.root"
	#event_ID="0" #43
	event_ID="78" #3424
	binFileOut="gamma_on_nsb_1x_ev"$event_ID"_out.bin"
	rndseed="1312312"	
	./runana 333 $inRootFiles $event_ID $binFileOut $rndseed
    elif [ "$1" = "-evMAP" ]; then
	inRootFiles="../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/gamma_on_nsb_1x/root/0000/corsika_0000ID.root"
	txtFileOut="../scratch/mono-lst-sipm-pmma-3ns-v1_triggerless/gamma_on_nsb_1x/root/gamma_on_nsb_1x_0000_event_map.txt"
	./runana 3333 $inRootFiles $txtFileOut
    elif [ "$1" = "-c" ]; then
	make clean; make -f Makefileana clean ; make -j; make -f Makefileana -j	
    elif [ "$1" = "-h" ]; then
        printHelp
    else
        printHelp
    fi
fi

#espeak    "I have done"



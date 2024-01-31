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
    echo " [0] -sg    : gamma"
    echo " [0] -sgd   : gamma diffuse"
    echo " [0] -se    : electron"
    echo " [0] -sp    : proton"
    echo " [0] -fg    : fast gamma"
    echo " [0] -fgd   : fast gamma diffuse"
    echo " [0] -fe    : fast electron"
    echo " [0] -fp    : fast proton"
    echo " [0] -fall  : fast all"
    echo " [0] -pg    : PCA gamma"
    echo " [0] -pg_d  : Draw PCA gamma"
    echo " [0] -pgd   : PCA gamma diffuse"
    echo " [0] -pgd_d : Draw PCA gamma diffuse"
    echo " [0] -pp    : PCA proton"
    echo " [0] -pp_d  : Draw PCA proton"
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
	
	#./runana 4 $inRootFile $outHistF 4027 gamma   # n_pe 50 energy  GeV
	#./runana 4 $inRootFile $outHistF 4727 gamma   # n_pe 50 energy  GeV
	#./runana 4 $inRootFile $outHistF 6332 gamma   # n_pe 50 energy  GeV
	./runana 4 $inRootFile $outHistF 7347 gamma   # n_pe 50 energy  GeV
	#./runana 4 $inRootFile $outHistF 7478 gamma   # n_pe 50 energy  GeV
	#./runana 4 $inRootFile $outHistF 8142 gamma   # n_pe 50 energy  GeV
	#./runana 4 $inRootFile $outHistF 9855 gamma   # n_pe 50 energy  GeV
	#./runana 4 $inRootFile $outHistF 10392 gamma   # n_pe 50 energy  GeV
	#./runana 4 $inRootFile $outHistF 12409 gamma   # n_pe 50 energy  GeV
	
    elif [ "$1" = "-th" ]; then
	./runana 5
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
	#./runana 7 proton_nsb_1x.list hist_fast_proton_nsb_1x.root anaFast_p.conf | tee hist_fast_proton_nsb_1x.log
	./runana 7 proton_nsb_1x.list hist_fast_proton_nsb_1x.root anaFast_p.conf | tee hist_fast_proton_nsb_1x_cut.log
	#./runana 7 proton_nsb_1x.list hist_fast_proton_nsb_1x.root anaFast_p.conf | tee hist_fast_proton_nsb_1x_cut.log.ev
    elif [ "$1" = "-fall" ]; then	
	./runana 7 gamma_on_nsb_1x.list hist_fast_gamma_on_nsb_1x.root anaFast_g.conf | tee hist_fast_gamma_on_nsb_1x.log
	./runana 7 gamma_diffuse_nsb_1x.list hist_fast_gamma_diffuse_nsb_1x.root anaFast_gd.conf | tee hist_fast_gamma_diffuse_nsb_1x.log
	./runana 7 electron_nsb_1x.list hist_fast_electron_nsb_1x.root anaFast_e.conf | tee hist_fast_electron_nsb_1x.log
	./runana 7 proton_nsb_1x.list hist_fast_proton_nsb_1x.root anaFast_p.conf | tee hist_fast_proton_nsb_1x.log
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
    elif [ "$1" = "-c" ]; then
	make clean; make -f Makefileana clean ; make -j; make -f Makefileana -j		
    elif [ "$1" = "-h" ]; then
        printHelp
    else
        printHelp
    fi
fi

#espeak    "I have done"



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
#inRootFile="gamma_on_nsb_1x.list"
#outHistF="./hist_gamma_on_tw.root"
#inRootFile="gamma_diffuse_nsb_1x.list"
#outHistF="./hist_gamma_diffuse_nsb_1x_tw.root"
inRootFile="proton_nsb_1x.list"
outHistF="./hist_proton_nsb_1x_tw.root"

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
    echo " [0] -c     : recompile"
    echo " [0] -h     : print help"
}

if [ $# -eq 0 ] 
then    
    printHelp
else
    if [ "$1" = "-d" ]; then
	./runana 1 $inRootFile $outHistF
    elif [ "$1" = "-ds" ]; then
	./runana 3 $inRootFile $outHistF
    elif [ "$1" = "-tw" ]; then

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
	./runana 4 $inRootFile $outHistF 663307 proton
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
	./runana 7 gamma_on_nsb_1x.list hist_fast_gamma_on_nsb_1x_cut.root anaFast_g.conf | tee hist_fast_gamma_on_nsb_1x_cut.log
    elif [ "$1" = "-fgd" ]; then
	#./runana 7 gamma_diffuse_nsb_1x.list hist_fast_gamma_diffuse_nsb_1x.root anaFast_gd.conf | tee hist_fast_gamma_diffuse_nsb_1x.log
	./runana 7 gamma_diffuse_nsb_1x.list hist_fast_gamma_diffuse_nsb_1x_cut.root anaFast_gd.conf | tee hist_fast_gamma_diffuse_nsb_1x_cut_3000pe.log
    elif [ "$1" = "-fe" ]; then
	./runana 7 electron_nsb_1x.list hist_fast_electron_nsb_1x.root anaFast_e.conf | tee hist_fast_electron_nsb_1x.log
    elif [ "$1" = "-fp" ]; then
	./runana 7 proton_nsb_1x.list hist_fast_proton_nsb_1x.root anaFast_p.conf | tee hist_fast_proton_nsb_1x.log
    elif [ "$1" = "-fall" ]; then	
	./runana 7 gamma_on_nsb_1x.list hist_fast_gamma_on_nsb_1x.root anaFast_g.conf | tee hist_fast_gamma_on_nsb_1x.log
	./runana 7 gamma_diffuse_nsb_1x.list hist_fast_gamma_diffuse_nsb_1x.root anaFast_gd.conf | tee hist_fast_gamma_diffuse_nsb_1x.log
	./runana 7 electron_nsb_1x.list hist_fast_electron_nsb_1x.root anaFast_e.conf | tee hist_fast_electron_nsb_1x.log
	./runana 7 proton_nsb_1x.list hist_fast_proton_nsb_1x.root anaFast_p.conf | tee hist_fast_proton_nsb_1x.log
    elif [ "$1" = "-c" ]; then
	make clean; make -f Makefileana clean ; make -j; make -f Makefileana -j		
    elif [ "$1" = "-h" ]; then
        printHelp
    else
        printHelp
    fi
fi

#espeak    "I have done"

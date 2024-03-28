#!/bin/bash

function printHelp {
    echo " --> ERROR in input arguments "
    echo " [0] -pgd_d_reco : Draw reco. showers (PCA) gamma diffuse"
    echo " [0] -h          : print help"
}

if [ $# -eq 0 ] 
then    
    printHelp
else
    if [ "$1" = "-pgd_d_reco" ]; then
	nEvMax=50
	#./runana 63 gamma_diffuse_nsb_1x.list hist_gamma_diffuse_nsb_1x_PCA_reco.root $nEvMax reco_100.cvs
	./runana 63 gamma_diffuse_nsb_1x.list hist_gamma_diffuse_nsb_1x_PCA_reco.root $nEvMax reco_10.cvs
	./runana 63 gamma_diffuse_nsb_1x.list hist_gamma_diffuse_nsb_1x_PCA_reco.root $nEvMax reco_20.cvs
	./runana 63 gamma_diffuse_nsb_1x.list hist_gamma_diffuse_nsb_1x_PCA_reco.root $nEvMax reco_30.cvs
	./runana 63 gamma_diffuse_nsb_1x.list hist_gamma_diffuse_nsb_1x_PCA_reco.root $nEvMax reco_40.cvs
	./runana 63 gamma_diffuse_nsb_1x.list hist_gamma_diffuse_nsb_1x_PCA_reco.root $nEvMax reco_50.cvs
	./runana 63 gamma_diffuse_nsb_1x.list hist_gamma_diffuse_nsb_1x_PCA_reco.root $nEvMax reco_60.cvs
	./runana 63 gamma_diffuse_nsb_1x.list hist_gamma_diffuse_nsb_1x_PCA_reco.root $nEvMax reco_70.cvs
	./runana 63 gamma_diffuse_nsb_1x.list hist_gamma_diffuse_nsb_1x_PCA_reco.root $nEvMax reco_80.cvs
	./runana 63 gamma_diffuse_nsb_1x.list hist_gamma_diffuse_nsb_1x_PCA_reco.root $nEvMax reco_90.cvs
	./runana 63 gamma_diffuse_nsb_1x.list hist_gamma_diffuse_nsb_1x_PCA_reco.root $nEvMax reco_100.cvs
	./runana 63 gamma_diffuse_nsb_1x.list hist_gamma_diffuse_nsb_1x_PCA_reco.root $nEvMax reco_110.cvs
	./runana 63 gamma_diffuse_nsb_1x.list hist_gamma_diffuse_nsb_1x_PCA_reco.root $nEvMax reco_120.cvs
	./runana 63 gamma_diffuse_nsb_1x.list hist_gamma_diffuse_nsb_1x_PCA_reco.root $nEvMax reco_130.cvs
	./runana 63 gamma_diffuse_nsb_1x.list hist_gamma_diffuse_nsb_1x_PCA_reco.root $nEvMax reco_140.cvs
	./runana 63 gamma_diffuse_nsb_1x.list hist_gamma_diffuse_nsb_1x_PCA_reco.root $nEvMax reco_150.cvs
	./runana 63 gamma_diffuse_nsb_1x.list hist_gamma_diffuse_nsb_1x_PCA_reco.root $nEvMax reco_160.cvs
	./runana 63 gamma_diffuse_nsb_1x.list hist_gamma_diffuse_nsb_1x_PCA_reco.root $nEvMax reco_170.cvs
	./runana 63 gamma_diffuse_nsb_1x.list hist_gamma_diffuse_nsb_1x_PCA_reco.root $nEvMax reco_180.cvs
	./runana 63 gamma_diffuse_nsb_1x.list hist_gamma_diffuse_nsb_1x_PCA_reco.root $nEvMax reco_190.cvs
	./runana 63 gamma_diffuse_nsb_1x.list hist_gamma_diffuse_nsb_1x_PCA_reco.root $nEvMax reco_200.cvs
	./runana 63 gamma_diffuse_nsb_1x.list hist_gamma_diffuse_nsb_1x_PCA_reco.root $nEvMax reco_210.cvs
	./runana 63 gamma_diffuse_nsb_1x.list hist_gamma_diffuse_nsb_1x_PCA_reco.root $nEvMax reco_220.cvs
	./runana 63 gamma_diffuse_nsb_1x.list hist_gamma_diffuse_nsb_1x_PCA_reco.root $nEvMax reco_230.cvs
	./runana 63 gamma_diffuse_nsb_1x.list hist_gamma_diffuse_nsb_1x_PCA_reco.root $nEvMax reco_240.cvs
	./runana 63 gamma_diffuse_nsb_1x.list hist_gamma_diffuse_nsb_1x_PCA_reco.root $nEvMax reco_250.cvs
	./runana 63 gamma_diffuse_nsb_1x.list hist_gamma_diffuse_nsb_1x_PCA_reco.root $nEvMax reco_260.cvs
	./runana 63 gamma_diffuse_nsb_1x.list hist_gamma_diffuse_nsb_1x_PCA_reco.root $nEvMax reco_270.cvs
	./runana 63 gamma_diffuse_nsb_1x.list hist_gamma_diffuse_nsb_1x_PCA_reco.root $nEvMax reco_280.cvs
	./runana 63 gamma_diffuse_nsb_1x.list hist_gamma_diffuse_nsb_1x_PCA_reco.root $nEvMax reco_290.cvs
	./runana 63 gamma_diffuse_nsb_1x.list hist_gamma_diffuse_nsb_1x_PCA_reco.root $nEvMax reco_300.cvs
	./runana 63 gamma_diffuse_nsb_1x.list hist_gamma_diffuse_nsb_1x_PCA_reco.root $nEvMax reco_310.cvs
	./runana 63 gamma_diffuse_nsb_1x.list hist_gamma_diffuse_nsb_1x_PCA_reco.root $nEvMax reco_320.cvs
	./runana 63 gamma_diffuse_nsb_1x.list hist_gamma_diffuse_nsb_1x_PCA_reco.root $nEvMax reco_330.cvs
	./runana 63 gamma_diffuse_nsb_1x.list hist_gamma_diffuse_nsb_1x_PCA_reco.root $nEvMax reco_340.cvs
	./runana 63 gamma_diffuse_nsb_1x.list hist_gamma_diffuse_nsb_1x_PCA_reco.root $nEvMax reco_350.cvs
	./runana 63 gamma_diffuse_nsb_1x.list hist_gamma_diffuse_nsb_1x_PCA_reco.root $nEvMax reco_360.cvs
	./runana 63 gamma_diffuse_nsb_1x.list hist_gamma_diffuse_nsb_1x_PCA_reco.root $nEvMax reco_370.cvs
	./runana 63 gamma_diffuse_nsb_1x.list hist_gamma_diffuse_nsb_1x_PCA_reco.root $nEvMax reco_380.cvs
	./runana 63 gamma_diffuse_nsb_1x.list hist_gamma_diffuse_nsb_1x_PCA_reco.root $nEvMax reco_390.cvs
	./runana 63 gamma_diffuse_nsb_1x.list hist_gamma_diffuse_nsb_1x_PCA_reco.root $nEvMax reco_400.cvs
	./runana 63 gamma_diffuse_nsb_1x.list hist_gamma_diffuse_nsb_1x_PCA_reco.root $nEvMax reco_410.cvs
	./runana 63 gamma_diffuse_nsb_1x.list hist_gamma_diffuse_nsb_1x_PCA_reco.root $nEvMax reco_420.cvs
	./runana 63 gamma_diffuse_nsb_1x.list hist_gamma_diffuse_nsb_1x_PCA_reco.root $nEvMax reco_430.cvs
	./runana 63 gamma_diffuse_nsb_1x.list hist_gamma_diffuse_nsb_1x_PCA_reco.root $nEvMax reco_440.cvs
	./runana 63 gamma_diffuse_nsb_1x.list hist_gamma_diffuse_nsb_1x_PCA_reco.root $nEvMax reco_450.cvs
	./runana 63 gamma_diffuse_nsb_1x.list hist_gamma_diffuse_nsb_1x_PCA_reco.root $nEvMax reco_460.cvs
	./runana 63 gamma_diffuse_nsb_1x.list hist_gamma_diffuse_nsb_1x_PCA_reco.root $nEvMax reco_470.cvs
	./runana 63 gamma_diffuse_nsb_1x.list hist_gamma_diffuse_nsb_1x_PCA_reco.root $nEvMax reco_480.cvs
	./runana 63 gamma_diffuse_nsb_1x.list hist_gamma_diffuse_nsb_1x_PCA_reco.root $nEvMax reco_490.cvs
	./runana 63 gamma_diffuse_nsb_1x.list hist_gamma_diffuse_nsb_1x_PCA_reco.root $nEvMax reco_500.cvs
	./runana 63 gamma_diffuse_nsb_1x.list hist_gamma_diffuse_nsb_1x_PCA_reco.root $nEvMax reco_510.cvs
	./runana 63 gamma_diffuse_nsb_1x.list hist_gamma_diffuse_nsb_1x_PCA_reco.root $nEvMax reco_520.cvs
	./runana 63 gamma_diffuse_nsb_1x.list hist_gamma_diffuse_nsb_1x_PCA_reco.root $nEvMax reco_530.cvs
	./runana 63 gamma_diffuse_nsb_1x.list hist_gamma_diffuse_nsb_1x_PCA_reco.root $nEvMax reco_540.cvs
	./runana 63 gamma_diffuse_nsb_1x.list hist_gamma_diffuse_nsb_1x_PCA_reco.root $nEvMax reco_550.cvs
	./runana 63 gamma_diffuse_nsb_1x.list hist_gamma_diffuse_nsb_1x_PCA_reco.root $nEvMax reco_560.cvs
	./runana 63 gamma_diffuse_nsb_1x.list hist_gamma_diffuse_nsb_1x_PCA_reco.root $nEvMax reco_570.cvs
	./runana 63 gamma_diffuse_nsb_1x.list hist_gamma_diffuse_nsb_1x_PCA_reco.root $nEvMax reco_580.cvs
	./runana 63 gamma_diffuse_nsb_1x.list hist_gamma_diffuse_nsb_1x_PCA_reco.root $nEvMax reco_590.cvs
	./runana 63 gamma_diffuse_nsb_1x.list hist_gamma_diffuse_nsb_1x_PCA_reco.root $nEvMax reco_600.cvs
	./runana 63 gamma_diffuse_nsb_1x.list hist_gamma_diffuse_nsb_1x_PCA_reco.root $nEvMax reco_610.cvs
	./runana 63 gamma_diffuse_nsb_1x.list hist_gamma_diffuse_nsb_1x_PCA_reco.root $nEvMax reco_620.cvs
	./runana 63 gamma_diffuse_nsb_1x.list hist_gamma_diffuse_nsb_1x_PCA_reco.root $nEvMax reco_630.cvs
	./runana 63 gamma_diffuse_nsb_1x.list hist_gamma_diffuse_nsb_1x_PCA_reco.root $nEvMax reco_640.cvs
	./runana 63 gamma_diffuse_nsb_1x.list hist_gamma_diffuse_nsb_1x_PCA_reco.root $nEvMax reco_650.cvs
	./runana 63 gamma_diffuse_nsb_1x.list hist_gamma_diffuse_nsb_1x_PCA_reco.root $nEvMax reco_660.cvs
	./runana 63 gamma_diffuse_nsb_1x.list hist_gamma_diffuse_nsb_1x_PCA_reco.root $nEvMax reco_670.cvs
	./runana 63 gamma_diffuse_nsb_1x.list hist_gamma_diffuse_nsb_1x_PCA_reco.root $nEvMax reco_680.cvs
	./runana 63 gamma_diffuse_nsb_1x.list hist_gamma_diffuse_nsb_1x_PCA_reco.root $nEvMax reco_690.cvs
	./runana 63 gamma_diffuse_nsb_1x.list hist_gamma_diffuse_nsb_1x_PCA_reco.root $nEvMax reco_700.cvs
	./runana 63 gamma_diffuse_nsb_1x.list hist_gamma_diffuse_nsb_1x_PCA_reco.root $nEvMax reco_710.cvs
	./runana 63 gamma_diffuse_nsb_1x.list hist_gamma_diffuse_nsb_1x_PCA_reco.root $nEvMax reco_720.cvs
	./runana 63 gamma_diffuse_nsb_1x.list hist_gamma_diffuse_nsb_1x_PCA_reco.root $nEvMax reco_730.cvs
	./runana 63 gamma_diffuse_nsb_1x.list hist_gamma_diffuse_nsb_1x_PCA_reco.root $nEvMax reco_740.cvs
	./runana 63 gamma_diffuse_nsb_1x.list hist_gamma_diffuse_nsb_1x_PCA_reco.root $nEvMax reco_750.cvs
	./runana 63 gamma_diffuse_nsb_1x.list hist_gamma_diffuse_nsb_1x_PCA_reco.root $nEvMax reco_760.cvs
	./runana 63 gamma_diffuse_nsb_1x.list hist_gamma_diffuse_nsb_1x_PCA_reco.root $nEvMax reco_770.cvs
	./runana 63 gamma_diffuse_nsb_1x.list hist_gamma_diffuse_nsb_1x_PCA_reco.root $nEvMax reco_780.cvs
	./runana 63 gamma_diffuse_nsb_1x.list hist_gamma_diffuse_nsb_1x_PCA_reco.root $nEvMax reco_790.cvs
    elif [ "$1" = "-h" ]; then
        printHelp
    else
        printHelp
    fi
fi

#espeak    "I have done"



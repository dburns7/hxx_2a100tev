#!/usr/bin/python

# Branching ratios from PDG and https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageBR2#Higgs_2_gauge_bosons
# Assumes mH=126 GeV
BR_H_AA = 2.28E-03
# Z->v v~  : 20.0%      = 0.200
BR_Z_VV = 0.200
# W->l v   : 10.8 * 2%  = 0.216
BR_W_LV = 0.216


# Higgs production cross sections in fb from https://twiki.cern.ch/twiki/bin/view/LHCPhysics/HiggsEuropeanStrategy#SM_Higgs_decay_branching_ratio_M
PB_TO_FB = 1E3
SIGMA_H  = (740.3 + 82) * PB_TO_FB
SIGMA_WH = 15.9         * PB_TO_FB
SIGMA_ZH = 11.26        * PB_TO_FB  

# Background production cross sections from madgraph
SIGMA_AA  = 0.3605      * PB_TO_FB
SIGMA_ZAA = 0.004908    * PB_TO_FB

# Signal cross section is taken to be 100 pb:
SIGMA_HXX = 100 * PB_TO_FB

print "Sigma x BR in fb at 100 TeV"
print "Backgrounds:  "
print "Zh, Z->vv     ........................", SIGMA_ZH  * BR_Z_VV
print "Wh, W->lv     ........................", SIGMA_WH  * BR_W_LV
print "h->aa         ........................", SIGMA_H   * BR_H_AA
print "aa            ........................", SIGMA_AA
print "Zaa, Z->vv    ........................", SIGMA_ZAA * BR_Z_VV
print "Signal:  "
print "Hxx           ........................", SIGMA_HXX 

# #!/bin/sh

# ##This will generate figures 1 and 2 - gen window 10, alpha 0
# python make_figs.py -f 1,2 -d sim1 
# ##This command will generate figure 3, S8
# python make_figs.py -f 3 -d sim5 
# ##This command will generate figure 4, S9, S10
# python make_figs.py -f 4 -d sim1 

# ##other supplemental
# ##Figure S3
# python make_figs.py -f 1,2 -d sim2 
# ##Figure S1 and S6
# python make_figs.py -f 1,2 -d sim3 
# ##Figure S2 and S7
# python make_figs.py -f 1,2 -d sim4 
# ##Figure S5
# python make_figs.py -f 1,2 -d sim6
# #Figure S4
# python make_figs.py -f 1,2 -d sim7
# ##Figure 3, S8
# python make_figs.py -f 3 -d sim5 
##Figure S11
python make_figs.py -f 5 -d GxAncestry
##Figure S12
python make_figs.py -f 5 -d GxAncestry.pc.corrected
##Figure S13
python make_figs.py -f 5 -d GxAncestry_noCASS
##Figure S14
python make_figs.py -f 5 -d GxAncestry_noCASS.pc.corrected
##Figure S15
python make_figs.py -f 5 -d GxE
##Figure S16
python make_figs.py -f 5 -d GxE_noCASS
##Figure S17
python make_figs.py -f 5 -d ld2

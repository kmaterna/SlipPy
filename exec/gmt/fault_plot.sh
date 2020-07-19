#!/bin/bash
# Script to plot fault patches with components of slip, 
# such as strike slip, dip slip, or tensile. 

slip_file=$1
title=$2
ll1=$3
ll2=$4
ur1=$5
ur2=$6
proj=$7
cmin=$8
cmax=$9
cntv=${10}
observed_gps_file=${11}
predicted_gps_file=${12}
outdir=${13}


range=$ll1/$ur1/$ll2/$ur2
output=$outdir/$title.ps

gmt makecpt -T$cmin/$cmax/$cntv -Cwysiwyg > mycpt.cpt

# Plotting the slip distribution
gmt pscoast -R$range -J$proj -Wthin,black -Di -N1 -N2 -B1.0 -Y4 -Gwhite -Slightgray -K > $output
gmt psxy $slip_file -R$range -J$proj -Wthinnest,gray -Cmycpt.cpt -L -K -O >> $output

# Plotting obs. vs. data
horiz_scale=100
awk '{print $1, $2, $3, $4, $6, $7, 0}'  $observed_gps_file | gmt psvelo  -A+e+gblue+pthickest -Wblue -Se$horiz_scale/0.68/10 -R$range -J$proj -O -K >> $output
awk '{print $1, $2, $3, $4, $6, $7, 0}' $predicted_gps_file | gmt psvelo -A+e+gred+pthickest -Wred -Se$horiz_scale/0.68/10 -R$range -J$proj -O -K >> $output

gmt psscale -Cmycpt.cpt -Dx7c/-1c+w12c/0.5c+jTC+h -Bxaf+l$title" slip" -By+l"m" -K -O >> $output

# Scale vectors. Their placement and length will be problem-dependent. Might need to adjust. 
gmt psvelo -A+e+gred+pthickest -Wred -Se$horiz_scale/0.68/10 -R$range -J$proj -O -K <<EOF>> $output
`echo $ll1 + 0.1 | bc` `echo $ll2 + 0.02 | bc` .005 0 0 0 0 5mm model
EOF
gmt psvelo -A+e+gblue+pthickest -Wblue -Se$horiz_scale/0.68/10 -R$range -J$proj -O -K <<EOF>> $output
`echo $ll1 + 0.1 | bc` `echo $ll2 + 0.04 | bc` .005 0 0.001 0.001 0 5mm obs
EOF


rm gmt.history 
rm mycpt.cpt
# open $output


#!/bin/bash
# GMT plot to show the GPS inversion's fit to the data. 

observed_gps_file=$1
predicted_gps_file=$2
ll1=$3
ll2=$4
ur1=$5
ur2=$6
proj=$7
faultfile=$8
outdir=$9

range=$ll1/$ur1/$ll2/$ur2
output=$outdir"gps_displacements.ps"
horiz_scale=100

# Color scales
gmt makecpt -T-10/10/0.5 -Cpolar -D > mycpt.cpt
gmt makecpt -T0/0.2/0.001 -Cwysiwyg > slip.cpt

gmt pscoast -R$range -J$proj -Wthin,black -Di -N1 -N2 -B1.0 -Y4 -Gwhite -Slightgray -K > $output

# Put the fault annotations for kicks
gmt psxy $faultfile -R$range -J$proj -Wthinnest,gray -Cslip.cpt -L -K -O >> $output # Putting the faults on there for kicks

# Vertical
gmt pscoast -R$range -J$proj -Wthin,black -Di -N1 -N2 -B1.0 -B+t"GPS displacements" -K -O >> $output
awk '{print $1, $2, $5*1000}' $observed_gps_file | gmt psxy -Sh0.35 -Cmycpt.cpt -W0.1,black -R$range -J$proj -O -K >> $output
awk '{print $1, $2, $5*1000}' $predicted_gps_file | gmt psxy -Sh0.2 -Cmycpt.cpt -W0.1,red -R$range -J$proj -O -K >> $output

# Horizontal Vectors
awk '{print $1, $2, $3, $4, $6, $7, 0}'  $observed_gps_file | gmt psvelo -A+e+gblack+pthickest -Wblack -Se$horiz_scale/0.68/10 -R$range -J$proj -O -K >> $output
awk '{print $1, $2, $3, $4, $6, $7, 0}' $predicted_gps_file | gmt psvelo -A+e+gred+pthickest -Wred -Se$horiz_scale/0.68/10 -R$range -J$proj -O -K >> $output

# Scale annotations
gmt psscale -Cmycpt.cpt -Dx14c/8c+w5c/0.5c+jTC -Bxaf+l"Vertical" -By+l"mm" -K -O >> $output
gmt psscale -Cslip.cpt -Dx17c/8c+w5c/0.5c+jTC -Bxaf+l"Slip" -By+l"m" -K -O >> $output

# Scale vectors. Their placement and length will be problem-dependent. Might need to adjust. 
gmt psvelo -A+e+gred+pthickest -Wred -Se$horiz_scale/0.68/10 -R$range -J$proj -O -K <<EOF>> $output
`echo $ll1 + 0.15 | bc` `echo $ll2 + 0.02 | bc` .005 0 0 0 0 5mm model
EOF
gmt psvelo -A+e+gblack+pthickest -Wblack -Se$horiz_scale/0.68/10 -R$range -J$proj -O -K <<EOF>> $output
`echo $ll1 + 0.15 | bc` `echo $ll2 + 0.04 | bc` .005 0 0.001 0.001 0 5mm obs
EOF

rm mycpt.cpt
rm slip.cpt
rm gmt.history 
open $output

#!/bin/bash
# GMT plot to show the GPS inversion's fit to the data. 

observed_gps_file=$1
predicted_gps_file=$2
ll1=$3
ll2=$4
ur1=$5
ur2=$6
proj=$7

range=$ll1/$ur1/$ll2/$ur2
output="displacements.ps"
horiz_scale=10

gmt pscoast -R$range -J$proj -Wthin,black -Di -N1 -N2 -B1.0 -B+t"GPS displacements" -Y4 -Gwhite -Slightgray -K > $output
awk '{print $1, $2, $3, $4, $6, $7, 0}' $predicted_gps_file | gmt psvelo -A+e+gblue+pthickest -Wblue -Se$horiz_scale/0.68/10 -R$range -J$proj -O -K >> $output
awk '{print $1, $2, $3, $4, $6, $7, 0}'  $observed_gps_file | gmt psvelo -A+e+gred+pthickest -Wred -Se$horiz_scale/0.68/10 -R$range -J$proj -O -K >> $output

# Scale vectors. Their placement and length will be problem-dependent. Might need to adjust. 
gmt psvelo -A+e+gblue+pthickest -Wblue -Se$horiz_scale/0.68/10 -R$range -J$proj -O -K <<EOF>> $output
`echo $ll1 + 1.2 | bc` `echo $ll2 + 0.2 | bc` .050 0 0 0 0 5cm model
EOF
gmt psvelo -A+e+gred+pthickest -Wred -Se$horiz_scale/0.68/10 -R$range -J$proj -O -K <<EOF>> $output
`echo $ll1 + 1.2 | bc` `echo $ll2 + 0.3 | bc` .050 0 0.003 0.003 0 5cm obs
EOF

rm gmt.history 
open $output

#!/bin/bash
# GMT plot to show the InSAR's fit to the data. 
# Observations
# Models
# Residuals

obs_insar_file=$1
model_insar_file=$2
ll1=$3
ll2=$4
ur1=$5
ur2=$6
proj=$7
scale_low=$8
scale_high=$9
scale_int=${10}
slip_file=${11}
plottype=${12}

range=$ll1/$ur1/$ll2/$ur2
output=$plottype"_displacements.ps"

# LOS scale
gmt makecpt -T$scale_low/$scale_high/$scale_int -Cjet -D > mycpt.cpt

# Observation
gmt pscoast -R$range -J$proj -Wthin,black -Di -N1 -N2 -B1.0 -B+t"Observation" -Y4 -Gwhite -Slightgray -K > $output
awk '{print $1, $2, $3*1000}' $obs_insar_file | gmt psxy -Sc0.09 -Cmycpt.cpt -R$range -J$proj -O -K >> $output

# Model
gmt pscoast -R$range -J$proj -Wthin,black -Di -N1 -N2 -B1.0 -B+t"Model" -X3.5i -Gwhite -Slightgray -O -K >> $output
awk '{print $1, $2, $3*1000}' $model_insar_file | gmt psxy -Sc0.09 -Cmycpt.cpt -R$range -J$proj -O -K >> $output

# Residual
gmt pscoast -R$range -J$proj -Wthin,black -Di -N1 -N2 -B1.0 -B+t"Residual" -X3.5i -Gwhite -Slightgray -O -K >> $output
paste $obs_insar_file $model_insar_file | awk '{print $1, $2, $3*1000, $10*1000}'  > temp_insar.txt
awk '{print $1, $2, $3-$4}' temp_insar.txt | gmt psxy -Sc0.09 -Cmycpt.cpt -R$range -J$proj -O -K >> $output
gmt psxy $slip_file -R$range -J$proj -Wthinnest,gray -L -K -O >> $output # Putting the faults on there for kicks

gmt psscale -Cmycpt.cpt -Dx-7c/-1c+w9c/0.5c+jTC+h -Bxaf+l"Disp" -By+l"mm" -K -O >> $output

rm mycpt.cpt
rm gmt.history 
rm temp_insar.txt
open $output

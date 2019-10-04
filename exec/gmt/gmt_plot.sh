#!/bin/bash

slip_file=$1
title=$2
ll1=$3
ll2=$4
ur1=$5
ur2=$6
proj=$7
cmin=$8
cmax=$9
cntv=$10


range=$ll1/$ur1/$ll2/$ur2
output=$title.ps

gmt makecpt -T$cmin/$cmax/$cntv -Cwysiwyg > mycpt.cpt

gmt pscoast -R$range -J$proj -Wthin,black -Di -N1 -N2 -B1.0 -Y4 -Gwhite -Slightgray -K > $output
gmt psxy $slip_file -R$range -J$proj -Wthinner,white -Cmycpt.cpt -L -K -O >> $output

gmt psscale -Cmycpt.cpt -Dx7c/-1c+w12c/0.5c+jTC+h -Bxaf+l$title" slip" -By+l"m" -K -O >> $output

rm gmt.history 
rm mycpt.cpt

open $output


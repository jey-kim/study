#!/bin/bash


S="-Sx.001"

input_ext="principal_ext.out"
input_con="principal_con.out"
GNSS_data="GNSS_horizontal_synthetic.gmt"
M="average_strain_pred.ps"
J="-Jm12"
B="-BneSW -Ba0.2f0.2"
R_larger="-R-117.05/-115.95/33.15/34.15"
R="-R-117/-116/33.2/34.1"
CPT="vert.cpt"
T="-T-20/20/4"


awk '{print $1,$2}' $GNSS_data > stations.dat
gmt gmtset FORMAT_GEO_MAP D

#gmt makecpt -Chaxby $T -Z > $CPT
#gmt blockmean $vert_pred $R -I0.6m -V > vert_pred.dat
#gmt surface vert_pred.dat -Gvert_pred.grd -I1m $R -T0.2

gmt psbasemap $R_larger  $J $B -K -V -Y2.0i -P > $M
#gmt grdimage vert_pred.grd -C$CPT $R_larger $J -V -O -K >>$M
gmt psxy fault.dat $R_larger $J -W0.2,black $B -O -K -V >>$M
gmt psxy san_andras_fault.dat $R_larger $J -W0.2,black -O -K -V >>$M
gmt psxy stations.dat $R_larger $J -St0.4c -Gwhite -W1,black -O -K -V >>$M

gmt pscoast -Slightblue $R_larger $J -O -K -Df >>$M

gmt psvelo $input_con $J $R_larger $S -W2p,133/133/133 -A0.1p+a30,30 -K -O -V >> $M
gmt psvelo $input_ext $J $R_larger $S -W2p,red -A0.1p+a30,30 -K -O -V >> $M

gmt psxy $R_larger $J -O -K -L -Wthicker -V >>$M <<END
-117 33.20
-117 34.10
-116 34.10
-116 33.20
END

#gmt psscale -C$CPT -Dx5.6i/0i+w5.4i/0.2i -By+l"[mm/yr]" -L -O -K >>$M

gmt psvelo $R_larger $J -W1.5p,red $S -L -A0.1p+a30,30 -O -K <<EOT>> $M
-116.96 33.18  250  0 0
EOT

gmt pstext  $R_larger $J -P -O -V <<EOT>> $M
-116.79  33.18    15   0.0   0       6     2.5x10@+-7@+/yr
EOT



gmt psconvert $M -Tf -A


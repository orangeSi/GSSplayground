. ../env.sh
set -vex
sh ../ClustersPloter.sh tracks.1.list out1 . main.1.conf
sh ../ClustersPloter.sh tracks.2.list out2 . main.2.conf
sh ../ClustersPloter.sh tracks.3.list out3 . main.3.conf
sh ../ClustersPloter.sh  tracks.4.list  out4 . main.4.conf
sh ../ClustersPloter.sh tracks.5.list out5 . main.5.conf
sh ../ClustersPloter.sh  tracks.6.list  out6 . main.6.conf
sh ../ClustersPloter.sh  tracks.7.list  out7 . main.7.conf

sh ../ClustersPloter.sh  tracks.8.list out8 ./ main.8.conf
ls -lhtr *error 
ls -lhtr *html
sh ../ClustersPloter.sh tracks.3.list out3.1 . main.3.1.conf

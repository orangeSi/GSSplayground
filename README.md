# ClustersPloter
install:<br>
&nbsp;&nbsp;&nbsp;&nbsp;way1: git clone

usage:<br>
&nbsp;&nbsp;&nbsp;&nbsp;cd example <br>
&nbsp;&nbsp;&nbsp;&nbsp;cat test.sh <br>
&nbsp;&nbsp;&nbsp;&nbsp;plot gene clusters of many samples, one track means one sample, one track contain more than one fragments. one fragment contain gene cluster. you can defined every gene(color,label,font size) in clusters. And add crossing link for any pair of genes.<br>

main feature:<br>
&nbsp;&nbsp;&nbsp;&nbsp;1. every track mean one sample , one sample can has more than one fragments. you can defind the feature color/lable font size/label color/label rotaion in feature.color.label.conf <br>
&nbsp;&nbsp;&nbsp;&nbsp;2. you can draw crosslink or sysnteny among features of different tracks<br>

bug:<br>
&nbsp;&nbsp;&nbsp;&nbsp;1. text font size of legend may be oversize or undersize ，so you can adjust legend_font_size of main.conf by hand.

todo:<br>
&nbsp;&nbsp;&nbsp;&nbsp;1. plot tracks by sort sample list, or you can adjust the track order by adjust list file <br>
&nbsp;&nbsp;&nbsp;&nbsp;2. feature only has arrow type, add other type<br>
&nbsp;&nbsp;&nbsp;&nbsp;3. sort by feature,so same feature of different tracks can align centre<br>
&nbsp;&nbsp;&nbsp;&nbsp;4. add ratio scale in the bottom to recognise the approximate length of every feature <br> <br> 

changelog:<br>
&nbsp;&nbsp;&nbsp;&nbsp;1. update *conf for more easily specific the start or end position , remove the border line of legend , freely defind the crosslink color and opacity one by one <br> <br> 

![gene cluster image](example/out.svg)
<br>
![gene cluster image](example/out2.svg)
<br>
![gene cluster image](example/out3.svg)
<br>
![gene cluster image](example/out6.svg)

contact:<br>
&nbsp;&nbsp;&nbsp;&nbsp;QQ: 1522051171<br>
&nbsp;&nbsp;&nbsp;&nbsp;mail: ilikeorangeapple@gmail.com

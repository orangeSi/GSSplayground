# ClustersPloter
#### install:<br>
```
git clone http://github.com/orangeSi/ClustersPloter.git
```

#### usage:<br>
```
  cd  example 
  cat test.sh 
```
&nbsp;&nbsp;&nbsp;&nbsp;plot gene clusters of many samples, one track means one sample, one track contain more than one fragments. one fragment contain gene cluster. you can defined every gene or feature(rotation,color,label,order depth,font size) in clusters. And add crossing link for any pair of genes.<br>

#### main feature:<br>
- every track mean one sample , one sample can has more than one fragments. you can defind the feature color/lable font size/label color/label rotaion in feature.color.label.conf <br>
- you can draw crosslink or sysnteny among features of different tracks<br>

### update feture:<br>
#### 2018-11-12:<br>
- rewrite feature.crossing.link.conf so that defined color and opactity or order or anchor position of every pair of links.<br><br>
#### 2018-11-07:<br>
- feature_shaple now support <b>circle_point</b>/rect/arrow, add <b>feature_shift_y</b>, as example/out9.svg. combine with <b>feature_height_ratio</b> or <b>feature_shift_y</b>, try to plot <b>line</b> or <b>scatter</b> or <b>heatmap</b> or <b>histogram</b> is possible but hard for much data points, I will try do this by embed svg in svg, the embedded svg is producted by python or R tool, abd and the linux convert cmd not support embed svg, but cairosvg support it (should use the absoluted embedded svg path)!<br>
#### 2018-10-29:<br>
- redesign the the main.conf and feature.crossing.link.conf to make more freely to reset every feature, reset every link color and opacity one by one, and reset order depth of every feature or crosslink or track, reset feature height one by one as in out2.svg<br>
- add rect for feature shape, now have arrow and rect, you can use arrow and rect at the same time<br>
- remove the legend border line<br>

#### bug:<br>
&nbsp;&nbsp;&nbsp;&nbsp; welcome to call me at the issue<br>

#### todo:<br>
&nbsp;&nbsp;**2018-11-09**:<br>
&nbsp;&nbsp;&nbsp;&nbsp;1. plot tracks by sort sample list, or you can adjust the track order by adjust list file <br>
&nbsp;&nbsp;&nbsp;&nbsp;3. sort by feature,so same feature of different tracks can align centre<br>
&nbsp;&nbsp;&nbsp;&nbsp;4. add ratio scale in the bottom to recognise the approximate length of every feature <br>
&nbsp;&nbsp;&nbsp;&nbsp;5. embed svg in svg to make heatmap or scatter or line or histogram more easy to product, like example/out.svg*svg or \<foreignObject\> <br> <br>


![gene cluster image](example/out.svg)
<br><br><br>
![gene cluster image](example/out2.svg)
<br><br><br>
![gene cluster image](example/out3.svg)
<br><br>
![gene cluster image](example/out3.2.svg)
<br><br>
![gene cluster image](example/out9.svg)
<br><br><br>
![gene cluster image](example/out6.svg)
<br><br><br>
![gene cluster image](example/out7.svg)

contact:<br>
&nbsp;&nbsp;&nbsp;&nbsp;QQ: 1522051171<br>
&nbsp;&nbsp;&nbsp;&nbsp;mail: ilikeorangeapple@gmail.com

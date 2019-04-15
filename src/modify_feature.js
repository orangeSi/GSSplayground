

// 每个feature或者文字都可以点击后出现一个窗口，可以在里面设置这个feature的移动transform="translate(-110 0)"(即整体x-10,y不变)，color，opacity，text font size
function feature_show(event) { 
				var title=this.getElementsByTagName("title");
				var feature_info;
				var childs=[];
				console.log(title.length);
				if(title.length > 0){
					feature_info = title[0].innerHTML;
					childs = this.childNodes;
				}else{
					feature_info = this.innerHTML;
					childs.push(this);
				}
				console.log('feature_info is '+feature_info); 
				//var childs = this.getElementsByTagName("*");
				//var childs = this.querySelectorAll("*");
				
				console.log(childs.length);
				for (var i=0;i < childs.length;i++){
					var tagname=childs[i].tagName;
					var type=typeof tagname;
					console.log('tagname is '+tagname+", typeof is "+type);
					if(type == undefined || tagname  == 'title' || tagname == 'TSPAN' || tagname == 'defs' || tagname == 'clipPath' || tagname == undefined) continue;
					
					var info_btn= "<tspan style='background-color:green;opacity:0.7;' onMouseOver=\"this.style.cursor='pointer'\" onclick=\"div2.style.display = 'none';div1.style.display = 'block'; div1.getElementsByTagName('tspan')[0].style.opacity=0.7;div1.getElementsByTagName('tspan')[1].style.opacity=0.5;\">Feature_Info</tspan>";
					var attr_btn= "<tspan style='background-color:green;opacity:0.5;' onMouseOver=\"this.style.cursor='pointer'\" onclick=\"div1.style.display = 'none';div2.style.display = 'block'; div2.getElementsByTagName('tspan')[0].style.opacity=0.5;div2.getElementsByTagName('tspan')[1].style.opacity=0.7;console.log('div1.width'+div1.offsetWidth)\">Setting</tspan>";
					 
					 // highlight this
					 var rawopacity=childs[i].style.opacity;
					 if(rawopacity>=0.8){
						childs[i].style.opacity=0.4;
					 }else{
						childs[i].style.opacity=1;
					}
					 setTimeout(function(child, opa){child.style.opacity=opa;}, 3000, childs[i], rawopacity);
					 
					 /// popup dialog box
					 var x_pos=event.clientX;
					 var y_pos=event.clientY;
					 
					 if(div1){div1.remove()}
					 if(div2){div2.remove()}
					 
					 console.log("clientX: " + event.clientX + " - clientY: " + event.clientY);
					 
					 // info of div
					 div1 = document.createElement('div');
					 div1.style.position = "absolute";
					 div1.style.left = x_pos+'px';
					 div1.style.top = y_pos+'px';
					 div1.style.width = "auto";
					 div1.style.height = "auto";
					 div1.style.background = "white";
					 div1.style.border = '2px solid rgba(0, 128, 0, .5)'; //0.6 is opacity, color is green
					 div1.style.color = "black";
					 div1.innerHTML= "<br>&nbsp;"+feature_info;
					 div1.innerHTML = div1.innerHTML.replace(/<\/tspan>/g, '&nbsp;&nbsp;</tspan><br/>&nbsp;').replace(/<tspan>/, '&nbsp;<tspan>'); 
					 div1.innerHTML = info_btn+"  "+attr_btn+"<br>" + div1.innerHTML + "<br><br>";
					
					 //div1.innerHTML += "<br>" +tagname;
					 
					 // property of feature
					 div2 = document.createElement('div');
					 div2.style.position = "absolute";
					 div2.style.left = x_pos+'px';
					 div2.style.top = y_pos+'px';
					 div2.style.width = "auto";
					 div2.style.display = "none";
					 div2.style.height = "auto";
					 div2.style.background = "white";
					 div2.style.border = '2px solid rgba(0, 128, 0, .5)'; //0.6 is opacity, color is green
					 div2.style.color = "black";
					 var shift_x=0;
					 var shift_y=0;
					 var zoom_scale=1;
					 target = childs[i];
					 console.log('target style is ');
					 console.log(target.style);
					 var translate_raw = get_transform('translate', target);
					 var scale_raw = get_transform('scale', target);
					 var rotate_raw = get_transform('rotate', target);
					 
					 var shift_xy="&nbsp;&nbsp;move by x y: &nbsp;&nbsp;<input type=button value='x+' onclick='caculate_translate(target, \"shift_xy\", \"xy_step\", \"x\", \"+\")'>&nbsp;<input type=button value='x-'  onclick='caculate_translate(target, \"shift_xy\", \"xy_step\", \"x\", \"-\")'>&nbsp;<input type=button value='y+' onclick='caculate_translate(target, \"shift_xy\", \"xy_step\", \"y\", \"+\")'>&nbsp;<input type=button value='y-' onclick='caculate_translate(target, \"shift_xy\", \"xy_step\", \"y\", \"-\")'>&nbsp;&nbsp;<input id='shift_xy' onkeyup='transform_it(target, this.value, \"translate\")' type=text size=5 title='"+translate_raw+"' value='"+translate_raw+"' >&nbsp;&nbsp;step:<input id='xy_step' type=text value=1 size=3 >&nbsp;&nbsp;<br> "
					 var scale_xy="&nbsp;&nbsp;scale x and y: &nbsp;&nbsp;<input id='scale_xy' type=text size=5 onkeyup='transform_it(target, this.value, \"scale\")' value='"+scale_raw+"' ><br>";
					 var rotate="&nbsp;&nbsp;rotate: &nbsp;&nbsp;<input id='rotate' type=text size=5 onkeyup='transform_it(target, this.value, \"rotate\")' value='"+rotate_raw+"' ><br>";
					 if(tagname == "text"){
						//"<input type=value >"
						// font size , font color
						//	transform=(translate(shift_x shift_y))
						//div2.innerHTML+="<br>ddddddddhhhhhhhhhhhhddddd text";
						//transform="translate(-110 0)
						var font_size=childs[i].style['font-size'];
						var font_content=childs[i].innerHTML;
						var font_color=childs[i].style['fill'];
						var font_opacity=childs[i].style['opacity'];
						var alignment_baseline=childs[i].style['alignment-baseline'];
						var text_anchor=childs[i].style['text-anchor'];
						if(alignment_baseline){
							alignment_baseline="&nbsp;&nbsp;alignment_baseline:&nbsp;&nbsp;<input type=text size=5 title='"+alignment_baseline+" , baseline or middle' value='"+alignment_baseline+"'  onkeyup='change_style(target, this.value, \"alignment-baseline\")' ><br>"
						}
						if(text_anchor){
							text_anchor="&nbsp;&nbsp;text_anchor:&nbsp;&nbsp;<input type=text size=5 title='"+text_anchor+", start/end/middle' value='"+text_anchor+"'  onkeyup='change_style(target, this.value, \"text-anchor\")' ><br>"
						}
						if(font_opacity == undefined){
							font_opacity = 1;
						}
						if(font_size){
							font_size="&nbsp;&nbsp;font size:&nbsp;&nbsp;<input type=text size=5 title='"+font_size+"' value='"+font_size+"'  onkeyup='change_style(target, this.value, \"font-size\")' ><br>";
						}
						var value = "&nbsp;&nbsp;font conent:&nbsp;&nbsp;<input type=text size=20 title='"+font_content+"' value='"+font_content+"' onkeyup='change_style(target, this.value, \"value\")'><br>";
						if(font_color){
							font_color ="&nbsp;&nbsp;font color:&nbsp;&nbsp;<input type=text size=5 title='"+font_color+"' value='"+font_color+"' onkeyup='change_style(target, this.value, \"fill\")'><br>";
						}
						if(font_opacity){
							font_opacity="&nbsp;&nbsp;font opacity:&nbsp;&nbsp;<input type=text size=5 title='"+font_opacity+"' value='"+font_opacity+"' onkeyup='change_style(target, this.value, \"opacity\")'><br>";
						}
						var display="&nbsp;&nbsp;<input type=button value=show onclick='if(target.style[\"display\"] == \"none\"){ change_style(target, \"block\", \"display\"); this.value=\"show\"}else{change_style(target, \"none\", \"display\");this.value=\"hide\" }'> <br> ";

						//if ('opacity' in childs[i].style){
							//font_opacity = childs[i].style['opacity'];
						//}else{
							//font_opacity = 1;
						//}
						//alert('font size'+font_size);
						div2.innerHTML+="<br>"+font_size+value+font_color+font_opacity+alignment_baseline+text_anchor+scale_xy+rotate+shift_xy+display+"<br><br>";
						
					}else if(tagname == "rect" || tagname == "polygon" || tagname == "path" || tagname == "polygon" || tagname == "circle"){
						// fill color, opacity, strike-border, strike color
						//style="stroke:black;stroke-width:0.2;fill:orange;opacity:0.5"
						var fill_color=childs[i].style['fill'];
						var fill_opacity=rawopacity;
						console.log('fill_opacity is '+fill_opacity);
						
						if(!fill_opacity && fill_opacity != 0){fill_opacity=1}
						var strike_width=childs[i].style['stroke-width'];
						if(strike_width == ""){strike_width = 0.1}
					
						strike_width = parseFloat(strike_width);
						var strike_color=childs[i].style['stroke'];
						var strike_opacity=childs[i].style['stroke-opacity'];
						if(strike_opacity == ""){strike_opacity = 1}
						strike_opacity= parseFloat(strike_opacity);
						
						console.log("strike_width "+strike_width);
						console.log('strike_opacity is ' +strike_opacity);
						console.log('strike_color is ' +strike_color);
						
						//if(strike_opacity && strike_width){
							console.log('here');
							strike_opacity = "&nbsp;&nbsp;strike color opacity:&nbsp;&nbsp;<input type=text size=5 title='"+strike_opacity+"' value='"+strike_opacity+"' onkeyup='target.style[\"stroke-opacity\"]=this.value'><br>";
							strike_width = "&nbsp;&nbsp;strike width:&nbsp;&nbsp;<input type=text size=5 title='"+strike_width+"' value='"+strike_width+"' onkeyup='target.style[\"stroke-width\"]=this.value'><br>";
							strike_color = "&nbsp;&nbsp;strike color:&nbsp;&nbsp;<input type=text size=5 title='"+strike_color+"' value='"+strike_color+"' onkeyup='target.style[\"stroke\"]=this.value'><br>";
						//}else{
						//	strike_opacity = "";
						//	strike_width = "";
						//	strike_color = "";
						//}
						var fill_color = "&nbsp;&nbsp;fill color:&nbsp;&nbsp;<input type=text size=5 title='"+fill_color+"' value='"+fill_color+"' onkeyup='target.style[\"fill\"]=this.value'><br>";
						
						if(fill_opacity){
							fill_opacity="&nbsp;&nbsp;fill color opacity:&nbsp;&nbsp;<input type=text size=5 title='"+fill_opacity+"' value='"+fill_opacity+"' onkeyup='target.style[\"opacity\"]=this.value'><br>";
						}
						var display="&nbsp;&nbsp;<input type=button value=show onclick='if(target.style[\"display\"] == \"none\"){ change_style(target, \"block\", \"display\"); this.value=\"show\"}else{change_style(target, \"none\", \"display\");this.value=\"hide\" }'> <br> ";

						div2.innerHTML+="<br>"+fill_color+fill_opacity+shift_xy+strike_color+strike_opacity+strike_width+display+"<br><br>";
						//div2.innerHTML+="<br>ddddddddhhhhhhhhhhhhddddd "+tagname;
					}else{
						alert("not support "+tagname+" yet~");
						return 1;
					}
					 
					 
					 div2.innerHTML = info_btn+"  "+attr_btn+"<br>" + div2.innerHTML;
					 
					 console.log("inner is "+div1.innerHTML);
					 document.body.appendChild(div1);
					 document.body.appendChild(div2);
					 
					 console.log('tag name is '+childs[i].tagName);
					//continue if(childs[i].tag);
					 console.log(childs[i].style);
					 //childs[i].style.fill='black';
					//this.style.color="black";
				}
			}
function change_features(){
		  
		  var myths = document.getElementsByClassName("myth");
		  for(var j=0;j < myths.length;j++){
			//console.log(j);
			myths[j].addEventListener('click', feature_show);
		}
		window.onclick = function(event) {
			var tagname=event.target.tagName;
			console.log('target is '+tagname);
			if (tagname == 'HTML' || tagname == 'H1' || tagname == "BODY" || tagname =="svg") {
				if(div1){
					//div1.style.display = "none";
					div1.remove();
				}
				if(div2){
					//div2.style.display = "none";
					div2.remove();
				}
			}
		}
		//alert("find feature -> click it -> modify attributions");
}
	
var div1, div2, target; // div1 and div2 are used as global variable in change_feautures;
change_features();
	
//var rotate_raw = get_transform('rotate', target);
// get transform attribute
function get_transform(attr, target){
	var transform = d3.select(target).attr('transform'); //rotate(30 200 314)  translate(-10 -10) scale(1 0.5) 
	var ret;
	if(transform){
		var myRe = new RegExp(attr);
		
		if(attr == "translate"){
			if(transform.match(myRe)){
				ret = transform.replace(/^.*translate\s*\(\s*/, "").replace(/\s*\).*$/g, "");
				console.log('ret is '+ret);
			}else{
				ret = "0 0";
			}
		}else if(attr == "scale"){
			if(transform.match(myRe)){
				ret = transform.replace(/^.*scale\s*\(\s*/, "").replace(/\s*\).*$/g, "");
			}else{
				ret = "1 1";
			}
		}else if(attr == "rotate"){
			if(transform.match(myRe)){
				ret = transform.replace(/^.*rotate\s*\(\s*/, "").replace(/\s*\).*$/g, "").replace(/\s+.*$/g, "");
			}else{
				ret = "0";
			}
		}else{
			alert("error not support "+attr);
			return 1;
		}
	}else{
		if(attr == "translate"){
			ret = "0 0";
		}else if(attr == "scale"){
			ret = "1 1";
		}else if(attr == "rotate"){
			ret = "0"
		}
	}
	return ret;
}


//onkeyup='translate_it(target, this.value)'
function transform_it(target, value, typee){
	//target.setAttribute('transform', );
	var trans="";
	var transform = d3.select(target).attr('transform');
	console.log('raw transform is '+transform);
	var transform_new="";
	value = value.replace(/^\s+/, "").replace(/\s+$/, "");
	console.log('value is '+value);
	
	if(typee == "rotate"){
		console.log('rotate it');
		var x = d3.select(target).attr("x");
		var y = d3.select(target).attr("y");
		console.log('x ' +x +"; y is "+y);
		if( x || y ){
			value = value + " " + x + " " + y;
		}
	}
	var transform_tmp = typee + "( " + value + " ) ";
	if(transform){
		var myRe = new RegExp(typee);
		var match = myRe.exec(transform);
		if(match){
			myRe= new RegExp(typee+"[^a-z]+", "g");		
			transform_new = transform.replace(myRe, transform_tmp);
			console.log('1');
		}else{
			transform_new = transform + " " + transform_tmp;
			console.log('2');
		}
	}else{
		transform_new = transform_tmp;
		console.log('3');
	}
	
	target.setAttribute('transform', transform_new);
	console.log('new transform is '+transform_new);
}

//change_style(target, this.value, \"font-size\")' 
function change_style(target, value, style_name){
	console.log('2target is '+target);
	console.log('value is '+value)
	console.log('style name is '+style_name);
	if(style_name == "value"){
		target.innerHTML=value;
	}else{
		target.style[style_name]=value;
	}
}
//onclick='move_feature(\"shift_xy\", \"xy_step\", \"x\", \"+\")'
function caculate_translate(target, shift_xy_id, xy_step_id, xy, change){
	//console.log(d3.select(target).attr("transform"));
	var transform = d3.select(target).attr("transform"); //rotate(30 200 314)  translate(-10 -10) scale(1 0.5) 
	console.log('transform is '+transform);
	if(transform == null){transform=""}
	//var x = d3.select(target).attr("x");
	//var y = d3.select(target).attr("y");
	//var scale_x=1;
	//var scale_y=1;
	//var shift_xy = document.getElementById(shift_xy_id).value;
	var step = document.getElementById(xy_step_id).value;
	step = parseFloat(step);
	console.log('step is ' +step);
	var shift_x=0;
	var shift_y=0;
	
	var myRe = /translate\s*\(\s*([-\.\d]+)\s+([-\.\d]+)\s*\)/g;
	var translate = myRe.exec(transform);
	var translate_str;
	if(translate){
		translate_str = "translate\s*\S+\s+\S+";
		shift_x = translate[1];
		shift_y = translate[2];
		console.log('1 shift_x is '+shift_x);
		console.log('1 shift_y is '+shift_y);
		shift_x = parseFloat(shift_x);
		shift_y = parseFloat(shift_y);
	}
	
	if(change == "+"){
		change = 1;
		console.log('change is +');
	}else if(change == "-"){
		change = -1;
		console.log('change is -');
	}else{
		alert('error, need + or - , not '+change);
		return 1;
	}
	console.log('change is '+change);
	if(xy == "x" ){
		console.log(shift_x+"+="+step+"*"+change);
		shift_x = shift_x + step * change;
		console.log(shift_x);
	}else if(xy == "y"){
		console.log(shift_y+"+="+step+"*"+change);
		shift_y = shift_y + step * change;
		console.log(shift_y);
	}else{
		alert('error in move_feature '+ xy);
		return 1;
	}
	document.getElementById(shift_xy_id).value = shift_x + " " + shift_y;
	console.log('shift_xy value is now '+ document.getElementById(shift_xy_id).value);
	
	//transform_it(target, value, type)
	transform_it(target, shift_x + " " + shift_y, 'translate');
		
}

function saveSVG(classname, width, height, bgcolor){
	console.log("xxxxxx");
	rect2.style.display="none";
	rect1.style.display="none";
	var myDate=new Date();
	myDate=myDate.toString().replace(/\s+/g,"_").replace(/\(.*\)/g, "");
	//var flag = document.getElementById(svg_id).getElementsByClassName("svg-pan-zoom_viewport");
	//var svg = document.getElementsByClassName(classname)[0];
	var textToWrite = document.getElementsByClassName(classname)[0].innerHTML;
	textToWrite = textToWrite.replace(/<g id="viewport->[^>]+>/,"").replace(/<g id="svg-pan-zoom-controls.*<\/g>/, "");
	textToWrite="<svg id=\"demo-tiger\" width=\""+width+"\" height=\""+height+"\" style=\"background-color:"+bgcolor+"\" ersion=\"1.1\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\" >"+textToWrite+"</svg>";

    var textFileAsBlob = new Blob([textToWrite], {type:'text/plain'});
    //var fileNameToSaveAs = document.getElementById("inputFileNameToSaveAs").value;
	var  fileNameToSaveAs="out."+myDate+".svg";
      var downloadLink = document.createElement("a");
    downloadLink.download = fileNameToSaveAs;
    downloadLink.innerHTML = "Download File";
    if (window.webkitURL != null)
    {
        // Chrome allows the link to be clicked
        // without actually adding it to the DOM.
        downloadLink.href = window.webkitURL.createObjectURL(textFileAsBlob);
    }
    else
    {
        // Firefox requires the link to be added to the DOM
        // before it can be clicked.
        downloadLink.href = window.URL.createObjectURL(textFileAsBlob);
        downloadLink.onclick = destroyClickedElement;
        downloadLink.style.display = "none";
        document.body.appendChild(downloadLink);
    }

    downloadLink.click();	
	rect2.style.display="block";
	rect1.style.display="block";
}




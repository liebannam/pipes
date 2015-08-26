#include "/usr/local/share/povray-3.7/include/colors.inc"
#include "/usr/local/share/povray-3.7/include/metals.inc"
 #include "/usr/local/share/povray-3.7/include/textures.inc"
#declare ro = 1.05;// inner and outer radii for pipes, if drawn
#declare ri = .99;
//coordinate arrows if you're lost...
#declare px = <1,1,0>;
#declare py = <0,155/255,0>;
#declare pz = <0,0,205/255>;
/*
union{
union { cylinder { 0*x,  1*x,  .05 open pigment{rgb px }}cone {1*x,  .1, 1.1*x, 0.0 pigment{rgb px}}}
union { cylinder { 0*y,  1*y,  .05 open pigment{ rgb py }}cone {1*y,  .1, 1.1*y, 0.0 pigment{rgb py}}}
union { cylinder { 0*z,  1*z,  .05 open pigment{ rgb pz }}cone {1*z,  .1, 1.1*z, 0.0 pigment{rgb pz}}}
}*/
///////////from Chris Rycroft's Voro++ plotting tools:
// Right-handed coordinate system in which the z-axis points upwards
camera { 
//	orthographic
	//location <6,-15,1.5>
	sky z
//	up y*80
//	right -80*x*image_width/image_height
	right -0.24*x*image_width/image_height
	up 0.24*z
	right -0.20*x*image_width/image_height
	up 0.20*z
//	location <0,15,0>
//	look_at <0,0,0>
//	location<15,-5,1>
//	right -80*x*image_width/image_height
//	up 80*z
// for 7deSeptiembbre
	location <55, 60, 500>
	look_at <55,60,1>
//for dumb 3pipe.inp
//	location <70, -32, 25>
//	look_at <15,0,5>
//for slightly less dumb 3pipe2.inp
//	location <5, -52, 40>
//	look_at <7,0,5>
//bees
//	location <6, -32, 25>
//	look_at <4,-5,5>

}

/*
camera{ orthographic angle 45
        location <0,0,10>
        sky z
	right -0.24*x*image_width/image_height
	look_at  <0,0,0>
//        right x*image_width/image_height
 //       translate <0,2.00,0>
      } 
*/
// Two lights with slightly different colors
light_source{<-8,20,40> color rgb <0.77,0.75,0.75>}
light_source{<25,-12,12> color rgb <0.38,0.40,0.40>}
/////////////////////
/*union {
cylinder {
    0*z,  1*z,  .05
    open
    texture{
       pigment{ Blue }
    }
}
cone {
  1*z,  .1,
  1.1*z, 0.0
   texture{pigment{ Blue}}
}
}*/
#declare H2Oblues_old =
color_map {
       	[0.0 color Black] 
//	[0.2 color rgb <48/255, 128/255, 20/255>]//sap green
	[0.7 color rgb <25/255, 25/255, 112/255>] //midnight blue 
	[0.8 color rgb <70/255, 130/255, 180/255>]//steel blue
 //      	[0.9 color rgb <255/255, 255/255, 240/255>] //ivory
	[1.0 color rgb <100/225, 149/255,247/255>] //cornflower blue 
  //     	[0.95 color Red] 
  //     	[1.0 color rgb <255/255, 255/255, 240/255>] //ivory

}

#declare H2Oblues=
color_map {
//    [0.0 color rgb <255/255, 255/255, 240/255>] //ivory
//	[0.00 color rgb <100/225, 149/255,247/255>] //cornflower blue 
//	[0.00 color rgb <70/255, 130/255, 180/255>]//steel blue
	[0.00 color rgb <25/255, 25/255, 112/255>] //midnight blue 	
	[0.01 color rgb <75./255.,0.,130./255.>]//indigo
	[0.05 color Black]     
	[0.7 color Red]
	[1.0 color Yellow]
//	[0.7 color rgb <48/255, 128/255, 20/255>]//sap green
}

#declare H2Oblues=
color_map {
  //  [0.0 color rgb <255/255, 255/255, 240/255>] //ivory
//	[0.000 color rgb <100/225, 149/255,247/255>] //cornflower blue 
//	[0.00 color rgb <70/255, 130/255, 180/255>]//steel blue
	[0.00 color rgb <25/255, 25/255, 112/255>] //midnight blue 	
//	[0.5 color Black]     
	[0.2 color rgb <75./255.,0.,130./255.>]//indigo
	[0.7 color Red]
	[1.0 color Yellow]}
#declare H2Ofalse_old = 
color_map {
	 [0.0 color rgb <255/255, 255/255, 240/255>] //ivory
	 [0.0125 color rgb <70/255, 130/255, 180/255>] //steel blue
	 [0.025 color rgb <25/255, 25/255, 112/255>] //midnight blue 
//	 [0.1 color rgb <100/225, 149/255,247/255>] //cornflower blue 
	 [0.03725 color rgb <48/255, 128/255, 20/255>]//sap green
	 [1.0 color Red]
}

#declare H2Ofalse = 
color_map {
	 [0.0 color rgb <255/255, 255/255, 240/255>] //ivory
     [0.1 color rgb <100/225, 149/255,247/255>] //cornflower blue 
	 [0.25 color rgb <25/255, 25/255, 112/255>] //midnight blue 
	 [0.6 color rgb <70/255, 130/255, 180/255>] //steel blue
	 [0.8 color rgb <48/255, 128/255, 20/255>]//sap green
     [0.9 color Yellow]
	 [1.0 color Red]
}


#include "plottmp.pov"

//plane { <0, 0, 1>, 0 pigment { checker color White, color Black }}
//plane { <0, 0, 1>, 0 pigment {color Black}}
/*
#declare cbarmax = 30;
union {
 box {<0,0,0>,<1,cbarmax,0.001>
   pigment {
      gradient y
      color_map {H2Ofalse}
      scale cbarmax}
}
box{<-.5,0-.2,0><1.5,0, 0.001> pigment{color Black}}
text{
	ttf "/usr/local/share/povray-3.7/include/timrom.ttf" "0" 2,0
	pigment{rgb 0}
	scale 2
	translate <2.,-.5, 0>
}
box{<-.5,cbarmax/2,0><1.5,cbarmax/2+.2, 0.001> pigment{color Black}}
text{
	ttf "/usr/local/share/povray-3.7/include/timrom.ttf" "0.5" 2,0
	pigment{rgb 0}
	scale 2.
	translate <2.,cbarmax/2-.5, 0>
}
box{<-.5,cbarmax,0><1.5,cbarmax+.2, 0.001> pigment{color Black}}
text{
	ttf "/usr/local/share/povray-3.7/include/timrom.ttf" "1.0" 2,0
	pigment{rgb 0}
	scale 2
	translate <2.,cbarmax-.5, 0>
}

translate<10,60,0>

}*/

//[0.0 color rgb <100/225, 149/255,247/255>] //cornflower blue 
// [0.0 color rgb <48/255, 128/255, 20/255>]//sap green
// [0.0 color rgb <139/255, 69/255, 19/255>]//chocolate 4
 //[0.0 color rgb<0,1,0>] //dark green
// [0.0 color rgb <25/255, 25/255, 112/255>] //midnight blue 
 //[0.0 color Black]  
// [0.8 color rgb <70/255, 130/255, 180/255>] //steel blue
// [1.0 color rgb <255/255, 255/255, 240/255>] //ivory
// [0.0 color rgb <25/255, 25/255, 112/255>] //midnight blue 


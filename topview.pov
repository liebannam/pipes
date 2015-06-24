#include "/usr/local/share/povray-3.7/include/colors.inc"
#include "/usr/local/share/povray-3.7/include/metals.inc"
 #include "/usr/local/share/povray-3.7/include/textures.inc"
#declare ro = 1.05;// inner and outer radii for pipes, if drawn
#declare ri = .99;
///////////from Chris Rycroft's Voro++ plotting tools:
// Right-handed coordinate system in which the z-axis points upwards
camera { 
	orthographic
	//location <6,-15,1.5>
//for simple3pipes
	sky z
	up y*25
	right -25*x*image_width/image_height
	location <15, 0, 25>
	look_at <15,0,0>


//	right -0.24*x*image_width/image_height
//	up 0.24*z
//	location <0,15,0>
//	look_at <0,0,0>
//	location<15,-5,1>
//	right -80*x*image_width/image_height
//	up 80*z
// for 7deSeptiembbre
	 up y*50
	right -50*x*image_width/image_height
	location <55, 50, 250>
	look_at <55,50,1>
//for 3pipes2
//	up y*10
//	right -10*x*image_width/image_height
//	location<5,0,20>
//	look_at<5,0,1>
}


camera{ orthographic angle 45
  
        sky z
	right -0.24*x*image_width/image_height
//for 7 de Septiembre
		look_at  <55,50,0>
  	location <55,50,180>
//for 3pipes2
//	look_at<7,0,0>
//	location<7,0,25>
        right -x*image_width/image_height
 //       translate <0,2.00,0>
      } 

// Two lights with slightly different colors
light_source{<-8,-20,30> color rgb <0.77,0.75,0.75>}
light_source{<25,-12,12> color rgb <0.38,0.40,0.40>}
/////////////////////
//arrow at origin if you get lost
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
       	[0.0 color rgb <255/255, 255/255, 240/255>] //ivory
//	[0.2 color rgb <48/255, 128/255, 20/255>]//sap green
	[0.5 color rgb <70/255, 130/255, 180/255>]//steel blue
	[1.0 color rgb <25/255, 25/255, 112/255>] //midnight blue 
}

#declare H2Oblues_old2=
color_map {
    [0.0 color rgb <255/255, 255/255, 240/255>] //ivory
	[0.1 color rgb <100/225, 149/255,247/255>] //cornflower blue 
	[0.5 color rgb <70/255, 130/255, 180/255>]//steel blue
	[0.9 color rgb <25/255, 25/255, 112/255>] //midnight blue 	
	[1.0 color rgb <48/255, 128/255, 20/255>]//sap green
}
#declare H2Oblues=
color_map {
  //  [0.0 color rgb <255/255, 255/255, 240/255>] //ivory
	[0.000 color rgb <100/225, 149/255,247/255>] //cornflower blue 
	[0.05 color rgb <70/255, 130/255, 180/255>]//steel blue
	[0.2 color rgb <25/255, 25/255, 112/255>] //midnight blue 	
//	[0.5 color Black]     
	[0.5 color rgb <75./255.,0.,130./255.>]//indigo
	[0.7 color Red]
	[1.0 color Yellow]}

#declare H2Ofalse = 
color_map {
	 [0.0 color rgb <255/255, 255/255, 240/255>] //ivory
	 [0.03 color rgb <70/255, 130/255, 180/255>] //steel blue
	 [0.4 color rgb <25/255, 25/255, 112/255>] //midnight blue 
//	 [0.1 color rgb <100/225, 149/255,247/255>] //cornflower blue 
	 [0.5 color rgb <48/255, 128/255, 20/255>]//sap green
     [0.6 color Yellow]
	 [1.0 color Red]
}

#declare negpres = 
color_map {
//	 [0.03 color rgb <70/255, 130/255, 180/255>] //steel blue
//	 [0.5 color rgb <48/255, 128/255, 20/255>]//sap green
	 [0. color Red]
         [0.6 color Yellow]
	 [0.99 color rgb <25/255, 25/255, 112/255>] //midnight blue 
	 [1.0 color rgb <100/225, 149/255,247/255>] //cornflower blue 
}
//7deseptiembre colormap
/*#declare H2Ofalse = 
color_map {
	 [0.0 color rgb <255/255, 255/255, 240/255>] //ivory
	 [0.0125 color rgb <70/255, 130/255, 180/255>] //steel blue
	 [0.025 color rgb <25/255, 25/255, 112/255>] //midnight blue 
//	 [0.1 color rgb <100/225, 149/255,247/255>] //cornflower blue 
	 [0.03725 color rgb <48/255, 128/255, 20/255>]//sap green
         [0.05 color Yellow]
	 [1.0 color Red]
}*/

#include "plottmp.pov"

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


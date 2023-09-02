//external.xm
/*
 * Original shader from: https://www.shadertoy.com/view/XtfGzj
 */

#ifdef GL_ES
precision highp float;
#endif

// glslsandbox uniforms
uniform float time;
uniform vec2 resolution;

// shadertoy globals
#define iTime time/1000
#define iResolution resolution
const vec4  iMouse = vec4(0.0);

// Emulate a black texture
#define texture(s, uv) vec4(0.0)

// --------[ Original ShaderToy begins here ]---------- //
// "Space Racing Lite"  

// Distance function and initial design for space car is by eiffie:
// https://www.shadertoy.com/view/ldjGRh
// the rest is by me but he also helped to optimize the code.

// I removed some features by default because the original was crashing the Shadertoy browser
// for win7 users - try commenting this lines to see if the full version compiles for you: 

//#define LOW_QUALITY // No reflections, no shadows, no planet, reduced ray steps & detail
//#define NO_HUD 

#define LOOP_BREAKS // Could speed up, speed down, or just make your browser crash!



#ifdef LOW_QUALITY
	#define RAY_STEPS 65
	#define SHADOW_STEPS 0
	#define ITERATIONS 5
	#define MAX_DIST 30.
#else
	#define RAY_STEPS 75
	#define SHADOW_STEPS 40
	#define ITERATIONS 6
	#define MAX_DIST 35.
#endif
#define LIGHT_COLOR vec3(1.,.85,.6)
#define AMBIENT_COLOR vec3(.7,.85,1.)
#define SUN_COLOR vec3(1.,.8,.5)
#define TUBE_COLOR vec3(1.,.6,.25)*1.2
#define CAR_COLOR vec3(.4,.7,1.)
#define TURBINES_COLOR vec3(.6,.75,1.)
#define HUD_COLOR vec3(0.6,1.,0.3)
#define PLANET_COLOR vec3(101., 153., 189.)/256.
	
#define CAR_SCALE 2.
#define SPECULAR 0.4
#define DIFFUSE  2.0
#define AMBIENT  0.4

#define BRIGHTNESS .9
#define GAMMA 1.1
#define SATURATION .85


#define DETAIL .004
#define SPEED 8.
#define t (mod(iTime,500.)+10.)*SPEED

#define LIGHTDIR normalize(vec3(0.6,-0.2,-1.))

// ------------------------------------------------------------------
//    Global Variables
// ------------------------------------------------------------------

float FOLD=2.; // controls fractal fold and track width
const vec3 planetpos=vec3(-3.5,1.,-5.); // planet position
const vec2 tubepos=vec2(0.35,0.); // light tubes pos relative to FOLD
mat2 trmx=mat2(0.);//the turbine spinner
float det=0.; // detail level (varies with distance)
float backcam=0.; // back cam flag
vec3 carpos=vec3(0.); // car position
vec3 carvec=vec3(0.); // car pointing vector
mat3 carrot=mat3(0.); // car rotation
float hitcar=0.; // ray-car hit flag
mat2 fractrot=mat2(0.); // rot mat for fractal (set in main)
mat2 cartilt=mat2(0.); // secondary car rotation
float minT=0., minL=0.; // min distance traps for glows of tube and turbines
float ref=0.; // reflection flag
float tubeinterval=0.; // tube tiling (for glow and lighting)


// ------------------------------------------------------------------
//    General functions 
// ------------------------------------------------------------------

mat2 rot(float a) {
	return mat2(cos(a),sin(a),-sin(a),cos(a));	
}

mat3 lookat(vec3 fw,vec3 up){
	fw=normalize(fw);vec3 rt=normalize(cross(fw,normalize(up)));return mat3(rt,cross(rt,fw),fw);
}

float smin(float a,float b,float k){return -log(exp(-k*a)+exp(-k*b))/k;}//from iq

float Sphere(vec3 p, vec3 rd, float r){
	float b = dot( -p, rd );
	float inner = b * b - dot( p, p ) + r * r;
	if( inner < 0.0 ) return -1.0;
	return b - sqrt( abs(inner) );
}


// ------------------------------------------------------------------
//    Track 
// ------------------------------------------------------------------

// the track function, just some curves
vec3 path(float ti) { 
	float freq=.5, amp=1.; // for trying different settings
	ti*=freq;
	float x=cos(cos(ti*.35682)+ti*.2823)*cos(ti*.1322)*1.5;
	float y=sin(ti*.166453)*4.+cos(cos(ti*.125465)+ti*.17354)*cos(ti*.05123)*2.;
	vec3  p=vec3(x,y,ti/freq);
	return p;
}

// see if we are in the tunnel, and used also by distance function
float tunnel(float z) {
return abs(100.-mod(z+15.,200.))-30.;
}


// ------------------------------------------------------------------
//    DE functions
// ------------------------------------------------------------------


// carcarspacecar by eiffie // (with some mods by Kali)
// a reconfig of the carcar's tires (someday I will have to do an animation from the original to this)
//the DE looks a bit long but there are actually fewer instructions than the car

float carcarDE(in vec3 p0){
	p0.xy*=cartilt;
	p0*=CAR_SCALE;
	vec3 p=p0;
	p.y+=1.24;
	float r=length(p.yz);
	float d=length(max(vec3(abs(p.x)-0.35,r-1.92,-p.y+1.4),0.0))-0.05;
	d=max(d,p.z-1.05);
	p=p0+vec3(0.0,-0.22,0.39);
	p.xz=abs(p.xz);
	p.xyz-=vec3(0.72,0.0,1.06);
	float w1=0.23,w2=0.24;
	if(p0.z<0.0){//this is discontinuous yet works unless you stand in front of the rear turbines
		w1=0.23,w2=0.05; //where you would get sucked into the blades anyways
		p=p.xzy; //so i'm comfortable with it :)
	} 
	r=length(p.xy);
	d=smin(d,length(vec2(max(abs(p.z)-w2,0.0),r-w1))-0.02,8.0);//done with the car shape, the rest is just turbines and could be replaced with lights or something
	d=min(d,(length(p*vec3(1.,1.,.6))-.08-p0.z*.03));
	p.xy=trmx*p.xy;//spin
	float d2=min(abs(p.x),abs(p.y))*.15;//4 blades
	//p.xy=mat2(0.707,-0.707,0.707,0.707)*p.xy;//45 degree turn
	//d2=min(d2,min(abs(p.x),abs(p.y))*.2);//8 blades
	d2=max(r-w1-.05,max(d2-0.005,abs(p.z)-w2+0.04));
	d2=min(d2,(length(p)-.05-p0.z*.035)*.07);
	d2=min(d2,max(d+.02,max(abs(p0.y-.07),abs(p0.x)-.4+min(0.,p0.z)))*.18);
	minL=min(minL,d2);//catch the minimum distance to the glowing parts
	// I used d2 only for glows (Kali)
	return d/CAR_SCALE;// *carScale
}


vec3 carcarCE(in vec3 p0){//coloring
	p0*=CAR_SCALE;
	vec4 trpc=vec4(0.);//color trap (taking samples when finding the norm)// not for now (Kali)

	//right here you should inv-transform p0 as it gets used later
	//p0=(p0-carPos)*carRot/carScale;//or something like that??
	p0.xy*=cartilt;
	vec3 p=p0;
	p.y+=1.24;
	float r=length(p.yz);
	float d=length(max(vec3(abs(p.x)-0.35,r-1.92,-p.y+1.4),0.0))-0.05;
	d=max(d,p.z-1.0);
	p=p0+vec3(0.0,-0.22,0.39);
	p.xz=abs(p.xz);
	p.xyz-=vec3(0.72,0.0,1.06);
	float w1=0.2,w2=0.24;
	if(p0.z<0.0){//this is discontinuous yet works unless you stand in front of the rear turbines
		w1=0.23,w2=0.05; //where you would get sucked into the blades anyways
		p=p.xzy; //so i'm comfortable with it :)
	}
	r=length(p.xy);
	d=smin(d,length(vec2(max(abs(p.z)-w2,0.0),r-w1))-0.02,8.0);//done with the car shape, the rest is just turbines and could be replaced with lights or something
	p.xy=trmx*p.xy;
	float d2=min(abs(p.x),abs(p.y));//4 blades
	p.xy=mat2(0.707,-0.707,0.707,0.707)*p.xy;//45 degrees
	d2=min(d2,min(abs(p.x),abs(p.y)));//8 blades
	d2=max(r-w1+0.02,max(d2-0.005,abs(p.z)-w2+0.04));
	//up to here it is the same as DE, now accumulate colors
	if(d2<d){d=d2;trpc+=vec4(1.,0.6,0.3,256.0);}//turbines
	else {//the car's body
		p0.x=abs(p0.x);
		if((abs(p0.y-0.58)>0.05-p0.z*0.09 || p0.z>0.25) && 
		   length(max(abs(p0.xz+vec2(-p0.z*.03,-0.5))-vec2(0.15+p0.z*.03,0.4),0.0))>0.1)
			trpc+=vec4(CAR_COLOR,16.0);
		else trpc+=vec4(CAR_COLOR*.4,2.0);//the windsheild
	}
	return trpc.xyz; // *carScale
}

//-------------------------------------------

// DE for tubes
float tubes(vec3 pos) {
	pos.x=abs(pos.x)-tubepos.x-FOLD;
	pos.y+=tubepos.y;
	return (length(pos.xy)-.05);
}

// ------------------------------------------------------------------
//    SCENE DE
// ------------------------------------------------------------------

float de(vec3 pos) {
	vec3 carp=pos-carpos; // scale car coordinates
	carp=carrot*carp; // rotate car
	pos.xy-=path(pos.z).xy; // transform coordinates to follow track
	FOLD=1.7+pow(abs(100.-mod(pos.z,200.))/100.,2.)*2.; //varies fractal fold & track width
	pos.x-=FOLD;
	float hid=0.;
	vec3 tpos=pos;
	tpos.z=abs(2.-mod(tpos.z,4.));
	vec4 p=vec4(tpos,1.);
	for (int i=0; i<ITERATIONS; i++) { // calculate fractal
		p.xz=clamp(p.xz,-vec2(FOLD,2.),vec2(FOLD,2.))*2.0-p.xz;
		p=p*2.5/clamp(dot(p.xyz,p.xyz),.5,1.)-vec4(1.2,0.5,0.,-0.5);
		p.xy*=fractrot;
	}
	pos.x+=FOLD;
	float fr=min(max(pos.y+.4,abs(pos.x)-.15*FOLD),(max(p.x,p.y)-1.)/p.w); // fractal+pit
	float tub=tubes(pos);  
	minT=min(minT,tub*.5); // trap min distance to tubes	
	float car=carcarDE(carp); 
	float d=tub;
	d=min(fr,d);
	d=min(d,max(abs(pos.y-1.35+cos(3.1416+pos.x*.8)*.5)-.1,tunnel(pos.z))); // tunnel DE
	if (ref<1.) d=min(d,car);
	d=max(d,abs(pos.x)-FOLD*2.);
	if (car<det) hitcar=1.; // ray hits the car!
	return d;
}


// ------------------------------------------------------------------
//    General Shading Functions
// ------------------------------------------------------------------


vec3 normal(vec3 p) {
	vec3 e = vec3(0.0,det,0.0);
	
	return normalize(vec3(
			de(p+e.yxx)-de(p-e.yxx),
			de(p+e.xyx)-de(p-e.xyx),
			de(p+e.xxy)-de(p-e.xxy)
			)
		);	
}

#ifndef LOW_QUALITY 

float shadow(vec3 pos, vec3 sdir) {
	float sh=1.0;
	float totdist = DETAIL*10.;
	float dist=1000.;
		for (int steps=0; steps<SHADOW_STEPS; steps++) {
			if (totdist<MAX_DIST && dist>DETAIL) {
				vec3 p = pos - totdist * sdir;
				dist = de(p);
				sh = min(sh, 10.*max(0.0,dist)/totdist);
				sh*= sign(max(0.,dist-DETAIL));
				totdist += max(0.02,dist);
			}
#ifdef LOOP_BREAKS		
		else break;
#endif
		}
	
	return clamp(sh,0.,1.0);
}

#endif

float calcAO(vec3 pos, vec3 nor) {
	float hr,dd,aoi=0.,sca=1.,totao=0.;
	hr = .075*aoi*aoi;dd = de(nor * hr + pos);totao += (hr-dd)*sca;sca*=.6;aoi++;
	hr = .075*aoi*aoi;dd = de(nor * hr + pos);totao += (hr-dd)*sca;sca*=.55;aoi++;
	hr = .075*aoi*aoi;dd = de(nor * hr + pos);totao += (hr-dd)*sca;sca*=.55;aoi++;
	hr = .075*aoi*aoi;dd = de(nor * hr + pos);totao += (hr-dd)*sca;sca*=.55;aoi++;
	return clamp( 1.0 - 4.*totao, 0., 1.0 );
}


// ------------------------------------------------------------------
//    Light and Coloring
// ------------------------------------------------------------------



vec3 shade(in vec3 p, in vec3 dir, in vec3 n) {

	float savehitcar=hitcar;

	vec3 trackoffset=-vec3(path(p.z).xy,0.);
	vec3 pos=p;
	vec3 col=vec3(.5); // main color
	vec3 carp=pos-carpos; //scaled coordinates for the car
	carp=carrot*carp; // rotate car
	pos+=trackoffset; // apply track transformation to the coordinates
	// track lines
	if (pos.y<.5) col+=pow(max(0.,.2-abs(pos.x))/.2*abs(sin(pos.z*2.)),8.)*TUBE_COLOR*2.;
	pos.x=abs(pos.x);
	// fake AO for the tunnel's upper corners
	if(tunnel(pos.z)<0.)
		col*=1.-pow(max(0.5,1.-length(pos.xy+vec2(-FOLD*1.5,-.85))),5.)*max(0.,1.+pos.y);
	if (tubes(pos)<det) col=TUBE_COLOR; // hit tubes
	if (carcarDE(carp)<det) col=carcarCE(carp); // hit car, get coloring

	float ao=calcAO(p,n); // calc AO
	float camlight=max(0.,dot(dir,-n)); // camlight used for ambient

	// --- Tube lights ---

	vec3 tpos1=vec3((tubepos+vec2(FOLD,0.)),0.)+trackoffset; // get tube positions
	vec3 tpos2=tpos1-vec3((tubepos.x+FOLD)*2.,0.,0.);
	// light direction
	vec3 tube1lightdir=normalize(vec3(p.xy,0.)+vec3(tpos1.xy,0.)); 
	vec3 tube2lightdir=normalize(vec3(p.xy,0.)+vec3(tpos2.xy,0.));
	// light falloffs
	float falloff1,falloff2;	
	if (savehitcar>0.) {
		falloff1=pow(max(0.,1.-.15*distance(vec3(p.xy,0.),vec3(-tpos1.xy,0.))),4.);
		falloff2=pow(max(0.,1.-.15*distance(vec3(p.xy,0.),vec3(-tpos2.xy,0.))),4.);
	} else {
		falloff1=pow(max(0.,1.-.2*distance(vec3(p.xy,0.),vec3(-tpos1.xy,0.))),4.);
		falloff2=pow(max(0.,1.-.2*distance(vec3(p.xy,0.),vec3(-tpos2.xy,0.))),4.);
	}
	
	float diff, spec;
	
	vec3 r=reflect(LIGHTDIR,n);
	
	// tube1 calcs
	diff=max(0.,dot(tube1lightdir,-n)); 
	diff+=max(0.,dot(normalize(tube1lightdir+vec3(0.,0.,.2)),-n))*.5; // add 2 more 
	diff+=max(0.,dot(normalize(tube1lightdir-vec3(0.,0.,.2)),-n))*.5; // with Z shifted
	spec=pow(max(0.,dot(tube1lightdir+vec3(0.,0.,.4),r)),15.)*.7;
	spec+=pow(max(0.,dot(tube1lightdir-vec3(0.,0.,.4),r)),15.)*.7;
	float tl1=(falloff1*ao+diff+spec)*falloff1;

	// tube2 calcs
	diff=max(0.,dot(tube2lightdir,-n));
	diff+=max(0.,dot(normalize(tube2lightdir+vec3(0.,0.,.2)),-n))*.5;
	diff+=max(0.,dot(normalize(tube2lightdir-vec3(0.,0.,.2)),-n))*.5;
	spec=pow(max(0.,dot(tube2lightdir+vec3(0.,0.,.4),r)),15.)*.7;
	spec+=pow(max(0.,dot(tube2lightdir-vec3(0.,0.,.4),r)),15.)*.7;
	float tl2=(falloff2*ao+diff+spec)*falloff2;

	// sum tube lights - add ambient - apply tube intervall
	vec3 tl=((tl1+tl2)*(.5+tubeinterval*.5))*TUBE_COLOR;//*(1.+tun*.5);


	// --- Car lights ---

	// get the car turbines direction (aproximate)
	vec3 carlightdir1=normalize(p-carpos+vec3(.2,0.06,.15));
	vec3 carlightdir2=normalize(p-carpos+vec3(-.2,0.06,.15));
	vec3 carlightdir3=normalize(p-carpos+vec3(.2,0.06,-.35));
	vec3 carlightdir4=normalize(p-carpos+vec3(-.2,0.06,-.35));

	float cfalloff=pow(max(0.,1.-.1*distance(p,carpos)),13.); // car light falloff

	// accumulate diffuse
	diff=max(0.,dot(carlightdir1,-n))*.5;
	diff+=max(0.,dot(carlightdir2,-n))*.5;
	diff+=max(0.,dot(carlightdir3,-n))*.5;
	diff+=max(0.,dot(carlightdir4,-n))*.5;

	if (savehitcar<1.) diff*=clamp(1.-carlightdir1.y,0.,1.);
	
	// add ambient and save car lighting
	vec3 cl=TURBINES_COLOR*((diff+spec*.0)*cfalloff+cfalloff*.3)*1.2;
 	
	// --- Main light ---
	
#ifdef LOW_QUALITY
	float sh=ao;
#else
	float sh=shadow(p, LIGHTDIR); // get shadow
#endif

	diff=max(0.,dot(LIGHTDIR,-n))*sh*1.3; // diffuse
	float amb=(.4+.6*camlight)*.6; // ambient+camlight
	spec=pow(max(0.,dot(dir,-r))*sh,20.)*SPECULAR; //specular
	if (savehitcar>0.) {diff*=.8;amb*=.3;}
	vec3 l=(amb*ao*AMBIENT_COLOR+diff*LIGHT_COLOR)+spec*LIGHT_COLOR;	

	if (col==TUBE_COLOR) l=.3+vec3(camlight)*.7; // special lighting for tubes

	return col*(l+cl+tl); // apply all lights to the color
}

// the planet shading...
// very simple and fast made, but for low res windowed mode it does the job :)
vec3 shadeplanet(vec3 pos, vec3 k) { 

	vec3 n=normalize(planetpos+pos+.2); // tweaked sphere normal
	float c=max(0.,dot(LIGHTDIR,normalize(k-n))); // surface value
	vec3 col=PLANET_COLOR+vec3(c,c*c,c*c*c)*.7; // surface coloring
	// add some noise
	float noi=max(0.,texture(iChannel1,n.yz*.5).x-.6);
	noi+=max(0.,texture(iChannel1,n.yz).x-.6);
	noi+=max(0.,texture(iChannel1,n.yz*2.).x-.6);
	col+=noi*(1.5-c)*.7;
	return col*max(0.,dot(LIGHTDIR,-n)); // diff light
}

// ------------------------------------------------------------------
//    Raymarching + FX rendering
// ------------------------------------------------------------------


vec3 raymarch(in vec3 from, in vec3 dir) 

{
	hitcar=0.;
	ref=0.;
	float totdist=0.;
	float glow=0.;
	float d=1000.;
	vec3 p=from, col=vec3(0.5);

	float deta=DETAIL*(1.+backcam); // lower detail for HUD cam
	vec3 carp=vec3(0.); // coordinates for car hit
	vec3 carn=vec3(0.); // normal for car
	float cardist=0.; // ray length for car
	vec3 odir=dir; // save original ray direction

	for (int i=0; i<RAY_STEPS; i++) {
		if (d>det && totdist<MAX_DIST) {
			d=de(p);
			p+=d*dir;
			det=max(deta,deta*totdist*.5*(1.+ref)); // scale detail with distance or reflect
			totdist+=d; 
			float gldist=det*8.; // background glow distance 
			if(d<gldist&&totdist<20.) glow+=max(0.,gldist-d)/gldist*exp(-.1*totdist); //accum glow
#ifndef LOW_QUALITY
			if (hitcar>0. && ref<1.) { // hit car, bounce ray (only once)
				p=p-abs(d-det)*dir; // backstep
				carn=normal(p); // save car normal
				carp=p; // save car hit pos
				dir=reflect(dir,carn); // reflect ray
				p+=det*dir*10.; // advance ray
				d=10.; cardist=totdist;
				ref=1.;
			}
#endif
		} 
#ifdef LOOP_BREAKS		
		else break;
#endif
	}

	tubeinterval=abs(1.+cos(p.z*3.14159*.5))*.5; // set light tubes interval
	float cglow=1./(1.0+minL*minL*5000.0); // car glow
	float tglow=1./(1.0+minT*minT*5000.0); // tubes glow
	float l=max(0.,dot(normalize(-dir),normalize(LIGHTDIR))); // lightdir gradient
	vec3 backg=AMBIENT_COLOR*.4*max(0.1,pow(l,5.)); // background
	float lglow=pow(l,50.)*.5+pow(l,200.)*.5; // sun glow

	if (d<.5) { // hit surface
		vec3 norm=normal(p); // get normal
		p=p-abs(d-det)*dir; // backstep
		col=shade(p, dir, norm); // get shading 
		col+=tglow*TUBE_COLOR*pow(tubeinterval,1.5)*2.; // add tube glow
		col = mix(backg, col, exp(-.015*pow(abs(totdist),1.5))); // distance fading

	} else { // hit background
		col=backg; // set color to background
		col+=lglow*SUN_COLOR; // add sun glow
		col+=glow*pow(l,5.)*.035*LIGHT_COLOR; // borders glow
		
#ifdef LOW_QUALITY
		vec3 st = (dir * 3.+ vec3(1.3,2.5,1.25)) * .3;
		for (int i = 0; i < 14; i++) st = abs(st) / dot(st,st) - .9;

		col+= min( 1., pow( min( 5., length(st) ), 3. ) * .0025 ); // add stars
#else
		float planet=Sphere(planetpos,dir, 2.); // raytrace planet

		// kaliset formula - used for stars and planet surface 
		float c;
		if (planet>0.) c=1.; else c=.9; // different params for planet and stars
		vec3 st = (dir * 3.+ vec3(1.3,2.5,1.25)) * .3;
		for (int i = 0; i < 14; i++) st = abs(st) / dot(st,st) - c;

		col+= min( 1., pow( min( 5., length(st) ), 3. ) * .0025 ); // add stars
		
		// planet atmosphere
		col+=PLANET_COLOR*pow(max(0.,dot(dir,normalize(-planetpos))),100.)*150.*(1.-dir.x);
		// planet shading
		if (planet>0.) col=shadeplanet(planet*dir,st);
#endif
		
	}
	// car shading

		// add turbine glows
	
#ifdef LOW_QUALITY
	cglow*=1.15;
#else
	if (ref>0.) {
		ref=0.;
		col=shade(carp,odir,carn)+col*.3; // car shade + reflection
		// I wanted a lighter background for backward reflection
		l=max(0.,dot(normalize(-odir),normalize(LIGHTDIR)));
		backg=AMBIENT_COLOR*.4*max(0.1,pow(l,5.)); 
		col = mix(backg, col,exp(-.015*pow(abs(cardist),1.5))); // distance fading
	}
#endif 

	
	col+=TURBINES_COLOR*pow(abs(cglow),2.)*.4;
	col+=TURBINES_COLOR*cglow*.15;


	return col; 
}

// simple HUD track graph with transparency
vec4 showtrack(vec2 p) {
	p.x+=.25;
	vec2 uv=p;
	float uvpos=uv.x+1.5;
	vec3 pth=path((uvpos-1.5)*30.+t)*.05;
	float curpos=path(t).x*.05;
	float curs=pow(max(0.,1.-length(uv+vec2(0.,curpos))*2.),10.)*max(0.5,sin(iTime*10.))*2.;
	uv.xy=uv.xy-(vec2(uvpos,0.))*rot(pth.x/uvpos);
	float dotline=pow(max(0.,1.-abs(uv.y)*5.),30.);
	float graph=(curs+dotline);
	return vec4((graph+.4)*HUD_COLOR,1.-.5*pow(abs(.025-mod(p.y*2.,.05))/.025,2.));
}


// ------------------------------------------------------------------
//    Main code - Camera - HUD 
// ------------------------------------------------------------------


void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
	minL=minT=1000.; // initialize min distance glows
	fractrot=rot(.5); // mat2D for the fractal formula
	vec3 carpos0=vec3(0.,-0.2,.0); // set relative car pos (unused now, only sets height)
	carpos=vec3(carpos0+vec3(0.,0.,t)); // real car pos
	vec3 carmove=path(carpos.z); carmove.x*=1.+FOLD*.1; // get pos, exagerate x pos based on track width.
	carvec=normalize((path(carpos.z+2.)-carmove)*vec3(FOLD*.25,1.,1.)); // car fwd vector
	carrot=lookat(-carvec,vec3(0.,1.,0.)); // car rotation
	cartilt=rot(-carvec.x*2.); // xy rotation based on distance from center
	carpos.xy+=carmove.xy-vec2(carvec.x,0.)*FOLD*.5; // move away from center when turning
	float tim=iTime*12.0;
	trmx=mat2(cos(tim),-sin(tim),sin(tim),cos(tim));//the turbine spinner

	// --- camera & mouse ---
	
	vec2 uv = fragCoord.xy / iResolution.xy*2.-1.;
	uv.y*=iResolution.y/iResolution.x;
	vec2 mouse=(iMouse.xy/iResolution.xy-.5)*vec2(7.,1.5);
	if (iMouse.z<1.) { // if no mouse, alternate rear and back cam
		mouse=vec2(sin(iTime)*.7,2.+sin(iTime*.2)*.22)
			*min(0.,sign(10.-mod(iTime,20.)));
	}
	vec3 dir=normalize(vec3(uv*.8,1.));

	vec3 campos=vec3(0.,0.2,-.6); // original relative camera position
	//rotate camera with mouse
	campos.yz=(campos.yz-carpos0.yz)*rot(mouse.y)+carpos0.yz;
	campos.xz=(campos.xz-carpos0.xz)*rot(mouse.x)+carpos0.xz;
	campos.x-=carvec.x*FOLD; // follow car x pos a bit when turning

	vec3 from;
	
	float fixcam=5.;
	float mt=mod(t/SPEED,fixcam*2.);
	//fixed cam every 15 seconds, random position, point at car position
	if ((mod(iTime,20.)>15. && iMouse.z<1.)) {
		fixcam*=SPEED;
		from=path(floor(t/fixcam)*fixcam+fixcam*.5);
		//from.x+=.05;// from.y+=.5;
		vec2 fixpos=(texture(iChannel1,vec2(from.z*.21325)).xy-.5)*vec2(FOLD*2.-.3,1.);
		fixpos.x+=sign(fixpos.x)*.3; fixpos.y+=.2;
		from.xy+=fixpos;
		dir=lookat(normalize(carpos-from),vec3(0.,1.,0.))*normalize(dir+vec3(0.,0.,0.5));

	} else { //normal cam
		from=path(t+campos.z)+campos;
		dir.y-=.3*campos.z;
		dir=lookat(normalize(carpos-from),vec3(0.,1.,0.))*normalize(dir);
	}

	vec4 hud=vec4(0.);
	
#ifndef NO_HUD	
	//HUD (hud camera was going to be transparent but won't compile)
	backcam=0.;
		vec2 huv=uv+vec2(.75,.44);
		if (length(huv*huv*huv*vec2(5.,50.))<.05) hud=showtrack(huv*2.); // track HUD
		uv+=vec2(-.75,.44);
		if (length(uv*uv*uv*vec2(5.,50.))<.05) { //set ray data for HUD cam
				backcam=1.;
				uv*=6.;
				dir=normalize(vec3(uv.xy*.6,-1.));
				from=vec3(carvec.x*.5,0.1,0.)+path(t-campos.z*1.7);
				dir=lookat(-normalize(carpos-from),normalize(vec3(0.,1.,0.)))*dir;
				//color+=HUD_COLOR*(vec3(HUDraymarch(from,dir))+.1);
		}

#endif	

	vec3 color=raymarch(from,dir); 	// Raymarch scene
	color=clamp(color,vec3(.0),vec3(1.));
	if (backcam>0.) { //if HUD cam, apply post effect
		color=(.2+pow(length(color),1.7)*.5)*HUD_COLOR
			*(1.-.5*pow(abs(.025-mod(uv.y*.9,.05))/.025,2.))*.9;
	}
	
	color=hud.rgb*hud.a+color*(1.-hud.a);//HUD transparency

	//color adjustments
	color=pow(abs(color),vec3(GAMMA))*BRIGHTNESS;
	color=mix(vec3(length(color)),color,SATURATION);
	fragColor = vec4(color,1.);
}
// --------[ Original ShaderToy ends here ]---------- //

// --------------------------------- //
// --------[ ICON Drawing]---------- //
// --------------------------------- //

const vec2 bitmap_size = vec2(96, 96);
const int[] palette = int[] (
0x00000000, 0x003f2161, 0x003e1e73, 0x003f2a83, 0x004c2b85, 0x005a288c, 0x005d2692, 0x005a2e7f, 0x00542a8f, 0x00562896, 0x004e2b92, 0x004b299b, 0x00523577, 0x00472d95, 0x003e377c, 0x0048377a, 0x003d3495, 0x00413591, 0x0049358f, 0x00484074, 0x00543691, 0x00583f76, 0x003a3897, 0x00654179,
0x00584a6d, 0x0042428d, 0x006a447f, 0x00414783, 0x006039ab, 0x005f37b9, 0x005a429d, 0x00584691, 0x004c43a0, 0x005343a1, 0x00594a8a, 0x005e3eb0, 0x00654a8a, 0x005f41a9, 0x006c4e7d, 0x005c41b6, 0x005a44ae, 0x005742b8, 0x00665283, 0x007547a4, 0x00554f96, 0x004949ad, 0x005747b3, 0x004e4f9c, 
0x005c5886, 0x00664fab, 0x0064539f, 0x005f4fae, 0x005a4eb3, 0x004756a4, 0x00645799, 0x005a52af, 0x004057a7, 0x00525a99, 0x00485d9d, 0x00425e9f, 0x00595aa8, 0x006552cb, 0x006255c3, 0x006157bd, 0x00744dd9, 0x006f58c5, 0x006d52da, 0x006661ac, 0x006e61aa, 0x00625eb9, 0x007157ce, 0x006b5fbb,
0x006563b4, 0x006164b6, 0x00605bd4, 0x007561bc, 0x005f5fcf, 0x006e5cd6, 0x005e62cb, 0x005f65c4, 0x006d61d0, 0x006c6fa6, 0x005266cc, 0x00846fa3, 0x006273a2, 0x006b65cd, 0x007871a3, 0x005a6cbd, 0x006a69c5, 0x00626cbe, 0x006770b0, 0x00726dba, 0x007d5ee8, 0x006a6dbe, 0x007c6eb8, 0x005972b5,
0x00747a95, 0x006172b4, 0x005a79aa, 0x008c7c94, 0x007169e5, 0x00527cb4, 0x007e6ce1, 0x006370de, 0x005b71e1, 0x007470de, 0x007279c5, 0x006d72e1, 0x007076d4, 0x008c6fe4, 0x008177d2, 0x006477de, 0x006b7dca, 0x006a81c5, 0x008677e0, 0x007b82c3, 0x00987bd2, 0x006c85be, 0x005e85c2, 0x007e86b9,
0x006586c1, 0x00687edc, 0x008381cb, 0x007787bd, 0x00607fe0, 0x009483c8, 0x00998bbb, 0x007c84e3, 0x009a93ad, 0x008594b5, 0x007194bb, 0x00698bdd, 0x009086e3, 0x00838cd8, 0x00698ae8, 0x00858fd2, 0x00788ddc, 0x007790d5, 0x007995ca, 0x007295cc, 0x006c94d1, 0x007393d3, 0x008c94c9, 0x007c93d1, 
0x006f92da, 0x00758cec, 0x008595cf, 0x008991ed, 0x00679dd6, 0x007d9be2, 0x007599ec, 0x00869ce2, 0x00a096f1, 0x00939de0, 0x00869beb, 0x00889fde, 0x00799fe9, 0x0080a2de, 0x0078a1e4, 0x0099a2d8, 0x0091a2db, 0x008ba4db, 0x0088a8d0, 0x0071a1ed, 0x0078a5e0, 0x00969fed, 0x00b2a1e9, 0x00a3a2ef,
0x00b4a9d8, 0x00a1afd0, 0x0088b1dc, 0x00a1b6c6, 0x0091aaf2, 0x0082acf0, 0x008cafea, 0x0090b8cd, 0x009cb3da, 0x0079b1eb, 0x0096aeee, 0x0081b1eb, 0x0083b4e1, 0x0097b1e6, 0x0087b1eb, 0x00a4adf1, 0x00bdb5d1, 0x00a3afed, 0x00a0b2e5, 0x0078b6e7, 0x00c3bac7, 0x0096b6e1, 0x008eb6fb, 0x00b3b3ff,
0x00a0b5ff, 0x008fb9f7, 0x0082bef2, 0x00a2bafb, 0x00a1bcf7, 0x00c8c1dd, 0x00b9b9fd, 0x0090c0f3, 0x00a4c2f0, 0x009bc2f6, 0x00aec1f7, 0x00a0cade, 0x0092c7ef, 0x009cc7ed, 0x00a7c7ec, 0x00aecbde, 0x00bac3f8, 0x008ac9f3, 0x00b2c7f2, 0x00d0cddc, 0x00d1c4fa, 0x00bacceb, 0x00bed0e0, 0x00b4cdf0,
0x00afcdfc, 0x00b6d7de, 0x009dd2f8, 0x00bdd8dc, 0x00a7d4f7, 0x00b2d5f8, 0x00a0d8f4, 0x00b2d9ee, 0x00bcd3ff, 0x00ddd3f9, 0x00aadaf3, 0x00c0dbfb, 0x00c1def4, 0x00a9e0fd, 0x00cadefd, 0x00b1e1fd, 0x00d3dffd, 0x00bce2fd, 0x00e8dffd, 0x00ade6fa, 0x00b6e8fc, 0x00bee9f9, 0x00cae7fe, 0x00eaebef,
0x00d2eaff, 0x00ccecfb, 0x00b9eefa, 0x00edebfb, 0x00d5effc, 0x00c3f1fe, 0x00cbf2fe, 0x00e0f2ff, 0x00cbf9fb, 0x00d5f8fd, 0x00e4fafd, 0x00ebfafc, 0x00defbfd, 0x00fbfaff, 0x00ffffff
);
const int longs_per_line = 24;
const int[] bitmap = int[] (
0xf3efeffe, 0xf3f3fdfe, 0xfefefefd, 0xfefefefe, 0xfefefefe, 0x9580bcf3, 0x5d8db1b1, 0x5d595d5d, 0x1412195d, 0x8e808036, 0x6c4d467a, 0x77b88eba, 0x612a7e12, 0x5db6b1b1, 0x5d5da1b0, 0x7a445d5d, 0x00005e84, 0x7e000000, 0x4d4d4d4b, 0x00004441, 0xb54d4d6a, 0x00000077, 0x00000000, 0x00000000,
0xf3efeffd, 0xfdf3eff3, 0xfefefefe, 0xfefefefe, 0xfefefefe, 0x9e5bbcd3, 0x5d495d95, 0x5d59595d, 0x44308081, 0x4d475e8e, 0xd76c4646, 0x36a87eba, 0x8b172a36, 0xa8c9b3b1, 0x5d4873d5, 0xa532495d, 0x9fa0a0c4, 0x00007e9f, 0x4d4d664b, 0x00008766, 0xb9644d73, 0x000000a1, 0x00000000, 0x00000000,
0xf3efefef, 0xfefefdfd, 0xfefefefe, 0xfefefefe, 0xeffefefe, 0x705d80bc, 0x59495d5d, 0x565d595d, 0x41447e80, 0x4d4d4d41, 0xeed56c46, 0x118ea88e, 0xa4221a24, 0xeae1c4b1, 0x5d73d4ea, 0xce36485d, 0xbababadd, 0x00a88e8e, 0x697f6e00, 0x7a8ea572, 0xb764646e, 0x000000ce, 0x00000000, 0x36000000,
0xfefdfdfe, 0xfefefefe, 0xfefefefe, 0xfefefefe, 0xbcf3fefe, 0x595944bc, 0x32434359, 0x141e2121, 0x3f27251e, 0x464d463d, 0xeef0e37a, 0x0c2ca87e, 0xb1431717, 0xeae1d4b3, 0x5dc6f3ea, 0xce77435d, 0x0000b0dd, 0x00b80000, 0xa5b58700, 0x73d5c4a7, 0xb964646c, 0x00009bba, 0x00000000, 0x47000000,
0xfefefefe, 0xfefefefe, 0xfefefefe, 0xfefefefe, 0xbcbcfdfe, 0x435d5d80, 0x1c1c1e21, 0x1c1c1c1c, 0x1c1c1c1c, 0x87271c1c, 0xeef0f0e3, 0x1712449f, 0xb1751717, 0xe1e1d0b3, 0x5d6ad0ea, 0xc484515d, 0x0000badd, 0x00000000, 0xe3ba9900, 0x9be4e0d7, 0xca6b6464, 0x000000ce, 0x00000000, 0x50310000,
0xfefefefe, 0xfefefefe, 0xfefefefe, 0xfdfefefe, 0xbcbcd3fd, 0x251e495e, 0x25251c23, 0x41413f28, 0x3d3e3e3e, 0x315b7334, 0xeef0d773, 0x17151251, 0xb18c1717, 0x9f99b6b1, 0x5d5d487a, 0x9a733148, 0x0000cedd, 0x00000000, 0xee857200, 0xa1eeeeee, 0xd86b6b69, 0x000000a0, 0x00000000, 0x505b0000,
0xfefefefe, 0xfefefefe, 0xfefefefe, 0xfdfefefe, 0x80bcbcd3, 0x28232543, 0x7f886a48, 0x4d696969, 0x6e694d69, 0xd7e6e6ba, 0xe49f5d9f, 0x17171551, 0xb18d1717, 0x48489eb1, 0x5d5d4849, 0x885b1e43, 0x0000cedd, 0x00000000, 0xe3694100, 0xc8eee9e9, 0xd77f796b, 0x0000008f, 0x00000000, 0x6b586300,
0xfefefefe, 0xfefefefe, 0xfefefefe, 0xd3fdfefe, 0x36bcbcbc, 0xc88a3125, 0x9aacc4c9, 0xc4a57f9a, 0xddb099a7, 0xeee6e3e3, 0x879fe6f0, 0x17171a56, 0xb18d1a17, 0x3c4870b1, 0x48584545, 0x6958141f, 0x0000ced8, 0x00000000, 0xd27f4655, 0xe0e3ddd8, 0xddb29191, 0x00000016, 0x00000000, 0x6b6b5b00,
0xfefefefe, 0xfefefefe, 0xfefefefe, 0xbcd3f3fd, 0x3356b8bc, 0xe9e9d78e, 0x99c8dde9, 0xbaddd8c8, 0xe3d7ceba, 0xeee4e3e3, 0xd5f0f0f1, 0x17172451, 0xb1a41a17, 0x3c48599e, 0x48454545, 0x4c55141e, 0x0000c89b, 0x00000000, 0xceb54d46, 0xe9d8c8c8, 0xdfe79196, 0x00000000, 0x00000000, 0x58796c56,
0xfefefef3, 0xfefefefe, 0xfdfefefe, 0xbcb8d3e1, 0xc87356b8, 0xe9e3dde9, 0xddcedde9, 0xa09bbadd, 0xd2ceceb5, 0xe3d7d2d2, 0xf0f0f0ee, 0x171724cf, 0xb1a42217, 0x3c434889, 0x36494548, 0x4c4e321e, 0x0000ba7f, 0x17000000, 0xb0d76c46, 0xd8c4c9b2, 0xcef1ae9c, 0x00000000, 0x34000000, 0x5b6c866c,
0xfefefdd3, 0xfefefefe, 0xd3effdfd, 0xbcb8c5c5, 0x485a4353, 0x73484545, 0xcee9d2a0, 0x6a8f9292, 0xd2b9977f, 0xa5a5c8ca, 0xf0d79fa0, 0x171722e4, 0xb1a42217, 0x433c4859, 0x1e484848, 0x4c4c431e, 0x0000a06c, 0x6e7a0000, 0xbad5c846, 0xddc7b6ae, 0x2de4e79c, 0x00000000, 0x7a000000, 0x605b8383,
0xfdefd3d3, 0xf3f3fdfe, 0xc5c5c5d3, 0x77a8b8c5, 0x28272525, 0x23232323, 0x6a6a3f23, 0x69696c6c, 0xa0a59369, 0xd7d2c8ba, 0xd5d7dde0, 0x171751e6, 0x9ea42217, 0x43484845, 0x1e434848, 0x4c4c3236, 0x0000a055, 0x464b0000, 0xd5cff585, 0xdcc1adc1, 0x00cff6cd, 0x00000000, 0x88920000, 0x00605d88,
0xd3c5a877, 0xc5c5d3ef, 0xc5c5c5c5, 0x23317eb8, 0x2727231c, 0x23232327, 0x34272723, 0x64674c3f, 0x9b8f8869, 0xe3e3ceb5, 0xf0eee6e6, 0x171a81f0, 0x8c902c17, 0x43484845, 0x252c4849, 0x4c4c3c36, 0x00009955, 0x55413100, 0xa0cff6d7, 0xc9b6b3b6, 0x00d5e4e4, 0x00000000, 0x958f7700, 0x00000073,
0x80364731, 0xc5c5b8a9, 0x7eb8b8c5, 0x2323285e, 0x33332323, 0x33483433, 0x2927272e, 0x503e3427, 0x9ba1a06a, 0xd7baa1a1, 0xeee4e3e3, 0x1717b0f1, 0x59702c17, 0x3c484848, 0x1e1e4348, 0x3e3e3756, 0x00008755, 0x9b4d4725, 0x7f87eef1, 0xaeb3adb6, 0x0000cfe4, 0x00000000, 0x8fae97b0, 0x00000077,
0x1a324747, 0x2a1a1a1a, 0x1e1a1a2a, 0x45232525, 0xdfdfba8e, 0x6cceae95, 0x2e3e4c4c, 0x473e3d3e, 0x9b9b9b9b, 0xb5a1a1a1, 0xe4e4e3d7, 0x1726d5e6, 0x45593617, 0x3c484848, 0x1e1e3248, 0x3e3f3f36, 0x31007a55, 0xce503d3f, 0x7f69b0f5, 0xaeadadad, 0x0000d6df, 0xb0000000, 0x80b0baba, 0x00000000,
0x1748314b, 0x1a1a1a1a, 0x28241a1a, 0xdf9f3125, 0x9ee7edf5, 0x6f70cd90, 0x67676767, 0x504f3f4c, 0x9b9b9b6c, 0xa0a1a1a1, 0xe3e3d7b5, 0x172ae3e4, 0x48432217, 0x3c484848, 0x32211f48, 0x3f343436, 0x41005347, 0xdd6c3d3e, 0x696969c8, 0xaaadadad, 0x000000cf, 0xd5a90000, 0x0000b0d7, 0x00000000,
0x32434447, 0x171a1a1a, 0x25311717, 0xe9e9cb1e, 0x7996dde9, 0x6f6f7097, 0x4c67676f, 0x89875b34, 0x9ba06c6c, 0xa1a1a1a1, 0xddd2b5a1, 0x1756e3e3, 0x483c2217, 0x3c484848, 0x1e211f43, 0x3434341e, 0x41335345, 0xb5ae3e3e, 0x6969696c, 0xb0adad9c, 0x000000ab, 0xd9dfcfb0, 0x000000a2, 0x00000000,
0x32324743, 0x17171a1a, 0x32331f17, 0xc9dde436, 0x4f4f95c8, 0x4c6c7c4f, 0x344c4e4c, 0x87878748, 0x9b7f7a8f, 0x929292a0, 0xbaa0a1a0, 0x178eddd7, 0x483c2217, 0x3c484848, 0x32211e32, 0x34343432, 0x3e341e48, 0x5595413e, 0x69695050, 0xa99d9c96, 0xcf9f0000, 0x4bcfedd5, 0x00000000, 0x00000000,
0x433c4324, 0x1a171a1a, 0x77ba3117, 0xb6aee277, 0x28283488, 0x3e343f4f, 0x3c343e4c, 0x8787878f, 0x7a899b8f, 0xa2a2a28a, 0x9b929292, 0x17c8c8a1, 0x48432217, 0x43484848, 0x32281e24, 0x34343432, 0x3e3f3432, 0x5055583e, 0x6c695050, 0x818f9c96, 0xf1d9a900, 0x0000baee, 0x00000000, 0x00000000,
0x36321a8e, 0x1a1a1a24, 0x56eda11f, 0x6c96b68e, 0x1e253434, 0x3f343333, 0x6a1f3e3f, 0x8f878787, 0xb58e8792, 0xd2ceceb5, 0x9fc8cece, 0x36a0a07b, 0x48481a17, 0x48484848, 0x3233281e, 0x34343432, 0x3f3f3f28, 0x55553e34, 0x6c555555, 0xb0779588, 0xe4f9edcf, 0x00009fba, 0x00000000, 0x00000000,
0x242456ed, 0x171a1a32, 0x56ece96a, 0x3455908f, 0x331e342e, 0x34343333, 0x87222134, 0x7a928787, 0xc8c8c892, 0xeee3cec8, 0xdde3e6ee, 0x448ebad2, 0x4348151a, 0x48484949, 0x32373732, 0x37373733, 0x3f3f3f28, 0x55554737, 0x6c555555, 0xdfbd8588, 0xb0e4f6f6, 0x00000000, 0x00000000, 0x00000000,
0x1a36dff8, 0x171a2436, 0x77ecedcb, 0x34345890, 0x33331e34, 0x34343333, 0x872a1f34, 0x9f7b8787, 0xc8c8c8c8, 0xf1eecec8, 0xeef0f0f0, 0xced7e0e6, 0x4843518e, 0x48484948, 0x1e37373c, 0x37373733, 0x45453f21, 0x55554f37, 0x6c555555, 0xf5e9a16c, 0x00a9e4f6, 0x00000000, 0x00000000, 0x00000000,
0x1aa9f8f6, 0x511a2417, 0x8adce7e9, 0x37343479, 0x33332537, 0x21141421, 0x7a241521, 0xbdbd777a, 0xc8c8cac8, 0xf1eecec8, 0xf1f1f0f0, 0xe3eef1f1, 0x4392cedd, 0x48484948, 0x2137373c, 0x45373737, 0x45454537, 0x58585845, 0x58585858, 0xf1e99288, 0x00009fe4, 0x00000000, 0x00000000, 0x00000000,
0x5eedf5f5, 0x81171a1a, 0x89cccde5, 0x37373745, 0x28213133, 0x37333321, 0x6a361a1f, 0xc8bdb592, 0xc8c8c8c8, 0xeeddcece, 0xf1f1eef0, 0xf1f1f1f1, 0xc8dde3ee, 0x4549438a, 0x1e21373c, 0x45454537, 0x45454537, 0x58585845, 0x6c585858, 0xe4a0a19d, 0x0000bcd3, 0x00000000, 0x00000000, 0x00000000,
0xcee4ecec, 0xbd171a2a, 0x89c2c7c7, 0x3737373f, 0x21342833, 0x33333333, 0x8e511715, 0xc8bab5ba, 0xc8c8c8c8, 0xd7cec8ce, 0xf1f1f1f0, 0xf1f1f1f1, 0xe3f0f0f1, 0x45438ed7, 0x2145453c, 0x45453321, 0x37454545, 0x58585845, 0x6c585858, 0x7355dcb5, 0x00c5d300, 0x00000000, 0x00000000, 0x00000000,
0xcecedcdc, 0xcc1a1a81, 0x95b1b1b3, 0x3f454545, 0x45213e28, 0x21453337, 0xbaa02a24, 0xc8c8baba, 0xcec8c8c8, 0xc8cecece, 0xf0f0f1d7, 0xf4f4f4f4, 0xeef0f4f4, 0x43bdddf0, 0x3c45453c, 0x3728213c, 0x3c484545, 0x58585847, 0x6a585858, 0x4b6ec8aa, 0xa87ed47d, 0x00000000, 0x00000000, 0x00000000,
0xa1bdc7bd, 0xc72a36c7, 0x9eb1b1b1, 0x28454545, 0x08080a21, 0x140a0a08, 0xbab5731f, 0xc8c8c8bd, 0xcececec8, 0xbac8c8ce, 0xf4f0e6ce, 0xf7f7f7f7, 0xf4f4f7f7, 0xb0e9f0f4, 0x3c45453c, 0x3433483c, 0x213c453c, 0x585d3321, 0x89585858, 0x425e8a95, 0x8ea8006a, 0x000000ef, 0x00000000, 0x00000000,
0x8f9b9d8a, 0xb1368a8a, 0x70b1a4b1, 0x05144545, 0x08080505, 0x2c080808, 0xa06a5b48, 0xc8c8c8bd, 0xcececec8, 0xb999b2c8, 0xf7f0e3d2, 0xf7f7f7f7, 0xf7f7f7f7, 0xeef0f4f4, 0x4949438e, 0x3448493c, 0x4949453f, 0x45333749, 0x89585858, 0x425e5b89, 0x7aa8004b, 0x0000009f, 0x00000000, 0x00000000,
0x70718f7b, 0xa443735d, 0x499eb18c, 0x05051448, 0x08050505, 0x922c371e, 0xbab5baba, 0xcec8c8b5, 0xd2d7cece, 0xd8bab5c8, 0xf7f4f0e6, 0xfaf7f7f7, 0xf7f7f7f7, 0xe6d5d5e4, 0x49485bdd, 0x45494949, 0x49373f3d, 0x34375d5d, 0x8f704f45, 0x4d55635b, 0x7dc53641, 0x00007e6e, 0x00000000, 0x00000000,
0x6a707b56, 0x5d5b5d5d, 0x4870bba4, 0x08050833, 0x49371405, 0xb5734549, 0xc8babdbd, 0xcacecec8, 0xd2e3d7ca, 0x6c7f97ba, 0xf4d2a085, 0xfaf7f7fa, 0xf7fafafa, 0xd5e8f7fa, 0x454785a0, 0x4949453c, 0x493f3d28, 0x3f5d5d5d, 0x6a6a583e, 0xa55e445b, 0x8e004466, 0x00007d6e, 0x00000000, 0x00000000,
0x6a6a5e63, 0x78435d6a, 0x43499eb1, 0x1e060908, 0x49494832, 0xb5a04848, 0xc8c8bdbd, 0xbab5b5b5, 0xc8eee6d7, 0x41504148, 0xf79f2e3e, 0xfaf7f7fa, 0xf7fafafa, 0xf7f7fafa, 0x4d6acef4, 0x37333f41, 0x453d0933, 0x495d5d5d, 0x435d5d6c, 0xca6a315e, 0xa9005e84, 0x00008e6e, 0x00000000, 0x00000000,
0x6a5d7e00, 0xa437456a, 0x144971b1, 0x99430908, 0x48484848, 0xbab54848, 0x7f9bbaba, 0xf0e6ba93, 0xe4f4f4f0, 0x5e473492, 0xf7f08e7e, 0xfafafafa, 0xfbfbfbfb, 0xf7fafafa, 0x37baf4f7, 0x3e3d3d3d, 0x6e251427, 0x5d5d5d59, 0x77436a6c, 0xca6e2531, 0xc50087b7, 0x00007d6e, 0x00000000, 0x00000000,
0x6a5d5a00, 0x7137345d, 0x083259b1, 0x92ca5b08, 0x48484848, 0xa0a06a48, 0x857f6c87, 0xf7f4f0d7, 0xf7f7f7f7, 0xd5cfd2f4, 0xfbfbf4d5, 0xfbfbfbfb, 0xfbfbfbfb, 0xfafafbfb, 0x7be3f7f7, 0x3d3d3f37, 0x3245463d, 0x5d5d5d5d, 0x3044b26a, 0xc4852325, 0xc500a5ca, 0x00007d6e, 0x00000000, 0x00000000,
0x6a5b5100, 0x61453448, 0x06084990, 0x48d8d85b, 0x48484848, 0x588f7343, 0xe3876969, 0xf7f7f4f0, 0xf7f7f7f7, 0xfbfafbfb, 0xfbfbfafb, 0xfbfbfbfb, 0xfbfbfbfb, 0xfafafafb, 0xb0f0f7fa, 0x3d3f3c48, 0x9b936946, 0x5d595d3c, 0x319fd773, 0xb9841d23, 0xa800b9e0, 0x00007d6e, 0x00000000, 0x00000000,
0x5d440000, 0x59483734, 0x32083770, 0x48b0ddca, 0x433c3c48, 0x503e6a32, 0xf4e48769, 0xf7f7f7f4, 0xfbfbfafa, 0xfafbfbfb, 0xfbfbfbfa, 0xfbfbfbfb, 0xfbfbfbfb, 0xf7fafafb, 0xd7f0f4f7, 0x3f3c4875, 0xd89a7f69, 0x5d5d48ba, 0x259fe373, 0xa57f3e23, 0x7d00b0e0, 0x00007d72, 0x00000000, 0x00000000,
0x43770000, 0x5949373e, 0xba140a49, 0x4892dddd, 0x32434343, 0x50463e2c, 0xf4eee36e, 0xfbfaf7f4, 0xfbfbfbfb, 0xf7fafbfb, 0xfafaf7f7, 0xfbfbfbfb, 0xfbfbfbfb, 0xf7f7f7fb, 0xe3f0f0f4, 0x43434392, 0xd8c4ac7a, 0x5d5d9fdd, 0x238edd6a, 0xa5693d1d, 0x7d008ee0, 0x00007d84, 0x00000000, 0x00000000,
0x31820000, 0x5d59453d, 0xdd870a32, 0x4375ced8, 0x1f324343, 0x55503d3f, 0xf0f1eed2, 0xfbf7f7f4, 0xfbfbfbfb, 0xeeeefafb, 0xfaf7f7f4, 0xfbfbfbfa, 0xfbfbfbfb, 0xf4f7f7fb, 0xe6eef0f0, 0x433232b5, 0xc9c39948, 0x5d6ac8c4, 0x2392e981, 0xa5693d1d, 0x879f77d7, 0x00008e98, 0x00000000, 0x00000000,
0x47800000, 0x5949483e, 0xbaca3221, 0x4343baac, 0x3e323243, 0x87273d3d, 0xe3e3d2c8, 0xfbfaf7f4, 0xfbfbfbfb, 0xd7d7e6fb, 0xf7f4f4d7, 0xfbfbfaf7, 0xfbfbfbfb, 0xf4f7f7f7, 0xe3e6eef0, 0x442c5ace, 0x9aac7543, 0x5ab59a9a, 0x2392e9a9, 0x9a694d1d, 0xba9f32ce, 0x00007d99, 0x00000000, 0x00000000,
0x41530000, 0x4959483f, 0x85b58714, 0x32439f6c, 0x4a342c32, 0xb548234a, 0xe4d7d2b5, 0xfbf7e4e3, 0xfbfbfbfb, 0xd7d7d7f7, 0xf0f1eeba, 0xfbfbf7f4, 0xfafbfbfb, 0xeef0f7f7, 0xe6eeeeee, 0x3636a0d2, 0x796c3244, 0x8e976b7f, 0x2787dd77, 0xa5696b27, 0xd79f00cf, 0x0000a9ba, 0x00000000, 0x00000000,
0x3d440000, 0x485d5934, 0x4155a737, 0x36368e41, 0x3e4a2536, 0xa0921e14, 0xd7d7d2ba, 0xf4e6e4d7, 0xfafbfbf7, 0xd7d7d7e6, 0xf0f1e4c8, 0xfbfbf7f4, 0xf7fbfbfb, 0xeef0f0f7, 0xe3e3e6e3, 0x2c51a0dd, 0x671e1e32, 0x73796767, 0x2773c86a, 0xb26b674a, 0xe6cfa8cb, 0x00007e9f, 0x00000000, 0x00000000,
0x3d474400, 0x3c5d5834, 0x3e3e886a, 0x222a733e, 0x143e3e1f, 0xb58a431e, 0xc8cacac8, 0xe4e4cfba, 0xf7faf1e4, 0xd7d7d7e4, 0xeef7d7d2, 0xfbfbf4f0, 0xf7fafbfb, 0xeef0f0f0, 0xe3e3e3e3, 0x3673a0dd, 0x34251f32, 0x9b676767, 0x276ac46a, 0xa56f6767, 0xe6e6a877, 0x0000a89f, 0x00000000, 0x00000000,
0x4a3f5600, 0x215d5d34, 0x3e3d5585, 0x24245e3e, 0x14213e34, 0xc4c49b5b, 0xe4b5a19b, 0xfaf7f7f0, 0xf4fafbfb, 0xd7d7d7dd, 0xf0f7e4d2, 0xfbfbf4f0, 0xf7f7fafb, 0xe3e6f0f0, 0xe3e3e3e3, 0x229ba0dd, 0x282e2c36, 0x886f6f67, 0x3e47a16a, 0x7a6f6f6f, 0xd5e6bab8, 0x0000007e, 0x00000000, 0x00000100,
0x4a475100, 0x34595937, 0x4c4c4c6c, 0x2c24474e, 0x5825234c, 0x55729388, 0xeeeeba97, 0xf7f7f0ee, 0xf4fbfbfa, 0xd7d7d7e3, 0xf0f0e3d2, 0xfbfbf7f0, 0xe4f7fafa, 0xcf9f9fba, 0xe3e0e0e3, 0x51a1a0d8, 0x284e3c24, 0x6c79794f, 0x6f34ae6c, 0x6a797986, 0xa9e6b973, 0x000000a8, 0x807e7e63, 0x00808063,
0x67583200, 0x3f45483c, 0x4c67674e, 0x34223f4e, 0x6958212e, 0x854c4e55, 0xe3e3e3ce, 0xf7f0eee6, 0xf4fafbf7, 0xd7d7d7e3, 0xf4f4eed7, 0xfbfbf7f4, 0xe1d6f4fa, 0xfdfdfdf3, 0xd7d5f4fb, 0x92a1a0d8, 0x2e673c22, 0x7986863f, 0x913f9b70, 0x6a838686, 0x8edda57f, 0x9f5e8e00, 0x7e7e8e9f, 0x0000abab,
0x674e3c00, 0x4e495b45, 0x676f673f, 0x3f2c344c, 0x4c674f28, 0xc86c673e, 0xe3e3e3e3, 0xf7f0eee3, 0xf7fafaf7, 0xd7d7e3f1, 0xf0f0f4e4, 0xfbfbf7f4, 0xf3fdfbfb, 0xd0a7b7d4, 0xb5e1fef3, 0xa1a1a0d2, 0x346f4536, 0x7c86862e, 0x866f7391, 0x85888686, 0x63ba9a91, 0x9f858480, 0x0000009f, 0x00000000,
0x676c5b00, 0x4e49515d, 0x6f6f6f3f, 0x28212867, 0x3e676734, 0xe0a0674c, 0x8eb0dde3, 0xe6d7b0a9, 0xf7fbfaf7, 0xd7dde4f7, 0xf0f0f0e3, 0xfbfbf7f4, 0xd0d4fdfb, 0xa5c6e0e6, 0x7efde185, 0xa0a0b5a0, 0x346f4f5a, 0x86866c28, 0x86915896, 0x91898686, 0x6a7a9a91, 0x80999385, 0x00000000, 0x00000000,
0x6f9b8100, 0x675d444f, 0x6f6f672e, 0x28342167, 0x674c6767, 0xe0e07367, 0xf3f3d6dd, 0xfdfdf3f3, 0xf7fbfdfd, 0xe3e3f4f7, 0xf4f0eeee, 0xfbfbf7f7, 0xf4b7d4fb, 0xbfbfd0f3, 0xbcf37a72, 0xa1925a36, 0x3f797092, 0x96967034, 0x96967995, 0x91968696, 0x866a8596, 0x00b8a89b, 0x00000000, 0x00000000,
0x6fb56000, 0x6748516c, 0x6f6f4e34, 0x3f2e2e67, 0x6f4e676f, 0xd7e3a067, 0xc6eafefa, 0xc6a6b7bf, 0xf7fbfdf3, 0xe4f0f7f7, 0xf4f0f1ee, 0xfbfbfaf7, 0xd4bf84f3, 0x663d5c98, 0xb8a85066, 0xa0a0871f, 0x3e8379a0, 0x9696703e, 0x9196969e, 0x96a38696, 0x97868696, 0x0000aba9, 0x00000000, 0x00000000,
0x79a98100, 0x6f483679, 0x7c7c3f3f, 0x4e342167, 0x676f676f, 0xe1d7b549, 0xe8a7e1fd, 0x98d0d0e6, 0xf7fbfd99, 0xf0f4f7f7, 0xf4f0f1f1, 0xfafaf7f7, 0x66726ef7, 0x210b233d, 0xaba89966, 0xa0a0a02c, 0x3e838392, 0x9c9c703e, 0x7996969d, 0x96a396a3, 0x9fae9686, 0x000000d3, 0x00000000, 0x00000000,
0x9ba98100, 0x7c5a5a79, 0x7c7c344f, 0x6f34286c, 0x6f7c6f7c, 0xef7ba15d, 0xeae084ea, 0xa7a7b7e1, 0xf7fda866, 0xf0f4f7f7, 0xf4f0f0f0, 0xf7f7f7f4, 0x235c47d5, 0x1e0b0b0b, 0x81fdfde1, 0xa0a0a05a, 0x4e969673, 0xb6b38955, 0x8691a39d, 0x96a3a3a3, 0x7eb0b596, 0x00000000, 0x00000000, 0x00000000,
0xae600000, 0x7c5b7186, 0x86552e5d, 0x7c4e2870, 0x527c7c7c, 0xef5a6a8d, 0xa6b766c5, 0x50253d66, 0xf0e8735c, 0xf0f4f7f7, 0xf0f0f0f0, 0xf4f0f4f0, 0x142522d5, 0x0c020303, 0x2c7ea8a8, 0x92a0a17b, 0x55969675, 0xc1c7896c, 0x9679969e, 0x9ea3a3a3, 0x00b8cfcd, 0x00000000, 0x00000000, 0x00000000,
0xa0180000, 0x86777183, 0x86342e89, 0x8667334e, 0x577c7c86, 0xd6433c8c, 0x29664d7d, 0x0d0a2127, 0xeed69f6d, 0xf0f0f4f4, 0xf0f0f0f0, 0xeef4f4f0, 0x010c36cb, 0x0c010100, 0x2f190c1f, 0x92929b5a, 0x5596968c, 0xaec9ae88, 0xad7f969c, 0xc7ada39c, 0x000001cb, 0x00000000, 0x00000000, 0x00000000,
0xa9000000, 0x86758996, 0x793434a1, 0x8679344e, 0x49868686, 0xab3c5ab1, 0x0a463f51, 0x53040a0a, 0xe6f3fdf3, 0xf0f0f4f1, 0xf0f0eef0, 0xe3eef0ee, 0x010f43bd, 0x0c010101, 0x77613901, 0x8a929b92, 0x6c959e90, 0xaee7c888, 0xad9a7f9c, 0xcfb6b3ad, 0x0000008e, 0x00000000, 0x00000000, 0x00000000,
0x8e000000, 0x9e7b9cb6, 0x4f345dae, 0x86833458, 0x61798686, 0x773c71b1, 0x02040253, 0x56020202, 0xeea9a8b8, 0xddeee6ee, 0xcef4e6ce, 0xb5bde3e4, 0x39437ab5, 0x0e0f0f19, 0xa19d782f, 0x8a8f89a1, 0x7f959e9d, 0xadecdd88, 0xc1ad7f9c, 0x80ceb3b3, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
0x18000000, 0xb692b3cb, 0x2e34a19e, 0x7c793745, 0x5d70797c, 0x39738db1, 0x01010136, 0x070c0101, 0xdd9f1e22, 0xb5d2ced7, 0xc8cef1ce, 0xb5b5b5c8, 0x9db6b5b5, 0x2f3c7175, 0x9bb3a45a, 0x8f8f9da1, 0x889cb3b3, 0xaddce995, 0xc1b39a93, 0x00cfceb6, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0xcdbdcda2, 0x3445bdc2, 0x7c794545, 0x5d707c7c, 0x3992a4b1, 0x0101012f, 0x0c0c0101, 0xc8b03c1f, 0xa1bdcebd, 0xc7b5d2e3, 0xb5babdb5, 0xb3b39d9b, 0x5a8aa48d, 0x8fb3b37b, 0x9d8f9d92, 0xad9cb6b6, 0x9ab5f5ae, 0xc7c1c193, 0x007adfc7, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0xcedcce80, 0x346acde5, 0x79795d5d, 0x5d707c79, 0x8a8a8da4, 0x0e192f5a, 0x19120f0e, 0xb5a0923c, 0xa09daeb5, 0xa49dbae4, 0x9db5b5b5, 0xb3b1a492, 0xa192a48d, 0x8fb3b38f, 0xaa8fb392, 0xb697c1c7, 0x9aaeedda, 0xc7c7c79a, 0x0000cfdd, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0xdfecab80, 0x33b0e7f2, 0x79798f48, 0x5d6c7979, 0xbd8aa4a4, 0x437592ba, 0x775b3c2f, 0xb5b5a092, 0xc88db39b, 0xa1b19dba, 0x95a1b5b5, 0xb1b1b19b, 0xa192b18d, 0x8fb3b38f, 0xcd8fb592, 0xb3b3c9dc, 0x97b6dde7, 0xdcc7c7c1, 0x000000cf, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0xdfcf8000, 0x5de9f5f5, 0x7979c85d, 0x6a857979, 0xd27ba49e, 0x9fc8d7e3, 0x8e8a8e5b, 0x92b5b5bd, 0xc8899d9d, 0xb59595b5, 0x928da1a1, 0x9eb1b18f, 0xa192a471, 0x8a9e9e8f, 0xdda2dc92, 0xc7b6dce7, 0xacaccdf1, 0xe4dcc9c9, 0x000000d9, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0xedb00000, 0xa2f5f5f5, 0x8388dd45, 0x70928383, 0xba8da49e, 0xe4f0f0f0, 0xdddde3e4, 0x8fa1b5d2, 0xb59b8d9e, 0xa092959b, 0x9d9592a1, 0x959e9e8d, 0xa192b38f, 0x8f9d9d9d, 0xedcde29b, 0xcdc7dce9, 0xda9bc7f5, 0xdbeddcdc, 0x000000ab, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0xce770000, 0xe9f5f5f5, 0x95b5d75d, 0x70aa9e95, 0xa09da49e, 0xf0f0f0ee, 0xe3eeeef0, 0x958db5d7, 0xa09b8d9e, 0xa0a08995, 0x929d8da1, 0x95a49e9d, 0xa192bdae, 0x8aaeaeb6, 0xf6dfe4a2, 0xe7c7dce9, 0xedcdaee7, 0xbcf1eeed, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0xa2000000, 0xf5f5f5ed, 0xc9e2c8b0, 0x89cdccc7, 0x8fa4a49e, 0xf0f0f0e0, 0xf1eeeef0, 0x9595a0d2, 0x929b9595, 0xa0ba9b95, 0xa18fb692, 0x97a4b6c7, 0x92a1cdc7, 0x8b9dcccd, 0xf6e4f5ce, 0xf6c7dae4, 0xf4edb5c8, 0x00f1f9f9, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0x81000000, 0xf8f5f5df, 0xe2eda1ed, 0x9de2e7e7, 0x8da49590, 0xf0f0f0d2, 0xf1f1eef0, 0x959592c8, 0xbda19595, 0xa1c8c8cd, 0xcd9bcdaa, 0xb595c7cc, 0xaa8abddc, 0x8abdcdcd, 0xf6e4f6e4, 0xeddccced, 0xfcfcddae, 0x00bcf0fc, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0x00000000, 0xf8f5eda9, 0xece9dff8, 0xaee2f2f2, 0x9da4b39d, 0xf0f0eeb0, 0xe0eeeef0, 0x959e8db5, 0xdfa19d95, 0x9bbad2dc, 0xe5bdcee2, 0xbd71bdda, 0xd77592cd, 0xa1cde7e7, 0xf6f1edf6, 0xdcedc7ed, 0xfcfcf9c8, 0x0000d6f9, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0x00000000, 0xf8f6ab81, 0xebe7f8f8, 0xaae2e5e7, 0xc7ccc7b6, 0xceb08a8a, 0xbabac8d7, 0xae9d95a1, 0xecbd9dae, 0xaaa0d7ce, 0xece2bdec, 0xcea2d7e7, 0xf5aaa2df, 0xcebdeded, 0xf9f6dff9, 0xc9f1dae2, 0xfcfcfcf1, 0x0000cfe4, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0x00000000, 0xf6e48e00, 0xdcdff2f8, 0x8fcdcccd, 0xdedccdde, 0xa0cebabd, 0x737b7b8a, 0x8f8f9575, 0xedd7aaaa, 0xdd92ced7, 0xececbde2, 0xf2cecdec, 0xf6d78af2, 0xf1b0f6f6, 0xf9fcdff9, 0xe4dde7da, 0xf4fcfcfc, 0x0000cfce, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0x00000000, 0xf8ab3000, 0xcce2e5f6, 0x89b4adb3, 0xe7e7e2cd, 0xf1e3a0dc, 0xbad5e6f0, 0xcdc7c792, 0xecdfaae5, 0xeca1bddd, 0xe7ecdfcd, 0xf5ecaeec, 0xf6f592dd, 0xf9cbf1f6, 0xf6fce4f9, 0xfcddf1cd, 0xe4fcfcfc, 0x000000ce, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0x00000000, 0xe4510000, 0xdcaaccf5, 0x719e909e, 0xecece98f, 0xe6c8bdec, 0xd7f7f7f0, 0xece2e2cd, 0xe9e9aaec, 0xecce8add, 0xecededcb, 0xf5edcde2, 0xf9f9d7bd, 0xfcf1edf9, 0xedfcf9f1, 0xfcf9e9ed, 0xe4f4fafa, 0x000000a9, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0x00000000, 0x81510000, 0xcd83c7f6, 0x71909090, 0xedecbd8d, 0xd7a1e7f2, 0xd7f7f7f4, 0xebe7e7cd, 0xe9edaae2, 0xe7e77bba, 0xececede2, 0xedecdeda, 0xf9f9f1a2, 0xfcf9e4f9, 0xecfcfcee, 0xfafaf4f1, 0xe4f1fcfa, 0x000000a9, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0x00000000, 0x32560000, 0x9b83e7cb, 0x719090ae, 0xf2e98a5d, 0xbaa1ecf2, 0xcef7f7ee, 0xdedce7dc, 0xdff2aace, 0xe2e7a192, 0xece2edec, 0xdcdedab4, 0xf9f9f9cd, 0xf9fcedf9, 0xf6f1faf1, 0xfafafaee, 0xe3edf1fa, 0x00009fb0, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0x00000000, 0x36530000, 0x909ddf48, 0x909090c7, 0xecaa4971, 0x92cdebec, 0xcef7f7ce, 0xc7ccdae5, 0xcbe5aac7, 0xcddabd7b, 0xcddeebeb, 0xcddaccb3, 0xf9f9f9ec, 0xf1faf9f1, 0xf4f9faf9, 0xfafafcfa, 0xb2aee4f7, 0x00008ed2, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0x00000000, 0x36630000, 0x90cd7333, 0x9e9e9e90, 0xce6a718c, 0x9ddedee5, 0xcff7e6a0, 0xb6c2c7da, 0xa2cc8db3, 0xcccccc6a, 0xb6e5deda, 0xbdccd1b3, 0xf9f9f9e7, 0xf4fafcee, 0xf4f9fcfa, 0xfafafafa, 0x9691dcee, 0x000087c8, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0x00000000, 0x22000000, 0x9db03733, 0x9e9e9e9e, 0x8b718c89, 0xb3dedede, 0xb0f1ba9d, 0xb1b3b3cc, 0x92b3a4b1, 0xc7b3d175, 0xa4deccd1, 0xdab4d1c2, 0xf5f9f6eb, 0xf7f4fcf1, 0xfaf4fcfa, 0xf1fbfafa, 0x7995b5e4, 0x000087a5, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0x00000000, 0x362a0000, 0xcd492821, 0xb1a4a49e, 0x71a4599e, 0xccdada9d, 0xbac8a0a4, 0xbbb1b1bb, 0x71bba4b1, 0xc2b3c28b, 0xb3b6ccd1, 0xccc2d1d1, 0xe2f5f2ec, 0xfae4fcf5, 0xfafaf4fa, 0xeef4fafa, 0x649ad8dd, 0x007e7a7f, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0x00000000, 0x44220000, 0x926a2121, 0xc2b1b3cd, 0xa47071b1, 0xd1d1b471, 0xa1a08fc2, 0xb1b1b1bb, 0x71bba4b1, 0xc2a4c2a4, 0xd1a4d1d1, 0xccd1c2d1, 0xdee5ece7, 0xfaf9edf5, 0xfafafaf4, 0xe9e4f7fa, 0x69b7dde9, 0x007d9966, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0x00000000, 0x5b560000, 0xae49211e, 0xd1ccdacd, 0xa45da4d1, 0xd1c278b1, 0x8a92b1d1, 0xb1b1b1b1, 0x71bbb1b1, 0xc2b1c2b3, 0xc2c2c2c2, 0xdad1b3d1, 0xdec7dee5, 0xf1fceded, 0xfbfafaf9, 0xdfdde4fb, 0x69b2e9e9, 0x007aa764, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0x00000000, 0x48240000, 0xa1373314, 0xdae5cdcd, 0x8d6ac2da, 0xc28cb3a4, 0x8aa4c2d1, 0xbbb1b1a4, 0x8cb1bbbb, 0xd1c2b3b3, 0xd1d1c2c2, 0xd1d1d1d1, 0xe5c7b6da, 0xf1f1f9e4, 0xfbfafafc, 0xe9cec8e4, 0x93e3eeee, 0x007d985c, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0x00000000, 0x43240000, 0xae5d1414, 0xe7cecdcd, 0x719de5eb, 0x8cb1b1a4, 0x8bb1b1b3, 0xb1b1b1a4, 0x8db1b1b1, 0xd1c7b3b1, 0xd1d1b6d1, 0xd1d1d1c2, 0xdfccb6c7, 0xf9dcedf6, 0xeefafafc, 0xe9d7cdbd, 0xb7eef1ee, 0x00877242, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0x00000000, 0x33360000, 0xa1371414, 0xdfc8cdcd, 0x89e2f2f2, 0xb1a4a4a4, 0xa4b1b18c, 0xb1b1b18d, 0xa4b1b1b1, 0xd1d1b3b1, 0xd1d1ccd1, 0xd1d1d1d1, 0xeddaccb3, 0xfcf4c8e4, 0xaaeefafc, 0xe9d7bdb5, 0x84e3f1ee, 0x00a56666, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0x00000000, 0x4b360000, 0xa1330a0a, 0xcddcdcdc, 0xcef5f5e9, 0xb1a4b170, 0xa4a48db1, 0xb3a4a4a4, 0xb6b6b6b3, 0xd1c7b3b3, 0xccccccc7, 0xb6ccd1d1, 0xe2dadac7, 0xf9fcf4bd, 0xaaaae4fc, 0xf5ddbdaa, 0x66e3f1f1, 0x00994d42, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0x00000000, 0x5e360000, 0xa121040a, 0xdadcdcda, 0xf5f5e9a1, 0xa4a470cd, 0x9e9ea49e, 0xb39e9e9e, 0xc7dadac7, 0xc7b3b3a4, 0xccccccc7, 0xb6b3cccc, 0xcdb4ccc7, 0xf6f6fcf8, 0xaaaaaad7, 0xeee3cebd, 0xb7e6f1f1, 0x00a64d42, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0x00000000, 0x73510000, 0x31120408, 0xcddcdcbd, 0xf5e9a2a1, 0xa471bdf5, 0x9e9e9e8d, 0xc79e9e9e, 0xdce5deda, 0xb6b3b3b6, 0xccccc7c7, 0x9d9db6cc, 0xf6cdb4b6, 0xcde9f5f9, 0xaaaaaaaa, 0xeeddbdbd, 0xbaf1f1f1, 0x00a54d4d, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0x00000000, 0x73240000, 0x320a0414, 0xcddcdcc8, 0xdf928a92, 0x5db0f6f6, 0xa4a48870, 0xdaa49e9e, 0xececece5, 0xadb6c7e7, 0xc7c7c7b6, 0xb69d979d, 0xf2ebdab6, 0xaea1dfdf, 0xaaaaaaaa, 0xf1ddc8bd, 0xcaf1f1ee, 0x00a54d40, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0x00000000, 0x5e360000, 0x731e081e, 0xdcdcdcdc, 0x8a8a92aa, 0xdcf8f9ce, 0xae705975, 0xc7c7c7b3, 0xf5edf2e2, 0x95e2edf5, 0x95979795, 0xccb6b695, 0xe5dadacc, 0xaa8fb2bd, 0xbdaaaaaa, 0xeee9dfc8, 0x99eef1ee, 0x00846640, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0x00000000, 0x5e2a0000, 0x73140847, 0xdcdcddce, 0x8a92aacd, 0xf9f1b08a, 0x7173b0f5, 0xe5dac79b, 0xf5ece2ce, 0xdaceedf5, 0xb6b4b6ae, 0xd1ccc7c7, 0x97cddad1, 0x9bc8dcc7, 0xbdaaaaaa, 0xeee9dfce, 0xa7e3eef1, 0x0084724d, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0x00000000, 0x22000000, 0x3208084b, 0xdadcddb5, 0x8b9baacd, 0xa98a8a8a, 0xd7f1f6f1, 0xc89d8aa2, 0xecece7dc, 0xb6cde4f5, 0xc7c7ccda, 0xc7b6c1c7, 0xc1919dcc, 0xc4c9c9c9, 0xaab6aa97, 0xeeeee9d2, 0x5cb7e6ee, 0x00739840, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0x00000000, 0x1b000000, 0x1408095b, 0xbddcc86a, 0x8a8a8f9d, 0x738a8a8a, 0xf5dfce7b, 0xdfe4f1f6, 0xedece7df, 0xb3dff6f8, 0xb6b6b3b6, 0xa19d9db6, 0xa0857382, 0xc4c4c4c4, 0xa17a99c4, 0xe4e9ceaa, 0x4066b5d7, 0x0006995c, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0x00000000, 0x00000000, 0x08051e43, 0xbdceb55b, 0x8a8a92a2, 0x618a8a8a, 0xced8daae, 0xdfd5cfcb, 0xe4e4e4e4, 0x8e92b0ce, 0x8a928a8a, 0x00777373, 0x7e000000, 0xc4b2997a, 0x6eb2b9c4, 0xd2c8a085, 0x404072c8, 0x00008e84, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0x00000000, 0x00000000, 0x08084744, 0xbdb56a14, 0x8b8a9bbd, 0x895d8a8a, 0xb0c4b2ac, 0xa2818e77, 0xa9a9a2a2, 0x008e8e8e, 0x00000000, 0x00000000, 0x00000000, 0x7a7a7700, 0xb7b9a585, 0x6e728498, 0x72404d50, 0x00008e7d, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0x00000000, 0x00000000, 0x05145b2a, 0xb5874414, 0x8a8a92aa, 0x93586a8a, 0x567a9a9a, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x856e736e, 0x72727a85, 0x7d747272, 0x000000eb, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0x00000000, 0x00000000, 0x05333200, 0x875b0805, 0x8a8a8a8f, 0x7f7f3f73, 0x00775b7f, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0x00000000, 0x00000000, 0x095e4400, 0x5b331406, 0x737b738a, 0x6e72693f, 0x0000515b, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0x00000000, 0x00000000, 0x4b445300, 0x14090808, 0x3f1e1e31, 0x5e5b5566, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0x00000000, 0x00000000, 0x6e770000, 0x0606094b, 0x69502809, 0x0000736e, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0x00000000, 0x00000000, 0x7e000000, 0x50506e5e, 0x6a5e6e66, 0x0000004b, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x727a7d63, 0x00a85e6e, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000
);

int getPaletteIndexXY(in ivec2 fetch_pos)
{
    int palette_index = 0;
    if(fetch_pos.x >= 0 && fetch_pos.y >= 0
        && fetch_pos.x < int(bitmap_size.x) && fetch_pos.y < int(bitmap_size.y))
    {
        int line_index = fetch_pos.y * longs_per_line;

        int long_index = line_index + (fetch_pos.x >> 2);
        int bitmap_long = bitmap[long_index];

        int byte_index = fetch_pos.x & 0x03;
        palette_index = (bitmap_long >> (byte_index << 3)) & 0xff;
    }
    return palette_index;
}


int getPaletteIndex(in vec2 uv) {
    int palette_index = 0;
    ivec2 fetch_pos = ivec2(uv * bitmap_size);
    palette_index = getPaletteIndexXY(fetch_pos);
    return palette_index;
}

vec4 getColorFromPalette(in int palette_index) {
    int int_color = palette[palette_index];
    return vec4(float(int_color & 0xff)/ 255.0,
                float((int_color >> 8)& 0xff)/ 255.0,
                float((int_color >> 16)& 0xff)/ 255.0,
                0);
}

vec4 getBitmapColor(in vec2 uv) {
    return getColorFromPalette(getPaletteIndex(uv));
}

void iconImage(out vec4 fragColor, in vec2 fragCoord) {
    vec2 uv = fragCoord / bitmap_size;
    fragColor = getBitmapColor(uv);
}

// --------------------------------- //
// --------[ Text Drawing]---------- //
// --------------------------------- //

vec2 uv;
const vec2 ch_size  = vec2(1.0, 2.0) * 0.6;              // character size (Y,X)
const vec2 ch_space = ch_size + vec2(1.0, 1.0);    // character distance Vector(X,Y)
const vec2 ch_start = vec2 (ch_space.x * -1.75, 3.25); // start position
      vec2 ch_pos   = vec2 (0.0, 0.0);             // character position(X,Y)
#define REPEAT_SIGN false // True/False; True=Multiple, False=Single

#define n0 ddigit(0x22FF);
#define n1 ddigit(0x0281);
#define n2 ddigit(0x1177);
#define n3 ddigit(0x11E7);
#define n4 ddigit(0x5508);
#define n5 ddigit(0x11EE);
#define n6 ddigit(0x11FE);
#define n7 ddigit(0x2206);
#define n8 ddigit(0x11FF);
#define n9 ddigit(0x11EF);

#define A ddigit(0x119F);
#define B ddigit(0x927E);
#define C ddigit(0x007E);
#define D ddigit(0x44E7);
#define E ddigit(0x107E);
#define F ddigit(0x101E);
#define G ddigit(0x807E);
#define H ddigit(0x1199);
#define I ddigit(0x4466);
#define J ddigit(0x4436);
#define K ddigit(0x9218);
#define L ddigit(0x0078);
#define M ddigit(0x0A99);
#define N ddigit(0x8899);
#define O ddigit(0x00FF);
#define P ddigit(0x111F);
#define Q ddigit(0x80FF);
#define R ddigit(0x911F);
#define S ddigit(0x8866);
#define T ddigit(0x4406);
#define U ddigit(0x00F9);
#define V ddigit(0x2218);
#define W ddigit(0xA099);
#define X ddigit(0xAA00);
#define Y ddigit(0x4A00);
#define Z ddigit(0x2266);
#define _ ch_pos.x += ch_space.x;
#define s_dot     ddigit(0);
#define s_minus   ddigit(0x1100);
#define s_plus    ddigit(0x5500);
#define s_greater ddigit(0x2800);
#define s_less    ddigit(0x8200);
#define s_sqrt    ddigit(0x0C02);
#define nl1 ch_pos = ch_start;  ch_pos.y -= 3.0;
#define nl2 ch_pos = ch_start;  ch_pos.y -= 6.0;
#define nl3 ch_pos = ch_start;	ch_pos.y -= 9.0;
#define nl4 ch_pos = ch_start;	ch_pos.y -= 12.0;

float dseg(vec2 p0, vec2 p1)
{
	vec2 dir = normalize(p1 - p0);
	vec2 cp = (uv - ch_pos - p0) * mat2(dir.x, dir.y,-dir.y, dir.x);
	return distance(cp, clamp(cp, vec2(0), vec2(distance(p0, p1), 0)));   
}

bool bit(int n, int b)
{
	return mod(floor(float(n) / exp2(floor(float(b)))), 2.0) != 0.0;
}

float d = 1e6;

void ddigit(int n)
{
	float v = 1e6;	
	vec2 cp = uv - ch_pos;
	if (n == 0)     v = min(v, dseg(vec2(-0.405, -1.000), vec2(-0.500, -1.000)));
	if (bit(n,  0)) v = min(v, dseg(vec2( 0.500,  0.063), vec2( 0.500,  0.937)));
	if (bit(n,  1)) v = min(v, dseg(vec2( 0.438,  1.000), vec2( 0.063,  1.000)));
	if (bit(n,  2)) v = min(v, dseg(vec2(-0.063,  1.000), vec2(-0.438,  1.000)));
	if (bit(n,  3)) v = min(v, dseg(vec2(-0.500,  0.937), vec2(-0.500,  0.062)));
	if (bit(n,  4)) v = min(v, dseg(vec2(-0.500, -0.063), vec2(-0.500, -0.938)));
	if (bit(n,  5)) v = min(v, dseg(vec2(-0.438, -1.000), vec2(-0.063, -1.000)));
	if (bit(n,  6)) v = min(v, dseg(vec2( 0.063, -1.000), vec2( 0.438, -1.000)));
	if (bit(n,  7)) v = min(v, dseg(vec2( 0.500, -0.938), vec2( 0.500, -0.063)));
	if (bit(n,  8)) v = min(v, dseg(vec2( 0.063,  0.000), vec2( 0.438, -0.000)));
	if (bit(n,  9)) v = min(v, dseg(vec2( 0.063,  0.063), vec2( 0.438,  0.938)));
	if (bit(n, 10)) v = min(v, dseg(vec2( 0.000,  0.063), vec2( 0.000,  0.937)));
	if (bit(n, 11)) v = min(v, dseg(vec2(-0.063,  0.063), vec2(-0.438,  0.938)));
	if (bit(n, 12)) v = min(v, dseg(vec2(-0.438,  0.000), vec2(-0.063, -0.000)));
	if (bit(n, 13)) v = min(v, dseg(vec2(-0.063, -0.063), vec2(-0.438, -0.938)));
	if (bit(n, 14)) v = min(v, dseg(vec2( 0.000, -0.938), vec2( 0.000, -0.063)));
	if (bit(n, 15)) v = min(v, dseg(vec2( 0.063, -0.063), vec2( 0.438, -0.938)));
	ch_pos.x += ch_space.x;
	d = min(d, v);
}

vec3 hsv2rgb_smooth( in vec3 c )
{
    vec3 rgb = clamp( abs(mod(c.x*6.0+vec3(0.0,4.0,2.0),6.0)-3.0)-1.0, 0.0, 1.0 );

	rgb = rgb*rgb*(3.0-2.0*rgb); // cubic smoothing	

	return c.z * mix( vec3(1.0), rgb, c.y);
}

void mainText( out vec4 fragColor, in vec2 fragCoord ) 
{
	
	vec2 aspect = (resolution.xy / resolution.y) + 0.2;
	uv = ( fragCoord.xy / resolution.y ) - aspect / 2.0;
	float _d =  1.0-length(uv);
	uv *= 18.0 ;
	uv -= vec2(-7., 1.);

	vec3 ch_color = hsv2rgb_smooth(vec3(iTime*0.4+uv.y*0.1,0.5,0.5));
	uv.x += 0.7+sin(iTime+uv.y*0.7)*0.5;
	ch_pos = ch_start;
	
	
	nl1
	nl2
	_ _ L I G H T A M P nl3


	vec3 color = mix(ch_color, vec3(0,0,0), 1.0- (0.09 / d*2.0));  // shading
	fragColor = vec4(color, 1.0);
}

void secText( out vec4 fragColor, in vec2 fragCoord ) 
{
	
	vec2 aspect = (resolution.xy / resolution.y) - 1.3;
	uv = ( fragCoord.xy / resolution.y ) - aspect / 2.0;
	float _d =  1.0-length(uv);
	uv *= 18.0 ;
	uv -= vec2(-7., 1.);

	vec3 ch_color = hsv2rgb_smooth(vec3(iTime*0.4+uv.y*0.1,0.5,0.5));
	uv.x += -60.0 + (iTime*3);
	ch_pos = ch_start;
	
	s_dot s_dot s_dot A _ N E W _ F R E S H _ R E L E A S E _ F O R _ U s_dot s_dot s_dot
	vec3 color = mix(ch_color, vec3(0,0,0), 1.0- (0.09 / d*1.0));  // shading
	fragColor = vec4(color, 1.0);
}

void thirdText( out vec4 fragColor, in vec2 fragCoord ) 
{
	
	vec2 aspect = (resolution.xy / resolution.y) + 0.3;
	uv = ( fragCoord.xy / resolution.y ) - aspect / 2.0;
	float _d =  1.0-length(uv);
	uv *= 18.0 ;
	uv -= vec2(-7., 1.);

	vec3 ch_color = hsv2rgb_smooth(vec3(iTime*0.4+uv.y*0.1,0.5,0.5));
	uv.x += -60.0 + ((iTime-14)*3);
	ch_pos = ch_start;
	
	A T _ F I R S T _ T H A N K _ Y O U _ A L L _ F O R _ B E I N G _ A W E S O M E s_dot _ N O W _ T H X _ T O _ M Y _ C O M P U T E R _ F O R _ N O T _ C R A S H I N G s_dot _ 
	N O _ T H X _ T O _ W I N D O W S _ F O R _ B E I N G _ C R A P _ S I N C E _ n8 n5 s_dot _ _
	R U N N I N G _ O U T _ O F _ T E X T _ D O N T _ K N O W _ W H A T _ T O _ S A Y s_dot _ N A H s_dot  _ E N J O Y _ T H E _ A N I M A T I O N _ A N D _ M U S I C s_dot
	vec3 color = mix(ch_color, vec3(0,0,0), 1.0- (0.09 / d*1.0));  // shading
	fragColor = vec4(color, 1.0);
}

void main(void)
{
    vec4 text_FragColor;
    vec4 scroltext_FragColor;
    vec4 IconImage_FragColor;
    

    if (time < 7801 && time > 0)
    {
		mainText(text_FragColor, gl_FragCoord.xy);
		gl_FragColor = text_FragColor / (time/12000);
	}
	else if (time < 22500 && time > 7800)
	{
		mainText(text_FragColor, gl_FragCoord.xy);
		gl_FragColor = text_FragColor;
		
		iconImage(IconImage_FragColor, gl_FragCoord.xy);
		gl_FragColor += IconImage_FragColor;

		secText(scroltext_FragColor, gl_FragCoord.xy);
		gl_FragColor += scroltext_FragColor;
	}
	else
	{
		mainImage(gl_FragColor, gl_FragCoord.xy);
		gl_FragColor.a = 1.0;
		thirdText(scroltext_FragColor, gl_FragCoord.xy);
		gl_FragColor += scroltext_FragColor;
	}
}
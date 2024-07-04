// Not working
#version 120
#ifdef GL_ES
precision mediump float;
#endif

uniform float time;
varying vec2 surfacePosition;
#define time (time/2000) + length(surfacePosition)
uniform vec2 mouse;
uniform vec2 resolution;

#define ENABLE_HARD_SHADOWS // turn off to enable faster AO soft shadows 
#define ENABLE_VIBRATION

//Text
#define CHAR_SIZE vec2(6, 7)
#define CHAR_SPACING vec2(6, 9)

#define DOWN_SCALE 4.0

vec2 res = resolution.xy / DOWN_SCALE;
vec2 start_pos = vec2(0);
vec2 print_pos = vec2(0);
vec2 print_pos_pre_move = vec2(0);
vec3 text_color = vec3(1);
float hash(float x);

float linenum = 0.;


// Parameters for a shitty computer
#define RAY_STEPS 60
#define SHADOW_STEPS 12
#define LIGHT_COLOR vec3(.85,.9,1.)
#define AMBIENT_COLOR vec3(.3,.9,1.)
#define FLOOR_COLOR vec3(1.,.7,.9)
#define ENERGY_COLOR vec3(1.,.7,.4)
#define BRIGHTNESS .9
#define GAMMA 2.1
#define SATURATION .9
#define ENABLE_POSTPROCESS


#define detail .00005
#define t time*.1


vec3 lightdir=normalize(vec3(0.8,-0.3,-1.));
vec3 ambdir=normalize(vec3(0.,0.,1.));
const vec3 origin=vec3(0.,3.11,0.);
vec3 energy=vec3(0.01);
#ifdef ENABLE_VIBRATION
float vibration=sin(mouse.x*1000.)*.0013;
#else
float vibration=0.;
#endif
float det=0.0;
vec3 pth1;

float fnoo() {
	return sin(time) * 0.5 + 0.5;
}

mat2 rot(float a) {
	return mat2(cos(a),sin(a),-sin(a),cos(a));	
}


vec3 path(float ti) {
return vec3(sin(ti),.3-sin(ti*.632)*.3,cos(ti*.5))*.5;
}

float Sphere(vec3 p, vec3 rd, float r){//A RAY TRACED SPHERE
	float b = dot( -p, rd );
	float inner = b * b - dot( p, p ) + r * r;
	if( inner < 0.0 ) return -1.0;
	return b - sqrt( inner );
}

vec2 de(vec3 pos) {
	float hid=0.;
	vec3 tpos=pos;
	tpos.xz=abs(.5-mod(tpos.xz,1.));
	vec4 p=vec4(tpos,1.);
	float y=max(0.,.35-abs(pos.y-3.35))/.35;
	for (int i=0; i<5; i++) {//LOWERED THE ITERS
		p.xyz = abs(p.xyz)-vec3(-0.02,1.98,-0.02);
		p=p*(2.0+vibration*y)/clamp(dot(p.xyz,p.xyz),.4,1.)-vec4(0.5,1.,0.4,0.);
		p.xz*=mat2(-0.416*fnoo(),-0.91,0.91,-0.416*fnoo());
	}
	float fl=pos.y-3.013;
	float fr=(length(max(abs(p.xyz)-vec3(0.1,5.0,0.1),vec3(0.0)))-0.05)/p.w;//RETURN A RRECT
	//float fr=length(p.xyz)/p.w;
	float d=min(fl,fr);
	d=min(d,-pos.y+3.95);
	if (abs(d-fl)<.001) hid=1.;
	return vec2(d,hid);
}


vec3 normal(vec3 p) {
	vec3 e = vec3(0.0,det,0.0);
	
	return normalize(vec3(
			de(p+e.yxx).x-de(p-e.yxx).x,
			de(p+e.xyx).x-de(p-e.xyx).x,
			de(p+e.xxy).x-de(p-e.xxy).x
			)
		);	
}

float shadow(vec3 pos, vec3 sdir) {//THIS ONLY RUNS WHEN WITH HARD SHADOWS
	float sh=1.0;
	float totdist =2.0*det;
	float dist=10.;
	float t1=Sphere((pos-.005*sdir)-pth1,-sdir,0.015);
	if (t1>0. && t1<.5) {
		vec3 sphglowNorm=normalize(pos-t1*sdir-pth1);
		sh=1.-pow(max(.0,dot(sphglowNorm,sdir))*1.2,3.);
	} 
		for (int steps=0; steps<SHADOW_STEPS; steps++) {
			if (totdist<.6 && dist>detail) {
				vec3 p = pos - totdist * sdir;
				dist = de(p).x;
				sh = min( sh, max(50.*dist/totdist,0.0) );
				totdist += max(.01,dist);
			}
		}
	
    return clamp(sh,0.1,1.0);
}


float calcAO( const vec3 pos, const vec3 nor ) {
	float aodet=detail*40.;
	float totao = 0.0;
    float sca = 14.0;
    for( int aoi=0; aoi<5; aoi++ ) {
        float hr = aodet*float(aoi*aoi);
        vec3 aopos =  nor * hr + pos;
        float dd = de( aopos ).x;
        totao += -(dd-hr)*sca;
        sca *= 0.7;
    }
    return clamp( 1.0 - 5.0*totao, 0., 1.0 );
}

float texture(vec3 p) {
	p=abs(.5-fract(p*10.));
	vec3 c=vec3(3.);
	float es, l=es=0.;
	for (int i = 0; i < 10; i++) { 
			p = abs(p + c) - abs(p - c) - p; 
			p/= clamp(dot(p, p), .0, 1.);
			p = p* -1.5 + c;
			if ( mod(float(i), 2.) < 1. ) { 
				float pl = l;
				l = length(p);
				es+= exp(-1. / abs(l - pl));
			}
	}
	return es;
}

vec3 light(in vec3 p, in vec3 dir, in vec3 n, in float hid) {//PASSING IN THE NORMAL
	#ifdef ENABLE_HARD_SHADOWS
		float sh=shadow(p, lightdir);
	#else
		float sh=calcAO(p,-2.5*lightdir);//USING AO TO MAKE VERY SOFT SHADOWS
	#endif
	float ao=calcAO(p,n);
	float diff=max(0.,dot(lightdir,-n))*sh;
	float y=3.35-p.y;
	vec3 amb=max(.5,dot(dir,-n))*.5*AMBIENT_COLOR;
	if (hid<.5) {
		amb+=max(0.2,dot(vec3(0.,1.,0.),-n))*FLOOR_COLOR*pow(max(0.,.2-abs(3.-p.y))/.2,1.5)*2.;
		amb+=energy*pow(max(0.,.4-abs(y))/.4,2.)*max(0.2,dot(vec3(0.,-sign(y),0.),-n))*2.;
	}
	vec3 r = reflect(lightdir,n);
	float spec=pow(max(0.,dot(dir,-r))*sh,10.);
	vec3 col;
	float energysource=pow(max(0.,.04-abs(y))/.04,4.)*2.;
	if (hid>1.5) {col=vec3(1.); spec=spec*spec;}
	else{
		float k=texture(p)*.23+.2; 
		k=min(k,1.5-energysource);
		col=mix(vec3(k,k*k,k*k*k),vec3(k),.3);
		if (abs(hid-1.)<.001) col*=FLOOR_COLOR*1.3;
	}
	col=col*(amb+diff*LIGHT_COLOR)+spec*LIGHT_COLOR;	
	if (hid<.5) { 
		col=max(col,energy*2.*energysource);
	}
	col*=min(1.,ao+length(energy)*.5*max(0.,.1-abs(y))/.1);
	return col;
}

vec3 raymarch(in vec3 from, in vec3 dir) 

{
	float ey=mod(t*.5,1.);
	float glow,eglow,ref,sphdist,totdist=glow=eglow=ref=sphdist=0.;
	vec2 d=vec2(1.,0.);
	vec3 p, col=vec3(0.);
	vec3 origdir=dir,origfrom=from,sphNorm;
	
	//FAKING THE SQUISHY BALL BY MOVING A RAY TRACED BALL
	vec3 wob=cos(dir*500.0*length(from-pth1)+(from-pth1)*250.+time*10.)*0.0005;
	float t1=Sphere(from-pth1+wob,dir,0.015);
	float tg=Sphere(from-pth1+wob,dir,0.02);
	if(t1>0.){
		ref=1.0;from+=t1*dir;sphdist=t1;
		sphNorm=normalize(from-pth1+wob);
		dir=reflect(dir,sphNorm);
	} 
	else if (tg>0.) { 
		vec3 sphglowNorm=normalize(from+tg*dir-pth1+wob);
		glow+=pow(max(0.,dot(sphglowNorm,-dir)),5.);
	};
	
	for (int i=0; i<RAY_STEPS; i++) {
		if (d.x>det && totdist<3.0) {
			p=from+totdist*dir;
			d=de(p);
			det=detail*(1.+totdist*60.)*(1.+ref*5.);
			totdist+=d.x; 
			energy=ENERGY_COLOR*(1.5+sin(time*20.+p.z*10.))*.25;
			if(d.x<0.015)glow+=max(0.,.015-d.x)*exp(-totdist);
			if (d.y<.5 && d.x<0.03){//ONLY DOING THE GLOW WHEN IT IS CLOSE ENOUGH
				float glw=min(abs(3.35-p.y-ey),abs(3.35-p.y+ey));//2 glows at once
				eglow+=max(0.,.03-d.x)/.03*
				(pow(max(0.,.05-glw)/.05,5.)
				+pow(max(0.,.15-abs(3.35-p.y))/.15,8.))*1.5;
			}
		}
	}
	float l=pow(max(0.,dot(normalize(-dir.xz),normalize(lightdir.xz))),2.);
	l*=max(0.2,dot(-dir,lightdir));
	vec3 backg=.5*(1.2-l)+LIGHT_COLOR*l*.7;
	backg*=AMBIENT_COLOR;
	if (d.x<=det) {
		vec3 norm=normal(p-abs(d.x-det)*dir);//DO THE NORMAL CALC OUTSIDE OF LIGHTING (since we already have the sphere normal)
		col=light(p-abs(d.x-det)*dir, dir, norm, d.y)*exp(-.2*totdist*totdist); 
		col = mix(col, backg, 1.0-exp(-1.*pow(totdist,1.5)));
	} else { 
		col=backg;
	}
	vec3 lglow=LIGHT_COLOR*pow(l,30.)*.5;
	col+=glow*(backg+lglow)*1.3;
	col+=pow(eglow,2.)*energy*.015;
	col+=lglow*min(1.,totdist*totdist*.3);
	if (ref>0.5) {
		vec3 sphlight=light(origfrom+sphdist*origdir,origdir,sphNorm,2.);
		col=mix(col*.3+sphlight*.7,backg,1.0-exp(-1.*pow(sphdist,1.5)));
	}
	return col; 
}

vec3 move(inout mat2 rotview1,inout mat2 rotview2) {
	vec3 go=path(t);
	vec3 adv=path(t+.7);
	vec3 advec=normalize(adv-go);
	float an=atan(advec.x,advec.z);
	rotview1=mat2(cos(an),sin(an),-sin(an),cos(an));
		  an=advec.y*1.7;
	rotview2=mat2(cos(an),sin(an),-sin(an),cos(an));
	return go;
}

//Textsection
//Text coloring
#define HEX(i) text_color = mod(vec3(i / 65536,i / 256,i),vec3(256.0))/255.0;
#define RGB(r,g,b) text_color = vec3(r,g,b);

#define STRWIDTH(c) (c * CHAR_SPACING.x)
#define STRHEIGHT(c) (c * CHAR_SPACING.y)
#define BEGIN_TEXT(x,y) print_pos = floor(vec2(x,y)); start_pos = floor(vec2(x,y)); linenum = 0.;

//Automatically generated from the sprite sheet here: http://uzebox.org/wiki/index.php?title=File:Font6x8.png
#define _ col+=char(vec2(0.0,0.0),uv);
#define _spc col+=char(vec2(0.0,0.0),uv)*text_color;
#define _exc col+=char(vec2(276705.0,32776.0),uv)*text_color;
#define _quo col+=char(vec2(1797408.0,0.0),uv)*text_color;
#define _hsh col+=char(vec2(10738.0,21134484.0*hash(time)),uv)*text_color;
#define _dol col+=char(vec2(538883.0,19976.0),uv)*text_color;
#define _pct col+=char(vec2(1664033.0,68006.0),uv)*text_color;
#define _amp col+=char(vec2(545090.0,174362.0),uv)*text_color;
#define _apo col+=char(vec2(798848.0,0.0),uv)*text_color;
#define _lbr col+=char(vec2(270466.0,66568.0),uv)*text_color;
#define _rbr col+=char(vec2(528449.0,33296.0),uv)*text_color;
#define _ast col+=char(vec2(10471.0,1688832.0),uv)*text_color;
#define _crs col+=char(vec2(4167.0,1606144.0),uv)*text_color;
#define _per col+=char(vec2(0.0,1560.0),uv)*text_color;
#define _dsh col+=char(vec2(7.0,1572864.0),uv)*text_color;
#define _com col+=char(vec2(0.0,1544.0),uv)*text_color;
#define _lsl col+=char(vec2(1057.0,67584.0),uv)*text_color;
#define _0 col+=char(vec2(935221.0,731292.0),uv)*text_color;
#define _1 col+=char(vec2(274497.0,33308.0),uv)*text_color;
#define _2 col+=char(vec2(934929.0,1116222.0),uv)*text_color;
#define _3 col+=char(vec2(934931.0,1058972.0),uv)*text_color;
#define _4 col+=char(vec2(137380.0,1302788.0),uv)*text_color;
#define _5 col+=char(vec2(2048263.0,1058972.0),uv)*text_color;
#define _6 col+=char(vec2(401671.0,1190044.0),uv)*text_color;
#define _7 col+=char(vec2(2032673.0,66576.0),uv)*text_color;
#define _8 col+=char(vec2(935187.0,1190044.0),uv)*text_color;
#define _9 col+=char(vec2(935187.0,1581336.0),uv)*text_color;
#define _col col+=char(vec2(195.0,1560.0),uv)*text_color;
#define _scl col+=char(vec2(195.0,1544.0),uv)*text_color;
#define _les col+=char(vec2(135300.0,66052.0),uv)*text_color;
#define _equ col+=char(vec2(496.0,3968.0),uv)*text_color;
#define _grt col+=char(vec2(528416.0,541200.0),uv)*text_color;
#define _que col+=char(vec2(934929.0,1081352.0),uv)*text_color;
#define _ats col+=char(vec2(935285.0,714780.0),uv)*text_color;
#define _A col+=char(vec2(935188.0,780450.0),uv)*text_color;
#define _B col+=char(vec2(1983767.0,1190076.0),uv)*text_color;
#define _C col+=char(vec2(935172.0,133276.0),uv)*text_color;
#define _D col+=char(vec2(1983764.0,665788.0),uv)*text_color;
#define _E col+=char(vec2(2048263.0,1181758.0),uv)*text_color;
#define _F col+=char(vec2(2048263.0,1181728.0),uv)*text_color;
#define _G col+=char(vec2(935173.0,1714334.0),uv)*text_color;
#define _H col+=char(vec2(1131799.0,1714338.0),uv)*text_color;
#define _I col+=char(vec2(921665.0,33308.0),uv)*text_color;
#define _J col+=char(vec2(66576.0,665756.0),uv)*text_color;
#define _K col+=char(vec2(1132870.0,166178.0),uv)*text_color;
#define _L col+=char(vec2(1065220.0,133182.0),uv)*text_color;
#define _M col+=char(vec2(1142100.0,665762.0),uv)*text_color;
#define _N col+=char(vec2(1140052.0,1714338.0),uv)*text_color;
#define _O col+=char(vec2(935188.0,665756.0),uv)*text_color;
#define _P col+=char(vec2(1983767.0,1181728.0),uv)*text_color;
#define _Q col+=char(vec2(935188.0,698650.0),uv)*text_color;
#define _R col+=char(vec2(1983767.0,1198242.0),uv)*text_color;
#define _S col+=char(vec2(935171.0,1058972.0),uv)*text_color;
#define _T col+=char(vec2(2035777.0,33288.0),uv)*text_color;
#define _U col+=char(vec2(1131796.0,665756.0),uv)*text_color;
#define _V col+=char(vec2(1131796.0,664840.0),uv)*text_color;
#define _W col+=char(vec2(1131861.0,699028.0),uv)*text_color;
#define _X col+=char(vec2(1131681.0,84130.0),uv)*text_color;
#define _Y col+=char(vec2(1131794.0,1081864.0),uv)*text_color;
#define _Z col+=char(vec2(1968194.0,133180.0),uv)*text_color;
#define _lsb col+=char(vec2(925826.0,66588.0),uv)*text_color;
#define _rsl col+=char(vec2(16513.0,16512.0),uv)*text_color;
#define _rsb col+=char(vec2(919584.0,1065244.0),uv)*text_color;
#define _pow col+=char(vec2(272656.0,0.0),uv)*text_color;
#define _usc col+=char(vec2(0.0,62.0),uv)*text_color;
#define _a col+=char(vec2(224.0,649374.0),uv)*text_color;
#define _b col+=char(vec2(1065444.0,665788.0),uv)*text_color;
#define _c col+=char(vec2(228.0,657564.0),uv)*text_color;
#define _d col+=char(vec2(66804.0,665758.0),uv)*text_color;
#define _e col+=char(vec2(228.0,772124.0),uv)*text_color;
#define _f col+=char(vec2(401543.0,1115152.0),uv)*text_color;
#define _g col+=char(vec2(244.0,665474.0),uv)*text_color;
#define _h col+=char(vec2(1065444.0,665762.0),uv)*text_color;
#define _i col+=char(vec2(262209.0,33292.0),uv)*text_color;
#define _j col+=char(vec2(131168.0,1066252.0),uv)*text_color;
#define _k col+=char(vec2(1065253.0,199204.0),uv)*text_color;
#define _l col+=char(vec2(266305.0,33292.0),uv)*text_color;
#define _m col+=char(vec2(421.0,698530.0),uv)*text_color;
#define _n col+=char(vec2(452.0,1198372.0),uv)*text_color;
#define _o col+=char(vec2(228.0,665756.0),uv)*text_color;
#define _p col+=char(vec2(484.0,667424.0),uv)*text_color;
#define _q col+=char(vec2(244.0,665474.0),uv)*text_color;
#define _r col+=char(vec2(354.0,590904.0),uv)*text_color;
#define _s col+=char(vec2(228.0,114844.0),uv)*text_color;
#define _t col+=char(vec2(8674.0,66824.0),uv)*text_color;
#define _u col+=char(vec2(292.0,1198868.0),uv)*text_color;
#define _v col+=char(vec2(276.0,664840.0),uv)*text_color;
#define _w col+=char(vec2(276.0,700308.0),uv)*text_color;
#define _x col+=char(vec2(292.0,1149220.0),uv)*text_color;
#define _y col+=char(vec2(292.0,1163824.0),uv)*text_color;
#define _z col+=char(vec2(480.0,1148988.0),uv)*text_color;
#define _lpa col+=char(vec2(401542.0,66572.0),uv)*text_color;
#define _bar col+=char(vec2(266304.0,33288.0),uv)*text_color;
#define _rpa col+=char(vec2(788512.0,1589528.0),uv)*text_color;
#define _tid col+=char(vec2(675840.0,0.0),uv)*text_color;
#define _lar col+=char(vec2(8387.0,1147904.0),uv)*text_color;
#define _gay col+=char(vec2(133120.0,0.0),uv)*text_color;
#define _nl print_pos = start_pos - (++linenum)*vec2(0,CHAR_SPACING.y);

//Extracts bit b from the given number.
float extract_bit(float n, float b)
{
	b = clamp(b,-1.0,22.0);
	return floor(mod(floor(n / pow(2.0,floor(b))),2.0));   
}

//Returns the pixel at uv in the given bit-packed sprite.
float sprite(vec2 spr, vec2 size, vec2 uv)
{
	uv = floor(uv);
	float bit = (size.x-uv.x-1.0) + uv.y * size.x;  
	bool bounds = all(greaterThanEqual(uv,vec2(0)))&& all(lessThan(uv,size)); 
	return bounds ? extract_bit(spr.x, bit - 21.0) + extract_bit(spr.y, bit) : 0.0;
}

//Prints a character and moves the print position forward by 1 character width.
vec3 char(vec2 ch, vec2 uv)
{
	float px = sprite(ch, CHAR_SIZE, uv - print_pos);
	print_pos.x += CHAR_SPACING.x;
	return vec3(px);
}


vec3 Text(vec2 uv)
{
    	vec3 col = vec3(0.0);
    	
    	vec2 center_pos = vec2(res.x/2.0 - STRWIDTH(20.0)/2.0,res.y/1.5 - STRHEIGHT(1.0)/2.0);
       	
    	BEGIN_TEXT(center_pos.x,center_pos.y)
	print_pos += vec2(20.+cos(time*6.)*7.,sin(time*3.)*7.);
	
	RGB(sin(2.*time+uv.x*0.03)/2. +0.5, -sin(4.*time-uv.y*0.03)/2. +0.5, cos(5.*time+dot(cos(uv),sin(uv))/3.)/2. +0.5)
		_B _M _P _ _1 _per _5 _per _7 _nl _t _h _a _n _k _s _ _t _o _ _a _l _l _ _m _e _m _b _e _r _s _ _nl _ _ _ _ _o _f _ _t _h _e _ _B _O _L _exc _nl
		_A _n _d _ _t _o _ _a _l _l _ _t _h _e _ _b _a _r _d _s _nl _ _ _a _r _o _u _n _d _ _t _h _e _ _w _o _r _l _d _exc
		
    
    	return col;
}

float hash(float x)
{
	return fract(fract(x) * 1234.56789 * fract(fract(-x) * 1234.56789));	
}


void main()
{
	pth1 = path(t+.3)+origin+vec3(0.,.01,0.);
	vec2 uv = gl_FragCoord.xy / resolution.xy*2.-1.;
	vec2 uv2=uv;
#ifdef ENABLE_POSTPROCESS
	uv*=1.+pow(length(uv2*uv2*uv2*uv2),4.)*.07;
#endif
	uv.y*=resolution.y/resolution.x;
	vec2 mouse=(mouse.xy/resolution.xy-.5)*3.;
	if (mouse.x<1.) mouse=vec2(0.);
	mat2 rotview1, rotview2;
	vec3 from=origin+move(rotview1,rotview2);
	vec3 dir=normalize(vec3(uv*.8,1.));
	dir.yz*=rot(mouse.y);
	dir.xz*=rot(mouse.x);
	dir.yz*=rotview2;
	dir.xz*=rotview1;
	vec3 color=raymarch(from,dir); 
	color=clamp(color,vec3(.0),vec3(1.));
	color=pow(color,vec3(GAMMA))*BRIGHTNESS;
	color=mix(vec3(length(color)),color,SATURATION);
#ifdef ENABLE_POSTPROCESS
	vec3 rain=vec3(0.);//pow(texture2D(iChannel0,uv2+iGlobalTime*7.25468).rgb,vec3(1.5));
	color=mix(rain,color,clamp(time*.5-.5,0.,1.));
	color*=1.-pow(length(uv2*uv2*uv2*uv2)*1.1,6.);
	uv2.y *= resolution.y / 360.0;
	color.r*=(.5+abs(.5-mod(uv2.y     ,.021)/.021)*.5)*1.5;
	color.g*=(.5+abs(.5-mod(uv2.y+.007,.021)/.021)*.5)*1.5;
	color.b*=(.5+abs(.5-mod(uv2.y+.014,.021)/.021)*.5)*1.5;
	color*=.9+rain*.35;
#endif
	gl_FragColor = vec4(color,1.);


	uv = gl_FragCoord.xy / DOWN_SCALE;
	vec2 duv = floor(gl_FragCoord.xy / DOWN_SCALE);
    duv.y = duv.y +50;
	vec3 pixel = Text(duv);
	vec3 col = pixel*1.9+0.1;
	//col *= (1.-distance(mod(uv,vec2(1.0)),vec2(0.65)))*1.2;
	
	gl_FragColor += vec4(vec3(col), 1.0);
}
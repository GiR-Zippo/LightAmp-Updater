//cocio_finsprit.xm
#version 330

#extension GL_OES_standard_derivatives : enable

precision highp float;

uniform float time;
uniform vec2 mouse;
uniform vec2 resolution;

#define iTime time/1000

// initially based on http://glslsandbox.com/e#18223.0 remixed by @aday_net_au
// NES Controller in GLSL - Let's call him Nessy :)

#ifdef GL_ES
precision mediump float;
#endif

uniform vec4 inputColour; // Using these for ShaderLoader.dll -- GLSLSandbox.com doesn't support these :\

varying vec2 surfacePosition;

#define PI 3.14159265359 

#define border 0.0007



float distanceToSegment(vec2 pos,vec2 p1,vec2 p2){
	vec2 v=p2-p1;
		pos-=p1;
  float t = dot(pos , v)/ dot(v,v);
        if (t < 0.0){ return length(pos);}       // Beyond the 'v' end of the segment
	else if (t > 1.0){ return distance(pos, v);}  // Beyond the 'w' end of the segment
  return distance(pos, t*v);
}

vec4 line(vec2 p,vec2 l1, vec2 l2, float d,vec4 color) {
	float l=distanceToSegment(p,l1,l2);
  return abs(l-d)<border?vec4(0,0,0,1):l<d?color:vec4(0);
}

vec4 circle(vec2 p,vec2 o,float radius,vec4 color){
	float d=distance(p,o);
	if(abs(d-radius)<border){
		return vec4(0,0,0,1);
	}else{
		return step(distance(p,o),radius)*color;
	}
}

vec4 arc(vec2 p,vec2 o,float r1,float r2,float startangle,float endangle,vec4 color){
	p=p-o;
	float a=atan(p.y,p.x)+PI;
	if(a>endangle){
		return step(distance(p,r1*vec2(cos(endangle-PI),sin(endangle-PI))),r2)*color;
	}else if(a<startangle){
		return step(distance(p,r1*vec2(cos(startangle-PI),sin(startangle-PI))),r2)*color;
	}else{
		return step(abs(length(p)-r1),r2)*color;
	}
}

float powLength(vec2 p,float e){
	return pow(pow(abs(p.x),e)+pow(abs(p.y),e),1.0/e);
}



vec4 box(vec2 p,vec2 size,vec2 o,vec4 color){
	p-=o;
	vec2 d=abs(p);
	if(d.x<size.x&&d.y<size.y){
		return (abs(d.x-size.x)<border*2.0||abs(d.y-size.y)<border*2.0)?vec4(0,0,0,1):color;
	}
	return vec4(0);
}
vec4 roundBox(vec2 p,vec2 size,vec2 o,vec4 color,float factor){
	p-=o;
	p.y*=size.x/size.y;
	float f=powLength(p,factor);
	return abs(f-size.x)<border?vec4(0,0,0,1):step(f,size.x)*color;
}

vec4 overlay(vec4 bg,vec4 fg){
	return mix(vec4(bg.xyz,1.0),vec4(fg.xyz,1.0),fg.w);
}

vec3 background(vec2 p){
	return mix(mix(vec3(mouse.x,mouse.y,0.1),vec3(1.0,0.0,1.5),length(p)),vec3(0.3,0.5,1.0),step(mouse.y,fract(10.0*(atan(p.x,p.y)-0.1*iTime)/(mouse.x*PI))));
}

//Text
vec3 color;

#define LAYERS 5
#define SPEED 100.
#define SIZE 4.
#define goTYPE vec2 p = ( gl_FragCoord.xy /resolution.xy ) * vec2(64,32);vec3 c = vec3(0);vec2 cpos = vec2(1,26);vec3 txColor = vec3(1);
#define goPRINT gl_FragColor = vec4(c, 1.0);

#define slashN cpos = vec2(1,cpos.y-6.);

#define inBLK txColor = vec3(0);
#define inWHT txColor = vec3(1);
#define inRED txColor = vec3(1,0,0);
#define inGRN txColor = vec3(0,1,0);
#define inBLU txColor = vec3(0,0,1);

// via http://www.dafont.com/pixelzim3x5.font

float screen_ratio = resolution.y / resolution.x;

const float font_width = 3.0;
const int font_height = 5;
const int text_length = 56;
int text[56];

vec3 font_spc[font_height];	// 0
vec3 font_a[font_height];	// 1
vec3 font_b[font_height];	// 2
vec3 font_c[font_height];	// 3
vec3 font_d[font_height];	// 4
vec3 font_e[font_height];	// 5
vec3 font_f[font_height];	// 6
vec3 font_g[font_height];	// 7
vec3 font_h[font_height];	// 8
vec3 font_i[font_height];	// 9
vec3 font_j[font_height];	// 10
vec3 font_k[font_height];	// 11
vec3 font_l[font_height];	// 12
vec3 font_m[font_height];	// 13
vec3 font_n[font_height];	// 14
vec3 font_o[font_height];	// 15
vec3 font_p[font_height];	// 16
vec3 font_q[font_height];	// 17
vec3 font_r[font_height];	// 18
vec3 font_s[font_height];	// 19
vec3 font_t[font_height];	// 20
vec3 font_u[font_height];	// 21
vec3 font_v[font_height];	// 22
vec3 font_w[font_height];	// 23
vec3 font_x[font_height];	// 24
vec3 font_y[font_height];	// 25
vec3 font_z[font_height];	// 26

void init_text()
{
	text[0] = 12;
	text[1] = 9;
	text[2] = 7;
	text[3] = 8;
	text[4] = 20;
	text[5] = 1;
	text[6] = 13;
	text[7] = 16;
	text[8] = 0;
	text[9] = 1;
	text[10] = 0;
	text[11] = 14;
	text[12] = 5;
	text[13] = 23;
	text[14] = 0;
	text[15] = 18;
	text[16] = 5;
	text[17] = 12;
	text[18] = 5;
	text[19] = 1;
	text[20] = 19;
	text[21] = 5;
	text[22] = 0;

	text[23] = 14;
	text[24] = 15;
	text[25] = 0;
	text[26] = 20;
	text[27] = 9;
	text[28] = 13;
	text[29] = 5;
	text[30] = 0;
	text[31] = 6;
	text[32] = 15;
	text[33] = 18;
	text[34] = 0;

	text[35] = 7;
	text[36] = 18;
	text[37] = 5;
	text[38] = 5;
	text[39] = 20;
	text[40] = 9;
	text[41] = 14;
	text[42] = 24;
	text[43] = 0;
	text[44] = 3;
	text[45] = 21;
	text[46] = 26;
	text[47] = 0;
	text[48] = 20;
	text[49] = 8;
	text[50] = 5;
	text[51] = 0;
	text[52] = 2;
	text[53] = 5;
	text[54] = 1;
	text[55] = 20;

}

void init_fonts()
{
	font_spc[0] = vec3(0.0, 0.0, 0.0);
	font_spc[1] = vec3(0.0, 0.0, 0.0);
	font_spc[2] = vec3(0.0, 0.0, 0.0);
	font_spc[3] = vec3(0.0, 0.0, 0.0);
	font_spc[4] = vec3(0.0, 0.0, 0.0);

	font_a[0] = vec3(0.1, 1.0, 0.2);
	font_a[1] = vec3(1.0, 0.9, 1.0);
	font_a[2] = vec3(1.0, 1.0, 1.0);
	font_a[3] = vec3(1.0, 0.5, 1.0);
	font_a[4] = vec3(1.0, 0.0, 1.0);

	font_b[0] = vec3(1.0, 1.0, 0.0);
	font_b[1] = vec3(1.0, 0.0, 1.0);
	font_b[2] = vec3(1.0, 1.0, 1.0);
	font_b[3] = vec3(1.0, 0.0, 1.0);
	font_b[4] = vec3(1.0, 1.0, 0.0);

	font_c[0] = vec3(0.0, 1.0, 0.0);
	font_c[1] = vec3(1.0, 0.0, 1.0);
	font_c[2] = vec3(1.0, 0.0, 0.0);
	font_c[3] = vec3(1.0, 0.0, 1.0);
	font_c[4] = vec3(0.0, 1.0, 0.0);
	
	font_d[0] = vec3(1.0, 1.0, 0.0);
	font_d[1] = vec3(1.0, 0.0, 1.0);
	font_d[2] = vec3(1.0, 0.0, 1.0);
	font_d[3] = vec3(1.0, 0.0, 1.0);
	font_d[4] = vec3(1.0, 1.0, 0.0);
	
	font_e[0] = vec3(1.0, 1.0, 1.0);
	font_e[1] = vec3(1.0, 0.7, 0.0);
	font_e[2] = vec3(1.0, 1.0, 0.0);
	font_e[3] = vec3(1.0, 0.7, 0.0);
	font_e[4] = vec3(1.0, 1.0, 1.0);
	
	font_f[0] = vec3(1.0, 1.0, 1.0);
	font_f[1] = vec3(1.0, 0.0, 0.0);
	font_f[2] = vec3(1.0, 1.0, 0.0);
	font_f[3] = vec3(1.0, 0.0, 0.0);
	font_f[4] = vec3(1.0, 0.0, 0.0);
	
	font_g[0] = vec3(0.1, 1.0, 0.2);
	font_g[1] = vec3(1.0, 0.0, 1.0);
	font_g[2] = vec3(1.0, 0.0, 0.0);
	font_g[3] = vec3(1.0, 0.0, 1.0);
	font_g[4] = vec3(0.4, 0.3, 1.0);
	
	font_h[0] = vec3(1.0, 0.0, 1.0);
	font_h[1] = vec3(1.0, 0.0, 1.0);
	font_h[2] = vec3(1.0, 1.0, 1.0);
	font_h[3] = vec3(1.0, 0.0, 1.0);
	font_h[4] = vec3(1.0, 0.0, 1.0);

	font_i[0] = vec3(0.0, 1.0, 0.0);
	font_i[1] = vec3(0.0, 1.0, 0.0);
	font_i[2] = vec3(0.0, 1.0, 0.0);
	font_i[3] = vec3(0.0, 1.0, 0.0);
	font_i[4] = vec3(0.0, 1.0, 0.0);

	font_j[0] = vec3(0.0, 0.0, 1.0);
	font_j[1] = vec3(0.0, 0.0, 1.0);
	font_j[2] = vec3(0.0, 0.0, 1.0);
	font_j[3] = vec3(1.0, 0.0, 1.0);
	font_j[4] = vec3(0.0, 1.0, 0.0);
	
	font_k[0] = vec3(1.0, 0.0, 1.0);
	font_k[1] = vec3(1.0, 0.0, 1.0);
	font_k[2] = vec3(1.0, 1.0, 0.0);
	font_k[3] = vec3(1.0, 0.0, 1.0);
	font_k[4] = vec3(1.0, 0.0, 1.0);
	
	font_l[0] = vec3(1.0, 0.0, 0.0);
	font_l[1] = vec3(1.0, 0.0, 0.0);
	font_l[2] = vec3(1.0, 0.0, 0.0);
	font_l[3] = vec3(1.0, 0.0, 0.0);
	font_l[4] = vec3(1.0, 1.0, 1.0);
	
	font_m[0] = vec3(0.2, 0.6, 0.1);
	font_m[1] = vec3(1.0, 1.0, 1.0);
	font_m[2] = vec3(1.0, 0.0, 1.0);
	font_m[3] = vec3(1.0, 0.0, 1.0);
	font_m[4] = vec3(1.0, 0.0, 1.0);

	font_n[0] = vec3(1.0, 0.0, 0.0);
	font_n[1] = vec3(1.0, 1.0, 0.0);
	font_n[2] = vec3(1.0, 0.0, 1.0);
	font_n[3] = vec3(1.0, 0.0, 1.0);
	font_n[4] = vec3(1.0, 0.0, 1.0);
	
	font_o[0] = vec3(0.1, 1.0, 0.2);
	font_o[1] = vec3(1.0, 0.5, 1.0);
	font_o[2] = vec3(1.0, 0.0, 1.0);
	font_o[3] = vec3(1.0, 0.6, 1.0);
	font_o[4] = vec3(0.4, 1.0, 0.3);
	
	font_p[0] = vec3(1.0, 1.0, 0.2);
	font_p[1] = vec3(1.0, 0.9, 1.0);
	font_p[2] = vec3(1.0, 1.0, 0.3);
	font_p[3] = vec3(1.0, 0.0, 0.0);
	font_p[4] = vec3(1.0, 0.0, 0.0);

	font_q[0] = vec3(0.0, 1.0, 0.0);
	font_q[1] = vec3(1.0, 0.0, 1.0);
	font_q[2] = vec3(1.0, 0.0, 1.0);
	font_q[3] = vec3(1.0, 1.0, 0.0);
	font_q[4] = vec3(0.0, 0.0, 1.0);
	
	font_r[0] = vec3(1.0, 1.0, 0.2);
	font_r[1] = vec3(1.0, 0.9, 0.3);
	font_r[2] = vec3(1.0, 1.0, 0.7);
	font_r[3] = vec3(1.0, 0.5, 0.2);
	font_r[4] = vec3(1.0, 0.0, 1.0);
	
	font_s[0] = vec3(0.1, 1.0, 1.0);
	font_s[1] = vec3(1.0, 0.7, 0.0);
	font_s[2] = vec3(1.0, 1.0, 1.0);
	font_s[3] = vec3(0.0, 0.8, 1.0);
	font_s[4] = vec3(1.0, 1.0, 0.3);
	
	font_t[0] = vec3(1.0, 1.0, 1.0);
	font_t[1] = vec3(0.0, 1.0, 0.0);
	font_t[2] = vec3(0.0, 1.0, 0.0);
	font_t[3] = vec3(0.0, 1.0, 0.0);
	font_t[4] = vec3(0.0, 1.0, 0.0);
	
	font_u[0] = vec3(1.0, 0.0, 1.0);
	font_u[1] = vec3(1.0, 0.0, 1.0);
	font_u[2] = vec3(1.0, 0.0, 1.0);
	font_u[3] = vec3(1.0, 0.6, 1.0);
	font_u[4] = vec3(0.4, 1.0, 0.3);

	font_v[0] = vec3(1.0, 0.0, 1.0);
	font_v[1] = vec3(1.0, 0.0, 1.0);
	font_v[2] = vec3(1.0, 0.0, 1.0);
	font_v[3] = vec3(1.0, 0.0, 1.0);
	font_v[4] = vec3(0.0, 1.0, 0.0);
	
	font_w[0] = vec3(1.0, 0.0, 1.0);
	font_w[1] = vec3(1.0, 0.0, 1.0);
	font_w[2] = vec3(1.0, 0.0, 1.0);
	font_w[3] = vec3(1.0, 1.0, 1.0);
	font_w[4] = vec3(1.0, 0.5, 1.0);
	
	font_x[0] = vec3(1.0, 0.0, 1.0);
	font_x[1] = vec3(1.0, 0.0, 1.0);
	font_x[2] = vec3(0.0, 1.0, 0.0);
	font_x[3] = vec3(1.0, 0.0, 1.0);
	font_x[4] = vec3(1.0, 0.0, 1.0);

	font_y[0] = vec3(1.0, 0.0, 1.0);
	font_y[1] = vec3(1.0, 0.0, 1.0);
	font_y[2] = vec3(0.0, 1.0, 0.0);
	font_y[3] = vec3(0.0, 1.0, 0.0);
	font_y[4] = vec3(0.0, 1.0, 0.0);
	
	font_z[0] = vec3(1.0, 1.0, 1.0);
	font_z[1] = vec3(0.0, 0.0, 1.0);
	font_z[2] = vec3(0.0, 1.0, 0.0);
	font_z[3] = vec3(1.0, 0.0, 0.0);
	font_z[4] = vec3(1.0, 1.0, 1.0);
}

float draw_font_block(vec2 pixel_position, vec2 font_position, vec2 size, float tile_type)
{
	float gradient = 0.0;

	vec2 centered_font_position = font_position;

	if (abs(pixel_position.x - centered_font_position.x) <= size.x / 2.0 && 
	    abs(pixel_position.y - centered_font_position.y) <= size.y / 2.0)
	{
		if (tile_type==0.0 || tile_type==1.0)
		{
			gradient = tile_type;
		}
		else if (tile_type==0.1)
		{
			gradient = floor(2.58 - length(pixel_position - centered_font_position + vec2(-0.01, 0.01)) * 64.0);
		}
		else if (tile_type==0.2)
		{
			gradient = floor(2.58 - length(pixel_position - centered_font_position + vec2(0.01, 0.01)) * 64.0);
		}
		else if (tile_type==0.3)
		{
			gradient = floor(2.58 - length(pixel_position - centered_font_position + vec2(0.01, -0.01)) * 64.0);
		}
		else if (tile_type==0.4)
		{
			gradient = floor(2.58 - length(pixel_position - centered_font_position + vec2(-0.01, -0.01)) * 64.0);
		}
		else if (tile_type==0.5)
		{
			gradient = floor(length(pixel_position - centered_font_position) * 64.0);
			if (pixel_position.y <= centered_font_position.y) gradient = 0.0;
		}
		else if (tile_type==0.6)
		{
			gradient = floor(length(pixel_position - centered_font_position) * 64.0);
			if (pixel_position.y >= centered_font_position.y) gradient = 0.0;
		}
		else if (tile_type==0.7)
		{
			gradient = floor(length(pixel_position - centered_font_position) * 64.0);
			if (pixel_position.x >= centered_font_position.x) gradient = 0.0;
		}
		else if (tile_type==0.8)
		{
			gradient = floor(length(pixel_position - centered_font_position) * 64.0);
			if (pixel_position.x <= centered_font_position.x) gradient = 0.0;
		}
		else if (tile_type==0.9)
		{
			gradient = floor(length(pixel_position - centered_font_position) * 64.0);
		}

		if (gradient > 1.0) gradient = 1.0;
	}

	return gradient;
}

float draw_font_render(vec2 pixel_position, vec2 font_position, vec2 size, vec3 fontdata[5])
{
	float gradient = 0.0;
	font_position = (font_position - vec2(size.x * font_width, size.y * float(font_height)) / 2.0) * vec2(1.0, screen_ratio);

// for (float y=0.0; y<font_height; y++)
// {
//  for (float x=0.0; x<font_width; x++)
//  {
	gradient += draw_font_block(pixel_position, font_position + size * vec2(0.0, 4.0), size, fontdata[0].x);
	gradient += draw_font_block(pixel_position, font_position + size * vec2(1.0, 4.0), size, fontdata[0].y);
	gradient += draw_font_block(pixel_position, font_position + size * vec2(2.0, 4.0), size, fontdata[0].z);

	gradient += draw_font_block(pixel_position, font_position + size * vec2(0.0, 3.0), size, fontdata[1].x);
	gradient += draw_font_block(pixel_position, font_position + size * vec2(1.0, 3.0), size, fontdata[1].y);
	gradient += draw_font_block(pixel_position, font_position + size * vec2(2.0, 3.0), size, fontdata[1].z);

	gradient += draw_font_block(pixel_position, font_position + size * vec2(0.0, 2.0), size, fontdata[2].x);
	gradient += draw_font_block(pixel_position, font_position + size * vec2(1.0, 2.0), size, fontdata[2].y);
	gradient += draw_font_block(pixel_position, font_position + size * vec2(2.0, 2.0), size, fontdata[2].z);

	gradient += draw_font_block(pixel_position, font_position + size * vec2(0.0, 1.0), size, fontdata[3].x);
	gradient += draw_font_block(pixel_position, font_position + size * vec2(1.0, 1.0), size, fontdata[3].y);
	gradient += draw_font_block(pixel_position, font_position + size * vec2(2.0, 1.0), size, fontdata[3].z);

	gradient += draw_font_block(pixel_position, font_position + size * vec2(0.0, 0.0), size, fontdata[4].x);
	gradient += draw_font_block(pixel_position, font_position + size * vec2(1.0, 0.0), size, fontdata[4].y);
	gradient += draw_font_block(pixel_position, font_position + size * vec2(2.0, 0.0), size, fontdata[4].z);
// } }
	return gradient;
}

float draw_font(vec2 pixel_position, vec2 font_position, vec2 size, int font_number)
{
	float gradient = 0.0;
	if      (font_number==0)	gradient = draw_font_render(pixel_position, font_position, size, font_spc);
	else if (font_number==1)	gradient = draw_font_render(pixel_position, font_position, size, font_a);
	else if (font_number==2)	gradient = draw_font_render(pixel_position, font_position, size, font_b);
	else if (font_number==3)	gradient = draw_font_render(pixel_position, font_position, size, font_c);
	else if (font_number==4)	gradient = draw_font_render(pixel_position, font_position, size, font_d);
	else if (font_number==5)	gradient = draw_font_render(pixel_position, font_position, size, font_e);
	else if (font_number==6)	gradient = draw_font_render(pixel_position, font_position, size, font_f);
	else if (font_number==7)	gradient = draw_font_render(pixel_position, font_position, size, font_g);
	else if (font_number==8)	gradient = draw_font_render(pixel_position, font_position, size, font_h);
	else if (font_number==9)	gradient = draw_font_render(pixel_position, font_position, size, font_i);
	else if (font_number==10)	gradient = draw_font_render(pixel_position, font_position, size, font_j);
	else if (font_number==11)	gradient = draw_font_render(pixel_position, font_position, size, font_k);
	else if (font_number==12)	gradient = draw_font_render(pixel_position, font_position, size, font_l);
	else if (font_number==13)	gradient = draw_font_render(pixel_position, font_position, size, font_m);
	else if (font_number==14)	gradient = draw_font_render(pixel_position, font_position, size, font_n);
	else if (font_number==15)	gradient = draw_font_render(pixel_position, font_position, size, font_o);
	else if (font_number==16)	gradient = draw_font_render(pixel_position, font_position, size, font_p);
	else if (font_number==17)	gradient = draw_font_render(pixel_position, font_position, size, font_q);
	else if (font_number==18)	gradient = draw_font_render(pixel_position, font_position, size, font_r);
	else if (font_number==19)	gradient = draw_font_render(pixel_position, font_position, size, font_s);
	else if (font_number==20)	gradient = draw_font_render(pixel_position, font_position, size, font_t);
	else if (font_number==21)	gradient = draw_font_render(pixel_position, font_position, size, font_u);
	else if (font_number==22)	gradient = draw_font_render(pixel_position, font_position, size, font_v);
	else if (font_number==23)	gradient = draw_font_render(pixel_position, font_position, size, font_w);
	else if (font_number==24)	gradient = draw_font_render(pixel_position, font_position, size, font_x);
	else if (font_number==25)	gradient = draw_font_render(pixel_position, font_position, size, font_y);
	else if (font_number==26)	gradient = draw_font_render(pixel_position, font_position, size, font_z);	return gradient;
}

float draw_text(vec2 pixel_position, vec2 font_position, vec2 size)
{
	float gradient = 0.0;
	float font_size = size.x * (font_width + 1.0);
	gradient += draw_font(pixel_position, font_position + vec2(0.0, 0.0) * font_size, size, text[0]);
	gradient += draw_font(pixel_position, font_position + vec2(1.0, 0.0) * font_size, size, text[1]);
	gradient += draw_font(pixel_position, font_position + vec2(2.0, 0.0) * font_size, size, text[2]);
	gradient += draw_font(pixel_position, font_position + vec2(3.0, 0.0) * font_size, size, text[3]);
	gradient += draw_font(pixel_position, font_position + vec2(4.0, 0.0) * font_size, size, text[4]);
	gradient += draw_font(pixel_position, font_position + vec2(5.0, 0.0) * font_size, size, text[5]);
	gradient += draw_font(pixel_position, font_position + vec2(6.0, 0.0) * font_size, size, text[6]);
	gradient += draw_font(pixel_position, font_position + vec2(7.0, 0.0) * font_size, size, text[7]);
	gradient += draw_font(pixel_position, font_position + vec2(8.0, 0.0) * font_size, size, text[8]);
	gradient += draw_font(pixel_position, font_position + vec2(9.0, 0.0) * font_size, size, text[9]);
	gradient += draw_font(pixel_position, font_position + vec2(10.0, 0.0) * font_size, size, text[10]);
	gradient += draw_font(pixel_position, font_position + vec2(11.0, 0.0) * font_size, size, text[11]);
	gradient += draw_font(pixel_position, font_position + vec2(12.0, 0.0) * font_size, size, text[12]);
	gradient += draw_font(pixel_position, font_position + vec2(13.0, 0.0) * font_size, size, text[13]);
	gradient += draw_font(pixel_position, font_position + vec2(14.0, 0.0) * font_size, size, text[14]);
	gradient += draw_font(pixel_position, font_position + vec2(15.0, 0.0) * font_size, size, text[15]);
	gradient += draw_font(pixel_position, font_position + vec2(16.0, 0.0) * font_size, size, text[16]);
	gradient += draw_font(pixel_position, font_position + vec2(17.0, 0.0) * font_size, size, text[17]);
	gradient += draw_font(pixel_position, font_position + vec2(18.0, 0.0) * font_size, size, text[18]);
	gradient += draw_font(pixel_position, font_position + vec2(19.0, 0.0) * font_size, size, text[19]);
	gradient += draw_font(pixel_position, font_position + vec2(20.0, 0.0) * font_size, size, text[20]);
	gradient += draw_font(pixel_position, font_position + vec2(21.0, 0.0) * font_size, size, text[21]);
	gradient += draw_font(pixel_position, font_position + vec2(22.0, 0.0) * font_size, size, text[22]);
gradient += draw_font(pixel_position, font_position + vec2(23.0, 0.0) * font_size, size, text[23]);
gradient += draw_font(pixel_position, font_position + vec2(24.0, 0.0) * font_size, size, text[24]);
gradient += draw_font(pixel_position, font_position + vec2(25.0, 0.0) * font_size, size, text[25]);
gradient += draw_font(pixel_position, font_position + vec2(26.0, 0.0) * font_size, size, text[26]);
gradient += draw_font(pixel_position, font_position + vec2(27.0, 0.0) * font_size, size, text[27]);
gradient += draw_font(pixel_position, font_position + vec2(28.0, 0.0) * font_size, size, text[28]);
gradient += draw_font(pixel_position, font_position + vec2(29.0, 0.0) * font_size, size, text[29]);
gradient += draw_font(pixel_position, font_position + vec2(30.0, 0.0) * font_size, size, text[30]);
gradient += draw_font(pixel_position, font_position + vec2(31.0, 0.0) * font_size, size, text[31]);
gradient += draw_font(pixel_position, font_position + vec2(32.0, 0.0) * font_size, size, text[32]);
gradient += draw_font(pixel_position, font_position + vec2(33.0, 0.0) * font_size, size, text[33]);
gradient += draw_font(pixel_position, font_position + vec2(34.0, 0.0) * font_size, size, text[34]);
gradient += draw_font(pixel_position, font_position + vec2(35.0, 0.0) * font_size, size, text[35]);
gradient += draw_font(pixel_position, font_position + vec2(36.0, 0.0) * font_size, size, text[36]);
gradient += draw_font(pixel_position, font_position + vec2(37.0, 0.0) * font_size, size, text[37]);
gradient += draw_font(pixel_position, font_position + vec2(38.0, 0.0) * font_size, size, text[38]);
gradient += draw_font(pixel_position, font_position + vec2(39.0, 0.0) * font_size, size, text[39]);
gradient += draw_font(pixel_position, font_position + vec2(40.0, 0.0) * font_size, size, text[40]);
gradient += draw_font(pixel_position, font_position + vec2(41.0, 0.0) * font_size, size, text[41]);
gradient += draw_font(pixel_position, font_position + vec2(42.0, 0.0) * font_size, size, text[42]);
gradient += draw_font(pixel_position, font_position + vec2(43.0, 0.0) * font_size, size, text[43]);

gradient += draw_font(pixel_position, font_position + vec2(44.0, 0.0) * font_size, size, text[44]);
gradient += draw_font(pixel_position, font_position + vec2(45.0, 0.0) * font_size, size, text[45]);
gradient += draw_font(pixel_position, font_position + vec2(46.0, 0.0) * font_size, size, text[46]);
gradient += draw_font(pixel_position, font_position + vec2(47.0, 0.0) * font_size, size, text[47]);
gradient += draw_font(pixel_position, font_position + vec2(48.0, 0.0) * font_size, size, text[48]);
gradient += draw_font(pixel_position, font_position + vec2(49.0, 0.0) * font_size, size, text[49]);
gradient += draw_font(pixel_position, font_position + vec2(50.0, 0.0) * font_size, size, text[50]);
gradient += draw_font(pixel_position, font_position + vec2(51.0, 0.0) * font_size, size, text[51]);
gradient += draw_font(pixel_position, font_position + vec2(52.0, 0.0) * font_size, size, text[52]);
gradient += draw_font(pixel_position, font_position + vec2(53.0, 0.0) * font_size, size, text[53]);
gradient += draw_font(pixel_position, font_position + vec2(54.0, 0.0) * font_size, size, text[54]);
gradient += draw_font(pixel_position, font_position + vec2(55.0, 0.0) * font_size, size, text[55]);


	return gradient;
}


void main( void ) {
	vec2 p = ( gl_FragCoord.xy / resolution.xy );
	p.x*=resolution.x/resolution.y;
	p.y-=0.1;
	vec4 color;
	
	//legs
	color=overlay(color,max(line(p,vec2(0.6,0.2),vec2(0.6,0.0),0.02,vec4(0.2,0.6,0.8,1.0)),
			       line(p,vec2(0.5,0.0),vec2(0.58,0.0),0.02,vec4(0.2,0.6,0.8,1.0))));
	color=overlay(color,max(line(p,vec2(0.7,0.2),vec2(0.8,0.0),0.02,vec4(0.2,0.6,0.8,1.0)),
			       line(p,vec2(0.9,0.0),vec2(0.82,0.0),0.02,vec4(0.2,0.6,0.8,1.0))));
	
	p.y+=sin(iTime*6.)*0.1; // DANCE NES (Mouse.x dance with me)
	

	//arms
	float a = sin(iTime*12.)*0.1; // FLAIL ARMS (Mouse.y dance me)
	color=overlay(color,line(p,vec2(-0.0009-a*0.5,0.45+a),vec2(0.36,0.5),0.02,vec4(0.2,0.6,0.8,1.0))); // LEFT ARM
	color=overlay(color,line(p,vec2(1.66,0.5),vec2(8.0-a*15.15,1.5-a),0.02,vec4(0.2,0.6,0.8,1.0))); // RIGHT ARM

	

	// Gray border center square
	color=overlay(color,roundBox(p,vec2(0.65,0.30),vec2(1.0,0.4),vec4(0.9,0.9,0.8,0.9),20.0)); // 

	// The black inner center square
	color=overlay(color,roundBox(p,vec2(0.6,0.25),vec2(1.0,0.4),vec4(0.0,0.0,0.0,1.0),20.0));

   	color=overlay(color,circle(p,vec2(0.71,0.7),0.075,vec4(1,1,1,1)));
	color=overlay(color,circle(p,vec2(0.60,0.7),0.075,vec4(1,1,1,3)));
	color=overlay(color,circle(p,vec2(0.71,0.7),0.025,vec4(0,0,0,3)));
	color=overlay(color,circle(p,vec2(0.60,0.7),0.025,vec4(0,0,0,3)));


	
	
	
		// middle
	color=overlay(color,roundBox(p,vec2(0.20,0.03),vec2(1.0,0.61),vec4(0.5,0.5,0.5,1.0),20.0)); // 
	
	color=overlay(color,roundBox(p,vec2(0.20,0.03),vec2(1.0,0.51),vec4(0.5,0.5,0.5,1.0),20.0)); // 
        
	color=overlay(color,roundBox(p,vec2(0.20,0.03),vec2(1.0,0.41),vec4(0.5,0.5,0.5,1.0),20.0)); // 
        
	color=overlay(color,roundBox(p,vec2(0.20,0.06),vec2(1.0,0.29),vec4(0.5,0.5,0.5,1.5),20.0)); // 
        
	color=overlay(color,roundBox(p,vec2(0.20,0.02),vec2(1.0,0.18),vec4(0.5,0.5,0.5,1.0),20.0)); // 
        
	
	color=overlay(color,roundBox(p,vec2(0.07,0.07),vec2(1.48,0.25),vec4(0.5,0.5,0.5,1.0),20.0));
        color=overlay(color,roundBox(p,vec2(0.07,0.07),vec2(1.32,0.25),vec4(0.5,0.5,0.5,1.0),20.0));
	color=overlay(color,circle(p,vec2(1.48,0.25),0.05,vec4(1,0,0.4,1)*(0.9+0.1*cos(iTime*5.4+33.))));
	color=overlay(color,circle(p,vec2(1.32,0.25),0.05,vec4(1,0,0.4,1)*(0.9+0.1*cos(iTime*5.4+33.))));
	
	
	// dpad
	color=overlay(color,max(box(p,vec2(0.09,0.045),vec2(0.62,0.32),vec4(1.1,1.1,1.1,1.1)),
				box(p,vec2(0.045,0.09),vec2(0.62,0.32),vec4(1.1,1.1,1.1,1.1))));
	
	
	color=overlay(color,max(box(p,vec2(0.08,0.025),vec2(0.62,0.32),vec4(0.0,0.0,0.0,1.0)),
				box(p,vec2(0.04,0.08),vec2(0.62,0.32),vec4(0.0,0.0,0.0,1.0))));

	
	//round buttons start select
	color=overlay(color,roundBox(p,vec2(0.07,0.02),vec2(1.09,0.3),vec4(0.0,0.0,0.2,1.0),3.0));
	color=overlay(color,roundBox(p,vec2(0.07,0.02),vec2(0.9,0.3),vec4(0.0,0.0,0.2,1.0),3.0));
	
	vec4 preCol = color;

	//text
	init_fonts();
    init_text();
        
	float zoom = 0.03; // + sin(iTime) * 0.00001;
	float sx = zoom;
	float sy = zoom;

	vec2 pixel_position = vec2(gl_FragCoord.x / resolution.x, gl_FragCoord.y / resolution.x);

	float fx = 0.2;
	float fy = 0.5;
	float dispx = 1.0 - mod(0.25* iTime, float(text_length) * (font_width + 1.0) * sx + 1.5);
	float dispy = sin(pixel_position.x * 16.0 + 0.5 * iTime) * 0.025 + sin(pixel_position.x * 48.0 + 1.5 * iTime) * 0.0125;

	vec2 font_position = (vec2(fx, fy) + vec2(dispx, dispy));
	vec2 font_size = vec2(sx, sy);

	float gradient = draw_text(pixel_position, font_position, font_size);

	float r = 0.2 + sin(iTime * 2.2);
	float g = 0.3 + sin(iTime * 1.5);
	float b = 0.7 + sin(iTime + pixel_position.y * 32.0);
	
    vec4 font_color = vec4(vec3(r,g,b) * gradient, 1.0);
	if (time > 27500)
		font_color = vec4(vec3(r,g,b), 1.0);
	vec4 background_color = vec4(vec3(abs(pixel_position.y - 0.5 * screen_ratio)) * 0.0, 1.0);
	
     	float x = gl_FragCoord.x;
		p = gl_FragCoord.xy / resolution.xy;

        //Mixing it all together

	
	color += mix(background_color, font_color, 0.5);	
	
	gl_FragColor=color;

}
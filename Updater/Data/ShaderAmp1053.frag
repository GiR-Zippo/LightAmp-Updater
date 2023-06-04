//symmetric4.xm
#version 330

#extension GL_OES_standard_derivatives : enable

#ifdef GL_ES
precision highp float;
#endif

// glslsandbox uniforms
uniform float time;
uniform vec2 resolution;

// shadertoy emulation
#define iTime time/1000
#define iResolution resolution

// Emulate a black texture
#define texture(s, uv) vec4(0.0)

// Emulate some GLSL ES 3.x
vec3 tanh(vec3 x) {
    vec3 ex = exp(2.0 * x);
    return ((ex - 1.) / (ex + 1.));
}

// --------[ Original ShaderToy begins here ]---------- //
// CC0 - Neonwave sunrise
//  Inspired by a tweet by I wanted to create something that looked
//  a bit like the tweet. This is the result.

#define RESOLUTION    iResolution
#define TIME          iTime
#define PI            3.141592654
#define TAU           (2.0*PI)

#define SHOW_FFT


// License: WTFPL, author: sam hocevar, found: https://stackoverflow.com/a/17897228/418488
const vec4 hsv2rgb_K = vec4(1.0, 2.0 / 3.0, 1.0 / 3.0, 3.0);
vec3 hsv2rgb(vec3 c) {
  vec3 p = abs(fract(c.xxx + hsv2rgb_K.xyz) * 6.0 - hsv2rgb_K.www);
  return c.z * mix(hsv2rgb_K.xxx, clamp(p - hsv2rgb_K.xxx, 0.0, 1.0), c.y);
}
// License: WTFPL, author: sam hocevar, found: https://stackoverflow.com/a/17897228/418488
//  Macro version of above to enable compile-time constants
#define HSV2RGB(c)  (c.z * mix(hsv2rgb_K.xxx, clamp(abs(fract(c.xxx + hsv2rgb_K.xyz) * 6.0 - hsv2rgb_K.www) - hsv2rgb_K.xxx, 0.0, 1.0), c.y))

// License: Unknown, author: Unknown, found: don't remember
vec4 alphaBlend(vec4 back, vec4 front) {
  float w = front.w + back.w*(1.0-front.w);
  vec3 xyz = (front.xyz*front.w + back.xyz*back.w*(1.0-front.w))/w;
  return w > 0.0 ? vec4(xyz, w) : vec4(0.0);
}

// License: Unknown, author: Unknown, found: don't remember
vec3 alphaBlend(vec3 back, vec4 front) {
  return mix(back, front.xyz, front.w);
}

// License: Unknown, author: Unknown, found: don't remember
float tanh_approx(float x) {
  //  Found this somewhere on the interwebs
  //  return tanh(x);
  float x2 = x*x;
  return clamp(x*(27.0 + x2)/(27.0+9.0*x2), -1.0, 1.0);
}

// License: Unknown, author: Unknown, found: don't remember
float hash(float co) {
  return fract(sin(co*12.9898) * 13758.5453);
}

// License: Unknown, author: Unknown, found: don't remember
float hash(vec2 p) {
  float a = dot (p, vec2 (127.1, 311.7));
  return fract(sin(a)*43758.5453123);
}

// Value noise: https://iquilezles.org/articles/morenoise
float vnoise(vec2 p) {
  vec2 i = floor(p);
  vec2 f = fract(p);

  vec2 u = f*f*(3.0-2.0*f);
//  vec2 u = f;

  float a = hash(i + vec2(0.0,0.0));
  float b = hash(i + vec2(1.0,0.0));
  float c = hash(i + vec2(0.0,1.0));
  float d = hash(i + vec2(1.0,1.0));

  float m0 = mix(a, b, u.x);
  float m1 = mix(c, d, u.x);
  float m2 = mix(m0, m1, u.y);

  return m2;
}

// License: MIT, author: Inigo Quilez, found: https://iquilezles.org/www/articles/spherefunctions/spherefunctions.htm
vec2 raySphere(vec3 ro, vec3 rd, vec4 sph) {
  vec3 oc = ro - sph.xyz;
  float b = dot( oc, rd );
  float c = dot( oc, oc ) - sph.w*sph.w;
  float h = b*b - c;
  if( h<0.0 ) return vec2(-1.0);
  h = sqrt( h );
  return vec2(-b - h, -b + h);
}

// License: MIT OR CC-BY-NC-4.0, author: mercury, found: https://mercury.sexy/hg_sdf/
float mod1(inout float p, float size) {
  float halfsize = size*0.5;
  float c = floor((p + halfsize)/size);
  p = mod(p + halfsize, size) - halfsize;
  return c;
}

// License: MIT OR CC-BY-NC-4.0, author: mercury, found: https://mercury.sexy/hg_sdf/
vec2 mod2(inout vec2 p, vec2 size) {
  vec2 c = floor((p + size*0.5)/size);
  p = mod(p + size*0.5,size) - size*0.5;
  return c;
}

// License: Unknown, author: Unknown, found: don't remember
vec2 hash2(vec2 p) {
  p = vec2(dot (p, vec2 (127.1, 311.7)), dot (p, vec2 (269.5, 183.3)));
  return fract(sin(p)*43758.5453123);
}

float hifbm(vec2 p) {
  const float aa = 0.5;
  const float pp = 2.0-0.;

  float sum = 0.0;
  float a   = 1.0;

  for (int i = 0; i < 5; ++i) {
    sum += a*vnoise(p);
    a *= aa;
    p *= pp;
  }

  return sum;
}

float lofbm(vec2 p) {
  const float aa = 0.5;
  const float pp = 2.0-0.;

  float sum = 0.0;
  float a   = 1.0;

  for (int i = 0; i < 2; ++i) {
    sum += a*vnoise(p);
    a *= aa;
    p *= pp;
  }

  return sum;
}

float hiheight(vec2 p) {
  return hifbm(p)-1.8;
}

float loheight(vec2 p) {
  return lofbm(p)-2.15;
}

vec4 plane(vec3 ro, vec3 rd, vec3 pp, vec3 npp, vec3 off, float n) {
  float h = hash(n);
  float s = mix(0.05, 0.25, h);

  vec3 hn;
  vec2 p = (pp-off*2.0*vec3(1.0, 1.0, 0.0)).xy;

  const vec2 stp = vec2(0.5, 0.33);
  float he    = hiheight(vec2(p.x, pp.z)*stp);
  float lohe  = loheight(vec2(p.x, pp.z)*stp);

  float d = p.y-he;
  float lod = p.y - lohe;

  float aa = distance(pp, npp)*sqrt(1.0/3.0);
  float t = smoothstep(aa, -aa, d);

  float df = exp(-0.1*(distance(ro, pp)-2.));
  vec3 acol = hsv2rgb(vec3(mix(0.9, 0.6, df), 0.9, mix(1.0, 0.0, df)));
  vec3 gcol = hsv2rgb(vec3(0.6, 0.5, tanh_approx(exp(-mix(2.0, 8.0, df)*lod))));

  vec3 col = vec3(0.0);
  col += acol;
  col += 0.5*gcol;

  return vec4(col, t);
}

vec3 stars(vec2 sp, float hh) {
  const vec3 scol0 = HSV2RGB(vec3(0.85, 0.8, 1.0));
  const vec3 scol1 = HSV2RGB(vec3(0.65, 0.5, 1.0));
  vec3 col = vec3(0.0);

  const float m = 6.0;

  for (float i = 0.0; i < m; ++i) {
    vec2 pp = sp+0.5*i;
    float s = i/(m-1.0);
    vec2 dim  = vec2(mix(0.05, 0.003, s)*PI);
    vec2 np = mod2(pp, dim);
    vec2 h = hash2(np+127.0+i);
    vec2 o = -1.0+2.0*h;
    float y = sin(sp.x);
    pp += o*dim*0.5;
    pp.y *= y;
    float l = length(pp);

    float h1 = fract(h.x*1667.0);
    float h2 = fract(h.x*1887.0);
    float h3 = fract(h.x*2997.0);

    vec3 scol = mix(8.0*h2, 0.25*h2*h2, s)*mix(scol0, scol1, h1*h1);

    vec3 ccol = col + exp(-(mix(6000.0, 2000.0, hh)/mix(2.0, 0.25, s))*max(l-0.001, 0.0))*scol;
    ccol *= mix(0.125, 1.0, smoothstep(1.0, 0.99, sin(0.25*TIME+TAU*h.y)));
    col = h3 < y ? ccol : col;
  }

  return col;
}

vec3 toSpherical(vec3 p) {
  float r   = length(p);
  float t   = acos(p.z/r);
  float ph  = atan(p.y, p.x);
  return vec3(r, t, ph);
}

const vec3 lpos   = 1E6*vec3(0., -0.15, 1.0);
const vec3 ldir   = normalize(lpos);

vec4 moon(vec3 ro, vec3 rd) {
  const vec4 mdim   = vec4(1E5*vec3(0., 0.4, 1.0), 20000.0);
  const vec3 mcol0  = HSV2RGB(vec3(0.75, 0.7, 1.0));
  const vec3 mcol3  = HSV2RGB(vec3(0.75, 0.55, 1.0));

  vec2 md     = raySphere(ro, rd, mdim);
  vec3 mpos   = ro + rd*md.x;
  vec3 mnor   = normalize(mpos-mdim.xyz);
  float mdif  = max(dot(ldir, mnor), 0.0);
  float mf    = smoothstep(0.0, 10000.0, md.y - md.x);
  float mfre  = 1.0+dot(rd, mnor);
  float imfre = 1.0-mfre;

  vec3 col = vec3(0.0);
  col += mdif*mcol0*4.0;

#if defined(SHOW_FFT)
  vec3 fcol = vec3(0.0);
  vec2 msp    = toSpherical(-mnor.zxy).yz;
  vec2 omsp   = msp;
  float msf   = sin(msp.x);
  msp.x       -= PI*0.5;
  const float mszy = (TAU/(4.0))*0.125;
  float msny  = mod1(msp.y, mszy);
  msp.y *= msf;

  const int limit = 1;
  for (int i = -limit; i <= limit; ++i) {
    vec2 pp     = msp+vec2(0.0, mszy*float(i));
    float d0    = abs(pp.y);
    vec2 cp     = vec2(0.055*abs(msny-float(i)), 0.25);
    float fft   = texture(iChannel0, cp).x;
    float d1    = length(pp)-0.05*fft;
    float h     =mix(0.66, 0.99, fft);
    vec3 mcol1  = hsv2rgb(vec3(h, 0.55, 1.0));
    vec3 mcol2  = hsv2rgb(vec3(h, 0.85, 1.0));
    fcol += mcol1*0.5*tanh_approx(0.0025/max(d0, 0.0))*imfre*pow(msf, mix(100.0, 10.0, fft));
    fcol += mcol2*5.0*tanh_approx(0.00025/(max(d1, 0.0)*max(d1, 0.0)))*imfre*msf;
  }
  float d0   = abs(msp.x);
  fcol += mcol3*0.5*tanh_approx(0.0025/max(d0, 0.0))*imfre;

  const float start = 18.0;
  col += fcol*smoothstep(start, start+6.0+2.0*abs(omsp.y), TIME);

#endif

  return vec4(col, mf);
}


vec3 skyColor(vec3 ro, vec3 rd) {
  const vec3 acol   = HSV2RGB(vec3(0.6, 0.9, 0.075));
  const vec3 lpos   = 1E6*vec3(0., -0.15, 1.0);
  const vec3 lcol   = HSV2RGB(vec3(0.75, 0.8, 1.0));

  vec2 sp     = toSpherical(rd.xzy).yz;

  float lf    = pow(max(dot(ldir, rd), 0.0), 80.0);
  float li    = 0.02*mix(1.0, 10.0, lf)/(abs((rd.y+0.055))+0.025);
  float lz    = step(-0.055, rd.y);

  vec4 mcol   = moon(ro, rd);

  vec3 col = vec3(0.0);
  col += stars(sp, 0.25)*smoothstep(0.5, 0.0, li)*lz;
  col  = mix(col, mcol.xyz, mcol.w);
  col += smoothstep(-0.4, 0.0, (sp.x-PI*0.5))*acol;
  col += tanh(lcol*li);
  return col;
}

vec3 color(vec3 ww, vec3 uu, vec3 vv, vec3 ro, vec2 p) {
  float lp = length(p);
  vec2 np = p + 2.0/RESOLUTION.y;
  float rdd = 2.0;
  vec3 rd = normalize(p.x*uu + p.y*vv + rdd*ww);
  vec3 nrd = normalize(np.x*uu + np.y*vv + rdd*ww);

  const float planeDist = 1.0;
  const int furthest = 12;
  const int fadeFrom = 10;

  const float fadeDist = planeDist*float(fadeFrom);
  const float maxDist  = planeDist*float(furthest);
  float nz = floor(ro.z / planeDist);

  vec3 skyCol = skyColor(ro, rd);


  vec4 acol = vec4(0.0);
  const float cutOff = 0.95;
  bool cutOut = false;

  // Steps from nearest to furthest plane and accumulates the color
  for (int i = 1; i <= furthest; ++i) {
    float pz = planeDist*nz + planeDist*float(i);

    float pd = (pz - ro.z)/rd.z;

    vec3 pp = ro + rd*pd;

    if (pp.y < 0. && pd > 0.0 && acol.w < cutOff) {
      vec3 npp = ro + nrd*pd;

      vec3 off = vec3(0.0);

      vec4 pcol = plane(ro, rd, pp, npp, off, nz+float(i));

      float nz = pp.z-ro.z;
      float fadeIn = smoothstep(maxDist, fadeDist, pd);
      pcol.xyz = mix(skyCol, pcol.xyz, fadeIn);
      pcol = clamp(pcol, 0.0, 1.0);

      acol = alphaBlend(pcol, acol);
    } else {
      cutOut = true;
      acol.w = acol.w > cutOff ? 1.0 : acol.w;
      break;
    }

  }

  vec3 col = alphaBlend(skyCol, acol);
  return col;
}

vec3 effect(vec2 p, vec2 q) {
  float tm= TIME*1.05;
  vec3 ro = vec3(0.0, 0.0, tm);
  vec3 dro= normalize(vec3(0.0, 0.09, 1.0));
  vec3 ww = normalize(dro);
  vec3 uu = normalize(cross(normalize(vec3(0.0,1.0,0.0)), ww));
  vec3 vv = normalize(cross(ww, uu));

  vec3 col = color(ww, uu, vv, ro, p);

  return col;
}

// License: Unknown, author: nmz (twitter: @stormoid), found: https://www.shadertoy.com/view/NdfyRM
float sRGB(float t) { return mix(1.055*pow(t, 1./2.4) - 0.055, 12.92*t, step(t, 0.0031308)); }
// License: Unknown, author: nmz (twitter: @stormoid), found: https://www.shadertoy.com/view/NdfyRM
vec3 sRGB(in vec3 c) { return vec3 (sRGB(c.x), sRGB(c.y), sRGB(c.z)); }

// License: Unknown, author: Matt Taylor (https://github.com/64), found: https://64.github.io/tonemapping/
vec3 aces_approx(vec3 v) {
  v = max(v, 0.0);
  v *= 0.6;
  float a = 2.51;
  float b = 0.03;
  float c = 2.43;
  float d = 0.59;
  float e = 0.14;
  return clamp((v*(a*v+b))/(v*(c*v+d)+e), 0.0, 1.0);
}

void mainImage(out vec4 fragColor, in vec2 fragCoord) {
  vec2 q = fragCoord/RESOLUTION.xy;
  vec2 p = -1. + 2. * q;
  p.x *= RESOLUTION.x/RESOLUTION.y;
  vec3 col = vec3(0.0);
  col = effect(p, q);
  // col *= smoothstep(0.0, 8.0, TIME-abs(q.y));
  col = aces_approx(col);
  col = sRGB(col);
  fragColor = vec4(col, 1.0);
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
	
A _ N E W _ R E L E A S E _ F O R _ U s_dot B I G _ T H A N K S _ T O _ A L L _ T H E _ P P L _ O U T _ T H E R E _ W H O _ M A D E _ T H I S _ P O S S I B L E s_dot
Y O U _ A R E _ A W E S O M E s_dot s_dot s_dot

	vec3 color = mix(ch_color, vec3(0,0,0), 1.0- (0.09 / d*1.0));  // shading
	fragColor = vec4(color, 1.0);
}

void main(void)
{
    vec4 text_FragColor;
    vec4 scroltext_FragColor;
    mainImage(gl_FragColor, gl_FragCoord.xy);
    mainText(text_FragColor, gl_FragCoord.xy);
    secText(scroltext_FragColor, gl_FragCoord.xy);

    gl_FragColor += text_FragColor;
    gl_FragColor += scroltext_FragColor;
}
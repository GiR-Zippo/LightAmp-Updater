/*
 * Original shader from: https://www.shadertoy.com/view/wd23zz
 */
// Tune: ko0x_-_nur_ein_wort.it

#version 330

#ifdef GL_ES
precision mediump float;
#endif

// glslsandbox uniforms
uniform float time;
uniform vec2 resolution;
uniform float iTimeDelta;            // render time (in seconds)
uniform float iFrameRate;            // shader frame rate
uniform int iFrame;                // shader playback frame
uniform float iChannelTime[4];       // channel playback time (in seconds)
uniform vec3 iChannelResolution[4]; // channel resolution (in pixels)
uniform vec4 iMouse;                // mouse pixel coords. xy: current (if MLB down), zw: click
uniform vec4 iDate;                 // (year, month, day, time in seconds)             

out vec4 color;

// shadertoy globals
float iTime = 0.0;
vec3  iResolution = vec3(0.0);

// --------[ Original ShaderToy begins here ]---------- //
// CC0: For the neon style enjoyers
//  Or is it synthwave style? Don't know!
//  Anyone been tinkering with this for awhile and now want to get on with other stuff
//  Hopefully someone enjoys it.

//#define THAT_CRT_FEELING

#define TIME        iTime
#define RESOLUTION  iResolution
#define PI          3.141592654
#define PI_2        (0.5*PI)
#define TAU         (2.0*PI)
#define SCA(a)      vec2(sin(a), cos(a))
#define ROT(a)      mat2(cos(a), sin(a), -sin(a), cos(a))

// License: WTFPL, author: sam hocevar, found: https://stackoverflow.com/a/17897228/418488
const vec4 hsv2rgb_K = vec4(1.0, 2.0 / 3.0, 1.0 / 3.0, 3.0);
vec3 hsv2rgb(vec3 c) {
  vec3 p = abs(fract(c.xxx + hsv2rgb_K.xyz) * 6.0 - hsv2rgb_K.www);
  return c.z * mix(hsv2rgb_K.xxx, clamp(p - hsv2rgb_K.xxx, 0.0, 1.0), c.y);
}
// License: WTFPL, author: sam hocevar, found: https://stackoverflow.com/a/17897228/418488
//  Macro version of above to enable compile-time constants
#define HSV2RGB(c)  (c.z * mix(hsv2rgb_K.xxx, clamp(abs(fract(c.xxx + hsv2rgb_K.xyz) * 6.0 - hsv2rgb_K.www) - hsv2rgb_K.xxx, 0.0, 1.0), c.y))
// License: WTFPL, author: sam hocevar, found: https://stackoverflow.com/a/17897228/418488
vec3 rgb2hsv(vec3 c) {
  const vec4 K = vec4(0.0, -1.0 / 3.0, 2.0 / 3.0, -1.0);
  vec4 p = mix(vec4(c.bg, K.wz), vec4(c.gb, K.xy), step(c.b, c.g));
  vec4 q = mix(vec4(p.xyw, c.r), vec4(c.r, p.yzx), step(p.x, c.r));

  float d = q.x - min(q.w, q.y);
  float e = 1.0e-10;
  return vec3(abs(q.z + (q.w - q.y) / (6.0 * d + e)), d / (q.x + e), q.x);
}

const vec3 skyCol       = HSV2RGB(vec3(0.58, 0.86, 1.0));
const vec3 speCol1      = HSV2RGB(vec3(0.60, 0.25, 1.0));
const vec3 speCol2      = HSV2RGB(vec3(0.55, 0.25, 1.0));
const vec3 diffCol1     = HSV2RGB(vec3(0.60, 0.90, 1.0));
const vec3 diffCol2     = HSV2RGB(vec3(0.55, 0.90, 1.0));
const vec3 sunCol1      = HSV2RGB(vec3(0.60, 0.50, 0.5));
const vec3 sunDir2      = normalize(vec3(0., 0.82, 1.0));
const vec3 sunDir       = normalize(vec3(0.0, 0.05, 1.0));
const vec3 sunCol       = HSV2RGB(vec3(0.58, 0.86, 0.0005));
const float mountainPos = -20.0;

// License: MIT, author: Pascal Gilcher, found: https://www.shadertoy.com/view/flSXRV
float atan_approx(float y, float x) {
  float cosatan2 = x / (abs(x) + abs(y));
  float t = PI_2 - cosatan2 * PI_2;
  return y < 0.0 ? -t : t;
}

// License: Unknown, author: Unknown, found: don't remember
float tanh_approx(float x) {
  //  Found this somewhere on the interwebs
  //  return tanh(x);
  float x2 = x*x;
  return clamp(x*(27.0 + x2)/(27.0+9.0*x2), -1.0, 1.0);
}

vec3 toSpherical(vec3 p) {
  float r   = length(p);
  float t   = acos(p.z/r);
  float ph  = atan_approx(p.y, p.x);
  return vec3(r, t, ph);
}

// License: Unknown, author: nmz (twitter: @stormoid), found: https://www.shadertoy.com/view/NdfyRM
vec3 sRGB(vec3 t) {
  return mix(1.055*pow(t, vec3(1./2.4)) - 0.055, 12.92*t, step(t, vec3(0.0031308)));
}

// License: Unknown, author: Matt Taylor (https://github.com/64), found: https://64.github.io/tonemapping/
vec3 aces_approx(vec3 v) {
  v = max(v, 0.0);
  v *= 0.6f;
  float a = 2.51f;
  float b = 0.03f;
  float c = 2.43f;
  float d = 0.59f;
  float e = 0.14f;
  return clamp((v*(a*v+b))/(v*(c*v+d)+e), 0.0f, 1.0f);
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

// License: MIT, author: Inigo Quilez, found: https://iquilezles.org/www/articles/intersectors/intersectors.htm
float rayPlane(vec3 ro, vec3 rd, vec4 p) {
  return -(dot(ro,p.xyz)+p.w)/dot(rd,p.xyz);
}


// License: MIT, author: Inigo Quilez, found: https://iquilezles.org/www/articles/distfunctions2d/distfunctions2d.htm
float equilateralTriangle(vec2 p) {
  const float k = sqrt(3.0);
  p.x = abs(p.x) - 1.0;
  p.y = p.y + 1.0/k;
  if( p.x+k*p.y>0.0 ) p = vec2(p.x-k*p.y,-k*p.x-p.y)/2.0;
  p.x -= clamp( p.x, -2.0, 0.0 );
  return -length(p)*sign(p.y);
}

// License: MIT, author: Inigo Quilez, found: https://iquilezles.org/www/articles/distfunctions2d/distfunctions2d.htm
float box(vec2 p, vec2 b) {
  vec2 d = abs(p)-b;
  return length(max(d,0.0)) + min(max(d.x,d.y),0.0);
}

// License: MIT, author: Inigo Quilez, found: https://iquilezles.org/www/articles/distfunctions2d/distfunctions2d.htm
float segment(vec2 p, vec2 a, vec2 b) {
  vec2 pa = p-a, ba = b-a;
  float h = clamp( dot(pa,ba)/dot(ba,ba), 0.0, 1.0 );
  return length(pa - ba*h);
}

// License: Unknown, author: Unknown, found: don't remember
float hash(vec2 co) {
  return fract(sin(dot(co.xy ,vec2(12.9898,58.233))) * 13758.5453);
}

// License: MIT, author: Inigo Quilez, found: https://www.shadertoy.com/view/XslGRr
float vnoise(vec2 p) {
  vec2 i = floor(p);
  vec2 f = fract(p);
    
  vec2 u = f*f*(3.0-2.0*f);

  float a = hash(i + vec2(0.0,0.0));
  float b = hash(i + vec2(1.0,0.0));
  float c = hash(i + vec2(0.0,1.0));
  float d = hash(i + vec2(1.0,1.0));
  
  float m0 = mix(a, b, u.x);
  float m1 = mix(c, d, u.x);
  float m2 = mix(m0, m1, u.y);
  
  return m2;
}

// License: MIT, author: Inigo Quilez, found: https://www.iquilezles.org/www/articles/spherefunctions/spherefunctions.htm
vec2 raySphere(vec3 ro, vec3 rd, vec4 dim) {
  vec3 ce = dim.xyz;
  float ra = dim.w;
  vec3 oc = ro - ce;
  float b = dot( oc, rd );
  float c = dot( oc, oc ) - ra*ra;
  float h = b*b - c;
  if( h<0.0 ) return vec2(-1.0); // no intersection
  h = sqrt( h );
  return vec2( -b-h, -b+h );
}

vec3 skyRender(vec3 ro, vec3 rd) {
  vec3 col = vec3(0.0);
  col += 0.025*skyCol;
  col += skyCol*0.0033/pow((1.001+((dot(sunDir2, rd)))), 2.0);

  float tp0  = rayPlane(ro, rd, vec4(vec3(0.0, 1.0, 0.0), 4.0));
  float tp1  = rayPlane(ro, rd, vec4(vec3(0.0, -1.0, 0.0), 6.0));
  float tp = tp1;
  tp = max(tp0,tp1);


  if (tp1 > 0.0) {
    vec3 pos  = ro + tp1*rd;
    vec2 pp = pos.xz;
    float db = box(pp, vec2(5.0, 9.0))-3.0;
    
    col += vec3(4.0)*skyCol*rd.y*rd.y*smoothstep(0.25, 0.0, db);
    col += vec3(0.8)*skyCol*exp(-0.5*max(db, 0.0));
    col += 0.25*sqrt(skyCol)*max(-db, 0.0);
  }

  if (tp0 > 0.0) {
    vec3 pos  = ro + tp0*rd;
    vec2 pp = pos.xz;
    float ds = length(pp) - 0.5;
    
    col += (0.25)*skyCol*exp(-.5*max(ds, 0.0));
  }

  return clamp(col, 0.0, 10.0);
}

vec4 sphere(vec3 ro, vec3 rd, vec4 sdim) {
  vec2 si = raySphere(ro, rd, sdim);
  
  vec3 nsp = ro + rd*si.x;

  const vec3 lightPos1   = vec3(0.0, 10.0, 10.0);
  const vec3 lightPos2   = vec3(0.0, -80.0, 10.0);
  
  vec3 nld1   = normalize(lightPos1-nsp); 
  vec3 nld2   = normalize(lightPos2-nsp); 
  
  vec3 nnor   = normalize(nsp - sdim.xyz);

  vec3 nref   = reflect(rd, nnor);

  const float sf = 4.0;
  float ndif1 = max(dot(nld1, nnor), 0.0);
  ndif1       *= ndif1;
  vec3 nspe1  = pow(speCol1*max(dot(nld1, nref), 0.0), sf*vec3(1.0, 0.8, 0.5));

  float ndif2 = max(dot(nld2, nnor), 0.0);
  ndif2       *= ndif2;
  vec3 nspe2  = pow(speCol2*max(dot(nld2, nref), 0.0), sf*vec3(0.9, 0.5, 0.5));

  vec3 nsky   = skyRender(nsp, nref);
  float nfre  = 1.0+dot(rd, nnor);
  nfre        *= nfre;

  vec3 scol = vec3(0.0); 
  scol += nsky*mix(vec3(0.25), vec3(0.5, 0.5, 1.0), nfre);
  scol += diffCol1*ndif1;
  scol += diffCol2*ndif2;
  scol += nspe1;
  scol += nspe2;
  
  float t = tanh_approx(2.0*(si.y-si.x)/sdim.w);
  
  return vec4(scol, t);
}

vec3 sphereRender(vec3 ro, vec3 rd) {
  vec3 skyCol = skyRender(ro, rd);
  vec3 col = skyCol;
  const vec4 sdim0 = vec4(vec3(0.0), 2.0);
  vec4 scol0 = sphere(ro, rd, sdim0);
  col = mix(col, scol0.xyz, scol0.w);
  return col;
}

vec3 sphereEffect(vec2 p) {
  const float fov = tan(TAU/6.0);
  const vec3 ro = 1.0*vec3(0.0, 2.0, 5.0);
  const vec3 la = vec3(0.0, 0.0, 0.0);
  const vec3 up = vec3(0.0, 1.0, 0.0);

  vec3 ww = normalize(la - ro);
  vec3 uu = normalize(cross(up, ww));
  vec3 vv = cross(ww,uu);
  vec3 rd = normalize(-p.x*uu + p.y*vv + fov*ww);

  vec3 col = sphereRender(ro, rd);
  
  return col;
}

vec3 cityOfKali(vec2 p) {
  vec2 c = -vec2(0.5, 0.5)*1.12;

  float s = 2.0;
  vec2 kp = p/s;
 
  const float a = PI/4.0;
  const vec2 n = vec2(cos(a), sin(a));

  float ot2 = 1E6;
  float ot3 = 1E6;
  float n2 = 0.0;
  float n3 = 0.0;

  const float mx = 12.0;
  for (float i = 0.0; i < mx; ++i) {
    float m = (dot(kp, kp));
    s *= m;
    kp = abs(kp)/m + c;
    float d2 = (abs(dot(kp,n)))*s;
    if (d2 < ot2) {
      n2 = i;
      ot2 = d2;
    }
    float d3 = (dot(kp, kp));
    if (d3 < ot3) {
      n3 = i;
      ot3 = d3;
    }
  }
  vec3 col = vec3(0.0);
  n2 /= mx;
  n3 /= mx;
  col += 0.25*(hsv2rgb(vec3(0.8-0.2*n2*n2, 0.90, 0.025))/(sqrt(ot2)+0.0025));
  col += hsv2rgb(vec3(0.55+0.8*n3, 0.85, 0.00000025))/(ot3*ot3+0.000000025);
  return col;
}

vec3 outerSkyRender(vec3 ro, vec3 rd) {
  vec3 center = ro+vec3(-100.0, 40.0, 100.0);
  vec4 sdim = vec4(center, 50);
  vec2 pi = raySphere(ro, rd, sdim);
  const vec3 pn = normalize(vec3(0., 1.0, -0.8));
  vec4 pdim = vec4(pn, -dot(pn, center)); 
  float ri = rayPlane(ro, rd, pdim);

  vec3 col = vec3(0.0);
  
  col += sunCol/pow((1.001-((dot(sunDir, rd)))), 2.0);

  if (pi.x != -1.0) {
    vec3 pp = ro + rd*pi.x;
    vec3 psp= pp-sdim.xyz;
    vec3 pn = normalize(pp-sdim.xyz);
    psp = psp.zxy;
    psp.yz *= ROT(-0.5);
    psp.xy *= ROT(0.025*TIME);
    vec3 pss= toSpherical(psp);
    vec3 pcol = vec3(0.0);
    float dif = max(dot(pn, sunDir), 0.0);
    vec3 sc = 2000.0*sunCol;
    pcol += sc*dif;
    pcol += (cityOfKali(pss.yz))*smoothstep(0.125, 0.0, dif);
    pcol += pow(max(dot(reflect(rd, pn), sunDir), 0.0), 9.0)*sc;
    col = mix(col, pcol, tanh_approx(0.125*(pi.y-pi.x)));
    
  }

  vec3 gcol = vec3(0.0);

  vec3 rp = ro + rd*ri;
  float rl = length(rp-center);
  float rb = 1.55*sdim.w;
  float re = 2.45*sdim.w;
  float rw = 0.1*sdim.w;
  vec3 rcol = hsv2rgb(vec3(clamp((0.005*(rl+32.0)), 0.6, 0.8), 0.9, 1.0));
  gcol = rcol*0.025;
  if (ri > 0.0 && (pi.x == -1.0 || ri < pi.x)) {
    float mrl = rl;
    float nrl = mod1(mrl, rw);
    float rfre = 1.0+dot(rd, pn);
    vec3 rrcol = (rcol/max(abs(mrl), 0.1+smoothstep(0.7, 1.0, rfre))); 
    rrcol *= smoothstep(1.0, 0.3, rfre);
    rrcol *= smoothstep(re, re-0.5*rw, rl);
    rrcol *= smoothstep(rb-0.5*rw, rb, rl);
    col += rrcol;;
  }

  col += gcol/max(abs(rd.y), 0.0033);

return col;
}

vec3 triRender(vec3 col, vec3 ro, vec3 rd, inout float maxt) {
  const vec3 tpn = normalize(vec3(0.0, 0.0, 1.0));
  const vec4 tpdim = vec4(tpn, -2.0);
  float tpd = rayPlane(ro, rd, tpdim);

  if (tpd < 0.0 || tpd > maxt) {
    return col;
  }

  vec3 pp = ro+rd*tpd;
  vec2 p = pp.xy;
  p *= 0.5;

  const float off = 1.2-0.02;
  vec2 op = p; 
  p.y -= off;
  const vec2 n = SCA(-PI/3.0);
  vec2 gp = p;
  float hoff = 0.15*dot(n, p);
  vec3 gcol = hsv2rgb(vec3(clamp(0.7+hoff, 0.6, 0.8), 0.90, 0.02));
  vec2 pt = p;
  pt.y = -pt.y;
  const float zt = 1.0;
  float dt = equilateralTriangle(pt/zt)*zt;
//  col += 2.0*gcol;
  col = dt < 0.0 ? sphereEffect(1.5*(p)) : col;
  col += (gcol/max(abs(dt), 0.001))*smoothstep(0.25, 0.0, dt);
  if (dt < 0.0) {
    maxt = tpd;
  }
  return col;  
}

float heightFactor(vec2 p) {
  return 4.0*smoothstep(7.0, 0.5, abs(p.x))+.5;
}

float hifbm(vec2 p) {
  p *= 0.25;
  float hf = heightFactor(p);
  const float aa = 0.5;
  const float pp = 2.0-0.;

  float sum = 0.0;
  float a   = 1.0;
  
  for (int i = 0; i < 5; ++i) {
    sum += a*vnoise(p);
    a *= aa;
    p *= pp;
  }
  
  return hf*sum;
}

float hiheight(vec2 p) {
  return hifbm(p);
}

float lofbm(vec2 p) {
  p *= 0.25;
  float hf = heightFactor(p);
  const float aa = 0.5;
  const float pp = 2.0-0.;

  float sum = 0.0;
  float a   = 1.0;
  
  for (int i = 0; i < 3; ++i) {
    sum += a*vnoise(p);
    a *= aa;
    p *= pp;
  }
  
  return hf*sum;
}

float loheight(vec2 p) {
  return lofbm(p)-0.5;
}

vec3 mountainRender(vec3 col, vec3 ro, vec3 rd, bool flip, inout float maxt) {
  const vec3 tpn = normalize(vec3(0.0, 0.0, 1.0));
  const vec4 tpdim = vec4(tpn, mountainPos);
  float tpd = rayPlane(ro, rd, tpdim);

  if (tpd < 0.0 || tpd > maxt) {
    return col;
  }

  vec3 pp = ro+rd*tpd;
  vec2 p = pp.xy;
  const float cw = 1.0-0.25;
  float hz = 0.0*TIME+1.0;
  float lo = loheight(vec2(p.x, hz));
  vec2 cp = p;
  float cn = mod1(cp.x, cw);


  const float reps = 1.0;

  float d = 1E3;

  for (float i = -reps; i <= reps; ++i) {
    float x0 = (cn -0.5 + (i))*cw;
    float x1 = (cn -0.5 + (i + 1.0))*cw;
  
    float y0 = hiheight(vec2(x0, hz));
    float y1 = hiheight(vec2(x1, hz));
    
    float dd = segment(cp, vec2(-cw*0.5 + cw * float(i), y0), vec2(cw*0.5 + cw * float(i), y1));
    
    d = min(d, dd);
  }

  vec3 rcol = hsv2rgb(vec3(clamp(0.7+(0.5*(rd.x)), 0.6, 0.8), 0.95, 0.125));

  float sd = 1.0001-((dot(sunDir, rd)));

  vec3 mcol = col;
  float aa = fwidth(p.y);
  if ((dFdy(d) < 0.0) == !flip) {
    mcol *= mix(0.0, 1.0, smoothstep(aa, -aa, d-aa));
    mcol += HSV2RGB(vec3(0.55, 0.85, 0.8))*smoothstep(0.0, 5.0, lo-p.y);
    col = mcol;
    maxt = tpd;
  }
  col += 3.*rcol/(abs(d)+0.005+800.*sd*sd*sd*sd);
  col += HSV2RGB(vec3(0.55, 0.96, 0.075))/(abs(p.y)+0.05);

  return col;  
}

vec3 groundRender(vec3 col, vec3 ro, vec3 rd, inout float maxt) {
  const vec3 gpn = normalize(vec3(0.0, 1.0, 0.0));
  const vec4 gpdim = vec4(gpn, 0.0);
  float gpd = rayPlane(ro, rd, gpdim);

  if (gpd < 0.0) {
    return col;
  }
  
  maxt = gpd;
  
  vec3 gp     = ro + rd*gpd;
  float gpfre = 1.0 + dot(rd, gpn);
  gpfre *= gpfre;
  gpfre *= gpfre;
  gpfre *= gpfre;
  
  vec3 grr = reflect(rd, gpn);
  
  vec2 ggp    = gp.xz;
  ggp.y += TIME;
  float dfy   = dFdy(ggp.y);
  float gcf = sin(ggp.x)*sin(ggp.y);
  vec2 ggn    = mod2(ggp, vec2(1.0));
  float ggd   = min(abs(ggp.x), abs(ggp.y));

  vec3 gcol = hsv2rgb(vec3(0.7+0.1*gcf, 0.90, 0.02));
  
  float rmaxt = 1E6;
  vec3 rcol = outerSkyRender(gp, grr);
  rcol = mountainRender(rcol, gp, grr, true, rmaxt);
  rcol = triRender(rcol, gp, grr, rmaxt);

  col = gcol/max(ggd, 0.0+0.25*dfy)*exp(-0.25*gpd);
  rcol += HSV2RGB(vec3(0.65, 0.85, 1.0))*gpfre;
  rcol = 4.0*tanh(rcol*0.25);
  col += rcol*gpfre;

  return col;
}

vec3 render(vec3 ro, vec3 rd) {
  float maxt = 1E6;  

  vec3 col = outerSkyRender(ro, rd);
  col = groundRender(col, ro, rd, maxt);
  col = mountainRender(col, ro, rd, false, maxt);
  col = triRender(col, ro, rd, maxt);

  return col;
}

vec3 effect(vec2 p, vec2 pp) {
  const float fov = tan(TAU/6.0);
  const vec3 ro = 1.0*vec3(0.0, 1.0, -4.);
  const vec3 la = vec3(0.0, 1.0, 0.0);
  const vec3 up = vec3(0.0, 1.0, 0.0);

  vec3 ww = normalize(la - ro);
  vec3 uu = normalize(cross(up, ww));
  vec3 vv = cross(ww,uu);
  vec3 rd = normalize(-p.x*uu + p.y*vv + fov*ww);

  float aa = 2.0/RESOLUTION.y;

  vec3 col = render(ro, rd);
#if defined(THAT_CRT_FEELING)  
  col *= smoothstep(1.5, 0.5, length(pp));
  col *= 1.25*mix(vec3(0.5), vec3(1.0),smoothstep(-0.9, 0.9, sin(0.25*TAU*p.y/aa+TAU*vec3(0.0, 1., 2.0)/3.0)));
#endif  
  col -= 0.05*vec3(.00, 1.0, 2.0).zyx;
  col = aces_approx(col); 
  col = sRGB(col);
  return col;
}

void mainImage( out vec4 fragColor, in vec2 fragCoord ) {
  vec2 q = fragCoord/RESOLUTION.xy;

  vec2 p = -1. + 2. * q;
  vec2 pp = p;
  p.x *= RESOLUTION.x/RESOLUTION.y;
  vec3 col = effect(p, pp);

  fragColor = vec4(col, 1.0);
}
// --------[ Original ShaderToy ends here ]---------- //


const vec2 bitmap_size = vec2(96, 96);
const int[] palette = int[] (
0x00000000,0x00813f2c,0x008d4817,0x00744238,0x004b4561,0x00965c2c,0x00af6514,0x00a45c35,0x0054546d,0x004f5373,0x00835c4b,0x00724d9c,0x00565a7a,0x00605c73,0x00714ea6,
0x00796466,0x0056617d,0x005c627d,0x00b06d45,0x00586184,0x009a7339,0x00bb762e,0x0065657f,0x00616587,0x007257b9,0x00ba8323,0x00686e89,0x006d7083,0x00947765,0x00646e91,
0x00697186,0x006f708a,0x005c7196,0x00c5853e,0x007a7d7b,0x007a65ca,0x0067749f,0x0067769b,0x00727894,0x00677997,0x006d7b91,0x006077a5,0x006e7d8e,0x00657aa0,0x00b1923e,
0x00c59c27,0x00bd8b62,0x00d39d24,0x008376b9,0x00cc9349,0x007080a5,0x00727fab,0x006481b1,0x006981b0,0x008272da,0x006e85a5,0x00ac9170,0x006f84ab,0x007487a0,0x006a84ad,
0x007a889f,0x009b8c8a,0x00888d9a,0x006487b9,0x00668f9f,0x007989af,0x006d8ea3,0x007086bc,0x007e909a,0x008f948b,0x008187b9,0x009b9587,0x00768bb2,0x008b9393,0x00718ab9,
0x006c8db7,0x00728db6,0x006d8cbc,0x00e0b227,0x008982e6,0x007893bc,0x00939aa0,0x00dfac55,0x007493c4,0x009a9f96,0x007095c7,0x007697c2,0x00a9a48e,0x007198c4,0x008398c0,
0x008898bf,0x006c97ca,0x00daaa6f,0x008e9cb4,0x0099a39c,0x008da3a0,0x0099a2a1,0x007d9db9,0x00a5a599,0x00809bc5,0x00699acf,0x00a29cbd,0x00948ef3,0x007f9ed0,0x008ca1ca,
0x0075a0d5,0x00ada9ab,0x00d3b96e,0x0091acab,0x00aeaba7,0x0094a0d5,0x0082a5ca,0x006ba4d4,0x0072a3d9,0x0087adb5,0x0078a5d3,0x008ea8c4,0x009aaabc,0x0078a7cf,0x00b1b0a4,
0x00a5aeae,0x007fa6d3,0x009aaeb2,0x0096a6ce,0x0091a5d5,0x00b2b39f,0x00d6b58f,0x0085aac9,0x00e7c261,0x0074a9d2,0x009aa1e5,0x0092a8d1,0x009aa0f9,0x007aabe2,0x007fafd7,
0x0098aed8,0x0086b0d7,0x009fabe1,0x0099adde,0x00a3b4c6,0x0077afe6,0x007cb2db,0x0094b2d6,0x00a0b2d4,0x0094b4d2,0x00b1baba,0x008fb5d4,0x0094bac3,0x0079b3e4,0x00bfc0ac,
0x0085b7d9,0x00a4bace,0x00bbbebc,0x00d6ca95,0x008abae3,0x00ecd27d,0x0086bce6,0x007bbbee,0x0086bee2,0x00d8c6b2,0x00a9c3c9,0x00b1c2d4,0x00a7bde7,0x008ec2de,0x00a1bde9,
0x009ac4d6,0x009ac1e0,0x00a6c1dd,0x00a8bcee,0x007ec0ed,0x00bec7c7,0x00b3c0e1,0x008ec7d5,0x009bc7d2,0x00b9c8ca,0x00b3c8cc,0x0089c6f2,0x008dccf1,0x0093ceeb,0x00b4c8fa,
0x009acfeb,0x00b9caf4,0x00a6d1e3,0x00cad3d3,0x00adceec,0x00eadea8,0x00a6d5d9,0x00a7d0e9,0x00b3d0e6,0x009bd2e7,0x009ccff1,0x0095d4e3,0x009bd4e2,0x00b0cdfa,0x0097d4eb,
0x00c1d6da,0x00c2d1eb,0x00d7e0c6,0x008fd7f7,0x00c0dbe1,0x0096d8f7,0x00b6d3ff,0x00c3daee,0x00e1e0d4,0x00a5dcf2,0x00a5dfed,0x00acdfec,0x00b8e1e3,0x00b2e1e5,0x009adef6,
0x00d7e0e0,0x00a0e0f0,0x00c0d9fe,0x00b9dfef,0x0095dff9,0x00cce1e5,0x00c1def4,0x00c3defc,0x009ee3fb,0x00b3e5f9,0x009be8fa,0x00aee8f6,0x00a2e9fa,0x00a9e9f9,0x00b6eaf7,
0x00d6e7fa,0x00c9e7fe,0x00eeefdb,0x00beedf1,0x00b7edf3,0x00d1e8fd,0x00a5eef8,0x00c4eef0,0x00b2eef6,0x00d8edf1,0x00d3eef4,0x00a1f2f6,0x00f2efea,0x00dfeff0,0x009bf3f9,
0x00cfeffe,0x00d7f1fe,0x00b2f5fe,0x00baf6fe,0x00e0f3ff,0x00d3f5fd,0x00c1f7fd,0x00dbf7fd,0x00d0fafc,0x00cafafe,0x00defbfa,0x00d7fbfc,0x00f0fafa,0x00ebfbfd,0x00e5fbff,
0x00ffffff
);

const int longs_per_line = 24;
const int[] bitmap = int[] (
0x3fddf6fc, 0x5b5b5b34, 0x91916f71, 0x787891b7, 0xe5d07878, 0xdef2f8f8, 0x71d3e8f8, 0x3a5b649c, 0xcbedcb2c, 0x3c455757, 0xbc83a175, 0xf5f5f5f5, 0xd9e2f0f5, 0x838ea4b8, 0x8f8f8f68, 0xffffffd2, 0x15389395, 0x6d220206, 0xa3f9aa98, 0x7973863b, 0x96ac3b56, 0x79797979, 0x8e4c3339, 0xfcfcfcae,
0x3fe5f6fc, 0x5b5b5b34, 0x91919071, 0x78789191, 0xf9e57878, 0xe9dcf9f8, 0x71cdd1f5, 0x445b55a9, 0xcbedaa91, 0x3a2a5477, 0xbc7b1b44, 0xf0f5f0f5, 0x87a2d8f0, 0x83838387, 0x978f8f87, 0xffffffb7, 0x022a5e95, 0x2a0f0202, 0x7fe5cf54, 0x7986883b, 0x76e04c56, 0x79797979, 0xb76f3332, 0xffffffff,
0x58f3f8fd, 0x5b5b5b34, 0x91b79785, 0x78917878, 0xf9e57a78, 0xfbdcf6f8, 0x81c2c6db, 0xa0555b8d, 0xd298afc3, 0x5d425177, 0xbc3a161e, 0xe2f0f0f0, 0x838387b8, 0x837b7b83, 0xa18f8787, 0x98edfcae, 0x0a2a4777, 0x42020202, 0x63adf96c, 0x79797950, 0x76e06156, 0x4b797979, 0xfc784835, 0xffffffff,
0x96f3f8ea, 0x5b5b5b34, 0x91aa8b8d, 0x6078aaaa, 0xf9e57a60, 0xf5dddff8, 0x9cb2d1d1, 0xc3635b69, 0xe70000c3, 0x3a7f92ac, 0xbc3a2a0d, 0xa4d9f0f0, 0x7b7b8383, 0x7b7b7b7b, 0xa18f8787, 0x7d9591a0, 0x2a545757, 0x22020202, 0x4872e4ba, 0x86798656, 0x73cea54c, 0x39797979, 0xed955139, 0xffffffff,
0xa3f3f5fe, 0x69695b34, 0x91b7ae8d, 0xb791b791, 0xf9e5afaa, 0xccf5dcf9, 0x9c9ed1c6, 0xc7975b5b, 0xece700a0, 0x1661acbf, 0xd83a3a1e, 0x8387b8f0, 0x7b7b8383, 0x7b7b7b7b, 0xa1978787, 0x6260758b, 0x2a47777d, 0x1c05020a, 0x395da0f8, 0x79797379, 0x86bdce4b, 0x39677979, 0xcb7d7d51, 0xffffffff,
0xbdf3fefd, 0x73715b3f, 0xb7d2b78d, 0xd2aab7aa, 0xf9f9c7d2, 0xc8fbdff6, 0x9c8ddac6, 0xa0c33f55, 0xecec0000, 0x1128c2dc, 0xd81e3c3c, 0x688383a7, 0x7b7b7b7b, 0x7b7b7b7b, 0x1e758f87, 0x54607851, 0x47579895, 0x9505052a, 0x395d42e4, 0x58797379, 0x73a3f64b, 0x3c4b7979, 0x957d7d7d, 0xfffffffc,
0xd3f6fdfe, 0x69715b3f, 0xd2d2ee86, 0xd2d2d2d2, 0xf8f9eaee, 0xc6dbf9e7, 0x818ddac6, 0x00374b5b, 0xecefe700, 0x271a61ec, 0x87271e3a, 0x83688383, 0x7b7b7b7b, 0x7b7b7b7b, 0x1e268f83, 0x6d78d23e, 0x5757cbff, 0xeb142249, 0x395022a0, 0x4b867979, 0x8696f67f, 0x47327979, 0x957d7757, 0xffffffcb,
0xe7f9fdcf, 0x69715b5b, 0xeefcfc86, 0xaeb7d2ee, 0xf8f8eaee, 0xc6d6f8e7, 0x5ba9d1d6, 0x00483455, 0xf2ecec00, 0x3a1d27bf, 0x83321b3c, 0x83837b7b, 0x7b7b7b7b, 0x7b7b7b7b, 0x1b1b8f7b, 0x6091910d, 0x775777cb, 0xf95e3e62, 0x354c4738, 0x3b797979, 0x8686e5a5, 0x57495679, 0x57575757, 0xffffcb7d,
0xf2f8fd00, 0x85855b81, 0xfffcee88, 0xffeed7fc, 0xf9f8e8ee, 0xc8d6e5e9, 0x5b9cc2d6, 0x00003455, 0xd3ecefd3, 0x32251d3b, 0x7b411b1b, 0x837b7b7b, 0x7b7b7b7b, 0x7b7b7b7b, 0x161b757b, 0x78d21b08, 0x777d5745, 0xd0ba5747, 0x48487a38, 0x4c587379, 0x8686ccce, 0x77772a79, 0x7d777777, 0xedcb9577,
0xf6f5d200, 0x948d5b8d, 0xfcffffa4, 0xfceeeeee, 0xf6f9d5eb, 0xc8c6c2f6, 0x5b9cc2c6, 0x0000293f, 0x58e7ece7, 0x26321d1d, 0x835a1f1f, 0x8383836e, 0x83838383, 0x7b7b7b7b, 0x0d167583, 0xb7b71608, 0x4762625e, 0xa0f86257, 0x4c48ad45, 0x4b4b7379, 0x8686a3f6, 0x47472256, 0x7d7d7762, 0xffcb957d,
0xf8e80000, 0x858d5bb2, 0xfcffffa7, 0xfffffffc, 0xf6f9ebff, 0xc6c6b1e0, 0x558dc2c6, 0xcf004834, 0x2496d6ec, 0x1f321d17, 0x83681f1f, 0x87a2a287, 0x83838368, 0x7b7b7b7b, 0x0d167583, 0xb73e160d, 0x57455eaa, 0x91e4af7d, 0x4c39a578, 0x7f3b7979, 0x868686e0, 0x477d7d44, 0x7d624545, 0xff957d7d,
0xf7e40000, 0x858d5bc2, 0xffffffb8, 0xfffffcff, 0xf9f9f0fc, 0xc6c6c8b4, 0x5581d1c6, 0xe7004834, 0x1d2986d3, 0x1f1a2517, 0xd9a23226, 0xe2e2e2e2, 0x835987d8, 0x7b7b7b83, 0x0d0d5d7b, 0xaa0d1616, 0x624978b7, 0x6dcfcf95, 0x503990ba, 0xc04c5679, 0x6f8686b4, 0x627d7d62, 0x54495e5e, 0xcb957d62,
0xfa000000, 0x8c8d5bcc, 0xfffffcb8, 0xfcffffff, 0xf6f9e8fc, 0xc6c6d18d, 0x3f5bded6, 0xe7000034, 0x171d3476, 0x321a261a, 0xf0e2d890, 0xf0f0e2e2, 0x68b8e2f0, 0x7b688368, 0x160d5d87, 0x1b011116, 0x495e91aa, 0x5ea0e85e, 0x504863d0, 0xe04b3b79, 0x5d868686, 0x547d7d62, 0x5e495e5e, 0x9562455e,
0xfb000000, 0x8c7169e4, 0xfdeeee8e, 0xcfebcfea, 0xe0f9e5a6, 0xc6c6d181, 0x3f5bd1d6, 0xbf000034, 0x171a2934, 0xe2972525, 0xf0e2e2f0, 0xf0f0f0f0, 0xd9e2f0f0, 0x83835a83, 0x1b165d87, 0x1e101e16, 0x5e78aa78, 0x6c6ce46c, 0x50394be5, 0xcd7f3b56, 0x60738886, 0x5e577d57, 0x60495e5e, 0x5460495e,
0xf9000000, 0x8c7173f7, 0xfb000076, 0xe0f9f9f9, 0xccf9e594, 0xc6c6d19e, 0x3f5bc2d6, 0x4bc00034, 0x1d1a1d29, 0xf1f0d874, 0xf0f0f0f0, 0xf0f0f0f0, 0xe2e2f0f0, 0x835a87d9, 0x1b165d87, 0x1e1e1e11, 0x5eaab72a, 0x937ad0a0, 0x50394bdd, 0xb4ac3550, 0x62567986, 0x5e545757, 0x91915e5e, 0x78784949,
0xe9000000, 0x8c6988fe, 0xe8ae0076, 0xccf9f9f6, 0x9ef9e581, 0xd6c6d1b1, 0x34559ed6, 0x296f0034, 0xb8271325, 0xf1f0f5f0, 0xf0f0f0f1, 0xe2e2f0f0, 0xf0f0f0e2, 0x5a87d9e2, 0x161f3c8f, 0x1e1e1e1e, 0x91b77811, 0xd05fa0d0, 0x4c394ace, 0x86e93b3b, 0x545d7986, 0x5e5e4522, 0x9191915e, 0x78787849,
0xd3000000, 0x8c73bbfa, 0xcfae0076, 0xb2f6f9f6, 0x81e4f681, 0xc6c6d1b1, 0x34558dd6, 0x293b0035, 0xf5d8321d, 0xf1f7f7f7, 0xf0f0f0f1, 0xd9e2f0f0, 0xf0f0f0e2, 0x87d9e6f1, 0x161f415a, 0x1e1e1e1e, 0xaf911611, 0xe55f6ce4, 0x39394cbd, 0x86cd613b, 0x5e607673, 0x5e5e5e49, 0x5e919178, 0x5e605e5e,
0xd3000000, 0x8569f0f8, 0xc7000073, 0x9ef6f9f6, 0x70ccf681, 0xc6c8c8d1, 0x343f81da, 0x24294a39, 0xf7f7e241, 0xf7f7f7f7, 0xf0f0f1f1, 0xc9c9e2f0, 0xf1f1f0d9, 0xc1e2f1f1, 0x161f2787, 0x1a1e1e1e, 0xa02a1616, 0xe9727ad0, 0x39483ba3, 0x73a3ac39, 0x5e5e6173, 0x5e5e5e5e, 0x60607878, 0xaa91605e,
0xd3000000, 0x8586faf3, 0xea000073, 0x94f6f9f6, 0x8da3f676, 0xc6c8c6d1, 0x353464da, 0x25243541, 0xf7f7f7d8, 0xf7f7f7f7, 0xf0f0f7f7, 0xf1f1e2f0, 0xf1f1f0e6, 0xe2f1f1f7, 0x1f1a83a4, 0x161e1e1e, 0x721a1e0d, 0xe9ad7aba, 0x394839a3, 0x7979c035, 0x605e4479, 0x5e5e5e60, 0xaa605e5e, 0xd2b7ae91,
0x00000000, 0x69bbfbe7, 0xf7000073, 0x94e0f9e5, 0x9e81f696, 0xd1d1c6da, 0x292964d1, 0xca203325, 0xf7fef7f7, 0xf7f7f7f7, 0xf0f1f7f7, 0xe1e6d4e2, 0xf1f1f1e6, 0xf0f7f7f7, 0x267ba4d9, 0x161a1e26, 0xad111a11, 0xddba7a93, 0x39393996, 0x7979bd4c, 0x60606056, 0x60789160, 0xb7aeaa91, 0xfceed2ae,
0x00000000, 0x69e2f6e7, 0xfe000081, 0x70cdf9e4, 0xb271e096, 0xdadad6da, 0x24295bd1, 0xf7972525, 0xfefefefe, 0xf7f7fefe, 0xb5f1f7f7, 0xe1eae1e1, 0xf1f1d9b5, 0xf7f7f7f7, 0x5a87b8f0, 0x0d161f26, 0xce10161e, 0xdfe9443a, 0x3935396f, 0x79797f96, 0x60aaae75, 0x78aeaeaa, 0xd2b7b7b7, 0xeefdeed7,
0x00000000, 0x86fbf2e7, 0xeb000073, 0x58c2f9f9, 0xd371cca3, 0xdededed6, 0x252555d3, 0xf7f14132, 0xfefefef4, 0xf7f7fefe, 0xb582e6f5, 0xf4e1c4b5, 0xf1b5b5fe, 0xf7f7f7f7, 0x87a4d8f5, 0x1a171f32, 0xe5280816, 0xddf23720, 0x39353956, 0x565656ac, 0x9191ae91, 0xaeaeb7b7, 0xfdd7c3d2, 0xfcfcfffc,
0x00000000, 0xbdf9ded1, 0xe8000081, 0x3f9ef9f6, 0xe770becc, 0xdededed6, 0x32324dd3, 0xf7f7a732, 0xfefefefe, 0xf7f7fefe, 0x664f6ef0, 0x84848484, 0xca82b582, 0xf7f7f7f7, 0xa4b8f0f7, 0x1a1f177b, 0xc042080d, 0xdddf7f20, 0x4c35394c, 0x42505092, 0xaeae9178, 0xc3eed2ae, 0xfceec7d2, 0xffffffff,
0x00000000, 0xe8e9ded3, 0xe8000081, 0x588df9f6, 0xe79e9add, 0xdededed6, 0x1d324bc2, 0xf7f7f041, 0xfefefefe, 0xf7fefefe, 0x4f3623ab, 0x66666666, 0x82666666, 0xf7f7f7f1, 0xa4e2f7f7, 0x161f3287, 0xad721604, 0xcddfac2b, 0x7435394c, 0x444c484c, 0xaeaeaa5e, 0xeaeefcee, 0xfffcfcea, 0x3eaad2fc,
0x00000000, 0xf9dedebf, 0x000000a3, 0x8176f6f6, 0xe7c281ce, 0xdee7e7d6, 0x5a413bcc, 0xf7f7f5a7, 0xfefefef7, 0xebf7fefe, 0x36362330, 0x664f4f4f, 0x364f4f66, 0xf7f7f7ca, 0xb8f0f7f7, 0x0d1a5aa4, 0x72ad1108, 0xbddfce20, 0x632b394c, 0x60374848, 0xc3b7915e, 0xfcfdfdee, 0x6aedfffc, 0x00263c3e,
0x00000000, 0xe7dede00, 0x000000cd, 0x8d58e5f6, 0xf2d371bd, 0xf2f2f2d6, 0x597b74bd, 0xf7f7f5d8, 0xfef7f7f7, 0xbcf7f7fe, 0x2323180e, 0x4f4f4f36, 0x36363636, 0xf7f7f797, 0xe2f5f7f7, 0x081a87a4, 0x61ce0d16, 0xbddff920, 0x4839483b, 0x5e3c4848, 0xb7785e5e, 0xfffcfdee, 0x263eaaed, 0x001e2626,
0x00a00000, 0xdedede00, 0x000000e9, 0x8158bdf6, 0xdfe771bd, 0xf2f2f2d6, 0x687b68b4, 0xf7f5f0f0, 0xf7f7f7f7, 0x7bf7f7f7, 0x18304630, 0x36232323, 0x23363636, 0xf7f7f77b, 0xf0f7f7f7, 0x0d4187b8, 0x3be50d16, 0xa3e7f942, 0x3939354c, 0x22084848, 0x5e5e5e49, 0x3c8baa60, 0x2626323c, 0x00001d27,
0x00000000, 0xdededf00, 0x0000c0d1, 0x3f738df3, 0xdef381a3, 0xf3f3f3da, 0xa25a83a3, 0xf5f0f0f0, 0xf7f7f7f7, 0x7bf7f7f7, 0x0e0b5aab, 0x18181818, 0x18231818, 0xf7f7f75a, 0xf5f7f7f7, 0x177ba2d8, 0x2be5100d, 0x96f2f693, 0x3539354c, 0x00042b39, 0x1e000000, 0x27373a00, 0x27272732, 0x2a000027,
0x00000000, 0xdeded300, 0x0000cdd1, 0x347373c2, 0xdef38d96, 0xf3f6f6da, 0xb85a4196, 0xf0f0f0e2, 0xf7f7f7f5, 0x8ff7f7fe, 0x0e65edd7, 0x0e0e0e0e, 0x300e0e0e, 0xf7f7f75a, 0xf5f7f7f7, 0x4187a4f0, 0x20ce3a0d, 0x58f2f6b6, 0x3333394c, 0x00041339, 0x00000000, 0x371d4828, 0x27272727, 0x00280028,
0x00000000, 0xdedecd00, 0x0000bdde, 0x34795b9c, 0xdef39c8d, 0xf6f9f6e7, 0xd96e327f, 0xf0e2e2e2, 0xfef7f1f0, 0xcaf7fefe, 0xd2fcfc8f, 0x0b0b3065, 0x890b0b0b, 0xf7f7f78f, 0xf5f7f7f7, 0x5a87b8f0, 0x20c07216, 0x56f3f3e8, 0x3539354c, 0x000d0448, 0x28220000, 0x32271d48, 0x32272727, 0x00002800,
0x00000000, 0xdee70000, 0x00ceb1e7, 0x346f7369, 0xe7f3b281, 0xf6f9f9f2, 0xd96e466f, 0xf0e2e2e2, 0xfef7f1f0, 0xfefefefe, 0xfcd2abea, 0xc4eefcfc, 0xc4c4c4c4, 0xf7f7f7ca, 0xf7f7f7f7, 0x7b87d9f0, 0x40ac931a, 0x56e9f3f8, 0x39393539, 0x0000042b, 0x4c282a00, 0x2b371127, 0x0041272b, 0x0000002a,
0x00000000, 0xdee90000, 0x00cdc2de, 0x34568855, 0xf2f3c271, 0xe0f9f9f6, 0xd98a6e50, 0xe2e2d9d9, 0xfef7f1f1, 0xfefefefe, 0xb5e1fefe, 0xc4ababb5, 0x8f8fabc4, 0xf7f7f7f1, 0xf7f7f7f7, 0x83a2e2f1, 0x9393ad41, 0x56cdf3f9, 0x48393539, 0x00000913, 0x37502700, 0x3737270c, 0x2a2b372b, 0x00000000,
0x00000000, 0xdee800d2, 0xce9edede, 0x344b965b, 0xe9f3cc71, 0xdbf9f9f6, 0xd4a26e59, 0xe2e2d9d9, 0xfef7f1e6, 0xfdfefefe, 0xfefefdfd, 0xc4e6e1f4, 0xf7f1cac4, 0xf7f7f7f7, 0xf7f7f7f7, 0x7ca2f0f1, 0xba7fce68, 0x56cdf3f3, 0x4833393b, 0x00000d04, 0x0c2b4c2b, 0x3737371d, 0x002a3737, 0x00000000,
0x00000000, 0xdff500d2, 0xbd94dede, 0x3434a35b, 0xf3f3cc71, 0xccf9f9f6, 0xd4b36e59, 0xe2e2d9d9, 0xfef7f1f0, 0xfdfdfefe, 0xfefefdfd, 0xfefefefe, 0xfefefefe, 0xf7f7f7f7, 0xf7f7f7f7, 0x8ab8f0f7, 0xe47fe587, 0x56bdf2f2, 0x25393948, 0x32000d08, 0x110c324c, 0x37373737, 0x00002a37, 0x00000000,
0x00000000, 0xe5eb0000, 0x9eb1dede, 0x343fb476, 0xf3f3dd81, 0xbdf9f9f9, 0xd4d46e6e, 0xe2e2d9d4, 0xfef7f1f1, 0xfdfdfefe, 0xf4fefefd, 0xfefef4f1, 0xfefefefe, 0xf7f7f7f7, 0xf7f5f5f7, 0x8ab8f1f5, 0xf992e987, 0x56a3f2e7, 0x0c483939, 0x3928000c, 0x320c0c2b, 0x4c323748, 0x0000003a, 0x00000000,
0x00000000, 0xf500d200, 0x70dededf, 0x343fcd96, 0xf3f3dd81, 0xbbf9f9f9, 0xd4d46e6e, 0xd8d9d9d9, 0xfef1f1e6, 0xfdfdfefe, 0xe6fefefe, 0xfef1e6e6, 0xfefefefe, 0xf7f7f7f7, 0xf1f0f1f1, 0x7cd9f1f0, 0xf6a5e987, 0x567ff2e7, 0x08393339, 0x2525000c, 0x4c270c10, 0x424c3748, 0x00000000, 0x00000000,
0x00000000, 0xfa00ee00, 0x8de7dedf, 0x3f3fe98d, 0xf6e9e081, 0xa6f9f9f9, 0xcad4896e, 0x685a8bc4, 0xfef7f1a7, 0xfefdfefe, 0xd9cafefe, 0xf1e6e6d9, 0xfefefefe, 0xf0f7f7fe, 0xf0f0f0f0, 0x7cd9e6f0, 0xf2b6e98e, 0x5676f2de, 0x0c273339, 0x131d1d08, 0x4c4c130c, 0x00614c4c, 0x00000000, 0x00000000,
0x00000000, 0xfab70000, 0xb4dfdee4, 0x3f58e564, 0xf6e9e094, 0x92f9f9f8, 0xfccaa26e, 0x95999ffc, 0xfeca757a, 0xfefdfefe, 0xd9a2f1fe, 0xf1e6e6d4, 0xfefefefe, 0xf0f0f7f7, 0xe2e2e2f0, 0x8ad9e2e2, 0xe7e4cd8e, 0x4c56f3de, 0x0c0c3939, 0x13132410, 0x4c4c3713, 0x00005650, 0x00000000, 0x00000000,
0x00000000, 0xebd20000, 0xcddedef5, 0x3f56e571, 0xf6e9f694, 0x8ef9f9f8, 0xffd25a6e, 0xb99b6b9f, 0xa07ac5b9, 0xfefdfefe, 0xd4a2e6fe, 0xf1e6d9d9, 0xfefefefe, 0xf0f0f1f7, 0xe2e2e2f0, 0x8ad4d9e2, 0xdef6cda6, 0x4856e9de, 0x08093939, 0x13101d1d, 0x4c4b4c1d, 0x00003758, 0x00000000, 0x00000000,
0x00000000, 0xeb000000, 0xe9dcd3fa, 0x3f7fe081, 0xf9dff6a9, 0x6ff6f9f9, 0xedff3e6e, 0x6b9b9b52, 0xcbe36b2c, 0xfefefdc5, 0xd4b8f4fe, 0xf1e6d9d9, 0xfefefefe, 0xf0f1f1f7, 0xd9e2e2e2, 0x89d4d9d9, 0xd1f6e9bb, 0x3956cdde, 0x170c2739, 0x1313131d, 0x564b4c39, 0x00000042, 0x00000000, 0x00000000,
0x00000000, 0xfd000000, 0xddd6e4fa, 0x3f92e4a9, 0xf9dff6b1, 0x7ce0f9f9, 0x9fff915a, 0x192d9b52, 0xed992d2d, 0xfefeeeb9, 0xd9e1fdf4, 0xf1e6e6d9, 0xfefefef4, 0xe6f1f1f1, 0xd9d9e2e6, 0x89c1d4d9, 0xb1f3f6bd, 0x4856bdd3, 0x1d091348, 0x25131313, 0x56504c3b, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0x00000000, 0xe7d6fbfe, 0x58a3e4dd, 0xf9dff9c8, 0x68e0f9f8, 0x9fffaa48, 0x2f15805c, 0xe32c4e4e, 0xfefdcbb9, 0xe6f4fdf4, 0xf1e6e6d9, 0xfefefefe, 0xe6e6f1f7, 0xd4d9e2e6, 0x89b3d4d4, 0xa9e7f6ce, 0x3550a3d3, 0x1d0c0939, 0x3b131313, 0x20634c4c, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0x00000000, 0xe7ddfbfd, 0x58bde5e9, 0xf9dff9d1, 0x48cef8f8, 0xffffaa48, 0x1506317e, 0xb92d4e2f, 0xfefdcb80, 0xe6f4f4f1, 0xf1f1e6e6, 0xeafefefe, 0xa2b8caca, 0xd4d9d9a2, 0x7cb3d4d4, 0x9ed3f6ce, 0x394c7fd3, 0x131d0827, 0x4c2b1313, 0x0061564c, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0x00000000, 0xdef8fbea, 0x64ccf6f6, 0xf8ddf9de, 0x61bdf8f8, 0xfffcae56, 0x06067eff, 0x80192115, 0xfefdcb52, 0xe6e6f1f1, 0xf1f1e6e6, 0x9ffcfdfe, 0x91959599, 0xd9875d7a, 0x6eb3d4d4, 0x9ec2f3e5, 0x484c56bd, 0x131d0c1d, 0x4b3b1d13, 0x00936f4b, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0x00000000, 0xdff8fae8, 0x70cef9f6, 0xf8e7f9f2, 0x83b4f9f9, 0x9fd28b56, 0x0612ede3, 0x2f310606, 0xfefee331, 0xf1f1f1f1, 0xf7f1f1f1, 0x6bcbfffe, 0xb9b9b99b, 0xa195e3e3, 0x6ea8d4d4, 0x96b2f3f3, 0x483b50a3, 0x13131310, 0x4b4b3913, 0x006cc056, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0x00000000, 0xf6f8f7d5, 0x70e0f9e7, 0xf8e9f9f3, 0x59a3f9f9, 0x6aae6156, 0x06152e12, 0x31311506, 0xfefefc2e, 0xf1f1f1f7, 0xfef1f1f1, 0x807efffc, 0x2c2c6b9b, 0xcbb9e399, 0x7ca8d4e1, 0x96b2e9f3, 0x483b4c96, 0x13131326, 0x564b3b2b, 0x0000e561, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0x00000000, 0xf8f8f700, 0x81f9f9df, 0xf8e5f9f9, 0x68a3f9f9, 0x98ae5656, 0x1512123d, 0x31312115, 0xfefefd77, 0xf4f4f1f7, 0xfefef1f1, 0x805cedfc, 0x4e4e1980, 0x9fc5c52c, 0x8ea8d4ee, 0x869edff3, 0x59354a88, 0x13131332, 0x564b3b3b, 0x0000baac, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0x00000000, 0xf8fbe800, 0x94f9e9f9, 0xf8e5f9f9, 0x8396f9f9, 0xa15d7990, 0x06070175, 0x1c2e1212, 0xfefefed2, 0xfefefefe, 0xfdfefefe, 0x525cedff, 0x4e2f1519, 0x99b9994e, 0xa689b5ed, 0x869ed3f3, 0x594c4a6f, 0x2b131332, 0x4b50353b, 0x000000ce, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0x00000000, 0xf9fae800, 0x9ef9e7f8, 0xfbf6f9f9, 0x8e96f9f9, 0xa75667d9, 0x030a6aa1, 0xae6d0a02, 0xfefefed7, 0xfefefefe, 0xfefefefe, 0x7eedfffc, 0x2f210615, 0x7e806b2f, 0xbd8a8ffc, 0x7696b2f3, 0x874c4a56, 0x391d1332, 0x584c4a3b, 0x000000ce, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0x00000000, 0xf8fad500, 0xbdf6f6f9, 0xfaf9f9f9, 0x8e92f9f8, 0x966f63bc, 0x6aaac397, 0xa1aea16a, 0xeaebf7ca, 0xfefefefe, 0xfefefefe, 0xedffedfd, 0x1506062e, 0x7e525215, 0xe97c68ed, 0x6f88bedd, 0xa2594a50, 0x3b2b1313, 0x963b4c35, 0x000000ba, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0x00000000, 0xfaf80000, 0xccf6f9f9, 0xfef9f9f9, 0x836fe4f8, 0x88bed8e2, 0x75746192, 0x8f417497, 0xf7bccae6, 0xfefefef7, 0xfefefefe, 0x9f7e98d7, 0x0606062e, 0x9f312f21, 0xf36359b7, 0x5688a3dd, 0xa2594a4c, 0x3b3b250c, 0xc04c3b3b, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0x00000000, 0xf7e00000, 0xcef6f9f8, 0xfdf8f9f8, 0x836fe4f8, 0x92dbb4d5, 0x505663b6, 0xb38774a6, 0xe2e88ea2, 0xfefefefe, 0xfefefefe, 0x0705aac3, 0x21121506, 0xcb313131, 0xf3a66375, 0x506fa3cd, 0x7ba24d4c, 0x3b3b2b13, 0xc0564c3b, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0x00000000, 0xf8cc0000, 0xe5f9f9fb, 0xfdfbf9e8, 0x836fe4fb, 0xbbccccbb, 0xbccacabc, 0xa28a83a4, 0xf5d8a27c, 0xfefefefe, 0xeafefefe, 0x013daea0, 0x21150607, 0xaa1c2121, 0xf3bd5961, 0x5076a3c8, 0x1dd4634c, 0x3b3b2b25, 0x00964a3b, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0x00000000, 0xe9bd0000, 0xe4f8f8fb, 0xfdfaf9e8, 0x836fe0fb, 0xb6cdccb4, 0xb46f8eb8, 0x83a2a2a8, 0xfed5bb83, 0xfefefefe, 0xf7fefefe, 0x0faea190, 0x05020101, 0x75910a05, 0xe9e55063, 0x4b6f9ed1, 0x137b8e4d, 0x2b3b2b2b, 0x00ad4c4c, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0x00000000, 0xd1000000, 0xf8f8f8f9, 0xfdfef9f8, 0x8e8ecefb, 0xa6e0dab2, 0xb4b6f0f0, 0x8383a2b8, 0xfef0cc8e, 0xfefefefe, 0xe6fefefe, 0xaa8b63ca, 0x0a030a3d, 0x6397983d, 0xe7f37f59, 0x4d88b4df, 0x250c4153, 0x25322b25, 0x00ad564a, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0x00000000, 0xb2000000, 0xf8f8f9e9, 0xfdfef8f8, 0x9a8ebdfa, 0x8ecedfb1, 0xa6ccf0f0, 0x8e7cb8f5, 0xf7f7bbbe, 0xf7fef7f7, 0xe2fefefe, 0x7563a7b8, 0xaaaaaaae, 0x636f75ae, 0xe7f3a683, 0x4d88b4f2, 0x2b13084b, 0x39252525, 0x0000924a, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0x00000000, 0xa3000000, 0xf8f9f6d1, 0xfdfdf8f5, 0xa38ebdfa, 0x8eb6e7b1, 0xbcc8d5f1, 0x9a6f83d8, 0xf7f7d5be, 0xf7f7f7f7, 0xf1fef7f7, 0x6fa7b5c1, 0x8b757568, 0x68689068, 0xe7f3bd83, 0x4d88b2f2, 0x24250c2b, 0x4c2b2525, 0x0000a54a, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0x00000000, 0x00000000, 0xf8f8e9b1, 0xfdfdfbfb, 0x868eb4fa, 0x8ea7e9a9, 0xd8a3b4f0, 0xb09a8ea7, 0xf7f7f5be, 0xf7f7f7f7, 0xf4f7f7f7, 0xa7a2a4c9, 0x6f6f6f83, 0x7c68906f, 0xf2e7ce83, 0x4d539ef2, 0x2525131d, 0x4a2b2525, 0x00000056, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0x00000000, 0x00000000, 0xf8f6dea9, 0xfdfcfbfb, 0x6996b4fa, 0x878ee08d, 0xd892c2e2, 0xb29c9a83, 0xf7f7f5bb, 0xf7f7f7f7, 0xfef7f7f7, 0xa4c1c1f1, 0xd8bcb8a2, 0x7c68a7d8, 0xf2e7f36f, 0x4a4a58c2, 0x2525250c, 0x4d3b251d, 0x000000a5, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0x00000000, 0x00000000, 0xf8e7b19a, 0xfdfdfafb, 0x719aacfa, 0x8783d08d, 0x878eb2bc, 0xb2a99ea6, 0xf7f7f7d8, 0xf7f7f7f7, 0xf7f1f7f7, 0xe2f0f5f7, 0xf0e2e2e2, 0x7c8aa4d8, 0xd3e7f392, 0x354d4d58, 0x202b2517, 0x4b4c2525, 0x000000ad, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0x00000000, 0x00000000, 0xf9d19d00, 0xfdfdfafa, 0x7396bdfb, 0x8783b69a, 0x8e87a3a3, 0xb2b0b0b2, 0xf5f5f7f0, 0xf5f5f5f5, 0xf1e6e6f5, 0xf5f5f7f7, 0xf0f0f0f5, 0x7c8a8ed8, 0x79d3e7bb, 0x254d4d4d, 0x252b2920, 0x924a2925, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0x00000000, 0x00000000, 0xddb19400, 0xfdfdfefa, 0x819ecdfb, 0x8383a5a3, 0xbe838896, 0xccc8b1c8, 0xf5f5f5f5, 0xf5f5f5f5, 0xe2e2e2f0, 0xf5f5f5f5, 0xf0f0f0f5, 0x7c8a8eb8, 0x4d58c2ce, 0x2b4d4a4d, 0x252b2b4a, 0xac4a3525, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0x00000000, 0x00000000, 0xd19d9a00, 0xfdfdfefb, 0x8db2e9fb, 0x6f83a6bd, 0xc28e6f86, 0xd5dfd1de, 0xf5f5f5f5, 0xf5f5f5f5, 0xd9d9d9f0, 0xf5f5f0e2, 0xd5f0f0f0, 0xa78e8eb8, 0x4d4d79dd, 0x4d564d4a, 0x252b2029, 0x007f4a2b, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0x00000000, 0x00000000, 0xb19d0000, 0xfdfdfee8, 0x94c2f3f8, 0x6f637fcc, 0xd3a36855, 0xf0d3f2f2, 0xf0f0f5f0, 0xf0f5f5f0, 0xd4c9c9f0, 0xf0f0f0d9, 0xd5f0f0f0, 0x92b88ea4, 0x4d4d4d96, 0x294d584d, 0x252b2b20, 0x007f4a39, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0x00000000, 0x00000000, 0xb08d0000, 0xfdfefedb, 0x94d3f6f8, 0x566f8ece, 0xe9bd503f, 0xf0bdf6f6, 0xf0f0f0f0, 0xf0f0f0f0, 0xc9c9c9e2, 0xf1f0d9d4, 0xbbf5f0f0, 0x7f87d5a6, 0x3f4d4d58, 0x20354d73, 0x2b253b25, 0x0000503b, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0x00000000, 0x00000000, 0x9d880000, 0xfefafec8, 0xa9e7f6f8, 0x558e8edd, 0xf6cd9253, 0xf0bbf3f9, 0xf0f0f0f0, 0xf5f0f0f0, 0xc9c9c1e2, 0xf0d9d4c9, 0xbbd5f0f0, 0x4b688eb6, 0x73554d4d, 0x2b25354a, 0x33254c2b, 0x00007f4c, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0x00000000, 0x00000000, 0x9d000000, 0xfefafeb1, 0xb1e9f8f9, 0x5b8e88dd, 0xf6e9a353, 0xf0bce9f9, 0xf0f0f0f0, 0xf0f0f0f0, 0xc9c1d4f0, 0xd9c9c9c9, 0xccccf0e2, 0x4d5a7b6f, 0x4a73584a, 0x2b2b2b4a, 0x3b2b294c, 0x0000007f, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0x00000000, 0x00000000, 0x8d000000, 0xfafbf7a9, 0xdef3fbf9, 0x698e9add, 0xf6f3b46f, 0xf0bccdf6, 0xf0f0f0f0, 0xf0f0f0f0, 0xc9b3d4f0, 0xc9c9c9c9, 0x8eccbbd9, 0x4d397b7c, 0x4d43733f, 0x2b2b2b29, 0x563b204b, 0x000000ad, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0x00000000, 0x00000000, 0x00000000, 0xfbfbf19d, 0xf2f6fef9, 0x718888e9, 0xf9f6c288, 0xf5b8bdf6, 0xf0f0f0f0, 0xf0f0f0f0, 0xc1a4d8f5, 0xc9c9c1c1, 0x67a3e0c9, 0x34354159, 0x354a4d73, 0x3b252b39, 0xa34b253b, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0x00000000, 0x00000000, 0x00000000, 0xf8f8d89c, 0xf2f9fef9, 0x948688e5, 0xf6f9cd86, 0xf0a496f6, 0xf0f0f0f0, 0xf0f0f0f0, 0x83bcf0f0, 0xc1a4a4a4, 0x5367ddcc, 0x73342589, 0x2b4a4a53, 0x4a252b35, 0xad584c29, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0x00000000, 0x00000000, 0x00000000, 0xf9f8b496, 0xf9f8fef9, 0x9d8d9ae4, 0xf6f6e58d, 0xf08e76f6, 0xf0f0f0f5, 0xf0f0f0f0, 0xf0f0f0f0, 0xa3a7bcbc, 0x5367cdce, 0x5573205a, 0x4d2b4a4d, 0x3439253b, 0x00bd583b, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0x00000000, 0x00000000, 0x00000000, 0xf3f8a300, 0xf8fbfef6, 0xa99d9ae4, 0xf9f9f69e, 0xd8876fe0, 0xf0f0f0f0, 0xf0f0f0f0, 0xf0f0f0f0, 0xe8f0f0f0, 0x6353a3e9, 0x4a585825, 0x4d4d354a, 0x3b4b2725, 0x00ad9658, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0x00000000, 0x00000000, 0x00000000, 0xdfe89600, 0xfbfafee7, 0xd1b08df5, 0xf9f9f9b2, 0xb8836fce, 0xf0f0f0f0, 0xf0f0f0f0, 0xf0f0f0f0, 0xe0e8f0f0, 0x395376e9, 0x4d4d793b, 0x294d4b4d, 0x58354a1d, 0x0000ce58, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0x00000000, 0x00000000, 0x00000000, 0xd1d50000, 0xfefefec8, 0xe9d18df7, 0xf9f9f9b4, 0xa2838ebd, 0xf0f0f0f0, 0xf0f0f0f0, 0xf0f0f0f0, 0xddf6f0f0, 0x254b58cd, 0x4d4a5379, 0x1d354d4d, 0x554b4a3b, 0x000000bd, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0x00000000, 0x00000000, 0x00000000, 0xddbc0000, 0xfdf7fbb1, 0xe5e594fe, 0xf9f9f9cc, 0x87878ea3, 0xf0f0f0e2, 0xf0f0f0f0, 0xf5f0f0f0, 0xf3cce0e2, 0x563b55a3, 0x534d4a53, 0x391d354d, 0x96584d4a, 0x000000c0, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0x00000000, 0x00000000, 0x00000000, 0xce000000, 0xfefbf8b1, 0xe5f99efe, 0xf9f9f9dd, 0x87878eb4, 0xf0f0f0b8, 0xf0f0f0f0, 0xd9f0f0f0, 0xf6cce0be, 0x694d559c, 0x4d564a53, 0x4a2b2034, 0xdf76584d, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0x00000000, 0x00000000, 0x00000000, 0xe8000000, 0xfbf8e0b1, 0xe4f5a9fd, 0xf8f9f9e0, 0x87878ecc, 0xf0f0f087, 0xf0f0f0f0, 0xa4c1d8f0, 0xddf6b4cc, 0x535b699c, 0x344b534a, 0x4d4d2b2b, 0x00cd5858, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0xf8f8d1b2, 0xf5f7b2fd, 0xfaf8f8f5, 0x878fa4dd, 0xf5f0b887, 0xd9e2f0f0, 0xa4a8a4c1, 0xdff9bdd5, 0x53696994, 0x29345853, 0x584d4b3b, 0x00bfd358, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0xf8f9b1bd, 0xf5ebccfa, 0xf7faf7eb, 0x878fa6e5, 0xe2d88787, 0x8ea4b8d9, 0xb48e8a8e, 0xd1f9e8bb, 0x5585718c, 0x4d344d55, 0x73584d4d, 0x0000cdc0, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0xf8e4b100, 0xf7e8e0f8, 0xfafdfdf7, 0x8f8796e5, 0x8787878f, 0x83878787, 0xd58e8e83, 0xc8e5fba7, 0x69868d9d, 0x4d3f344d, 0xbf86584b, 0x000000e9, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0xf9c8be00, 0xfee4e5f9, 0xf7fcfdf5, 0x596896f9, 0x83838783, 0x83838383, 0xa6a66868, 0xa9e9faeb, 0x559d9c94, 0x4d4d555b, 0xf2bf9658, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0xe9a9b600, 0xfee0f9f3, 0xfafdfaf5, 0x8368a3f9, 0x83838383, 0x83838387, 0xa0cc796f, 0x9dd3e8fe, 0x718cb094, 0x4d4d4d69, 0x00e9bf96, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0xc8b40000, 0xfbe9f6e7, 0xfafefbf9, 0x7c63acf6, 0x837b8368, 0x687b7b68, 0xeea59268, 0x9db1ddfe, 0x818d94a9, 0x764d4d4d, 0x00acf3c0, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0xc8000000, 0xf8d3d1d1, 0xf7fbf8f6, 0x6369bde5, 0x7b596f59, 0x415a6883, 0xfda0b641, 0xb1a9dffb, 0x53738da9, 0xc0584d4d, 0x0000ade5, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0xce000000, 0xf6c8b1a9, 0xf7f8f9f6, 0x3958bde5, 0x686f6767, 0x67676763, 0xf7eaac67, 0xa9b2b1e9, 0x4d556994, 0xe5c0564d, 0x000000ad, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0xddb1a9cd, 0xf8f9f9f6, 0x4d3bbddd, 0x59a37653, 0x534a5056, 0xf8fdb6a3, 0x94a9c8df, 0x4d535581, 0xade5c058, 0x00000000, 0x00000000, 0x00a00000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0xd1a9b4b6, 0xf6f9f9f9, 0x557692d3, 0x39ceb24d, 0x4c354d53, 0xe9fbfe9e, 0x8d94b1d1, 0x6f535373, 0x0000ced0, 0x00000000, 0x00000000, 0x0000a500, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0xc6b2d500, 0xf6f9f9f8, 0x349676c2, 0x61e5dd73, 0x76354d35, 0xdff9fad5, 0x7394b0b0, 0xd06f5555, 0x000000ba, 0x00000000, 0x96000000, 0x00000096, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0xc8d00000, 0xf3f9f9e4, 0x349658b1, 0xade5f6a9, 0xc25b4d34, 0xb1ddf9f8, 0x698c9db0, 0x6cd09673, 0x00000000, 0x00000000, 0x868d0000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0xd0000000, 0xf2f2f8e5, 0x8d865ba9, 0xd0e4e4cd, 0xf68d3f4b, 0xa9a9d1f9, 0x8d949db0, 0x0000babd, 0x00000000, 0x00000000, 0x00969486, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0xc8dff0ba, 0xa37355a9, 0xcf000000, 0xdddf58b6, 0xa9b09db1, 0xe5cca99d, 0x00000000, 0x00000000, 0x94000000, 0x00007694, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0xb1e8cf00, 0x008c8cb0, 0x00000000, 0xa9b0c2fb, 0xb0b0a99d, 0x00bae5c8, 0x00000000, 0x00000000, 0x8d947600, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0xe0d70000, 0x008da9b0, 0x00000000, 0x8dbdf7d7, 0x76817181, 0x00000000, 0x00000000, 0x814b0000, 0x00889494, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0xc7000000, 0x0000bdcc, 0x00000000, 0xd0e8d200, 0x0072a6ce, 0x00000000, 0x00000000, 0x9d9d8d76, 0x00000081, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x000000c7, 0x00000000, 0x00000000, 0xe4e8cf00, 0xe8cfd0e4, 0xa9b1e9e8, 0x5870949d, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0xb6000000, 0xf8e8e8d5, 0x9ecdf9f8, 0x00000081, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000
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

void Picture(out vec4 fragColor, in vec2 fragCoord) {
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
	
	vec2 aspect = (resolution.xy / (resolution.y-0.4)) + 0.2;
	uv = ( fragCoord.xy / (resolution.y-0.4) ) - aspect / 2.0;
	float _d =  1.0-length(uv);
	uv *= 18.0 ;
	uv -= vec2(-7., 1.);

	vec3 ch_color = hsv2rgb_smooth(vec3(iTime*0.4+uv.y*0.1,0.5,0.5));
    uv.y -=7.0;
	uv.x += 0.7+sin(iTime+uv.y*0.7)*0.5;
	ch_pos = ch_start;
	
	
	nl1
	nl2
	_ _ L I G H T A M P nl3


	vec3 color = mix(ch_color, vec3(0,0,0), 1.0- (0.09 / d*2.0));  // shading
	fragColor = vec4(color, 1.0);
}


//fadeStartTime and FadeTime absolute time
void Blend(out vec4 fragColor, in vec4 inColor, in int fadeStartTime, in int fadeTime) 
{
    float col = (time-fadeStartTime) / fadeTime;
    if (col < 0.0)
        fragColor = vec4(0.0, 0.0, 0.0, 0.0);
    else
    {
        inColor *= vec4(col,col,col,col);
        fragColor = inColor;
    }
}

void main(void)
{
    vec4 main_FragColor;
    vec4 text_FragColor;
	vec4 scroltext_FragColor;
    vec4 picture_FragColor;

    iTime = time /1000;
    iResolution = vec3(resolution, 0.0);
    mainImage(color, gl_FragCoord.xy);

    //Fade the main image in
    if (time < 7600)
        Blend(color, color, 1, 7600);


	//Text
    mainText(text_FragColor, gl_FragCoord.xy);
    if (time < 22000)
        Blend(text_FragColor, text_FragColor, 7600, 22000);

    Picture(picture_FragColor, gl_FragCoord.xy);

	//Compose
	color += text_FragColor;
	color += scroltext_FragColor;
    color += picture_FragColor;
}
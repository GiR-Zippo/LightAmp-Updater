// Tune: summer_memories_3.xm
#version 330

/*
 * Original shader from: https://www.shadertoy.com/view/lt2Bz3
 */

#ifdef GL_ES
precision mediump float;
#endif

// glslsandbox uniforms
uniform float time;
uniform vec2 resolution;

// shadertoy emulation
#define iTime time/1000
#define iResolution resolution
const vec4 iMouse = vec4(0.);

// --------[ Original ShaderToy begins here ]---------- //
// ⊂(◉‿◉)つ Kirby !
// Leon 2018/01/24
// Using code from IQ, Mercury, LJ, Duke, Koltes

#define STEPS 30.
#define VOLUME .01
#define FAR 10.
#define PI 3.14159
#define TAU 2.*PI

const vec3 pink = vec3(0.917,0.482,0.663);
const vec3 red = vec3(0.825,0.142,0.111);
const vec3 beige = vec3(0.905, 0.670, 0.235);
const vec3 blue = vec3(0.058, 0.074, 0.560);
const vec3 blueSky = vec3(0.741, 0.941, 1);
const vec3 green1 = vec3(0.298,0.830,0.153);
const vec3 green2 = vec3(0.038,0.260,0.047);
const vec3 gold = vec3(1, 0.858, 0.058);

// sdf toolbox
float rng (vec2 seed) { return fract(sin(dot(seed*.1684,vec2(54.649,321.547)))*450315.); }
mat2 rot (float a) { float c=cos(a),s=sin(a); return mat2(c,-s,s,c); }
float sdSphere (vec3 p, float r) { return length(p)-r; }
float sdCylinder (vec2 p, float r) { return length(p)-r; }
float sdIso(vec3 p, float r) { return max(0.,dot(p,normalize(sign(p))))-r; }
float sdBox( vec3 p, vec3 b ) { vec3 d = abs(p) - b; return min(max(d.x,max(d.y,d.z)),0.0) + length(max(d,0.0)); }
float sdTorus( vec3 p, vec2 t ) { vec2 q = vec2(length(p.xz)-t.x,p.y); return length(q)-t.y; }
float amod (inout vec2 p, float count) { float an = TAU/count; float a = atan(p.y,p.x)+an/2.; float c = floor(a/an); c = mix(c,abs(c),step(count*.5,abs(c))); a = mod(a,an)-an/2.; p.xy = vec2(cos(a),sin(a))*length(p); return c; }
float repeat (float v, float c) { return mod(v,c)-c/2.; }
vec2 repeat (vec2 v, float c) { return mod(v,c)-c/2.; }
vec3 repeat (vec3 v, float c) { return mod(v,c)-c/2.; }
float smoo (float a, float b, float r) { return clamp(.5+.5*(b-a)/r, 0., 1.); }
float smin (float a, float b, float r) { float h = smoo(a,b,r); return mix(b,a,h)-r*h*(1.-h); }
float smax (float a, float b, float r) { float h = smoo(a,b,r); return mix(a,b,h)+r*h*(1.-h); }
vec2 toroidal (vec2 p, float r) { return vec2(length(p.xy)-r, atan(p.y,p.x)); }
vec3 lookAt (vec3 eye, vec3 target, vec2 uv) {
  vec3 forward = normalize(target-eye);
  vec3 right = normalize(cross(vec3(0,1,0), forward));
  vec3 up = normalize(cross(forward, right));
  return normalize(forward + uv.x * right + uv.y * up);
}

struct Shape {
    float dist;
    vec3 color;
    float spec;
    float glow;
};
Shape newShape () { Shape shape; shape.dist = 1000.; shape.color = vec3(1.); shape.spec = 0.; shape.glow = 0.; return shape; }
Shape add (Shape a, Shape b) { Shape c = newShape(); c.dist = min(a.dist, b.dist); float op = step(b.dist, a.dist); c.color = mix(a.color, b.color, op); c.spec = mix(a.spec, b.spec, op); c.glow = mix(a.glow, b.glow, op); return c; }
Shape map (vec3 p);

vec3 getNormal (vec3 p) { vec2 e = vec2(.01,0); return normalize(vec3(map(p+e.xyy).dist-map(p-e.xyy).dist,map(p+e.yxy).dist-map(p-e.yxy).dist,map(p+e.yyx).dist-map(p-e.yyx).dist)); }
float getShadow (vec3 pos, vec3 at, float k) {
    vec3 dir = normalize(at - pos);
    float maxt = length(at - pos);
    float f = 1.;
    float t = VOLUME*10.;
    for (float i = 0.; i <= 1.; i += 1./STEPS) {
        float dist = map(pos + dir * t).dist;
        if (dist < VOLUME) return 0.;
        f = min(f, k * dist / t);
        t += dist;
        if (t >= maxt) break;
    }
    return f;
}

Shape sdKirby (vec3 pos) {
    Shape kirby = newShape();
    vec3 p;

    // foot
    p = pos;
    p.x = abs(p.x)-.4;
    p.y += 1.;
    p.z += .35;
    p.z *= .65;
    float foot = sdSphere(p, .4);
    foot = smax(foot, -p.y, .2);

    // breath animation
    float wave = .5+.5*sin(iTime*5.);
    pos.y += 1.;
    pos.y *= 1.+.1*wave;
    pos.y -= 1.;
    pos.xz *= 1.-.1*wave;

    // body
    p = pos;
    float body = sdSphere(p, 1.);

    // hand
    p = pos;
    p.x = abs(p.x)-1.;
    p.xy *= rot(PI/3.);
    p.x *= .75;
    p.y *= 1.5;
    float hand = sdSphere(p, .4);

    // body compo
    kirby.dist = min(min(body, foot), hand);
    kirby.color = mix(pink, red, step(foot, body));
    // kirby.spec = 1.;
    kirby.glow = 1.;

    // eyes
    p = pos;
    p.y -= .3;
    p.y *= 1./smoothstep(.0,.1,abs(sin(iTime)));
    p.x = abs(p.x) - .2;
    p.x *= 2.;
    p.y *= .75;
    kirby.color = mix(kirby.color, vec3(0), step(length(p.xy), .2));
    p.y -= .1;
    p.y *= 1.6;
    kirby.color = mix(kirby.color, vec3(1), step(length(p.xy), .1));
    p.y += .2;
    p.x *= .75;
    p.y *= .4;
    kirby.color = mix(kirby.color, blue, clamp(-p.y*10.,0.,1.)*step(length(p.xy), .1));

    // open mouth
    p = pos;
    p.x *= .5;
    p.y += .2-abs(p.x*.2);
    float d = length(p.xy);
    float mouth = step(d, .1);
    kirby.color = mix(kirby.color, red, mouth);
    kirby.color = mix(kirby.color, red*.1, mouth*(1.-clamp(-p.y*10.+.5,0.,1.)));

    // cheeks
    p = pos;
    p.x = abs(p.x) - .5;
    p *= 6.;
    p.x *= .75;
    kirby.color = mix(kirby.color, red, .75*(1.-smoothstep(0.5,1.,length(p.xy))));

    return kirby;
}

Shape sdGround (vec3 pos) {
    Shape ground = newShape();
    vec3 p;
    p = pos;
    p.y += 1.;
    float cell = .5;
    float height = .2;
    float padding = .45;
    p.xz = repeat(p.xz, cell);
    p.y += 1. + height;
    ground.dist = smin(p.y, sdBox(p, vec3(cell*padding, height, cell*padding)), .2);
    p.y -= 1.;
    ground.dist = smin(ground.dist, max(sdBox(pos,vec3(1,3,1)),sdBox(p, vec3(cell*padding, height, cell*padding))), .2);
    ground.color = beige;
    return ground;
}

Shape sdPlant (vec3 pos) {
    Shape plant = newShape();
    plant.spec = .5;
    plant.glow = 1.;
    float radius = 2.;
    pos.y += 1.;
    pos.xyz = pos.zxy;
    vec3 p = pos;
    p.xy = toroidal(p.xy, radius);
    p.y *= 2.;
    p.xz *= rot(p.y * 2.+sin(p.y+iTime));
    float id = amod(p.xz,2.);
    p.x -= .2;
    p.xz *= rot(-p.y+iTime+sin(p.y-iTime*2.)*5.);
    id += amod(p.xz, 4.);
    p.x -= .1;
    plant.dist = sdCylinder(p.xz, .04);
    plant.color = mix(green1, green2, mod(id,2.));
    return plant;
}


Shape sdStar (vec3 pos) {
    Shape star = newShape();
    star.spec = 1.;
    star.glow = 1.;
    float radius = 5.;
    float size = .2;
    vec3 p = pos;
    p.y -= radius;
    p.xy *= rot(-iTime*.5);
    float index = amod(p.xy, 16.);
    p.x -= radius-1.5;
    p.xy *= rot(iTime+index);
    amod(p.xy, 5.);
    star.dist = sdIso(p, size);
    star.color = gold;
    return star;
}

Shape map (vec3 pos) {
    Shape scene = newShape();
    vec3 p = pos;
    scene = add(scene, sdKirby(p));
    scene = add(scene, sdGround(p));
    scene = add(scene, sdPlant(p));
    scene = add(scene, sdStar(p));
    return scene;
}


vec3 camera (vec3 p) {
    float click = clamp(iMouse.z, 0., 1.);
    p.yz *= rot(click*(-.25*PI*(iMouse.y/iResolution.y-.5)));
    p.xz *= rot(click*(-.5*PI*(iMouse.x/iResolution.x-.5)));
    p.xz *= rot((1.-click)*(.5*PI*(.4*sin(iTime*.1))));
    return p;
}

vec3 raymarch (vec2 uv) {
    vec3 eye = camera(vec3(0,1,-5.5));
    vec3 ray = lookAt(eye, vec3(0), uv);
    float shade = 0., dist = 0.;
    vec3 pos = eye;
    Shape shape;
    for (float i = 0.; i <= 1.; i += 1./STEPS) {
        shape = map(pos);
        if (shape.dist < VOLUME || dist > FAR) { shade = 1.-i; break; }
        shape.dist *= .9 + .1 * rng(uv+fract(iTime));
        dist += shape.dist;
        pos = eye + ray * dist;
    }
	   vec3 color = shape.color;
    vec3 normal = getNormal(pos);
    vec3 view = normalize(eye-pos);
    vec3 lightPos = vec3(1.,0,-2.);
    lightPos.xy *= rot(iTime*.5);
    lightPos.y += 2.;
    vec3 lightDir = normalize(lightPos-pos);
    float light = clamp(dot(lightDir, normal),0.,1.);
    float ambient = .8;
    color = mix(shape.color * ambient, shape.color, light);
    color = mix(color, vec3(1), shape.spec*clamp(pow(light,16.),0.,1.));
    color *= .3+.7*getShadow(pos, lightPos, 64.)*clamp(dot(lightDir,normal),0.,1.);
    color = mix(color, shape.color, (1.-abs(dot(normal, view))) * shape.glow);
    color *= shade;
    float far = 1.-smoothstep(FAR/2.,FAR,dist);
    float y = uv.y + sin(uv.x*3.-iTime+sin(uv.x*6.+iTime))*.05;
    color = mix(beige*abs(y), color, far);
    color = pow(color, vec3(1./2.));
    color = mix(color, mix(blueSky, blue, y), (1.-far)*step(0.,y));
    return color;
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
	vec2 uv = (fragCoord.xy-.5*iResolution.xy)/iResolution.y;
	fragColor = vec4(raymarch(uv),1.0);
}


// --------[ Intro Text Begins here ]---------- //

#define SC 0.07 // line scale
#define LW 0.03 // line width
#define C1 (vec2 p){return min(l(p,
#define C2 ),min(l(p,

// line segment
float l(vec2 p, float ax, float ay, float bx, float by)
{
    p = floor(p*500.+0.5)/500.;

    vec2 a = vec2(ax,ay)*SC;
    vec2 b = vec2(bx,by)*SC;
    vec2 ab = b-a;return length(p-a-ab*clamp(dot(p-a,ab)/dot(ab,ab),0.0,1.0))-LW;
}

float A_I C1 1.,-8.,1.,-1.5 C2 1.,-1.5,5.,-1.5 C2 5.,-1.5,5.,-5. C2 5.,-5.,1.,-5. C2 1.,-5.,5.,-5.),l(p,5.,-5.,5.,-8.))))));}
float G_I C1 5.,-2.5,5.,-1.5 C2 5.,-1.5,1.,-1.5 C2 1.,-1.5,1.,-8. C2 1.,-8.,5.,-8. C2 5.,-8.,5.,-5.),l(p,5.,-5.,3.5,-5.))))));}
float H_I C1 1.,-1.5,1.,-8. C2 1.,-8.,1.,-5. C2 1.,-5.,5.,-5. C2 5.,-5.,5.,-1.5),l(p,5.,-1.5,5.,-8.)))));}
float I_I C1 1.5,-1.5,4.5,-1.5 C2 4.5,-1.5,3.,-1.5 C2 3.,-1.5,3.,-8. C2 3.,-8.,1.5,-8.),l(p,1.5,-8.,4.5,-8.)))));}
float L_I C1 1.,-1.5,1.,-8.),l(p,1.,-8.,5.,-8.));}
float M_I C1 1.,-8.,1.,-1.5 C2 1.,-1.5,3.,-4. C2 3.,-4.,5.,-1.5),l(p,5.,-1.5,5.,-8.))));}
float P_I C1 1.,-8.,1.,-1.5 C2 1.,-1.5,5.,-1.5 C2 5.,-1.5,5.,-5.),l(p,5.,-5.,1.,-5.))));}
float T_I C1 3.,-8.,3.,-1.5 C2 3.,-1.5,1.,-1.5),l(p,1.,-1.5,5.,-1.5)));}

vec4 mainImage2(in vec3 fragColor, in vec2 fragCoord )
{   
    float t = iTime;  
    vec2 uv = (2.0 * fragCoord.xy - iResolution.xy) / iResolution.yy;

    uv.x += 0.13;
    uv.x *= abs(sin(uv.x+t*2.0)*0.5+1.0)+1.0;
    uv.y *= abs(sin(uv.x+t*2.0)+1.0)+1.0;

    vec3 col = fragColor;
    float d = 1.0;

    d = min(d,L_I(uv-vec2(-2.0,0.5)));
    d = min(d,I_I(uv-vec2(-1.5,0.5)));
    d = min(d,G_I(uv-vec2(-1.0,0.5)));
    d = min(d,H_I(uv-vec2(-0.5,0.5)));
    d = min(d,T_I(uv-vec2( 0.0,0.5)));

    d = min(d,A_I(uv-vec2( 1.0,0.5)));
    d = min(d,M_I(uv-vec2( 1.5,0.5)));
    d = min(d,P_I(uv-vec2( 2.0,0.5)));

    col = mix(vec3(0),col,smoothstep(d,d+0.01,0.0));

    return -vec4(vec3(0.-40.1*d), -1.0);
}


///////////////////
// Text
///////////////////

float horizon;

float shift_right(float v, float amt) {
	v = floor(v) + 0.5;
	return floor(v * exp2(-amt));
}

float shift_left(float v, float amt) {
	return floor(v * exp2(amt) + 0.5);
}

float mask_last(float v, float bits) {
	return mod(v, shift_left(1.0, bits));
}
float extract_bits(float v, float from, float n) {
	return mask_last(shift_right(v, from), n);
}

float tex(float val, float j) {
	float sign = val > 0.0 ? 0.0 : 1.0;
	val = abs(val);
	float exponent = floor(log2(val));
	float biased_exponent = exponent + 63.0;
	float mantissa = ((val * exp2(-exponent)) - 1.0) * 8388608.0;

	float index = j * 2.0;
	if (j == 0.0) {
		return sign + mod(biased_exponent, 2.0) * 2.0;
	} else if (j < 4.){
		index -= 1.;
		return extract_bits(biased_exponent, index, 2.0);
	} else {
		index -= 7.;
		return extract_bits(mantissa, index, 2.0);
	}
}

float ch(mat4 m, float c, vec2 uv) {
	if (max(abs(uv.x), abs(uv.y)) < .5) {
		float i = floor(uv.x * 16.) + 8.;
		float j = floor(uv.y * 16.) + 8.;
		vec4 v = m[3];
		if (i < 8.) {
			v = (i < 4.) ? m[0] : m[1];
		} else if (i < 12.) {
			v = m[2];
		}
		i = mod(i, 4.);
		if (i < 2.) {
			if (i==0.) return tex(v.x, j);
			return tex(v.y, j);
		} else {
			if (i==2.) return tex(v.z, j);
			return tex(v.w, j);
		}
	}
	return c;
}

const mat4 A = mat4(1.19218981354e-07, 2050.70800781, -527189.25, -537941.25, -578901.25, -2261.33300781, -6.75373598824e-07, 6.74907482789e-07, 1.53595161289e-19, 1.31626876509e-07, 2261.33300781, -537941.25, -527189.25, -2050.70800781, -4.76875925415e-07, 4.76839545627e-07);
const mat4 B = mat4(1.21692806943e-07, 2229.33300781, -709973.25, -742741.25, -742741.25, -742741.25, -740050.75, -740180.75, -742741.25, -742741.25, -742741.25, -578869.25, -2261.20800781, -4.8910192163e-07, 4.79313371216e-07, 1.19354950812e-07);
const mat4 C = mat4(1.08984897264e-19, 4.75968708891e-10, 1.29764259782e-07, 2261.33300781, 2773.33300781, -742741.25, -742709.25,-742573.25, -742573.25, -742701.25, -2773.20800781, -2261.20800781, 5.18911519976e-07, 4.86773501507e-07, 1.19830161793e-07, 4.65670724203e-10);
const mat4 D = mat4(1.29143387539e-07, 2773.33300781, -742741.25, -742741.25, -742741.25, -742741.25, -740011.25, -742061.25, -742741.25, -2773.33300781, -2261.33300781, 5.19057039128e-07, 4.87391957904e-07, 1.19868815318e-07, 4.66267580101e-10, 1.08429006043e-19);
const mat4 E = mat4(4.75362527119e-10, 557.333251953, 709973.25, -742741.25, -742741.25, -742741.25, -742741.25, -742227.25, -742227.25, -742227.25, -742227.25, -742091.25, -709291.25, -2218.66699219, -4.86616613671e-07, 4.7683727189e-07);
const mat4 F = mat4(1.29143387539e-07, 2773.33300781, -742741.25, -742741.25, -742741.25, -2901.33300781, -6.75051182952e-07, 6.75051182952e-07, 1.53488842731e-19, 1.53488842731e-19, 1.5346071859e-19, 1.53453997747e-19, 1.53446966712e-19, 1.19140477987e-19, 1.10670148515e-19, 1.08949612841e-19);
const mat4 G = mat4(1.08984871414e-19, 4.75968708891e-10, 1.29764259782e-07, 2261.33300781, 2773.33300781, -742741.25, -742573.25,-742059.25, -742099.25, -742101.25, -2898.83300781, -2898.83300781, 5.25925543116e-07, 4.86780777464e-07, 1.19830161793e-07, 4.65663618776e-10);
const mat4 H = mat4(1.29143387539e-07, 2773.33300781, -742741.25, -742741.25, -742741.25, -709973.25, -2219.33300781, -5.1671906931e-07, 6.45707245894e-07, 2901.33300781, -742741.25, -742741.25, -742741.25, -2773.33300781, -5.16573550158e-07, 4.8677122777e-07);
const mat4 I = mat4(1.08420217249e-19, 1.08420217249e-19, 1.08420217249e-19, 1.08420217249e-19, 1.29104734015e-07, 2770.66699219, -742741.25, -742741.25, -742741.25, -742741.25, -742741.25, -2770.66699219, -5.16573550158e-07, 4.86732574245e-07, 1.08420217249e-19, 1.08420217249e-19);
const mat4 J = mat4(1.08420217249e-19, 4.65663618776e-10, 1.19211847505e-07, 2048.05175781, 2048.05175781, -524301.25, -524291.25,-524299.25, -567981.25, -709973.25, -2901.33300781, -2901.33300781, 6.75519231663e-07, 6.45716795589e-07, 1.29143387539e-07, 4.75362527119e-10);
const mat4 K = mat4(1.29143387539e-07, 2773.33300781, -742741.25, -742741.25, -742741.25, -709973.25, -2219.33300781, -4.87391957904e-07, 5.19057039128e-07, 2261.33300781, -742741.25, -742709.25, -742573.25, -2770.66894531, -5.165662742e-07, 4.86616613671e-07);
const mat4 L = mat4(1.29143387539e-07, 2773.33300781, -742741.25, -742741.25, -742741.25, -709973.25, -567981.25, -535213.25, -524301.25, -524301.25, -524301.25, -2048.05175781, -4.76847390019e-07, 4.76839545627e-07, 1.19209431659e-07, 1.08420217249e-19);
const mat4 M = mat4(1.29143387539e-07, 2773.33300781, -742741.25, -742741.25, -742741.25, -2773.33300781, -5.19047375747e-07, 4.87389684167e-07, 1.29761843937e-07, 2773.33300781, -742741.25, -742741.25, -742741.25, -2773.33300781, -5.16573550158e-07, 4.8677122777e-07);
const mat4 N = mat4(1.29143387539e-07, 2773.33300781, -742741.25, -742741.25, -742741.25, -2773.33300781, -5.19056470694e-07, 4.87391957904e-07, 1.29182183173e-07, 2773.33300781, -742741.25, -742741.25, -742741.25, -2773.33300781, -5.16573550158e-07, 4.8677122777e-07);
const mat4 O = mat4(1.08984897264e-19, 4.75968708891e-10, 1.29764259782e-07, 2261.33300781, 2773.33300781, -742741.25, -742581.25,-742061.25, -742581.25, -742741.25, -2773.33300781, -2261.33300781, 5.19057039128e-07, 4.87391957904e-07, 1.19830161793e-07, 4.66267580101e-10);
const mat4 P = mat4(1.29143387539e-07, 2773.33300781, -742741.25, -742741.25, -742741.25, -742741.25, -2901.33300781, -6.73190356792e-07, 6.75053001942e-07, 1.53595161289e-19, 1.53595161289e-19, 1.19713429808e-19, 1.19713429808e-19, 1.11241420125e-19, 1.08984871414e-19, 1.08561251543e-19);
const mat4 Q = mat4(1.08984897264e-19, 4.75968708891e-10, 1.29764259782e-07, 2261.33300781, 2773.33300781, -742741.25, -742581.25,-742061.25, -2900.70800781, -2773.33300781, -578901.25, -570709.25, -535891.25, -2058.66699219, 4.77457547277e-07, 1.19209289551e-07);
const mat4 R = mat4(1.21692806943e-07, 2229.33300781, -709973.25, -742741.25, -742741.25, -742741.25, -2901.33300781, -6.7319081154e-07, 6.75053570376e-07, 2901.33300781, 742741.25, -578901.25, -578893.25, -2101.16894531, -4.79320647173e-07, 4.77455728287e-07);
const mat4 S = mat4(4.75326999982e-10, 1.29725606257e-07, 2261.17675781, 2773.30175781, -742733.25, -742739.25, -742739.25, -742741.25, -742229.25, -742229.25, -742613.25, -2772.83300781, 2260.70800781, 5.18438582731e-07, 1.21690987953e-07, 4.66230276608e-10);
const mat4 T = mat4(1.17419942313e-19, 1.46669048773e-19, 1.53557808914e-19, 1.53557808914e-19, 1.68841012282e-07, 2901.33300781, -742741.25, -742741.25, -742741.25, -2901.33300781, -6.7536404913e-07, 6.7536404913e-07, 1.53557808914e-19, 1.46675666218e-19, 1.17453029538e-19, 1.10670148515e-19);
const mat4 U = mat4(1.17455226736e-19, 6.30582808192e-10, 1.68879807916e-07, 2901.33300781, 2901.33300781, -709973.25, -567989.25,-535213.25, -567989.25, -709973.25, -2901.33300781, -2901.33300781, 6.75519231663e-07, 6.45716795589e-07, 1.29143387539e-07, 4.75362527119e-10);
const mat4 V = mat4(1.17454683899e-19, 1.4681680391e-19, 1.53595161289e-19, 6.59686638649e-10, 1.68879807916e-07, 2773.33300781, -568021.25, -535221.25, -568021.25, -2773.33300781, -6.75519231663e-07, 6.75519117976e-07, 1.68879665807e-07, 6.30573482319e-10, 1.17455226736e-19, 1.10678833911e-19);
const mat4 W = mat4(1.29143387539e-07, 2773.33300781, -742741.25, -742741.25, -742741.25, -2773.33300781, -5.16583213539e-07, 4.86809881295e-07, 1.29145803385e-07, 2773.33300781, -742741.25, -742741.25, -742741.25, -2773.33300781, -5.16573550158e-07, 4.8677122777e-07);
const mat4 X = mat4(1.29105188762e-07, 2770.54394531, -742573.25, -742709.25, -742741.25, -2773.33300781, -5.19057039128e-07, 4.87391957904e-07, 1.29764259782e-07, 2773.33300781, -742741.25, -742709.25, -742573.25, -2770.66894531, -5.165662742e-07, 4.86733028993e-07);
const mat4 Y = mat4(1.17446412093e-19, 1.46781545336e-19, 1.53585933055e-19, 1.5359296409e-19, 1.6887921106e-07, 2773.33300781, -570709.25, -535893.25, -570709.25, -2773.33300781, -6.7551684424e-07, 6.75509568282e-07, 1.53586450043e-19, 1.46783613288e-19, 1.17454683899e-19, 1.1067676596e-19);
const mat4 Z = mat4(5.03860619894e-10, 692.503173828, 741901.25, -741941.25, -742101.25, -742229.25, -742741.25, -742741.25, -742741.25, -742739.25, -742731.25, -742699.25, -578731.25, -2090.62792969, -4.79281993648e-07, 4.76837158203e-07);
const mat4 _0 = mat4(1.08455372425e-19, 1.0899316907e-19, 4.75968708891e-10, 1.29764259782e-07, 2773.33300781, -2773.33300781, -742069.25, -740013.25, -742069.25, -2773.33300781, -2773.33300781, -5.19057039128e-07, 4.87391957904e-07, 1.1983925674e-07, 4.66267580101e-10, 1.08429006043e-19);
const mat4 _1 = mat4(1.08420217249e-19, 1.08420217249e-19, 1.08420217249e-19, 1.10670148515e-19, 1.29725464149e-07, 2773.33300781, -742741.25, -742741.25, -742741.25, -2773.33300781, -5.16573550158e-07, 4.8677122777e-07, 1.08420217249e-19, 1.08420217249e-19, 1.08420217249e-19, 1.08420217249e-19);
const mat4 _2 = mat4(1.21654608165e-07, 2226.55175781, -578741.25, -709845.25, -742613.25, -742229.25, -742229.25, -742741.25, -742739.25, -742739.25, -742731.25, -709963.25, -578859.25, -2228.66699219, -4.86769408781e-07, 4.79311097479e-07);
const mat4 _3 = mat4(4.75213313145e-10, 1.29601914978e-07, 2260.67675781, 2772.67675781, -742605.25, -742227.25, -742227.25, -742741.25, -742741.25, -742741.25, -742741.25, -2773.33300781, 2261.20800781, 5.18904244018e-07, 1.21692806943e-07, 4.68051042368e-10);
const mat4 _4 = mat4(1.08455501672e-19, 4.66305327684e-10, 4.68238114948e-10, 4.75968819913e-10, 5.06891639773e-10, 5.14167597387e-10, 1.68872531958e-07, 2900.70800781, -742741.25, -742741.25, -742741.25, -2773.33300781, -5.16583213539e-07, 4.86773501507e-07, 4.65812832751e-10, 1.08420760086e-19);
const mat4 _5 = mat4(5.04459252149e-10, 1.61419677625e-07, 2901.30175781, 2901.30175781, -742733.25, -742731.25, -742731.25, -742739.25, -742229.25, -742229.25, -742229.25, -2898.83300781, -2898.83300781, 6.43243083687e-07, 5.16544446327e-07, 1.21654608165e-07);
const mat4 _6 = mat4(1.08984897264e-19, 4.75968708891e-10, 1.29764259782e-07, 2261.33300781, 2773.33300781, -742741.25, -742731.25,-742219.25, -742227.25, -742229.25, -2899.33300781, -2763.33300781, 5.16612317369e-07, 4.86780777464e-07, 1.19218981354e-07, 4.65670724203e-10);
const mat4 _7 = mat4(1.17419942313e-19, 1.46669048773e-19, 1.53557808914e-19, 1.68838624859e-07, 2900.63574219, -742573.25, -742613.25, -742741.25, -742741.25, -2901.33300781, -6.75519117976e-07, 6.7551684424e-07, 6.59648891066e-10, 1.46677837567e-19, 1.17454683899e-19, 1.10670148515e-19);
const mat4 _8 = mat4(4.75326999982e-10, 1.29726061004e-07, 2261.20800781, -709973.25, -742741.25, -742741.25, -742738.75, -740178.75, -740181.25, -742741.25, -742741.25, -709973.25, -2261.20800781, 5.18904244018e-07, 1.21692806943e-07, 4.68051042368e-10);
const mat4 _9 = mat4(1.08982700065e-19, 1.10810769219e-19, 1.29761843937e-07, 2261.32324219, -709971.25, -742739.25, -742227.25, -740179.25, -742229.25, -742741.25, -2773.33300781, -2261.33300781, 5.19057039128e-07, 4.87391957904e-07, 1.19830161793e-07, 4.66267580101e-10);
const mat4 dot = mat4(1.08420217249e-19, 1.08420217249e-19, 1.08420217249e-19, 1.08420217249e-19, 1.19209431659e-07, 2048.01074219, -524291.25, -524291.25, -524291.25, -2048.01074219, -4.76837726637e-07, 4.7683727189e-07, 1.08420217249e-19, 1.08420217249e-19, 1.08420217249e-19, 1.08420217249e-19);
const mat4 comma = mat4(1.08420217249e-19, 1.08420217249e-19, 1.08420217249e-19, 1.08420217249e-19, 1.08420760086e-19, 1.19211819083e-07, 0.000122073397506, -524301.25, -2048.05175781, -4.76847390019e-07, 4.76839545627e-07, 1.19209431659e-07, 1.08420217249e-19, 1.08420217249e-19, 1.08420217249e-19, 1.08420217249e-19);
const mat4 bang = mat4(1.08420217249e-19, 1.08420217249e-19, 1.08420217249e-19, 1.08420217249e-19, 1.29142932792e-07, 2773.29199219, -742739.25, -742739.25, -742739.25, -2773.29199219, -5.16573550158e-07, 4.86770773023e-07, 1.08420217249e-19, 1.08420217249e-19, 1.08420217249e-19, 1.08420217249e-19);
const mat4 qm = mat4(1.08420217249e-19, 1.10643678735e-19, 1.17869928566e-19, 1.19676491024e-19, 1.68838738546e-07, 2900.66699219,-742219.25, -742739.25, -742739.25, -2901.32324219, -6.75509568282e-07, 6.45678142064e-07, 1.17984493081e-19, 1.10678420321e-19, 1.08982700065e-19, 1.08420217249e-19);

#define CHAR_W (17./16.)
#define _(cc) col = ch(cc, col, uv); uv.x -= CHAR_W;
#define SP uv.x -= CHAR_W;

const float text_width = 65. * CHAR_W;

float scroll(vec2 uv) {
	float zoom = 4.0;
	uv *= zoom;

	uv.x = mod(uv.x +4.*iTime, text_width) - .5;
	uv.y -= (2.0*horizon)+1.5;
	uv.y += sin(-.65*uv.x*(uv.x+cos(iTime*2.+uv.y*10.)) + iTime * 24.) * -.04;

	float col = 0.;

	_(bang)_(G)_(R)_(E)_(E)_(T)_(I)_(N)_(X)_(dot)_(T)_(O)_(dot)_(A)_(L)_(L)_(dot)_(T)_(H)_(E)_(dot)_(B)_(A)_(R)_(D)_(S)_(dot)
    _(O)_(U)_(T)_(dot)_(T)_(H)_(E)_(R)_(E)_(bang)_(dot)_(H)_(A)_(P)_(P)_(Y)_(dot)
    _(B)_(A)_(R)_(D)_(I)_(N)_(G)_(dot)

	return pow(col,3.)*1.66*(sin(iTime)*0.5+0.7);
}

vec3 copper(in vec3 c) {
    vec3 rgb = clamp( abs(mod(c.x*6.0+vec3(0.0,4.0,2.0),6.0)-3.0)-1.0, 0.0, 1.0 );
	rgb = rgb*rgb*(3.0-2.0*rgb);
	return c.z * mix( vec3(1.0), rgb, c.y) * vec3(1.0,1.0,2.5);
}



float blend = 0;
void main(void)
{

    if (iTime == 0)
    {
        blend = 0;
    }
        

        blend = iTime/11;
        if (blend <= 1)
        {
            vec4 fragColor;
            mainImage(fragColor, gl_FragCoord.xy);
            gl_FragColor = mix(mainImage2(fragColor.xyz, gl_FragCoord.xy), fragColor, blend );
        }
        else
        {
            vec4 fragColor;
            mainImage(fragColor, gl_FragCoord.xy);
	        vec2 uv = ( gl_FragCoord.xy/resolution.xy ) -.5;
	        uv.x*=resolution.x/resolution.y;
	        vec3 col;
            if (iTime > 15)
            {
	            if (uv.y > -0.4) 
	            {
                    vec3 copper = copper(vec3(uv.y*4.-time*.4,.3,.4));
		            float text = scroll(uv)+scroll(uv*vec2(2.,-1.)-vec2(0.,.68))*.5;
		            col = text > 0. ? vec3(.5)*text*copper : col;
	            }
            }
            gl_FragColor = fragColor;
            gl_FragColor += vec4(col,0.5);
        }
    
}
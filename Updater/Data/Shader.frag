// Tune: hotdogs.xm

#version 330

#ifdef GL_ES
precision highp float;
#endif

uniform vec2 resolution;
uniform float time;

float iTime = 0.0;
vec3  iResolution = vec3(0.0);

// --------------------------------- //
// --------[ Discord Logo ]--------- //
// --------------------------------- //
float sdEllipse(vec2 p, vec2 r)
{
    float f = length(p / r),
          g = length(p / r / r);
    
    return f * (f - 1.) / g;
}

float sdDiscord(vec2 p)
{
    p.x = abs(p.x);
    
    float d = length(p + vec2(0, .52)) - .91;
          d = max(d, length(p - vec2(0, .41)) - .83);
          d = max(d, length(p - vec2(.82, .09)) - .74);
          d = max(d, min(.54 - length(p - vec2(0, .21)), 
                         -(.81 * p.x + p.y + .01) / length(vec2(.81, 1))));
          d = min(d, max(length(p - vec2(0, .3)) - .59,
                         length(p + vec2(0, .36)) - .7));
          d = max(d, length(p + vec2(.34, .16)) - .84);
          d = max(d, -sdEllipse(p - vec2(.165, -.038), vec2(.09, .1)));
    
    return d;
}

float helper(float x, float a, float b, float c, float d, float m, float n)
{
    float x2 = x * x, x3 = x2 * x;
    return (a * x3 + b * x2 + c * x + d) / (x2 + m * x + n);
}

vec2 rotate(vec2 p, float a) 
{
    return mat2(cos(a), sin(a), -sin(a), cos(a)) * p;
}

void Discord(out vec4 O, vec2 I)
{
    vec2 p = (I + I - iResolution.xy) / iResolution.y;
    float r1 = .5, r2 = .351;
    
    float t = 0;
    if (iTime < 10)
        t= fract(iTime);
    else
        t= iTime -10;

    float k = t * (t - 1.), k2 = k * k;
    
    float s1 = k2 * helper(t, 18.69441, -38.86151, 23.25180, -3.53363, -1.27120, 0.535271) + 1.;
    float s2 = k2 * helper(t, 10.50167, -22.56208, 13.91412, -2.31023, -1.07310, 0.364973) + 1.;
    float s3 = k2 * helper(t, 21.72620, -40.39696, 22.28075, -3.43598, -1.22944, 0.605636) + 1.;
    
    float a1 = t * (k * helper(t, -0.668041,  1.870240, -2.58381, 0.798835, -1.029650, 0.321038) + t);
    float a2 = t * (k * helper(t,  0.274534, -0.241024, -1.22258, 0.473114, -0.881859, 0.239837) + t);
    float a3 = t * (k * helper(t, -2.536390,  5.595520, -4.82078, 1.054856, -0.926966, 0.280363) + t);
    
    a1 = mix(a1, a2, smoothstep(.15, .25, abs(t - .47)));
    a3 = mix(a3, a2, smoothstep(.15, .25, abs(t - .47)));
    
    vec2 p1 = rotate(p, radians(360. * a1)) / s1;
    vec2 p2 = rotate(p, radians(360. * a2)) / s2;
    vec2 p3 = rotate(p, radians(360. * a3)) / s3;
    
    float d1 = max(sdDiscord(p1), r1 - length(p1)) * s1;
    float d2 = max(sdDiscord(p2), max(length(p2) - r1, r2 - length(p2))) * s2;
    float d3 = max(sdDiscord(p3), length(p3) - r2) * s3;
    
    float d = (abs(t - .46) < .39 ? min(min(d1, d2), d3) : sdDiscord(p2) * s2) * iResolution.y;
    
    O = vec4(mix(vec3(1), vec3(28, 29, 35) / 255., smoothstep(-1., 1., d)), 1);
}


// --------------------------------- //
// ----------[ Main Image ]--------- //
// --------------------------------- //
#define ZERO 0
#define MIN_FLOAT 1e-6
#define MAX_FLOAT 1e6
#define EPSILON 1e-3
#define saturate(x) clamp(x, 0., 1.)
#define UP vec3(0., 1., 0.)
const float PI = acos(-1.);
float utime = 0.;

struct Ray{ vec3 origin, dir; };
struct Sphere{vec3 origin; float rad;};

vec3 rayDirection(float fieldOfView, vec2 size, vec2 fragCoord) {
    vec2 xy = fragCoord - size / 2.;
    float z = size.y / tan(radians(fieldOfView) / 2.);
    return normalize(vec3(xy, -z));
}

mat4 viewMatrix(vec3 eye, vec3 center, vec3 up) {
    vec3 f = normalize(center - eye),
         s = normalize(cross(f, up)),
         u = cross(s, f);
    return mat4(vec4(s, 0.), vec4(u, 0.), vec4(-f, 0.), vec4(vec3(0.), 1.));
}

Ray makeViewRay(vec2 coord, vec2 res, float a){
    //vec3 lookAt = vec3(0., 0., 0.);
    //vec3 origin = vec3(15. * cos(a), 1., 15. * sin(a));
    
    vec3 lookAt = vec3(0., 1.5 ,0.);
    vec3 origin = vec3(0., 2., 18.);
    vec3 viewDir = rayDirection(60., res, coord);
    mat4 viewToWorld = viewMatrix(origin, lookAt, vec3(0., 1., 0.));
    vec3 rd = (viewToWorld * vec4(viewDir, 1.0)).xyz;
    
    return Ray(origin, rd);
}

struct Box{ vec3 origin; vec3 size; };
#define MIN x
#define MAX y
bool box_hit(const in Box inbox, const in Ray inray){
    vec2 tx, ty, tz;
    vec3 maxbounds = inbox.origin + vec3(inbox.size);
    vec3 minbounds = inbox.origin + vec3(-inbox.size);
    tx = ((inray.dir.x >= 0.?vec2(minbounds.x, maxbounds.x):vec2(maxbounds.x, minbounds.x)) - inray.origin.x) / inray.dir.x;
	ty = ((inray.dir.y >= 0.?vec2(minbounds.y, maxbounds.y):vec2(maxbounds.y, minbounds.y)) - inray.origin.y) / inray.dir.y;
    if ((tx.MIN > ty.MAX) || (ty.MIN > tx.MAX))
        return false;
    tx = vec2(max(tx.MIN, ty.MIN), min(tx.MAX, ty.MAX));
	tz = ((inray.dir.z >= 0.?vec2(minbounds.z, maxbounds.z):vec2(maxbounds.z, minbounds.z)) - inray.origin.z) / inray.dir.z;
    if ((tx.MIN > tz.MAX) || (tz.MIN > tx.MAX))
        return false;
    tx = vec2(max(tx.MIN, tz.MIN), min(tx.MAX, tz.MAX));
    
    if(tx.MIN >= 0.){
    	return true;
    }
        
    return false;
}

bool sphere_hit(const in Sphere sphere, const in Ray inray) {
    vec3 oc = inray.origin - sphere.origin;
    float a = dot(inray.dir, inray.dir);
    float b = dot(oc, inray.dir);
    float c = dot(oc, oc) - sphere.rad*sphere.rad;
    float discriminant = b*b - a*c;
    if (discriminant > 0.) {
        return true;
        //return (-b - sqrt(discriminant))/a;
    }
    return false;
}

float Hash21(vec2 uv){
    float f = uv.x + uv.y * 37.0;
    return fract(sin(f)*104003.9);
}

vec2 Hash22(vec2 uv){
    float f = uv.x + uv.y * 37.0;
    return fract(cos(f)*vec2(10003.579, 37049.7));
}

float sdBox(vec3 p, vec3 radius){
  vec3 dist = abs(p) - radius;
  return min(max(dist.x, max(dist.y, dist.z)), 0.0) + length(max(dist, 0.0));
}

// capped cylinder distance field
float cylCap(vec3 p, float r, float lenRad){
    float a = length(p.xy) - r;
    a = max(a, abs(p.z) - lenRad);
    return a;
}

vec3 hsv2rgb(vec3 c) {
  const vec4 K = vec4(1.0, 2.0 / 3.0, 1.0 / 3.0, 3.0);
  vec3 p = abs(fract(c.xxx + K.xyz) * 6.0 - K.www);
  return c.z * mix(K.xxx, clamp(p - K.xxx, 0.0, 1.0), c.y);
}

float sun(vec2 p) {
  vec2 sp = p + vec2(0., .55);
  float clr = length(sp) - 1.;
  return clr * .5;
}

vec3 sunEffect(vec2 p) {
  vec3 res = vec3(0.1);
  vec3 skyCol1 = hsv2rgb(vec3(283.0/360.0, 0.83, 0.16));
  vec3 skyCol2 = hsv2rgb(vec3(297.0/360.0, 0.79, 0.43));
  res = mix(skyCol1, skyCol2, pow(clamp(0.5*(1.0+p.y+0.1*sin(4.0*p.x)), 0.0, 1.0), 4.0));
  
  p.y -= .375;
  float ds = sun(p);
  vec3 sunCol = mix(vec3(1.0, 0.0, 1.0), vec3(1.0, 1.0, 0.0), clamp(0.5 - .5 * p.y, 0.0, 1.0));
  vec3 glareCol = sqrt(sunCol);
  
  res += glareCol*(exp(-30.0*ds))*step(0.0, ds);
  res = mix(res, sunCol, smoothstep(-.01, 0., -ds));

  return res;
}

float sdEllipsoid(vec3 p, vec3 r){
  float k0 = length(p/r);
  float k1 = length(p/(r*r));
  return k0*(k0-1.0)/k1;
}

float sdRoundBox(vec3 p, vec3 b, float r){
  vec3 q = abs(p) - b;
  return length(max(q,0.0)) + min(max(q.x,max(q.y,q.z)),0.0) - r;
}

float sdPlane(vec3 p, vec3 n, float h ) {
    return dot(p,n) + h;
}

vec2 opMin(vec2 a, vec2 b){
    if(a.x <= b.x) return a; else return b;
}

//by iq

vec2 hash( in vec2 x )  // replace this by something better
{
    const vec2 k = vec2( 0.3183099, 0.3678794 );
    x = x*k + k.yx;
    return -1.0 + 2.0*fract( 16.0 * k*fract( x.x*x.y*(x.x+x.y)) );
}

#define quintic 0

// return gradient noise (in x) and its derivatives (in yz)
vec3 noised( in vec2 p )
{
    vec2 i = floor( p );
    vec2 f = fract( p );

#if quintic
    // quintic interpolation
    vec2 u = f*f*f*(f*(f*6.0-15.0)+10.0);
    vec2 du = 30.0*f*f*(f*(f-2.0)+1.0);
#else
    // cubic interpolation
    vec2 u = f*f*(3.0-2.0*f);
    vec2 du = 6.0*f*(1.0-f);
#endif    
    
    vec2 ga = hash( i + vec2(0.0,0.0) );
    vec2 gb = hash( i + vec2(1.0,0.0) );
    vec2 gc = hash( i + vec2(0.0,1.0) );
    vec2 gd = hash( i + vec2(1.0,1.0) );
    
    float va = dot( ga, f - vec2(0.0,0.0) );
    float vb = dot( gb, f - vec2(1.0,0.0) );
    float vc = dot( gc, f - vec2(0.0,1.0) );
    float vd = dot( gd, f - vec2(1.0,1.0) );

    return vec3( va + u.x*(vb-va) + u.y*(vc-va) + u.x*u.y*(va-vb-vc+vd),   // value
                 ga + u.x*(gb-ga) + u.y*(gc-ga) + u.x*u.y*(ga-gb-gc+gd) +  // derivatives
                 du * (u.yx*(va-vb-vc+vd) + vec2(vb,vc) - va));
}

float trail(vec2 uv){
    return smoothstep(distance(uv.y, 9. + noised(uv + vec2(utime)).x), .5 - uv.x * .5 - .5, - uv.x * .5 - .5);
}

vec3 sea(vec2 uv){
    
    uv *= vec2(1., 3.);
    
    const vec2 size = vec2(10., 0.);
    const vec3 off = vec3(-20.,0, 50.);
    
    float s11 = noised(uv).x;
    float s01 = noised(uv + off.xy).x;
    float s21 = noised(uv + off.zy).x;
    float s10 = noised(uv + off.yx).x;
    float s12 = noised(uv + off.yz).x;
    vec3 va = normalize(vec3(size.xy, s21-s01));
    vec3 vb = normalize(vec3(size.yx, s12-s10));
    vec4 bump = vec4( cross(va,vb), s11 );
    return bump.xzy;
}

#define MOD_POS(poz) poz.x += iTime * 4. + mod(floor(poz.z), 2.) * .35;
//#define MOD_POS(poz) poz.x += mod(floor(poz.z), 2.) * .35;

vec2 CityBlock(vec3 p, vec2 pint){
    // Get random numbers for this block by hashing the city block variable
    vec4 rand;
    rand.xy = Hash22(pint);
    rand.zw = Hash22(rand.xy);
    vec2 rand2 = Hash22(rand.zw);

    // Radius of the building
    float baseRad = 0.2 + (rand.x) * 0.1;
    baseRad = floor(baseRad * 20.0+0.5)/20.0;   // try to snap this for window texture

    // make position relative to the middle of the block
    vec3 baseCenter = p - vec3(0.5, 0.0, 0.5);
    float height = .75 * rand.w*rand.z + 0.3; // height of first building block
    // Make the city skyline higher in the middle of the city.
    height *= 1.5+(baseRad-0.15)*20.0;
    height += 0.1;  // minimum building height
    //height += sin(iTime + pint.x);    // animate the building heights if you're feeling silly
    height = floor(height*20.0)*0.05;   // height is in floor units - each floor is 0.05 high.
    float d = sdBox(baseCenter, vec3(baseRad, height, baseRad)); // large building piece

    // road
    d = min(d, p.y);

    //if (length(pint.xy) > 8.0) return vec2(d, mat);   // Hack to LOD in the distance

    // height of second building section
    float height2 = max(0.0, rand.y * 2.0 - 1.0);
    height2 = floor(height2*20.0)*0.05; // floor units
    rand2 = floor(rand2*20.0)*0.05; // floor units
    // size pieces of building
    d = min(d, sdBox(baseCenter - vec3(0.0, height, 0.0), vec3(baseRad, height2 - rand2.y, baseRad*0.4)));
    d = min(d, sdBox(baseCenter - vec3(0.0, height, 0.0), vec3(baseRad*0.4, height2 - rand2.x, baseRad)));
    // second building section
    if (rand2.y > 0.25)
    {
        d = min(d, sdBox(baseCenter - vec3(0.0, height, 0.0), vec3(baseRad*0.8, height2, baseRad*0.8)));
        // subtract off piece from top so it looks like there's a wall around the roof.
        float topWidth = baseRad;
        if (height2 > 0.0) topWidth = baseRad * 0.8;
        d = max(d, -sdBox(baseCenter - vec3(0.0, height+height2, 0.0), vec3(topWidth-0.0125, 0.015, topWidth-0.0125)));
    }
    else
    {
        // Cylinder top section of building
        if (height2 > 0.0) d = min(d, cylCap((baseCenter - vec3(0.0, height, 0.0)).xzy, baseRad*0.8, height2));
    }
    // mini elevator shaft boxes on top of building
    d = min(d, sdBox(baseCenter - vec3((rand.x-0.5)*baseRad, height+height2, (rand.y-0.5)*baseRad),
                     vec3(baseRad*0.3*rand.z, 0.1*rand2.y, baseRad*0.3*rand2.x+0.025)));
    // mirror another box (and scale it) so we get 2 boxes for the price of 1.
    vec3 boxPos = baseCenter - vec3((rand2.x-0.5)*baseRad, height+height2, (rand2.y-0.5)*baseRad);
    float big = sign(boxPos.x);
    boxPos.x = abs(boxPos.x)-0.02 - baseRad*0.3*rand.w;
    d = min(d, sdBox(boxPos,
    vec3(baseRad*0.3*rand.w, 0.07*rand.y, baseRad*0.2*rand.x + big*0.025)));

    // Put domes on some building tops for variety
    if (rand.y < 0.04)
    {
        d = min(d, length(baseCenter - vec3(0.0, height, 0.0)) - baseRad*0.8);
    }

    return vec2(d, 0.0);
}

float city(vec3 p){
    MOD_POS(p)
    vec3 rep = p;
    rep.xz = fract(p.xz);
    return CityBlock(rep, floor(p.xz)).x;
}

vec3 estimateCityNormal(vec3 p) {
    return normalize(vec3(
        city(vec3(p.x + EPSILON, p.y, p.z)) - city(vec3(p.x - EPSILON, p.y, p.z)),
        city(vec3(p.x, p.y + EPSILON, p.z)) - city(vec3(p.x, p.y - EPSILON, p.z)),
        city(vec3(p.x, p.y, p.z  + EPSILON)) - city(vec3(p.x, p.y, p.z - EPSILON))
    ));
}

const float voxelPad = .2;
float marchCity(vec3 eye, vec3 dir, float start, float end) {
    float t = start;
    vec3 pos = vec3(0.0);
    for (int i = ZERO; i < 128; i++) {
        pos = eye + dir * t;
        float dist = city(pos);
        if(abs(dist) < EPSILON){
            return t;
        }
        
        float walk = dist;
        float dx = -fract(pos.x);
        float dz = -fract(pos.z);
        float nearestVoxel = min(fract(dx/dir.x), fract(dz/dir.z)) + voxelPad;
        nearestVoxel = max(voxelPad, nearestVoxel);
        walk = min(walk, nearestVoxel);


        t += walk;
        if ((t > end)) break;
    }
    return -1.;
}

vec3 cityColor(vec3 pos, vec3 view){
    vec3 normal = estimateCityNormal(pos);
    MOD_POS(pos);
    float winid = Hash21(floor(pos.xy * vec2(2., 6.)));
    float winStencil = step(winid, Hash21(floor(pos.xz)) * .75);
    float win = step(distance(.125, abs(fract(pos.x) - .5)), .075)
              * step(distance(fract(pos.y * 6.), .5), .4);
    return hsv2rgb(vec3(.525 + (winid * .2 - .1), 1., 1.)) * winStencil * win * step(.000001, dot(normal, vec3(0., 0., 1.)))
         + vec3(1.0, 1.0, 0.0) * abs(dot(normal, vec3(1., 0., 0.))) * abs(dot(view, normal)) * .2;
}

const vec3 boat_center = vec3(0., -.3, -9.);
vec2 boat(vec3 pos){
    float scale = .65;
    float ang = noised(vec2(iTime)).x * .2 - .05;
    pos.xy *= mat2(cos(ang), -sin(ang), sin(ang), cos(ang));
    pos += boat_center + vec3(0., ang, 0.);
    
    vec2 res = vec2(0.);
    
    float c = length(pos + vec3(0., 7., 0.) * scale) - 7. * scale;
    res.x = max(res.x, c);
    float d = length(pos + vec3(-.75, 7., 0.) * scale) - 7. * scale;
    res.x = max(res.x, d);
    
    {
        float a = sdRoundBox(pos + vec3(.65, .35, 0.) * scale, vec3(.9, .25, .65 + pos.x * .2) * scale, .02 * scale);
        res.x = min(res.x, a);
    }
    
    {
        float a = sdEllipsoid(pos + vec3(1.8, .38, 0.) * scale, vec3(3., .75, .86) * scale);
        a = max(a, sdPlane(pos, normalize(vec3(-1., 1., 0.)), -.05 * scale));
        a = max(a, sdPlane(pos, normalize(vec3(0., -1., 0.)), -0.1));
        res = opMin(res, vec2(a, 1.));
    }
    
    {
        vec3 mpos = pos + vec3(0., .2 - pos.x * .1, 0.) * scale;
        float a = sdEllipsoid(mpos + vec3(0., 0., .25) * scale, vec3(2., .5, 1.) * scale);
        float b = sdEllipsoid(mpos - vec3(0., 0., .25) * scale, vec3(2., .5, 1.) * scale);
        res.x = max(max(a, b), res.x);
    }
    
    {
        float a = sdRoundBox(pos + vec3(.3, 0.2, 0.) * scale, vec3(.9, .25 + pos.y * .1, .42 + pos.x * .2) * scale, .25 * scale);
        float b = sdRoundBox(pos + vec3(.3, 0.2, 0.) * scale, vec3(.85, .24, .4 + pos.x * .2) * scale, .25 * scale);
        float comb = max(a, -b);
        comb = max(comb, sdPlane(pos, normalize(vec3(1., -.2, 0.)), .4 * scale));
        comb = max(comb, -sdPlane(pos, normalize(vec3(1., -.8, 0.)), .6 * scale));
        comb = max(comb, -pos.y - .2 * scale);
        res.x = min(res.x, comb);
    }

    float e = sdPlane(pos, -normalize(vec3(-1., .5, 0.) * scale), 1.2 * scale);
    res.x = max(res.x, -e);
    
    return res;
}

vec3 boatNormals(vec3 pos){
    vec2 eps = vec2(0.0, EPSILON);
    vec3 n = normalize(vec3(
        boat(pos + eps.yxx).x - boat(pos - eps.yxx).x,
        boat(pos + eps.xyx).x - boat(pos - eps.xyx).x,
        boat(pos + eps.xxy).x - boat(pos - eps.xxy).x));
    return n;
}

vec2 marchBoat(in Ray r){
    float t = .01;
    for(int i = ZERO; i <= 64; i++){
        vec3 p = r.origin + r.dir * t;
        vec2 dst = boat(p);
        if(dst.x < .01)
            return vec2(t, dst.y);
        t += dst.x;
    }
    return vec2(-1.);
}

vec3 boatColor(vec3 p, float matid){
    vec3 bn = boatNormals(p);
    vec3 albedo = mix(vec3(1.), vec3(1., 0., 0.), step(abs(-.4 - dot(bn, vec3(0., 1., 0.))), .2));
    
    vec3 sunPos = vec3(0., 3., -5.);
    vec3 sun = max(dot(bn, normalize(sunPos - p)), .15) * mix(vec3(1.0, 0.0, 1.0), vec3(1.0, 1.0, 0.0), .4);
    
    if(matid == 0.)
        return albedo * sun + smoothstep(.15, 1., dot(bn, normalize(vec3(1., -1., 0.)))) * vec3(1.0, .0, 1.0) * (noised(p.xz * 3. + vec2(iTime * 4., 0.)).x * .5 + .5) * .25;
    else if(matid == 1.)
        return vec3(0.);
    else
        return vec3(0.);
}

vec3 geometry(Ray r){
    vec3 color = vec3(0.);

    
    float start = (-r.origin.z)/r.dir.z;
    float end = (-3.-r.origin.z)/r.dir.z;
    float cityDist = -1.;
    if(box_hit(Box(vec3(0., 2.25, -2.), vec3(10., 2.25, 2.)), r)){
        cityDist = marchCity(r.origin, r.dir, start, end);
    }
    if (cityDist >= 0.) {
        color = cityColor(r.origin + cityDist * r.dir, r.dir);
    }else{
        vec3 backPlane = r.origin + end * r.dir;
        color = sunEffect(backPlane.xy * .125);
    }
    
    if(box_hit(Box(vec3(0., .4, 9.), vec3(1.2, .35, .5)), r)){
        vec2 boatDist = marchBoat(r);
        vec3 boatPoint = r.origin + boatDist.x * r.dir;
        if (boatDist.x >= 0. && boatPoint.y >= 0.) {
            color = boatColor(boatPoint, boatDist.y);
        }
    }
    
    return color;
}

void mainImage( out vec4 fragColor, in vec2 fragCoord) {
    vec2 uv = (2. * fragCoord.xy - iResolution.xy)/iResolution.y;
    utime = iTime;
    Ray originalRay = makeViewRay(fragCoord, iResolution.xy, iResolution.x * 5.);
    Ray r = originalRay;
    
    vec3 color = vec3(0.);
    
    float groundDst = (-r.origin.y)/r.dir.y;
    vec3 groundPos = r.origin + groundDst * r.dir;
    vec3 nor = vec3(0., 1., 0.);
    if(groundPos.z > 0. && r.dir.y < 0.) {
        vec3 nor;
        {
            vec2 uv = mod(groundPos.xz + vec2(iTime * 2., 0.), vec2(100.)) * vec2(5., 10.);
    
            const vec2 size = vec2(10., 0.);
            const vec3 off = vec3(-20.,0, 50.);

            float s11 = noised(uv).x + trail(groundPos.xz);
            float s01 = noised(uv + off.xy).x + trail(groundPos.xz + vec2(.4, 0.));
            float s21 = noised(uv + off.zy).x + trail(groundPos.xz + vec2(-.4, 0.));
            float s10 = noised(uv + off.yx).x + trail(groundPos.xz + vec2(.0, .4));
            float s12 = noised(uv + off.yz).x + trail(groundPos.xz + vec2(.0, -.4));
            vec3 va = normalize(vec3(size.xy, s21-s01));
            vec3 vb = normalize(vec3(size.yx, s12-s10));
            vec4 bump = vec4( cross(va,vb), s11 );
            nor =  bump.xzy;
        }
        
        vec3 reflected = reflect(r.dir, nor);
        r = Ray(groundPos, reflected);
    }
    color = geometry(r);
    
    if(box_hit(Box(vec3(0., .3, 9.), vec3(1.2, .3, .5)), originalRay)){
        vec2 boatDist = marchBoat(originalRay);
        if (boatDist.x >= 0. && boatDist.x < groundDst) {
            vec3 boatPoint = originalRay.origin + boatDist.x * originalRay.dir;
            color = boatColor(boatPoint, boatDist.y);
        }
    }
    
    fragColor = vec4(color, 1.);
}


// --------------------------------- //
// --------[ Text Drawing]---------- //
// --------------------------------- //
float d_neon = 1e6;
vec2 g_uv;
vec2 ch_pos;

const vec2 ch_space = vec2(1.8, 0.0);
const vec2 ch_start = vec2(-15.0, 0.0);

float neonSeg(vec2 p, vec2 a, vec2 b) {
    vec2 pa = p - a, ba = b - a;
    float h = clamp(dot(pa, ba) / dot(ba, ba), 0.0, 1.0);
    float dist = length(pa - ba * h);
    return abs(dist - 0.05); 
}

bool bit(int n, int b) {
    return mod(floor(float(n) / exp2(float(b))), 2.0) != 0.0;
}

void ddigit(int n) {
    float v = 1e6;
    vec2 p = g_uv - ch_pos;
    
    if (bit(n,  0)) v = min(v, neonSeg(p, vec2( 0.5,  0.06), vec2( 0.5,  0.94)));
    if (bit(n,  1)) v = min(v, neonSeg(p, vec2( 0.44, 1.0),  vec2( 0.06, 1.0)));
    if (bit(n,  2)) v = min(v, neonSeg(p, vec2(-0.06, 1.0),  vec2(-0.44, 1.0)));
    if (bit(n,  3)) v = min(v, neonSeg(p, vec2(-0.5,  0.94), vec2(-0.5,  0.06)));
    if (bit(n,  4)) v = min(v, neonSeg(p, vec2(-0.5, -0.06), vec2(-0.5, -0.94)));
    if (bit(n,  5)) v = min(v, neonSeg(p, vec2(-0.44,-1.0),  vec2(-0.06,-1.0)));
    if (bit(n,  6)) v = min(v, neonSeg(p, vec2( 0.06,-1.0),  vec2( 0.44,-1.0)));
    if (bit(n,  7)) v = min(v, neonSeg(p, vec2( 0.5, -0.94), vec2( 0.5, -0.06)));
    if (bit(n,  8)) v = min(v, neonSeg(p, vec2( 0.06, 0.0),  vec2( 0.44, 0.0)));
    if (bit(n,  9)) v = min(v, neonSeg(p, vec2( 0.06, 0.06), vec2( 0.44, 0.94)));
    if (bit(n, 10)) v = min(v, neonSeg(p, vec2( 0.0,  0.06), vec2( 0.0,  0.94)));
    if (bit(n, 11)) v = min(v, neonSeg(p, vec2(-0.06, 0.06), vec2(-0.44, 0.94)));
    if (bit(n, 12)) v = min(v, neonSeg(p, vec2(-0.44, 0.0),  vec2(-0.06, 0.0)));
    if (bit(n, 13)) v = min(v, neonSeg(p, vec2(-0.06,-0.06), vec2(-0.44,-0.94)));
    if (bit(n, 14)) v = min(v, neonSeg(p, vec2( 0.0, -0.94), vec2( 0.0, -0.06)));
    if (bit(n, 15)) v = min(v, neonSeg(p, vec2( 0.06,-0.06), vec2( 0.44,-0.94)));
    if (bit(n, 16)) v = min(v, length(p - vec2(0.5, -1.0)) - 0.05);
    if (bit(n, 17)) v = min(v, neonSeg(p, vec2(0.5, -0.9), vec2(0.4, -1.2)));

    d_neon = min(d_neon, v);
    ch_pos.x += 1.5;
}

#define _A ddigit(0x119F);
#define _B ddigit(0x927E);
#define _C ddigit(0x007E);
#define _D ddigit(0x44E7);
#define _E ddigit(0x107E);
#define _F ddigit(0x101E);
#define _G ddigit(0x807E);
#define _H ddigit(0x1199);
#define _I ddigit(0x4466);
#define _J ddigit(0x4436);
#define _K ddigit(0x9218);
#define _L ddigit(0x0078);
#define _M ddigit(0x0A99);
#define _N ddigit(0x8899);
#define _O ddigit(0x00FF);
#define _P ddigit(0x111F);
#define _Q ddigit(0x80FF);
#define _R ddigit(0x911F);
#define _S ddigit(0x8866);
#define _T ddigit(0x4406);
#define _U ddigit(0x00F9);
#define _V ddigit(0x2218);
#define _W ddigit(0xA099);
#define _X ddigit(0xAA00);
#define _Y ddigit(0x4A00);
#define _Z ddigit(0x2266);
#define _0 ddigit(0x22FF);
#define _1 ddigit(0x0281);
#define _2 ddigit(0x1177);
#define _3 ddigit(0x11E7);
#define _4 ddigit(0x5508);
#define _5 ddigit(0x11EE);
#define _6 ddigit(0x11FE);
#define _7 ddigit(0x2206);
#define _8 ddigit(0x11FF);
#define _9 ddigit(0x11EF);
#define _EXCL  ddigit(0x4002); // !
#define _QUES  ddigit(0x5102); // ?
#define _DASH  ddigit(0x1100); // -
#define _PLUS  ddigit(0x5500); // +
#define _SPC ch_pos.x += 1.5;
#define _DOT     ddigit(0x10000);

void mainText( out vec4 fragColor, in vec2 fragCoord, float tTime )
{
    vec2 res = iResolution.xy;
    g_uv = (fragCoord - 0.8 * res) / res.y;
    g_uv *= 10.0;
    
    g_uv.y += sin(tTime * 1.5) * 0.2;
    g_uv.x += (tTime * 3.0)-20;
    
    d_neon = 1e6;
    ch_pos = ch_start;
    
    _H _E _R _E _SPC _W _E _SPC _G _O _SPC _A _G _A _I _N _DASH _A _SPC _N _E _W _SPC _L _I _G _H _T _A _M _P _SPC _R _E _L _E _A _S _E _DOT _DOT _DOT
    _SPC _W _E _SPC _N _O _W _SPC _H _A _V _E _SPC _A _SPC _S _I _M _P _L _E _SPC _T _R _A _C _K _DASH _E _D _I _T _O _R _DOT _DOT _DOT
    _SPC _J _U _S _T _SPC _F _O _R _SPC _F _U _N _DOT _SPC _A _N _D _SPC _S _T _U _F _F _SPC _I _SPC _F _O _R _G _O _T _DOT _DOT _DOT
    _SPC _G _R _E _E _T _I _N _X _SPC _T _O _SPC _A _L _L _SPC _B _A _R _D _S _SPC _O _U _T _SPC _T _H _E _R _E _DOT _DOT _DOT

    float core = smoothstep(0.04, 0.01, d_neon);
    vec3 neonColor = vec3(1.0, 0.05, 0.6); // Hot Pink
    neonColor = mix(neonColor, vec3(0.0, 0.9, 1.0), step(0.5, fract(tTime * 0.2)));
    float glow = 0.08 / (d_neon + 0.02);
    float flicker = 0.95 + 0.05 * sin(tTime * 40.0) * sin(tTime * 12.0);
    vec3 finalCol = (neonColor * glow * 1.5 + core * 1.2) * flicker;
    fragColor = vec4(finalCol, 1.0);
}


void main() {
    iTime = time /1000;
    iResolution = vec3(resolution, 0.0);
    vec4 scroltext_FragColor;

    if (time <= 12500)
        Discord(gl_FragColor, gl_FragCoord.xy);
    else
    {
        mainImage(gl_FragColor, gl_FragCoord.xy);
        if (time >= 17500)
            mainText(scroltext_FragColor, gl_FragCoord.xy, iTime - 17.5);
        gl_FragColor += scroltext_FragColor;
    }
}

#version 330

#ifdef GL_ES
precision mediump float;
#endif

uniform float time;
uniform vec2 mouse;
uniform vec2 resolution;

#define iTime time /1000
// A recreation of the extra little shoreline scene from stage 2 of 
// Gimmick! (Or in PAL regions, Mr. Gimmick).
//     
// Original game by Sunsoft: https://en.wikipedia.org/wiki/Gimmick!
// 
// Original graphics design: Hiroyuki Kagoya
// Shader graphics: Gerard Geer (https://github.com/gerard-geer)
// Original music composition: Masashi Kageyama
// Shader sound and ShaderTracker sound engine: Michael Moffitt (https://github.com/mikejmoffitt)
// This shader on github: https://github.com/gerard-geer/GimmickShader/

// Nah we don't need precision.


// Quality definitions. If your browser doesn't set a limit on compilation time, give
// these a shot! (Cough, cough, this shader works better on FF.)
// #define DRAW_ALL_BIRDS
// #define DRAW_WAVES

// A 2,4,8,16, or 32 element array implemented as a binary search.
int ARR2(in int x, in int a, in int b)
{
    if(x<1) return a;
    else return b;
}

vec4 ARR2(in int x, in vec4 a, in vec4 b)
{
    if(x<1) return a;
    else return b;
}

int ARR4(in int x, in int a, in int b, in int c, in int d)
{
    if(x<2) return ARR2(x,a,b);
    else return ARR2(x-2,c,d);
}

vec4 ARR4(in int x, in vec4 a, in vec4 b, in vec4 c, in vec4 d)
{
    if(x<2) return ARR2(x,a,b);
    else return ARR2(x-2,c,d);
}

int ARR8(in int x, in int a, in int b, in int c, in int d,
			       in int e, in int f, in int g, in int h)
{
    if(x<4) return ARR4(x, a,b,c,d);
    else return ARR4(x-4, e,f,g,h);
}

vec4 ARR8(in int x, in vec4 a, in vec4 b, in vec4 c, in vec4 d,
			   	    in vec4 e, in vec4 f, in vec4 g, in vec4 h)
{
    if(x<4) return ARR4(x, a,b,c,d);
    else return ARR4(x-4, e,f,g,h);
}

int ARR16(in int x, in int a, in int b, in int c, in int d,
			        in int e, in int f, in int g, in int h,
         			in int i, in int j, in int k, in int l,
         			in int m, in int n, in int o, in int p)
{
    if(x<8) return ARR8(x, a,b,c,d, e,f,g,h);
    else return ARR8(x-8, i,j,k,l, m,n,o,p);
}

vec4 ARR16(in int x, in vec4 a, in vec4 b, in vec4 c, in vec4 d,
			   	     in vec4 e, in vec4 f, in vec4 g, in vec4 h,
         			 in vec4 i, in vec4 j, in vec4 k, in vec4 l,
         			 in vec4 m, in vec4 n, in vec4 o, in vec4 p)
{
    if(x<8) return ARR8(x, a,b,c,d, e,f,g,h);
    else return ARR8(x-8, i,j,k,l, m,n,o,p);
}

int ARR32(in int _x, in int a, in int b, in int c, in int d, in int e, in int f, in int g, in int h,
          			in int i, in int j, in int k, in int l, in int m, in int n, in int o, in int p,
          			in int q, in int r, in int s, in int t, in int u, in int v, in int w, in int x,
          			in int y, in int z, in int aa,in int ab,in int ac,in int ad,in int ae,in int af)
{
    if(_x<16) return ARR16(_x, a,b,c,d, e,f,g,h, i,j,k,l, m,n,o,p);
    else return ARR16(_x-16, q,r,s,t, u,v,w,x, y,z,aa,ab, ac,ad,ae,af); 
}

vec4 ARR32(in int _x, in vec4 a, in vec4 b, in vec4 c, in vec4 d, in vec4 e, in vec4 f, in vec4 g, in vec4 h,
          			 in vec4 i, in vec4 j, in vec4 k, in vec4 l, in vec4 m, in vec4 n, in vec4 o, in vec4 p,
          			 in vec4 q, in vec4 r, in vec4 s, in vec4 t, in vec4 u, in vec4 v, in vec4 w, in vec4 x,
          			 in vec4 y, in vec4 z, in vec4 aa,in vec4 ab,in vec4 ac,in vec4 ad,in vec4 ae,in vec4 af)
{
    if(_x<16) return ARR16(_x, a,b,c,d, e,f,g,h, i,j,k,l, m,n,o,p);
    else return ARR16(_x-16, q,r,s,t, u,v,w,x, y,z,aa,ab, ac,ad,ae,af); 
}

// Constant color vectors so palette functions don't continually have to initialize new stuff.
const vec4 D_BLUE  = vec4(.235, .737, .988, 1.00);
const vec4 L_BLUE  = vec4(.659, .894, .988, 1.00);
const vec4 WHITE   = vec4(.988, .988, .988, 1.00);
const vec4 BLACK   = vec4(.000, .000, .000, 1.00);
const vec4 GRAY    = vec4(.455, .455, .455, 1.00);
const vec4 GRASS   = vec4(.502, .816, .063, 1.00);
const vec4 D_GREEN = vec4(.000, .235, .078, 1.00);
const vec4 L_GREEN = vec4(.298, .863, .282, 1.00);
const vec4 D_GOLD  = vec4(.486, .031, .000, 1.00);
const vec4 L_GOLD  = vec4(.988, .596, .219, 1.00);
const vec4 BROWN   = vec4(.486, .031, .000, 1.00);
const vec4 TRANS   = vec4(.000, .000, .000, .000);

// Define out stuff so we don't have to pass the values as parameters.
const int YUMETAROU_X = 52;
const int YUMETAROU_Y = 117;
const int SHORE_Y = 136;
const int SHORE_END = 79;
const int FAR_CLOUD_Y = 128;
const int WAVES_Y = 168;
const int BIRD_A_Y = 20;
const int BIRD_B_Y = 32;
const int BIRD_C_Y = 45;
const int BIRD_D_Y = 53;
const int BIRD_E_Y = 62;
const int BIRD_F_Y = 69;
const int BIRD_G_Y = 72;
const float BIRD_FLIP_FREQUENCY = .23438;

// The big cloud takes a lot of orchestration. These are the coordinates
// of the individual tiles.
// The cloud tiles represent only the detailed upper portions of it.
// Anything below them is drawn in as white.
const int CLOUD_A_X = 97;
const int CLOUD_A_Y = 160;
const int CLOUD_B_X = 105;
const int CLOUD_B_Y = 152;
const int CLOUD_C_X = 113;
const int CLOUD_C_Y = 153;
const int CLOUD_D_X = 129;
const int CLOUD_D_Y = 144;
const int CLOUD_E_X = 137;
const int CLOUD_E_Y = 136;
const int CLOUD_F_X = 145;
const int CLOUD_F_Y = 128;
const int CLOUD_G_X = 161;
const int CLOUD_G_Y = 128;
const int CLOUD_H_X = 169;
const int CLOUD_H_Y = 128;
const int CLOUD_I_X = 177;
const int CLOUD_I_Y = 136;
const int CLOUD_J_X = 185;
const int CLOUD_J_Y = 144;
const int CLOUD_K_X = 193;
const int CLOUD_K_Y = 153;
const int CLOUD_L_X = 201;
const int CLOUD_L_Y = 153;
const int CLOUD_M_X = 217;
const int CLOUD_M_Y = 154;
const int CLOUD_N_X = 225;
const int CLOUD_N_Y = 152;

// The positioning of the smaller cloud.
const int S_CLOUD_A_X = 184;
const int S_CLOUD_A_Y = 115;
const int S_CLOUD_B_X = 192;
const int S_CLOUD_B_Y = 112;
const int S_CLOUD_C_X = 216;
const int S_CLOUD_C_Y = 115;
    
/*
*	Yumetarou's palette.
*   
*	Returns a color given a palette index.
*
*	c: The color index to look up.
*
*	Returns: The corresponding color.
*/
vec4 yumetarouPalette(in int c)
{
    if(c < 4)
    {
        return ARR4(c,  WHITE,  	// The slightly not white white.
                        L_GREEN,	// Light green.
                        D_GREEN,	// Dark green.
                        L_GOLD); 	// Light gold.
    }
    else
    {
        c-=4;
        return ARR2(c, 	D_GOLD, 	// Dark gold.
                		TRANS);  	// Transparency.
    }
}

/*
*	Yumetarou's eyes-open sprite frame.
*   
*	Returns a palette index given a position.
*
*	x: The x position of the current fragment.
*	y: The y position of the current fragment.
*
*	Returns: The corresponding palette index.
*/
int yumetarouEyesOpen(in int x, in int y)
{
    if(y<16)
    {
        return ARR16(y,  ARR16(x,5,5,5,5,4,4,5,5,5,5,5,5,5,5,5,5),
                         ARR16(x,5,5,5,5,4,3,4,5,5,5,5,5,5,5,5,5),
                         ARR16(x,5,5,5,5,4,3,3,4,5,5,5,5,5,5,5,5),
                         ARR16(x,5,5,5,5,2,2,2,2,2,2,5,5,5,5,5,5),
                         ARR16(x,5,5,5,5,2,2,1,1,1,1,2,2,5,5,5,5),
                         ARR16(x,5,5,5,2,1,1,1,1,1,1,1,1,2,5,5,5),
                         ARR16(x,5,5,2,1,1,1,1,1,1,1,1,1,2,5,5,5),
                         ARR16(x,5,5,2,1,1,1,0,0,0,1,1,0,0,2,5,5),
                         ARR16(x,5,2,1,1,1,0,0,2,2,2,0,0,2,2,5,5),
                         ARR16(x,5,2,1,1,1,0,0,0,2,2,0,0,0,2,5,5),
                         ARR16(x,2,1,1,1,1,0,0,2,2,2,0,0,2,2,2,5),
                         ARR16(x,2,1,1,1,1,1,0,0,0,1,1,0,0,1,2,5),
                         ARR16(x,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2),
                         ARR16(x,2,1,1,2,1,1,1,1,1,1,1,1,1,1,1,2),
                         ARR16(x,5,2,2,1,1,1,1,1,1,2,2,2,1,1,1,2),
                         ARR16(x,5,2,1,1,1,1,1,1,1,1,1,1,1,1,1,2));
    } 
    else if(y==16) return 	 ARR16(x,5,2,1,1,1,1,1,1,1,1,1,1,1,1,2,5);
    else if(y==17) return 	 ARR16(x,5,5,2,2,1,1,1,1,1,1,1,1,1,2,2,5);
    else		   return 	 ARR16(x,1,0,0,0,2,2,2,2,2,2,2,2,2,0,0,0);
}

/*
*	Yumetarou's eyes-closed sprite frame.
*   
*	Returns a palette index given a position.
*
*	x: The x position of the current fragment.
*	y: The y position of the current fragment.
*
*	Returns: The corresponding palette index.
*/
int yumetarouEyesClosed(in int x, in int y)
{
    if(y<16){
        return ARR16(y,  ARR16(x,5,5,5,5,4,4,5,5,5,5,5,5,5,5,5,5),
                         ARR16(x,5,5,5,5,4,3,4,5,5,5,5,5,5,5,5,5),
                         ARR16(x,5,5,5,5,4,3,3,4,5,5,5,5,5,5,5,5),
                         ARR16(x,5,5,5,5,2,2,2,2,2,2,5,5,5,5,5,5),
                         ARR16(x,5,5,5,5,2,2,1,1,1,1,2,2,5,5,5,5),
                         ARR16(x,5,5,5,2,1,1,1,1,1,1,1,1,2,5,5,5),
                         ARR16(x,5,5,2,1,1,1,1,1,1,1,1,1,2,5,5,5),
						 ARR16(x,5,5,2,1,1,1,0,0,0,1,1,0,0,2,5,5),
						 ARR16(x,5,2,1,1,1,0,0,0,0,0,0,0,0,0,5,5),
						 ARR16(x,5,2,1,1,1,0,2,2,2,2,0,2,2,2,5,5),
						 ARR16(x,2,1,1,1,1,0,0,0,0,0,0,0,0,0,2,5),
						 ARR16(x,2,1,1,1,1,1,0,0,0,1,1,0,0,1,2,5),
                         ARR16(x,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2),
                         ARR16(x,2,1,1,2,1,1,1,1,1,1,1,1,1,1,1,2),
                         ARR16(x,5,2,2,1,1,1,1,1,1,2,2,2,1,1,1,2),
                         ARR16(x,5,2,1,1,1,1,1,1,1,1,1,1,1,1,1,2));
    } 
    else if(y==16) return 	 ARR16(x,5,2,1,1,1,1,1,1,1,1,1,1,1,1,2,5);
    else if(y==17) return 	 ARR16(x,5,5,2,2,1,1,1,1,1,1,1,1,1,2,2,5);
    else		   return 	 ARR16(x,1,0,0,0,2,2,2,2,2,2,2,2,2,0,0,0);
}

/*
*	Yumetarou's draw function.
*   
*	Draws Yumetarou to the screen.
*
*	x: The x position of the current fragment.
*	y: The y position of the current fragment.
*
*	Returns: The color of Yumetarou from under the current texel.
*/
vec4 drawYumetarou(in int x, in int y)
{
    if(x < YUMETAROU_X || x > YUMETAROU_X + 15) return TRANS;
    else if(y < YUMETAROU_Y || y > YUMETAROU_Y + 18) return TRANS;
    else
    {
        x -= YUMETAROU_X;
        y -= YUMETAROU_Y;

        // Yummy yummy frame counting.
        float t = mod(iTime, 3.67);
        if( t > .066 && (t < .533 || t >.600) )
            return yumetarouPalette(yumetarouEyesOpen(x,y));
        else
            return yumetarouPalette(yumetarouEyesClosed(x,y));
	}
}

/*
*	The birds' palette.
*   
*	Returns a color given a palette index.
*
*	c: The color index to look up.
*
*	Returns: The corresponding color.
*/
vec4 birdPalette(in int c)
{
    return ARR4(c,  WHITE,  // The slightly not white white.
                	GRAY,  	// Gray
                	BLACK,  // Black
                    TRANS); // Transparency.
}

/*
*	A bird's wing-level tile frame.
*   
*	Returns a palette index given a position.
*
*	x: The x position of the current fragment.
*	y: The y position of the current fragment.
*
*	Returns: The corresponding palette index.
*/
int birdWingsLevel(in int x,in int y)
{
    if(y<4){
		return ARR4(y,	 3,
						 ARR8(x,3,0,0,0,1,0,0,3),
						 ARR8(x,2,3,3,0,0,3,3,2),
						 ARR8(x,3,3,3,1,0,3,3,3));
    }
    else return 3;
}

/*
*	The frame of the bird with its wings up.
*   
*	Returns a palette index given a position.
*
*	x: The x position of the current fragment.
*	y: The y position of the current fragment.
*
*	Returns: The corresponding palette index.
*/
int birdWingsUp(in int x,in int y)
{
    if(y<4){
		return ARR4(y,	 ARR8(x,3,2,0,3,3,3,2,3),
						 ARR8(x,3,3,0,0,3,0,3,3),
						 ARR8(x,3,3,1,0,1,0,3,3),
						 ARR8(x,3,3,3,0,0,3,3,3));
    }
    else return		 ARR8(x,3,3,3,1,0,3,3,3);
}

/*
*	The frame of the bird with its wings down.
*   
*	Returns a palette index given a position.
*
*	x: The x position of the current fragment.
*	y: The y position of the current fragment.
*
*	Returns: The corresponding palette index.
*/
int birdWingsDown(in int x,in int y)
{
    if(y<4){
		return ARR4(y,	 ARR8(x,3,3,3,0,1,3,3,3),
						 ARR8(x,3,3,3,0,0,0,3,3),
						 ARR8(x,3,3,0,1,0,1,3,3),
						 ARR8(x,3,3,0,3,3,3,0,3));
    }
    else return		 ARR8(x,3,3,2,3,3,3,2,3);
}

/*
*	The bird draw function.
*   
*	Draws a single bird to the screen.
*
*	x: The x position of the current fragment.
*	y: The y position of the current fragment.
*	atx: The x position at which to draw the bird.
*	aty: The y position at which to draw the bird.
*	flip: Whether or not to flip the bird. (along the x axis.)
*
*	Returns: The color of the bird from under the current texel.
*/
vec4 drawBird(in int x, in int y, in int atx, in int aty, bool flip)
{
    // Bounds checking.
    if(x < atx || x > atx + 7) return TRANS;
    if(y < aty || y > aty + 4) return TRANS;
    
    // Transform coordinates to bird space.
    x -= atx;
    y -= aty;
    
    // Flip the bird if necessary.
    if(flip) x = 7-x;
    
    // This animation is less framecounting and more dividing an amount
    // of time by four.
    float t = mod(iTime, .533);
    if(t < .133)	return birdPalette(birdWingsLevel(x,y));
    else if(t < .266)	return birdPalette(birdWingsUp(x,y));
    else if(t < .400)	return birdPalette(birdWingsLevel(x,y));
    else return birdPalette(birdWingsDown(x,y));
}



/*
*	The birds' animation function.
*   
*	Returns a modulated value by adding a triangle wave to the
*	starting value s.
*		
*	s: The starting position.
*	t: The current time within the function.
*	a: The amplitude of the triangle wave.
*	d: The boolean first derivative of the triangle function.
*
*	Returns: The modulated position.
*/
int anim(in int s, in float t, in float a, out bool d)
{
    // Triangle wave = |saw wave|
    
    // Let's get the derivative first.
    d = 2.0*(mod((t+1.0)*BIRD_FLIP_FREQUENCY, 1.0))-1.0 < 0.0;
    
    // Now that we've stored the direction let's go back and 
    // calculate the position.
    float val = abs( (mod((t)*BIRD_FLIP_FREQUENCY, 1.0)*2.0)-1.0 )*a;
    
    // Return the animated position.
	return s + int(val);
}

/*
*	The whole flock's draw function.
*
*	Draws all the birds to the screen, animated.
*
*	x: The x position of the current fragment.
*	y: The y position of the current fragment.
*	
*	Returns: The color of the birds from under the current fragment.
*/
vec4 drawBirds(in int x, in int y)
{
    // Since birds never cross we can use additive blending.
    // And as we've learned from the sound let's divvy up addition.
	
	// Getting the positioning and timing accurate to the actual game
	// was not fun. Frame-counting and screen-shooting dominated an
	// evening of mine. Should have bit the bullet and looked at a
	// disassembly.
	// Each bird's flight path lasts 128 frames each way. However
	// those path start times differ, as well as the length of the path.
    
    bool f; // For directional awareness.
	
	// Bird 1.
    int a = anim(110, iTime, 32.0, f);
    vec4 result = drawBird(x,y,a,BIRD_A_Y,f);
	
    #ifdef DRAW_ALL_BIRDS
	// Bird 2.
    a = anim(140, iTime+3.267, 24.0, f);
    result += drawBird(x,y,a,BIRD_B_Y,f);
	
	// Bird 3.
    a = anim(77, iTime+1.533, 40.0, f);
    result += drawBird(x,y,a,BIRD_C_Y,f);
	
	// Bird 4.
    a = anim(198, iTime+.1667, 32.0, f);
    result += drawBird(x,y,a,BIRD_D_Y,f);
	
	// Bird 5.
    a = anim(141, iTime+.5667, 32.0, f);
    result += drawBird(x,y,a,BIRD_E_Y,f);
	
	// Bird 6.
    a = anim(85, iTime+1.067, 24.0, f);
    result += drawBird(x,y,a,BIRD_F_Y,f);
	#endif
    
	// Bird 7.
    a = anim(165, iTime+1.167, 24.0, f);
    result += drawBird(x,y,a,BIRD_G_Y,f);
    return result;
    
}

/*
*	The rocky shoreline's palette.
*   
*	Returns a color given a palette index.
*
*	c: The color index to look up.
*
*	Returns: The corresponding color.
*/
vec4 shorePalette(in int x)
{
    return ARR4(x, WHITE,
                   GRASS,
                   GRAY,
                   BLACK);
}

/*
*	The repeated interior portion of the shore.
*   
*	Returns a palette index given a position.
*
*	x: The x position of the current fragment.
*	y: The y position of the current fragment.
*
*	Returns: The corresponding palette index.
*/
int shoreInterior(in int x, in int y)
{
    // This element repeats in the X direction.
    x = int(mod(float(x),32.0));
    return ARR32(y, 
             3,
			 0,
			 1,
			 1,
			 1,
			 1,
			 1,
			 1,
			 ARR32(x,1,1,1,1,1,1,2,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),
			 ARR32(x,1,1,1,1,1,2,2,2,1,1,1,1,2,2,1,1,1,1,1,1,1,1,1,1,2,2,2,1,1,1,1,1),
			 ARR32(x,1,1,1,1,1,2,2,3,2,1,2,3,3,2,2,1,1,1,1,1,1,1,2,2,3,3,3,2,1,1,1,1),
			 ARR32(x,1,1,1,1,2,2,3,3,3,3,3,3,3,3,2,1,1,1,1,1,1,2,2,3,2,2,2,2,2,1,1,1),
			 ARR32(x,1,1,1,1,2,2,3,3,2,0,0,0,0,2,2,1,1,1,1,1,2,2,2,0,0,0,0,0,2,2,1,1),
			 ARR32(x,1,1,2,2,2,2,3,0,0,0,0,0,0,0,2,2,1,1,1,2,2,0,0,0,0,0,0,0,0,2,2,1),
			 ARR32(x,2,2,2,2,2,3,0,0,0,0,0,0,0,0,2,3,2,2,2,2,0,0,0,0,0,0,0,0,0,2,3,2),
			 ARR32(x,2,2,2,2,3,3,0,0,0,0,0,0,0,0,2,2,3,2,2,0,0,0,0,0,0,0,0,0,0,2,3,2),
			 ARR32(x,2,2,2,3,3,2,0,0,0,0,0,0,0,0,2,2,3,2,0,0,0,0,0,0,0,0,0,0,0,2,3,3),
			 ARR32(x,3,3,3,3,3,2,0,0,0,0,0,0,0,2,2,2,3,2,0,0,0,0,0,0,0,0,0,0,0,2,3,3),
			 ARR32(x,2,2,2,3,3,2,0,0,0,0,0,0,2,2,2,3,3,2,0,0,0,0,0,0,0,0,0,0,2,2,3,3),
			 ARR32(x,0,0,0,0,2,2,2,0,0,0,0,2,2,2,3,3,3,2,2,0,0,0,0,0,0,0,0,2,2,2,3,2),
			 ARR32(x,0,0,0,0,0,2,2,2,2,2,2,2,2,2,2,0,0,0,2,2,0,0,0,0,0,0,2,2,2,3,2,0),
			 ARR32(x,0,0,0,0,0,0,2,3,3,2,2,2,2,2,0,0,0,0,0,2,3,2,2,2,2,2,2,2,2,3,0,0),
			 ARR32(x,0,0,0,0,0,0,0,2,3,3,3,3,2,0,0,0,0,0,0,0,2,3,2,2,2,2,2,2,3,2,0,0),
			 ARR32(x,0,0,0,0,0,0,0,0,2,3,3,2,0,0,0,0,0,0,0,0,2,3,3,2,2,2,3,3,3,2,0,0),
			 ARR32(x,0,0,0,0,0,0,0,0,2,3,2,0,0,0,0,0,0,0,0,0,2,3,0,0,0,0,2,3,3,2,0,0),
			 ARR32(x,0,0,0,0,0,0,0,0,2,3,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,2,3,2,0,0),
			 ARR32(x,0,0,0,0,0,0,0,2,2,3,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,2,3,0,0),
			 ARR32(x,0,0,0,0,0,0,2,2,3,2,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,2,3,2,0),
			 ARR32(x,0,0,0,0,0,2,2,2,3,2,0,0,0,0,0,0,0,0,0,2,2,0,0,0,0,0,0,0,2,3,2,2),
			 ARR32(x,2,0,0,0,2,2,2,2,3,2,0,0,0,0,0,0,0,2,3,3,2,0,0,0,0,0,0,0,2,3,2,2),
			 ARR32(x,2,2,2,2,2,2,2,3,3,2,2,0,0,0,0,0,2,3,2,2,2,2,0,0,0,0,0,0,2,3,2,2),
			 ARR32(x,2,2,2,2,2,2,3,3,3,2,2,2,2,2,2,2,2,0,0,0,0,0,2,0,0,0,0,2,2,3,3,2));
}

/*
*	The non-repeated exterior portion of the shore.
*   
*	Returns a palette index given a position.
*
*	x: The x position of the current fragment.
*	y: The y position of the current fragment.
*
*	Returns: The corresponding palette index.
*/
int shoreExterior(in int x, in int y)
{
    return ARR32(y,
            ARR16(x,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,2),
            ARR16(x,0,0,0,0,0,0,0,0,0,0,0,0,0,1,2,3),
            ARR16(x,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,3),
            ARR16(x,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,3),
            ARR16(x,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,3),
            ARR16(x,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,3),
            ARR16(x,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,3),
            ARR16(x,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,3),
            ARR16(x,1,1,1,1,1,1,2,2,1,1,1,1,2,2,2,3),
            ARR16(x,1,1,1,1,1,2,2,2,2,2,2,2,2,2,3,2),
            ARR16(x,1,1,1,1,1,2,2,3,2,2,2,2,3,3,3,2),
            ARR16(x,1,1,1,1,2,2,3,3,3,3,2,2,2,3,3,2),
            ARR16(x,1,1,1,1,2,2,3,3,0,0,0,0,0,2,3,2),
            ARR16(x,1,1,2,2,2,2,3,0,0,0,0,0,0,0,2,3),
            ARR16(x,2,2,2,2,2,3,0,0,0,0,0,0,0,0,2,3),
            ARR16(x,2,2,2,2,3,3,0,0,0,0,0,0,0,0,2,3),
            ARR16(x,2,2,2,3,3,2,0,0,0,0,0,0,0,0,2,3),
            ARR16(x,3,3,3,3,3,2,0,0,0,0,0,0,0,0,2,3),
            ARR16(x,2,2,2,3,3,2,0,0,0,0,0,0,0,2,2,3),
            ARR16(x,0,0,0,0,2,2,2,0,0,0,0,0,2,2,3,2),
            ARR16(x,0,0,0,0,0,2,2,2,0,0,0,2,2,3,3,2),
            ARR16(x,0,0,0,0,0,0,2,3,2,2,2,2,0,0,3,2),
            ARR16(x,0,0,0,0,0,0,0,2,3,3,2,0,0,0,2,3),
            ARR16(x,0,0,0,0,0,0,0,0,3,2,0,0,0,0,2,3),
            ARR16(x,0,0,0,0,0,0,0,0,3,0,0,0,0,0,2,3),
            ARR16(x,0,0,0,0,0,0,0,0,2,0,0,0,0,0,2,3),
            ARR16(x,0,0,0,0,0,0,0,2,2,0,0,0,0,0,2,3),
            ARR16(x,0,0,0,0,0,0,2,2,2,0,0,0,0,0,2,3),
            ARR16(x,0,0,0,0,0,2,2,2,3,0,0,0,0,0,2,3),
            ARR16(x,2,0,0,0,2,2,2,2,3,0,0,0,0,2,3,3),
            ARR16(x,2,2,2,2,2,2,2,3,3,0,0,0,2,2,2,3),
            ARR16(x,2,2,2,2,2,2,3,3,3,2,2,2,0,0,2,3));
}

/*
*	The shoreline's draw function.
*
*	Draws the two interior segments of the shore, then the endcap.
*
*	x: The x position of the current fragment.
*	y: The y position of the current fragment.
*	
*	Returns: The color of the shore from under the current fragment.
*/
vec4 drawShore(in int x, in int y)
{
    // Bounds checking.
    if(x > SHORE_END) return TRANS;
    else if(y < SHORE_Y || y > SHORE_Y + 31) return TRANS;
    else
    {
        // Transform to be relative to the shore tiles.
        y -= SHORE_Y;

        // Draw the interior of the shore.
        if(x < 64) return shorePalette(shoreInterior(x,y));
        // Draw the endcap exterior.
        else
        {
            x -= 64;
            return shorePalette(shoreExterior(x,y));
        }
    }
}

/*
*	The palette of those distant clouds.
*   
*	Returns a color given a palette index.
*
*	c: The color index to look up.
*
*	Returns: The corresponding color.
*/
vec4 farCloudsPalette(in int x)
{
    return ARR2(x, TRANS,
                   L_BLUE);
}

/*
*	The tile function of those clouds in the distance.
*   
*	Returns a palette index given a position.
*
*	x: The x position of the current fragment.
*	y: The y position of the current fragment.
*
*	Returns: The corresponding palette index.
*/
int farClouds(in int x, in int y)
{
    // The clouds repeat along the X axis across the entire screen.
    x = int(mod(float(x),32.0));
    if(y < 4)
    {
        return ARR4(y, 
        ARR32(x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0),
        ARR32(x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,0,0),
        ARR32(x,0,0,0,1,1,0,0,0,0,0,1,1,1,1,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,0,0,0),
        ARR32(x,1,1,1,1,0,1,1,0,1,1,1,1,1,1,1,1,0,0,0,1,1,1,1,1,1,1,1,1,1,1,0,0));
    }
    else return ARR32(x,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1);
}

/*
*	The distant cloud draw function.
*   
*	Draws the dark, distant clouds to the screen.
*
*	x: The x position of the current fragment.
*	y: The y position of the current fragment.
*
*	Returns: The color of the far clouds at the given position.
*/
vec4 drawFarClouds(in int x, in int y)
{
    // Above? Nada.
    if(y < FAR_CLOUD_Y) return TRANS;
    // Below? Fill'er'in.
    else if(y > FAR_CLOUD_Y+5) return L_BLUE;
    // Within the narrow band designated for the clouds?
    else return farCloudsPalette(farClouds(x,y-FAR_CLOUD_Y));
}

/*
*	The palette of the waves when under the shoreline.
*   
*	Returns a color given a palette index.
*
*	c: The color index to look up.
*
*	Returns: The corresponding color.
*/
vec4 wavesShadowPalette(in int x)
{
    if(x<4)
    {
        return ARR4(x,  D_BLUE,
			   			WHITE,
			   			L_BLUE,
			   			WHITE);
    }
    else return ARR2(x-4, D_BLUE, L_BLUE);
}

/*
*	The palette of the waves in the sun.
*   
*	Returns a color given a palette index.
*
*	c: The color index to look up.
*
*	Returns: The corresponding color.
*/
vec4 wavesSunnyPalette(in int x)
{
    if(x<4)
    {
        return ARR4(x, L_BLUE,
					   L_BLUE,
					   L_BLUE,
					   WHITE);
    }
    else return ARR2(x-4, WHITE, WHITE);
}

/*
*	One frame of the waves.
*	Note: The palette of the sunny and shadowed waves are
*	consolidated into a single map with a larger palette.
*   
*	Returns a palette index given a position.
*
*	x: The x position of the current fragment.
*	y: The y position of the current fragment.
*
*	Returns: The corresponding palette index.
*/
int wavesA(in int x, in int y)
{
    if(x < 32) // ARR64 would be a really long line.
    {
        return ARR8(y,
       	ARR32(x,3,3,3,3,3,3,3,3,3,5,5,5,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,2,3,3,3),
		ARR32(x,3,3,3,2,2,0,0,0,0,2,5,5,5,5,5,5,3,3,3,3,3,5,5,2,2,2,0,0,0,0,0,0),
		ARR32(x,3,2,0,0,0,0,0,0,0,0,0,0,2,3,3,3,3,2,0,0,0,0,5,5,5,5,5,5,2,2,2,2),
		ARR32(x,0,0,2,2,2,2,2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,5,5,5,5,5,5),
		ARR32(x,2,2,2,2,2,2,2,2,2,2,2,2,2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,5,5),
		ARR32(x,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,0,0,0,0,0,0,0,0,0,0,0,0),
		2,
        2);
    }
    else
    {
        x -= 32;
        return ARR8(y,
		ARR32(x,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,2,3,3,3,0,0,0,0,0,0,0,0,0,3,3,3),
        ARR32(x,0,0,0,0,2,5,5,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3,3,3,3,3,3,3,3),
        ARR32(x,2,2,1,3,3,3,2,2,2,5,0,0,0,0,0,0,0,0,5,5,3,3,3,3,3,3,3,3,2,3,3,3),
        ARR32(x,3,3,3,2,0,0,0,0,0,0,0,0,0,0,5,5,5,3,3,3,3,3,3,2,2,2,3,3,3,3,2,0),
        ARR32(x,5,5,5,5,5,5,2,0,0,0,0,5,5,3,3,3,3,3,3,3,2,2,2,3,3,3,3,2,0,0,0,2),
        ARR32(x,0,0,0,0,0,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,2,2,2,2,2,2,2,2,2),
        2,
        2);
    }
}

/*
*	Another frame of the waves.
*	Note: The palette of the sunny and shadowed waves are
*	consolidated into a single map with a larger palette.
*   
*	Returns a palette index given a position.
*
*	x: The x position of the current fragment.
*	y: The y position of the current fragment.
*
*	Returns: The corresponding palette index.
*/
int wavesB(in int x, in int y)
{
    if(x < 32) // ARR64 would be a really long line.
    {
        return ARR8(y,
		ARR32(x,2,2,2,2,2,0,0,0,0,0,0,0,0,2,2,3,3,3,3,3,3,3,3,3,5,5,2,0,0,0,0,0),
        ARR32(x,0,0,0,0,0,0,0,2,2,2,3,3,3,3,3,3,3,3,3,3,3,0,0,0,0,0,5,5,2,2,2,0),
        ARR32(x,0,0,0,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,0,0,0,0,0,0,0,0,0,0,4,4,4,2),
        ARR32(x,2,2,3,3,3,3,3,3,3,3,3,2,2,3,3,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4),
        ARR32(x,3,3,3,3,3,3,2,2,2,2,3,3,3,2,0,0,2,2,2,2,2,2,2,2,2,2,2,0,0,0,0,0),
        ARR32(x,3,3,3,2,2,2,3,3,3,3,0,0,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2),
		2,
        2);
	}
    else
    {
        x -= 32;
        return ARR8(y,
        ARR32(x,0,0,2,2,2,3,3,3,3,3,3,3,5,5,5,5,5,5,5,2,2,2,2,0,0,0,0,0,0,0,0,0),
        ARR32(x,0,0,0,0,0,0,0,2,2,3,3,5,5,5,5,2,0,0,0,0,0,0,0,0,0,0,2,2,2,2,2,2),
        ARR32(x,2,2,2,2,2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,2,2,3,3,3,3,2,3,2,0,0,0,0),
        ARR32(x,4,4,4,4,4,2,2,2,2,2,2,0,0,0,0,0,3,3,3,3,3,3,0,0,0,0,0,0,0,0,0,2),
        ARR32(x,0,0,0,4,4,4,4,4,5,5,5,5,5,5,2,2,0,0,0,0,0,0,0,0,0,0,0,2,3,3,3,3),
        ARR32(x,0,0,0,0,0,0,0,0,0,0,0,0,0,5,5,5,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3),
        ARR32(x,2,2,2,2,2,2,2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3,3,3,3,3,3,3,3,3,2),
        2);
    }
}

/*
*	A third frame of sweet wave action.
*	Note: The palette of the sunny and shadowed waves are
*	consolidated into a single map with a larger palette.
*   
*	Returns a palette index given a position.
*
*	x: The x position of the current fragment.
*	y: The y position of the current fragment.
*
*	Returns: The corresponding palette index.
*/
int wavesC(in int x, in int y)
{
    if(x < 32) // ARR64 would be a really long line.
    {
        return ARR8(y,
        ARR32(x,3,3,3,3,0,0,0,0,0,0,0,3,3,3,0,0,0,0,0,0,2,3,3,3,2,0,0,0,0,0,0,0),
        ARR32(x,2,0,0,0,5,3,3,2,0,0,0,0,0,0,0,0,2,2,3,3,3,2,2,2,3,2,2,0,0,3,3,3),
        ARR32(x,0,0,0,0,0,0,0,4,4,4,2,2,2,2,3,3,3,3,3,5,0,0,0,0,0,0,5,5,3,3,3,3),
        ARR32(x,0,0,0,0,0,0,0,0,0,4,4,4,5,5,5,5,5,5,0,0,0,0,0,5,5,5,3,3,3,1,3,4),
        ARR32(x,2,2,2,2,2,0,0,0,0,0,0,0,4,4,5,5,5,5,5,5,5,5,5,5,3,3,1,1,3,3,4,0),
        ARR32(x,2,2,2,2,2,2,2,2,2,2,2,0,0,0,0,0,2,3,3,3,3,3,3,3,3,3,3,4,0,0,0,0),
        ARR32(x,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,0,0,0,0,0,2,2,2),
        2);

    }
    else
    {
        x -= 32;
        return ARR8(y,
		ARR32(x,0,2,2,2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4,4,5,5,5,3,3,3,3,3),
		ARR32(x,3,3,3,3,5,2,2,2,0,0,0,0,0,0,4,4,4,4,4,4,5,5,5,5,5,3,3,3,3,3,3,3),
		ARR32(x,3,4,4,4,4,5,5,5,5,5,5,5,4,4,4,4,5,5,5,5,5,5,5,1,1,1,3,3,3,2,0,0),
		ARR32(x,4,4,0,0,0,0,0,0,0,0,2,5,5,5,5,5,3,1,2,2,2,1,3,3,3,3,4,0,0,0,2,2),
		ARR32(x,0,2,2,2,2,0,0,0,0,0,0,0,0,0,2,2,2,2,1,3,3,5,5,0,0,0,2,2,2,2,2,2),
		ARR32(x,2,2,2,2,2,2,2,2,2,2,2,2,0,0,0,0,0,0,0,0,0,2,2,2,2,2,2,2,2,2,2,2),
		2,
		2);
    }
}

/*
*	The forth frame of waves.
*	Note: The palette of the sunny and shadowed waves are
*	consolidated into a single map with a larger palette.
*   
*	Returns a palette index given a position.
*
*	x: The x position of the current fragment.
*	y: The y position of the current fragment.
*
*	Returns: The corresponding palette index.
*/
int wavesD(in int x, in int y)
{
    if(x < 32) // ARR64 would be a really long line.
    {
        return ARR8(y,
		ARR32(x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,2,2,3,3,2,0,0,0,0,0,0,0,0,0,2),
		ARR32(x,2,2,2,2,2,3,3,3,3,1,2,2,0,0,0,0,0,0,0,0,0,0,0,0,2,2,2,5,3,3,3,3),
		ARR32(x,3,3,3,3,3,1,0,0,0,0,0,5,5,2,2,2,2,0,0,0,0,2,5,5,5,3,3,3,3,3,3,3),
		ARR32(x,5,5,5,2,0,0,0,0,0,0,0,0,0,0,5,5,5,5,5,5,5,5,3,3,3,3,3,3,3,3,5,0),
		ARR32(x,0,4,4,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,3,3,3,3,3,3,3,3,5,0,0,0,0,0),
		ARR32(x,0,0,0,0,0,0,0,0,0,0,0,4,3,3,3,3,3,3,3,3,0,0,0,0,0,0,2,2,2,2,2,2),
		ARR32(x,0,0,0,0,0,0,0,0,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2),
		2);
    }
    else
    {
        x -= 32;
        return ARR8(y,
		ARR32(x,5,5,3,3,3,3,3,3,3,3,3,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,0,0,0,0,0),
		ARR32(x,3,3,3,5,0,0,0,0,0,0,0,4,4,4,4,5,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2),
		ARR32(x,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5),
		ARR32(x,0,0,0,0,0,2,2,2,2,2,2,2,0,0,0,0,0,0,0,0,0,0,0,0,4,4,4,4,4,4,4,5),
		ARR32(x,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
		ARR32(x,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,0,0,0,0,0,0),
		2,
		2);
    }
}

/*
*	The wave draw function.
*   
*	Draws the waves using the appropriate palette given the position.
*	This also animates the waves.
*
*	x: The x position of the current fragment.
*	y: The y position of the current fragment.
*
*	Returns: The color of the waves at the given position.
*/
vec4 drawWaves(in int x, in int y)
{
    // Bounds checking.
    if(y < WAVES_Y) return TRANS;
    else if(y > WAVES_Y+7) return L_BLUE;
    else
    {
        // Modulo the time and cast it to an int so the value returned
        // can be used as an index for which frame of animation to use.
        int t = int(mod(iTime*6.0,4.0));

        // We need to do the usual transform here as well.
        y -= WAVES_Y;
        
		// If we are under the shoreline, we need to use the palette
        // that reflects the shore.
        if(x > SHORE_END)
        {
            // The prior comparison required x to be pristine, so
            // we have to perform this modulo in here.
            x = int(mod(float(x),64.0));
            return ARR4(t,
                        wavesSunnyPalette(wavesC(x,y)),
                        wavesSunnyPalette(wavesA(x,y)),
                        wavesSunnyPalette(wavesB(x,y)),
                        wavesSunnyPalette(wavesD(x,y)));
        }
        // otherwise we use the palette that reflects the clouds.
        else
        {
            x = int(mod(float(x),64.0));
            return ARR4(t,
                        wavesShadowPalette(wavesC(x,y)),
                        wavesShadowPalette(wavesA(x,y)),
                        wavesShadowPalette(wavesB(x,y)),
                        wavesShadowPalette(wavesD(x,y)));
        }
    }
}

/*
*	The palette of the white clouds.
*   
*	Returns a color given a palette index.
*
*	c: The color index to look up.
*
*	Returns: The corresponding color.
*/
vec4 nearCloudsPalette(in int x)
{
	return ARR2(x, TRANS, WHITE);
}

/*
*	Cloud tile functions.
*	
*	What follows are the tile functions for the large white cloud.
*	Only the topmost sections with actual features are encoded.
*	The solid white interior is assumed.
*   
*	Returns a palette index given a position.
*
*	x: The x position of the current fragment.
*	y: The y position of the current fragment.
*
*	Returns: The corresponding palette index.
*/
int cloudA(in int x, in int y)
{
    // Do some bounds checking.
    // To the left or right? TRANSPARENT FOR YOU!
	if(x < CLOUD_A_X || x >= CLOUD_B_X) return 0;
    // Above this cloud tile? TRANSPARENT YOU AS WELL!
	else if(y < CLOUD_A_Y) return 0;
    // Below the tile? OH YOU ARE MORE CLOUD HAVE CLOUD COLOR.
	else if(y > CLOUD_A_Y+7) return 1;
	
    else
    {
        // Transform the coordinates to cloud space.
        x -= CLOUD_A_X;
        y -= CLOUD_A_Y;

        // Finally do the 2D binary lookup to get the actual color.
        return
        ARR8(y,
          0,
          ARR8(x,0,0,0,0,0,0,1,1),
          ARR8(x,0,0,0,0,1,1,1,1),
          ARR8(x,0,0,0,1,1,1,1,1),
          ARR8(x,0,0,1,1,1,1,1,1),
          ARR8(x,0,0,1,1,1,1,1,1),
          ARR8(x,0,1,1,1,1,1,1,1),
          ARR8(x,0,1,1,1,1,1,1,1)
        );
    }
}
// Cloud tile B.
int cloudB(in int x, in int y)
{
	if(x < CLOUD_B_X || x >= CLOUD_C_X) return 0;
	else if(y < CLOUD_B_Y) return 0;
	else if(y > CLOUD_B_Y+7) return 1;
	else
    {
        x -= CLOUD_B_X;
        y -= CLOUD_B_Y;

        return
        ARR8(y,
          0,
          0,
          ARR8(x,0,0,0,0,0,1,1,1),
          ARR8(x,0,1,1,1,1,0,1,0),
          ARR8(x,1,1,1,1,0,0,1,1),
          ARR8(x,1,1,1,0,0,1,1,1),
          ARR8(x,0,1,1,0,1,1,1,1),
          ARR8(x,1,0,0,0,1,1,1,1)
        );
    }
}
// Cloud tile C.
int cloudC(in int x, in int y)
{
	if(x < CLOUD_C_X || x >= CLOUD_D_X) return 0;
	else if(y < CLOUD_C_Y) return 0;
	else if(y > CLOUD_C_Y+1) return 1;
	else
    {
        x -= CLOUD_C_X;
        y -= CLOUD_C_Y;

        return
        ARR2(y,
          ARR16(x,0,0,1,0,0,0,1,1,1,1,0,0,0,0,0,0),
          ARR16(x,1,1,1,1,0,1,1,1,1,1,1,0,0,1,1,0)
        );
    }
}
// Cloud tile D.
int cloudD(in int x, in int y)
{
	if(x < CLOUD_D_X || x >= CLOUD_E_X) return 0;
	else if(y < CLOUD_D_Y) return 0;
	else if(y > CLOUD_D_Y+7) return 1;
	else
    {
        x -= CLOUD_D_X;
        y -= CLOUD_D_Y;

        return
        ARR8(y,
          0,
          ARR8(x,0,0,0,0,0,0,1,1),
          ARR8(x,0,0,0,0,1,1,1,1),
          ARR8(x,0,0,0,1,1,1,1,1),
          ARR8(x,0,0,1,1,1,1,1,1),
          ARR8(x,0,0,1,1,1,1,1,1),
          ARR8(x,0,1,1,1,1,1,1,1),
          ARR8(x,0,1,1,1,1,1,1,1)
        );
    }
}
// Cloud tile E.
int cloudE(in int x, in int y)
{
	if(x < CLOUD_E_X || x >= CLOUD_F_X) return 0;
	else if(y < CLOUD_E_Y) return 0;
	else if(y > CLOUD_E_Y+7) return 1;
	else
    {
        x -= CLOUD_E_X;
        y -= CLOUD_E_Y;

        return
        ARR8(y,
          0,
          0,
          ARR8(x,0,0,0,0,0,1,1,1),
          ARR8(x,0,1,1,1,1,0,1,0),
          ARR8(x,1,1,1,1,0,0,1,1),
          ARR8(x,1,1,1,0,0,1,1,1),
          ARR8(x,0,1,1,0,1,1,1,1),
          ARR8(x,1,0,0,0,1,1,1,1)
        );
    }
}
// Cloud tile F.
int cloudF(in int x, in int y)
{
	if(x < CLOUD_F_X || x >= CLOUD_G_X) return 0;
	else if(y < CLOUD_F_Y) return 0;
	else if(y > CLOUD_F_Y+15) return 1;
	else
    {
        x -= CLOUD_F_X;
        y -= CLOUD_F_Y;

        return
        ARR16(y,
          ARR16(x,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0),
          ARR16(x,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0),
          ARR16(x,1,1,0,1,1,0,1,1,1,1,0,1,1,0,1,0),
          ARR16(x,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0),
          0,
          ARR16(x,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0),
          ARR16(x,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0),
          0,
          ARR16(x,0,0,0,0,1,1,1,1,0,0,0,0,0,0,1,1),
          ARR16(x,0,0,1,1,1,1,1,1,1,1,0,1,1,1,1,1),
          ARR16(x,0,0,1,1,1,1,1,1,0,1,1,1,1,1,1,1),
          ARR16(x,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1),
          ARR16(x,1,1,1,1,1,0,0,0,1,1,1,1,1,1,1,1),
          ARR16(x,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1),
          1,
          ARR16(x,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1)
        );
    }
}
// Cloud tile G.
int cloudG(in int x, in int y)
{
	if(x < CLOUD_G_X || x >= CLOUD_H_X) return 0;
	else if(y < CLOUD_G_Y) return 0;
	else if(y > CLOUD_G_Y+7) return 1;
	else
    {
        x -= CLOUD_G_X;
        y -= CLOUD_G_Y;

        return
        ARR8(y,
          ARR8(x,0,0,0,0,0,0,1,1),
          ARR8(x,0,0,0,0,1,1,1,1),
          ARR8(x,1,1,0,1,1,1,1,1),
          ARR8(x,1,0,1,1,1,1,1,1),
          ARR8(x,0,1,1,1,1,1,1,1),
          ARR8(x,0,1,1,1,1,1,1,1),
          1,
          1
        );
    }
}
// Cloud tile H.
int cloudH(in int x, in int y)
{
	if(x < CLOUD_H_X || x >= CLOUD_I_X) return 0;
	else if(y < CLOUD_H_Y) return 0;
	else if(y > CLOUD_H_Y+7) return 1;
	else
    {
        x -= CLOUD_H_X;
        y -= CLOUD_H_Y;

        return
        ARR8(y,
          ARR8(x,1,1,1,0,0,0,0,0),
          ARR8(x,1,1,1,1,1,0,0,0),
          ARR8(x,1,1,1,1,1,1,0,0),
          ARR8(x,1,1,1,1,1,1,1,0),
          ARR8(x,1,1,1,1,1,1,1,0),
          1,
          ARR8(x,1,1,1,1,1,0,1,0),
          1
        );
    }
}
// Cloud tile I.
int cloudI(in int x, in int y)
{
	if(x < CLOUD_I_X || x >= CLOUD_J_X) return 0;
	else if(y < CLOUD_I_Y) return 0;
	else if(y > CLOUD_I_Y+7) return 1;
	else
    {
        x -= CLOUD_I_X;
        y -= CLOUD_I_Y;

        return
        ARR8(y,
          ARR8(x,1,1,0,0,0,0,0,0),
          ARR8(x,1,1,1,0,0,1,1,0),
          ARR8(x,1,1,1,0,0,0,1,0),
          ARR8(x,1,1,1,1,0,0,0,0),
          ARR8(x,1,1,1,1,0,0,0,0),
          ARR8(x,1,1,1,1,0,0,0,0),
          ARR8(x,1,1,1,0,1,0,0,0),
          1
        );
    }
}
// Cloud tile J.
int cloudJ(in int x, in int y)
{
	if(x < CLOUD_J_X || x >= CLOUD_K_X) return 0;
	else if(y < CLOUD_J_Y) return 0;
	else if(y > CLOUD_J_Y+7) return 1;
	else
    {
        x -= CLOUD_J_X;
        y -= CLOUD_J_Y;

        return
        ARR8(y,
          ARR8(x,1,1,0,0,1,1,0,0),
          ARR8(x,1,1,1,0,1,1,0,0),
          ARR8(x,1,1,1,1,0,0,0,0),
          ARR8(x,1,1,1,1,0,1,0,0),
          ARR8(x,1,1,1,1,0,0,0,0),
          ARR8(x,1,1,1,1,1,0,0,0),
          ARR8(x,1,1,1,1,1,0,0,0),
          ARR8(x,1,1,1,1,1,0,0,0)
        );
    }
}
// Cloud tile K.
int cloudK(in int x, in int y)
{
	if(x < CLOUD_K_X || x >= CLOUD_L_X) return 0;
	else if(y < CLOUD_K_Y) return 0;
	else if(y > CLOUD_K_Y+1) return 1;
	else
    {
        x -= CLOUD_K_X;
        y -= CLOUD_K_Y;

        return
        ARR2(y,
          ARR8(x,0,0,1,0,0,0,1,1),
          ARR8(x,1,1,1,1,0,1,1,1)
        );
    }
}
// Cloud tile L. This one is repeated twice along X.
int cloudL(in int x, in int y)
{
	if(x < CLOUD_L_X || x >= CLOUD_M_X) return 0;
	else if(y < CLOUD_L_Y) return 0;
	else if(y > CLOUD_L_Y+1) return 1;
	else
    {
        x -= CLOUD_L_X;
        y -= CLOUD_L_Y;

        x = int(mod(float(x),8.0));

        return
        ARR2(y,
          ARR8(x,1,1,0,0,0,0,0,0),
          ARR8(x,1,1,1,0,0,1,1,0)
        );
    }
}
// CLoud tile M.
int cloudM(in int x, in int y)
{
	if(x < CLOUD_M_X || x >= CLOUD_N_X) return 0;
	else if(y < CLOUD_M_Y) return 0;
	else if(y > CLOUD_M_Y+3) return 1;
	else
    {
        x -= CLOUD_M_X;
        y -= CLOUD_M_Y;

        return
        ARR4(y,
          ARR8(x,0,0,0,1,1,0,0,0),
          ARR8(x,0,1,0,1,1,0,0,0),
          ARR8(x,0,0,0,0,0,0,0,1),
          ARR8(x,1,1,1,1,0,0,1,1)
        );
    }
}
// Cloud tile N. This is repeated to coda.
int cloudN(in int x, in int y)
{
	if(x < CLOUD_N_X) return 0;
	else if(y < CLOUD_N_Y) return 0;
	else if(y > CLOUD_N_Y+1) return 1;
	else
    {
        x -= CLOUD_N_X;
        y -= CLOUD_N_Y;

        x = int(mod(float(x),32.0));

        return
        ARR2(y,
          ARR32(x,0,0,1,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,1,0,0,0,1,1,1,1,0,0,0,0,0,0),
          ARR32(x,1,1,1,1,0,1,1,1,1,1,1,0,0,1,1,0,1,1,1,1,0,1,1,1,1,1,1,0,0,1,1,0)
        );
    }
}

/*
*	The large cloud draw function.
*
*	Composites all the large cloud tiles together and
*	draws them to the screen.
*
*	x: The x position of the current fragment.
*	y: The y position of the current fragment.
*
*	Returns: The color of the cloud from under the current fragment.
*/
vec4 drawNearClouds(in int x, in int y)
{
    if(y > WAVES_Y) return TRANS;
    else if(y < FAR_CLOUD_Y) return TRANS;
    else if(x < CLOUD_A_X) return TRANS;
    else
    {
        // The usual broken-apart additive blending.
        vec4 result = nearCloudsPalette(cloudA(x,y));
        result += nearCloudsPalette(cloudB(x,y));
        result += nearCloudsPalette(cloudC(x,y));
        result += nearCloudsPalette(cloudD(x,y));
        result += nearCloudsPalette(cloudE(x,y));
        result += nearCloudsPalette(cloudF(x,y));
        result += nearCloudsPalette(cloudG(x,y));
        result += nearCloudsPalette(cloudH(x,y));
        result += nearCloudsPalette(cloudI(x,y));
        result += nearCloudsPalette(cloudJ(x,y));
        result += nearCloudsPalette(cloudK(x,y));
        result += nearCloudsPalette(cloudL(x,y));
        result += nearCloudsPalette(cloudM(x,y));
        result += nearCloudsPalette(cloudN(x,y));
		return result;
    }
}

/*
*	The palette of the smaller clouds floating above.
*   
*	Returns a color given a palette index.
*
*	c: The color index to look up.
*
*	Returns: The corresponding color.
*/
vec4 smallCloudPalette(in int x)
{
	return ARR4(x, TRANS, WHITE, L_BLUE, TRANS);
}

/*
*	The tile function of the smaller part of the small clouds.
*   
*	Returns a palette index given a position.
*	Since this tile is repeated within the cloud, we have to
*	be able to specify where to draw it.
*
*	x: The x position of the current fragment.
*	y: The y position of the current fragment.
*	atx: The x position at which to draw the cloud.
*	aty: The y position at which to draw the cloud.
*
*	Returns: The corresponding palette index.
*/
int smallCloudA(in int x, in int y, in int atx, in int aty)
{
	if(x < atx || x > atx+7) return 0;
	else if(y < aty || y > aty+3) return 0;
	else
    {
        x -= atx;
        y -= aty;

        return
        ARR4(y,
          ARR8(x,0,0,0,0,2,2,0,0),
          ARR8(x,0,0,2,1,1,0,2,0),
          ARR8(x,1,0,0,1,2,0,0,0),
          ARR8(x,2,0,0,2,0,0,0,0)
        );
    }
}

/*
*	The tile representing the large part of the small cloud.
*   
*	Returns a palette index given a position.
*
*	x: The x position of the current fragment.
*	y: The y position of the current fragment.
*
*	Returns: The corresponding palette index.
*/
int smallCloudB(in int x, in int y)
{
	if(x < S_CLOUD_B_X || x > S_CLOUD_B_X+15) return 0;
	else if(y < S_CLOUD_B_Y || y > S_CLOUD_B_Y+7) return 0;
	else
    {
        x -= S_CLOUD_B_X;
        y -= S_CLOUD_B_Y;

        return
        ARR8(y,
          ARR16(x,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
          ARR16(x,1,2,0,0,1,1,1,0,0,0,0,0,0,0,0,0),
          ARR16(x,0,0,0,1,1,1,1,2,0,1,1,0,0,0,0,0),
          ARR16(x,0,0,2,1,1,1,1,2,1,1,1,2,0,1,0,0),
          ARR16(x,0,0,2,1,1,1,2,2,2,0,2,1,0,0,0,1),
          ARR16(x,2,1,0,2,2,2,0,2,0,0,0,0,0,0,0,0),
          0,
          0
        );
    }
}

/*
*	The small cloud's draw function.
*
*	Draws the smaller cloud to the screen.
*
*	x: The x position of the current fragment.
*	y: The y position of the current fragment.
*
*	Returns: The color of the cloud from under the current fragment.
*/
vec4 drawSmallCloud(in int x, in int y)
{
    // smallCloudA actually appears twice.
	vec4 result = smallCloudPalette(smallCloudA(x,y,S_CLOUD_A_X,S_CLOUD_A_Y));
	result += smallCloudPalette(smallCloudB(x,y));
	result += smallCloudPalette(smallCloudA(x,y,S_CLOUD_C_X,S_CLOUD_C_Y));
	return result;
}

/*
*	The boat draw function.
*
*	Draws the boat in the corner. Unlike all the other art,
*	this doesn't use the LUT approach.
*
*	x: The x position of the current fragment.
*	y: The y position of the current fragment.
*
*	Returns: The color of the cloud from under the current fragment.
*/
vec4 drawBoat(in int x, in int y)
{
	// Oh look the boat looks just like elementary inequality graphs...
    x = -x; // save time, negate x.
    // Most common case is the first checked.
	if (y > 2*x + 71 || y > x + 40) return TRANS;
	else if(y > 2*x + 24) return BLACK;
    else return BROWN;
	
}

/*
*	The global draw function.
*
*	Calculates the contribution of all scene elements to 
*	the current fragment.
*
*	x: The x position of the current fragment.
*	y: The y position of the current fragment.
*
*	Returns: The color of the entire scene under the current fragment.
*/
vec4 drawElements(in int x, in int y)
{
    // Reuse some variables.
    vec4 result = drawFarClouds(x,y);
    vec4 element = drawNearClouds(x,y);
    
    result = mix(result, element, element.a);
    element = drawSmallCloud(x,y);
    result = mix(result, element, element.a);
    element = drawBirds(x,y);
    result = mix(result, element, element.a);
    element = drawBoat(x,y);
    result = mix(result, element, element.a);
    element = drawShore(x,y);
    result = mix(result, element, element.a);
    element = drawYumetarou(x,y);
    result = mix(result, element, element.a);
    #ifdef DRAW_WAVES
    element = drawWaves(x,y);
    result = mix(result, element, element.a);
    #endif
    return result;
}

/*
*	The main draw function.
*
*	Computes the color of the current fragment.
*
*	fragColor: The computed fragment color.
*	fragCoord: The coordinate of this fragment.
*/
vec4 mainImage(in vec2 fragCoord)
{
    // Normalize coordinates.
    fragCoord = (fragCoord.xy / resolution.xy);
    
    // Invert the Y axis.
    fragCoord.y = 1.0-fragCoord.y;
    
    // Account for aspect ratio.
    fragCoord.x *= resolution.x/resolution.y;
    
    // Let's get NES sized pixels. This is the Y-resolution of Gimmick's screen sans-HUD.
    // We also have to account for the fact that the NES didn't have square pixels.
    fragCoord *= vec2(184.0*0.85736,184.0);
    
    // Determine and store the texel of the scene elements this pixel occupies.
    vec4 imageElements = drawElements(int(fragCoord.x), int(fragCoord.y));
    
	// Mix this texel with the background.
	return mix(D_BLUE, imageElements, imageElements.a);
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
	uv.y -= (2.0*horizon)+1.0;
	uv.y += sin(-.65*uv.x*(uv.x+cos(iTime*2.+uv.y*10.)) + iTime * 24.) * -.04;

	float col = 0.;

	_(L)_(I)_(G)_(H)_(T)_(A)_(M)_(P)_(dot)_(B)_(R)_(O)_(U)_(G)_(T)_(dot)_(T)_(O)_(dot)_(U)_(dot)_(B)_(Y)_(dot)_(T)_(H)_(E)_(dot)_(B)_(O)_(L)_(dot)_(dot)_(dot)_(H)_(A)_(P)_(P)_(Y)_(dot)
    _(B)_(A)_(R)_(D)_(I)_(N)_(G)_(dot)

	return pow(col,3.)*1.66*(sin(iTime)*0.5+0.7);
}

vec3 copper(in vec3 c) {
    vec3 rgb = clamp( abs(mod(c.x*6.0+vec3(0.0,4.0,2.0),6.0)-3.0)-1.0, 0.0, 1.0 );
	rgb = rgb*rgb*(3.0-2.0*rgb);
	return c.z * mix( vec3(1.0), rgb, c.y) * vec3(1.0,1.0,2.5);
}

void main(void)
{
	vec2 uv = ( gl_FragCoord.xy/resolution.xy ) -.5;
	uv.x*=resolution.x/resolution.y;
	vec3 col;
	if (uv.y > -0.4) 
	{
        vec3 copper = copper(vec3(uv.y*4.-time*.4,.3,.4));
		float text = scroll(uv)+scroll(uv*vec2(2.,-1.)-vec2(0.,.68))*.5;
		col = text > 0. ? vec3(.5)*text*copper : col;
	}
    gl_FragColor = mainImage(gl_FragCoord.xy);
    gl_FragColor += vec4(col,0.5);

}
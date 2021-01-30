#include "texture.h"
#include "color.h"

#include <assert.h>
#include <iostream>
#include <algorithm>

using namespace std;

namespace CS248 {

inline void uint8_to_float( float dst[4], unsigned char* src ) {
  uint8_t* src_uint8 = (uint8_t *)src;
  dst[0] = src_uint8[0] / 255.f;
  dst[1] = src_uint8[1] / 255.f;
  dst[2] = src_uint8[2] / 255.f;
  dst[3] = src_uint8[3] / 255.f;
}

inline void float_to_uint8( unsigned char* dst, float src[4] ) {
  uint8_t* dst_uint8 = (uint8_t *)dst;
  dst_uint8[0] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[0])));
  dst_uint8[1] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[1])));
  dst_uint8[2] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[2])));
  dst_uint8[3] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[3])));
}

void Sampler2DImp::generate_mips(Texture& tex, int startLevel) {

  // NOTE: 
  // This starter code allocates the mip levels and generates a level 
  // map by filling each level with placeholder data in the form of a 
  // color that differs from its neighbours'. You should instead fill
  // with the correct data!

  // Extra credit: Implement this

  // check start level
  if ( startLevel >= tex.mipmap.size() ) {
    std::cerr << "Invalid start level"; 
  }

  // allocate sublevels
  int baseWidth  = tex.mipmap[startLevel].width;
  int baseHeight = tex.mipmap[startLevel].height;
  int numSubLevels = (int)(log2f( (float)max(baseWidth, baseHeight)));

  numSubLevels = min(numSubLevels, kMaxMipLevels - startLevel - 1);
  tex.mipmap.resize(startLevel + numSubLevels + 1);

  int width  = baseWidth;
  int height = baseHeight;
  for (int i = 1; i <= numSubLevels; i++) {

    MipLevel& level = tex.mipmap[startLevel + i];

    // handle odd size texture by rounding down
    width  = max( 1, width  / 2); assert(width  > 0);
    height = max( 1, height / 2); assert(height > 0);

    level.width = width;
    level.height = height;
    level.texels = vector<unsigned char>(4 * width * height);

  }

  // fill all 0 sub levels with interchanging colors
  Color colors[3] = { Color(1,0,0,1), Color(0,1,0,1), Color(0,0,1,1) };
  for(size_t i = 1; i < tex.mipmap.size(); ++i) {

    Color c = colors[i % 3];
    MipLevel& mip = tex.mipmap[i];

    for(size_t i = 0; i < 4 * mip.width * mip.height; i += 4) {
      float_to_uint8( &mip.texels[i], &c.r );
    }
  }

}

Color Sampler2DImp::sample_nearest(Texture& tex, 
                                   float u, float v, 
                                   int level) {

  // Task 4: Implement nearest neighbour interpolation

  // to return magenta in an invalid level:
  if (level > tex.mipmap.size() || level < 0) {
    // return magenta for invalid level
    return Color(1,0,1,1);
  }
  

  // CHECK! getting quite a bit of pixel offset using this

  MipLevel& thismip = tex.mipmap[level];
  int mip_width = thismip.width;
  int mip_height = thismip.height;
  // std::vector<unsigned char> mip_texels = thismip.texels;



  // Note: the texels here are stored in the 4 pixels frame buffer structure
  int nearest_x = round(mip_width*u - 0.5f); // offset the added 0.5 that was added in software_renderer to get correct coordinate
  int nearest_y = round(mip_height*v - 0.5f);

  Color color_to_return;

  color_to_return.r = thismip.texels[4*(nearest_x + nearest_y*mip_width)] / 255.f; // cast to float
  color_to_return.g = thismip.texels[4*(nearest_x + nearest_y*mip_width) + 1] / 255.f; // cast to float
  color_to_return.b = thismip.texels[4*(nearest_x + nearest_y*mip_width) + 2] / 255.f; // cast to float
  color_to_return.a = thismip.texels[4*(nearest_x + nearest_y*mip_width) + 3] / 255.f; // cast to float

  return color_to_return;

}

inline Color lerp(float x, Color c0, Color c1) {
        return (1 - x) * c0 + x * c1;
    }

Color Sampler2DImp::sample_bilinear(Texture& tex, 
                                    float u, float v, 
                                    int level) {


  
  // Task 4: Implement bilinear filtering

  // to return magenta in an invalid level:
  if (level > tex.mipmap.size() || level < 0) {
    // return magenta for invalid level
    return Color(1,0,1,1);
  }

  MipLevel& thismip = tex.mipmap[level];
  int mip_width = thismip.width;
  int mip_height = thismip.height;
  // std::vector<unsigned char> mip_texels = thismip.texels;
  float pixel_u = mip_width*u - 0.5f;
  float pixel_v = mip_height*v - 0.5f;

  // CLAMPING IMPLEMENTATION:
  // if any of the bounds are outside of the render square, set them to the bound value\
  // Technically, these statements are contradictory. Not an issue for textures larger than 1x1...
  if (pixel_u < 0) {
    pixel_u = 0;
  }
  if (pixel_v < 0) {
    pixel_v = 0;
  }
  if (pixel_u > mip_width - 1) {
    pixel_u = mip_width - 1;
  }
  if (pixel_v > mip_height - 1) {
    pixel_v = mip_height - 1;
  }



  // Note: the texels here are stored in the 4 pixels frame buffer structure

  int bounds_u0 = floor(pixel_u);
  int bounds_u1 = bounds_u0 + 1;
  int bounds_v0 = floor(pixel_v);
  int bounds_v1 = bounds_v0 + 1;

 
  // Bilinear takes 3 interpolations:
  // Get colors for first 2 based on outter bounds


  // Lower horizontal line
  Color lh_0, lh_1, lh_comb;

  int base_indx = 4*(bounds_u0 + bounds_v0*mip_width);
  lh_0.r = thismip.texels[base_indx] / 255.f; // cast to float
  lh_0.g = thismip.texels[base_indx + 1] / 255.f; // cast to float
  lh_0.b = thismip.texels[base_indx + 2] / 255.f; // cast to float
  lh_0.a = thismip.texels[base_indx + 3] / 255.f; // cast to float

  base_indx = 4*(bounds_u1 + bounds_v0*mip_width);
  lh_1.r = thismip.texels[base_indx] / 255.f; // cast to float
  lh_1.g = thismip.texels[base_indx + 1] / 255.f; // cast to float
  lh_1.b = thismip.texels[base_indx + 2] / 255.f; // cast to float
  lh_1.a = thismip.texels[base_indx + 3] / 255.f; // cast to float

  lh_comb.r = lh_0.r + ( pixel_u-bounds_u0 ) * (lh_1.r - lh_0.r);
  lh_comb.g = lh_0.g + ( pixel_u-bounds_u0 ) * (lh_1.g - lh_0.g);
  lh_comb.b = lh_0.b + ( pixel_u-bounds_u0 ) * (lh_1.b - lh_0.b);
  lh_comb.a = lh_0.a + ( pixel_u-bounds_u0 ) * (lh_1.a - lh_0.a);


  // upper horizontal line
  Color uh_0, uh_1, uh_comb;

  base_indx = 4*(bounds_u0 + bounds_v1*mip_width);
  uh_0.r = thismip.texels[base_indx] / 255.f; // cast to float
  uh_0.g = thismip.texels[base_indx + 1] / 255.f; // cast to float
  uh_0.b = thismip.texels[base_indx + 2] / 255.f; // cast to float
  uh_0.a = thismip.texels[base_indx + 3] / 255.f; // cast to float
  
  base_indx = 4*(bounds_u1 + bounds_v1*mip_width);
  uh_1.r = thismip.texels[base_indx] / 255.f; // cast to float
  uh_1.g = thismip.texels[base_indx + 1] / 255.f; // cast to float
  uh_1.b = thismip.texels[base_indx + 2] / 255.f; // cast to float
  uh_1.a = thismip.texels[base_indx + 3] / 255.f; // cast to float

  uh_comb.r = uh_0.r + ( pixel_u-bounds_u0 ) * (uh_1.r - uh_0.r);
  uh_comb.g = uh_0.g + ( pixel_u-bounds_u0 ) * (uh_1.g - uh_0.g);
  uh_comb.b = uh_0.b + ( pixel_u-bounds_u0 ) * (uh_1.b - uh_0.b);
  uh_comb.a = uh_0.a + ( pixel_u-bounds_u0 ) * (uh_1.a - uh_0.a);

  // Combination interpolation (along y axis between 2 horizontal interps)
  Color color_to_return;

  color_to_return.r = lh_comb.r + ( pixel_v-bounds_v0 ) * (uh_comb.r - lh_comb.r);
  color_to_return.g = lh_comb.g + ( pixel_v-bounds_v0 ) * (uh_comb.g - lh_comb.g);
  color_to_return.b = lh_comb.b + ( pixel_v-bounds_v0 ) * (uh_comb.b - lh_comb.b);
  color_to_return.a = lh_comb.a + ( pixel_v-bounds_v0 ) * (uh_comb.a - lh_comb.a);


  // printf("color @ u0,v0: red- %f, green- %f, blue- %f, alpha- %f\n", lh_0.r, lh_0.g, lh_0.b, lh_0.a);
  // printf("color @ u1,v0: red- %f, green- %f, blue- %f, alpha- %f\n", lh_1.r, lh_1.g, lh_1.b, lh_1.a);
  // printf("color @ u0,v1: red- %f, green- %f, blue- %f, alpha- %f\n", uh_0.r, uh_0.g, uh_0.b, uh_0.a);
  // printf("color @ u1,v1: red- %f, green- %f, blue- %f, alpha- %f\n", uh_1.r, uh_1.g, uh_1.b, uh_1.a);
  // printf("combo (%f, %f): red- %f, green- %f, blue- %f, alpha- %f\n", pixel_u, pixel_v, color_to_return.r, color_to_return.g, color_to_return.b, color_to_return.a);

  // printf("Width, height: %d, %d", mip_width, mip_width);
  // printf("Input  (u,v) = (%f, %f)", u,v);
  // printf("Output pixelation (u,v) = (%f, %f)", pixel_u,pixel_v);

  // printf("\n\n");


  return color_to_return;

}

Color Sampler2DImp::sample_trilinear(Texture& tex, 
                                     float u, float v, 
                                     float u_scale, float v_scale) {

  // Extra credit: Implement trilinear filtering

  // return magenta for invalid level
  return Color(1,0,1,1);

}

} // namespace CS248

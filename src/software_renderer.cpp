#include "software_renderer.h"

#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>

#include "triangulation.h"

#include <iostream>
#include <numeric>

using namespace std;

namespace CS248 {

// Global variable definitions!
std::vector<unsigned char> supersample_render_target;
size_t supersample_width = 0;
size_t supersample_height = 0;


// Implements SoftwareRenderer //

// fill a sample location with color
void SoftwareRendererImp::fill_sample(int sx, int sy, const Color &color) {
  // Here, we assume sx and sy are in the rendering boundaries

  // First, fill frame buffer at given sx,sy points with appropriate color

  supersample_render_target[4 * (sx + sy * supersample_width)] = (uint8_t)(color.r * 255);
  supersample_render_target[4 * (sx + sy * supersample_width) + 1] = (uint8_t)(color.g * 255);
  supersample_render_target[4 * (sx + sy * supersample_width) + 2] = (uint8_t)(color.b * 255);
  supersample_render_target[4 * (sx + sy * supersample_width) + 3] = (uint8_t)(color.a * 255);



  // Next, apply appropriate transform to give the sample alpha
  // CHECK! Is this correct? is this alpha blending needed? alpha_blending_helper doesnt seem to be implemented...
  Color pixel_color;
	float inv255 = 1.0 / 255.0;

  pixel_color.r = supersample_render_target[4 * (sx + sy * supersample_width)] * inv255;
	pixel_color.g = supersample_render_target[4 * (sx + sy * supersample_width) + 1] * inv255;
	pixel_color.b = supersample_render_target[4 * (sx + sy * supersample_width) + 2] * inv255;
	pixel_color.a = supersample_render_target[4 * (sx + sy * supersample_width) + 3] * inv255;

	pixel_color = ref->alpha_blending_helper(pixel_color, color);

  supersample_render_target[4 * (sx + sy * supersample_width)] = (uint8_t)(pixel_color.r * 255);
	supersample_render_target[4 * (sx + sy * supersample_width) + 1] = (uint8_t)(pixel_color.g * 255);
	supersample_render_target[4 * (sx + sy * supersample_width) + 2] = (uint8_t)(pixel_color.b * 255);
	supersample_render_target[4 * (sx + sy * supersample_width) + 3] = (uint8_t)(pixel_color.a * 255);


}

// fill samples in the entire pixel specified by pixel coordinates
void SoftwareRendererImp::fill_pixel(int x, int y, const Color &color) {

	// Task 2: Re-implement this function

	// check bounds
	if (x < 0 || x >= target_w) return;
	if (y < 0 || y >= target_h) return;

	Color pixel_color;
	float inv255 = 1.0 / 255.0;
	pixel_color.r = render_target[4 * (x + y * target_w)] * inv255;
	pixel_color.g = render_target[4 * (x + y * target_w) + 1] * inv255;
	pixel_color.b = render_target[4 * (x + y * target_w) + 2] * inv255;
	pixel_color.a = render_target[4 * (x + y * target_w) + 3] * inv255;

	pixel_color = ref->alpha_blending_helper(pixel_color, color);

	render_target[4 * (x + y * target_w)] = (uint8_t)(pixel_color.r * 255);
	render_target[4 * (x + y * target_w) + 1] = (uint8_t)(pixel_color.g * 255);
	render_target[4 * (x + y * target_w) + 2] = (uint8_t)(pixel_color.b * 255);
	render_target[4 * (x + y * target_w) + 3] = (uint8_t)(pixel_color.a * 255);


  // supersample_render_target[4 * (x + y * target_w)] = (uint8_t)(pixel_color.r * 255);
	// supersample_render_target[4 * (x + y * target_w) + 1] = (uint8_t)(pixel_color.g * 255);
	// supersample_render_target[4 * (x + y * target_w) + 2] = (uint8_t)(pixel_color.b * 255);
	// supersample_render_target[4 * (x + y * target_w) + 3] = (uint8_t)(pixel_color.a * 255);

}

void SoftwareRendererImp::draw_svg( SVG& svg ) {

  // set top level transformation
  transformation = canvas_to_screen;

  // draw all elements
  for ( size_t i = 0; i < svg.elements.size(); ++i ) {
    draw_element(svg.elements[i]);
  }

  // draw canvas outline
  Vector2D a = transform(Vector2D(    0    ,     0    )); a.x--; a.y--;
  Vector2D b = transform(Vector2D(svg.width,     0    )); b.x++; b.y--;
  Vector2D c = transform(Vector2D(    0    ,svg.height)); c.x--; c.y++;
  Vector2D d = transform(Vector2D(svg.width,svg.height)); d.x++; d.y++;

  rasterize_line(a.x, a.y, b.x, b.y, Color::Black);
  rasterize_line(a.x, a.y, c.x, c.y, Color::Black);
  rasterize_line(d.x, d.y, b.x, b.y, Color::Black);
  rasterize_line(d.x, d.y, c.x, c.y, Color::Black);

  // resolve and send to render target
  resolve();

}

void SoftwareRendererImp::set_sample_rate( size_t sample_rate ) {

  // Task 2: 
  // You may want to modify this for supersampling support
  this->sample_rate = sample_rate;

  // CHECK! what are we supposed to modify here? handled in drawsvg
  supersample_render_target.resize( target_h*target_w*4*sample_rate*sample_rate );
  supersample_width = this->target_w*sample_rate;
  supersample_height = this->target_h*sample_rate;

  printf("Sample rate ois %d\n", (int)sample_rate);


}

void SoftwareRendererImp::set_render_target( unsigned char* render_target,
                                             size_t width, size_t height ) {

  // Task 2: 
  // You may want to modify this for supersampling support
  this->render_target = render_target;
  this->target_w = width;
  this->target_h = height;

  // CHECK! what are we supposed to modify here? made the supersample buffer global
  cout << "width is " << width;
  cout << "height is " << height;
  supersample_render_target.resize(width*height*4);
  cout << "total array length is " << supersample_render_target.size();

  supersample_width = width*this->sample_rate;
  supersample_height = height*this->sample_rate;


  // // Supersample support...create the supersample_render_target
  // int arrSize = render_target.size();
  // cout << "Full array size is: " << arrSize;


}

void SoftwareRendererImp::draw_element( SVGElement* element ) {

	// Task 3 (part 1):
	// Modify this to implement the transformation stack

	switch (element->type) {
	case POINT:
		draw_point(static_cast<Point&>(*element));
		break;
	case LINE:
		draw_line(static_cast<Line&>(*element));
		break;
	case POLYLINE:
		draw_polyline(static_cast<Polyline&>(*element));
		break;
	case RECT:
		draw_rect(static_cast<Rect&>(*element));
		break;
	case POLYGON:
		draw_polygon(static_cast<Polygon&>(*element));
		break;
	case ELLIPSE:
		draw_ellipse(static_cast<Ellipse&>(*element));
		break;
	case IMAGE:
		draw_image(static_cast<Image&>(*element));
		break;
	case GROUP:
		draw_group(static_cast<Group&>(*element));
		break;
	default:
		break;
	}

}


// Primitive Drawing //

void SoftwareRendererImp::draw_point( Point& point ) {

  Vector2D p = transform(point.position);
  rasterize_point( p.x, p.y, point.style.fillColor );

}

void SoftwareRendererImp::draw_line( Line& line ) { 

  Vector2D p0 = transform(line.from);
  Vector2D p1 = transform(line.to);
  rasterize_line( p0.x, p0.y, p1.x, p1.y, line.style.strokeColor );

}

void SoftwareRendererImp::draw_polyline( Polyline& polyline ) {

  Color c = polyline.style.strokeColor;

  if( c.a != 0 ) {
    int nPoints = polyline.points.size();
    for( int i = 0; i < nPoints - 1; i++ ) {
      Vector2D p0 = transform(polyline.points[(i+0) % nPoints]);
      Vector2D p1 = transform(polyline.points[(i+1) % nPoints]);
      rasterize_line( p0.x, p0.y, p1.x, p1.y, c );
    }
  }
}

void SoftwareRendererImp::draw_rect( Rect& rect ) {

  Color c;
  
  // draw as two triangles
  float x = rect.position.x;
  float y = rect.position.y;
  float w = rect.dimension.x;
  float h = rect.dimension.y;

  Vector2D p0 = transform(Vector2D(   x   ,   y   ));
  Vector2D p1 = transform(Vector2D( x + w ,   y   ));
  Vector2D p2 = transform(Vector2D(   x   , y + h ));
  Vector2D p3 = transform(Vector2D( x + w , y + h ));
  
  // draw fill
  c = rect.style.fillColor;
  if (c.a != 0 ) {
    rasterize_triangle( p0.x, p0.y, p1.x, p1.y, p2.x, p2.y, c );
    rasterize_triangle( p2.x, p2.y, p1.x, p1.y, p3.x, p3.y, c );
  }

  // draw outline
  c = rect.style.strokeColor;
  if( c.a != 0 ) {
    rasterize_line( p0.x, p0.y, p1.x, p1.y, c );
    rasterize_line( p1.x, p1.y, p3.x, p3.y, c );
    rasterize_line( p3.x, p3.y, p2.x, p2.y, c );
    rasterize_line( p2.x, p2.y, p0.x, p0.y, c );
  }

}

void SoftwareRendererImp::draw_polygon( Polygon& polygon ) {

  Color c;

  // draw fill
  c = polygon.style.fillColor;
  if( c.a != 0 ) {

    // triangulate
    vector<Vector2D> triangles;
    triangulate( polygon, triangles );

    // draw as triangles
    for (size_t i = 0; i < triangles.size(); i += 3) {
      Vector2D p0 = transform(triangles[i + 0]);
      Vector2D p1 = transform(triangles[i + 1]);
      Vector2D p2 = transform(triangles[i + 2]);
      rasterize_triangle( p0.x, p0.y, p1.x, p1.y, p2.x, p2.y, c );
    }
  }

  // draw outline
  c = polygon.style.strokeColor;
  if( c.a != 0 ) {
    int nPoints = polygon.points.size();
    for( int i = 0; i < nPoints; i++ ) {
      Vector2D p0 = transform(polygon.points[(i+0) % nPoints]);
      Vector2D p1 = transform(polygon.points[(i+1) % nPoints]);
      rasterize_line( p0.x, p0.y, p1.x, p1.y, c );
    }
  }
}

void SoftwareRendererImp::draw_ellipse( Ellipse& ellipse ) {

  // Extra credit 

}

void SoftwareRendererImp::draw_image( Image& image ) {

  Vector2D p0 = transform(image.position);
  Vector2D p1 = transform(image.position + image.dimension);

  rasterize_image( p0.x, p0.y, p1.x, p1.y, image.tex );
}

void SoftwareRendererImp::draw_group( Group& group ) {

  for ( size_t i = 0; i < group.elements.size(); ++i ) {
    draw_element(group.elements[i]);
  }

}

// Rasterization //

// The input arguments in the rasterization functions 
// below are all defined in screen space coordinates

void SoftwareRendererImp::rasterize_point( float x, float y, Color color ) {

  // fill in the nearest pixel
  int sx = (int)floor(x);
  int sy = (int)floor(y);

  // check bounds
  if (sx < 0 || sx >= target_w) return;
  if (sy < 0 || sy >= target_h) return;

  // fill sample - NOT doing alpha blending!
  // TODO: Call fill_pixel here to run alpha blending
  render_target[4 * (sx + sy * target_w)] = (uint8_t)(color.r * 255);
  render_target[4 * (sx + sy * target_w) + 1] = (uint8_t)(color.g * 255);
  render_target[4 * (sx + sy * target_w) + 2] = (uint8_t)(color.b * 255);
  render_target[4 * (sx + sy * target_w) + 3] = (uint8_t)(color.a * 255);

  supersample_render_target[4 * (sx + sy * target_w)] = (uint8_t)(color.r * 255);
  supersample_render_target[4 * (sx + sy * target_w) + 1] = (uint8_t)(color.g * 255);
  supersample_render_target[4 * (sx + sy * target_w) + 2] = (uint8_t)(color.b * 255);
  supersample_render_target[4 * (sx + sy * target_w) + 3] = (uint8_t)(color.a * 255);

  // Alpha blending function
  // Want this to be after we fill the frame buffer because fill_pixel refers to what's at the frame buffer already and applies alpha to that
  fill_pixel(sx, sy, color); 


}

void SoftwareRendererImp::rasterize_line( float x0, float y0,
                                          float x1, float y1,
                                          Color color) {

  // Extra credit (delete the line below and implement your own)
  ref->rasterize_line_helper(x0, y0, x1, y1, target_w, target_h, color, this);

}

void SoftwareRendererImp::rasterize_triangle( float x0, float y0,
                                              float x1, float y1,
                                              float x2, float y2,
                                              Color color ) {
  // Task 1: 
  // Implement triangle rasterization (you may want to call fill_sample here)

  // How to rasterize triangle:
  //    1) Check that point is in the canvas boundary
  //    2) Implement any speed improvements for reducing the search space >not implemented<
  //    3) Iterate through points and check that the triple bounding condition is met
  //      --> If met, call fill_sample

  printf("sample rate = %d", sample_rate);

  float ss_x0 = x0*sample_rate;
  float ss_y0 = y0*sample_rate;
  float ss_x1 = x1*sample_rate;
  float ss_y1 = y1*sample_rate;
  float ss_x2 = x2*sample_rate;
  float ss_y2 = y2*sample_rate;

  // 1) check bounds
  // Point 0
  if (ss_x0 < 0 || ss_x0 >= supersample_width) return;
  if (ss_y0 < 0 || ss_y0 >= supersample_height) return;

  // Point 1
  if (ss_x1 < 0 || ss_x1 >= supersample_width) return;
  if (ss_y1 < 0 || ss_y1 >= supersample_height) return;

  // Point 2
  if (ss_x2 < 0 || ss_x2 >= supersample_width) return;
  if (ss_y2 < 0 || ss_y2 >= supersample_height) return;



  // 2)
  // Get minimum X and Y coords from P0,P1,P2
  int row_iter_Lo_bound = std::min({ ss_y0, ss_y1, ss_y2 });
  int col_iter_Lo_bound = std::min({ ss_x0, ss_x1, ss_x2 });

  // Get maximum X and Y coords from P0,P1,P2
  int row_iter_Hi_bound = std::max({ ss_y0, ss_y1, ss_y2 });
  int col_iter_Hi_bound = std::max({ ss_x0, ss_x1, ss_x2 });



  // 3) Go through all sample points
  // Iterate in row-dominant form
  // First, set up all of the 'inside line' based dot product equations
  float sample_x, sample_y;
  float L_i[3];
  bool in_triangle;

  // CHECK! If there's an error with task 1 output, check the triangle windings thing... might be the reason
  float A[3] = {ss_y1 - ss_y0, 
              ss_y2 - ss_y1,
              ss_y0 - ss_y2}; // windings line
  float B[3] = {ss_x1 - ss_x0,
              ss_x2 - ss_x1,
              ss_x0 - ss_x2}; // windings line
  float C[3] = {ss_y0*B[0] - ss_x0*A[0],
              ss_y1*B[1] - ss_x1*A[1],
              ss_y2*B[2] - ss_x2*A[2]};



  // // Attempting fast incremental update
  // // Instead of calculating all L_i every loop, we just calculate the first one and add whatever we need 
  // // CHECK! Attempted to make it faster...not sure if it is with all the extra assigns and stuffsample_rate
  // float prev_L_i[3];
  // float og_sample_x = col_iter_Lo_bound + 0.5;
  // float og_sample_y = row_iter_Lo_bound + 0.5;
  // int delta_iters = 0;
  // prev_L_i[0] = A[0]*og_sample_x - B[0]*og_sample_y + C[0];
  // prev_L_i[1] = A[1]*og_sample_x - B[1]*og_sample_y + C[1];
  // prev_L_i[2] = A[2]*og_sample_x - B[2]*og_sample_y + C[2];

  // float bound_to_check = 0.01;


  // for (int count_row = row_iter_Lo_bound; count_row < row_iter_Hi_bound; count_row++) {
  //   for (int count_col = col_iter_Lo_bound; count_col < col_iter_Hi_bound; count_col++) {

  //     sample_x = count_col + 0.5;
  //     sample_y = count_row + 0.5;

  //     if ((count_col == col_iter_Lo_bound+1) && (count_row == row_iter_Lo_bound+1)) {
  //       printf("sample_x = %f", sample_x);
  //       printf("count_col = %d", count_col);
  //     }

  //     in_triangle = (prev_L_i[0] <= bound_to_check) && (prev_L_i[1] <= bound_to_check) && (prev_L_i[2] <= bound_to_check);

  //     if ( in_triangle ) {
  //       fill_sample(count_col, count_row, color); // If we made it here, we're in the triangle, so Color it!
  //     }

  //     prev_L_i[0] += A[0];
  //     prev_L_i[1] += A[1];
  //     prev_L_i[2] += A[2];

  //   }

  //   prev_L_i[0] -= B[0];
  //   prev_L_i[1] -= B[1];
  //   prev_L_i[2] -= B[2];

  //   delta_iters = (col_iter_Hi_bound - col_iter_Lo_bound);

  //   prev_L_i[0] -= delta_iters*A[0];
  //   prev_L_i[1] -= delta_iters*A[1];
  //   prev_L_i[2] -= delta_iters*A[2];

  // }



  // Original method!
  for (int count_row = row_iter_Lo_bound; count_row < row_iter_Hi_bound; count_row++) {
    for (int count_col = col_iter_Lo_bound; count_col < col_iter_Hi_bound; count_col++) {

      sample_x = count_col + 0.5;
      sample_y = count_row + 0.5;


      for(int ixx = 0; ixx < 3; ixx++) {
        L_i[ixx] =  A[ixx]*sample_x - B[ixx]*sample_y + C[ixx];
        // if (L_i[ixx] == 0) {
        //   // We have encountered a literal edge case!
        //   // Let's consistently count it as 'inside'
        //   L_i[ixx] = -1;
        //   cout << "This happened";
        //   // CHECK! Is this how we're supposed to handle the literal 'edge cases'??
        // }
      }

      // in_triangle = (L_i[0] <= 0.001) && (L_i[1] <= 0.001) && (L_i[2] <= 0.001);
      in_triangle = (L_i[0] <= 0) && (L_i[1] <= 0) && (L_i[2] <= 0);

      if ( in_triangle ) {
        fill_sample(count_col, count_row, color); // If we made it here, we're in the triangle, so Color it!
      }

    }

  }




}

void SoftwareRendererImp::rasterize_image( float x0, float y0,
                                           float x1, float y1,
                                           Texture& tex ) {
  // Task 4: 
  // Implement image rasterization (you may want to call fill_sample here)

}

// resolve samples to render target
void SoftwareRendererImp::resolve( void ) {

  // Task 2: 
  // Implement supersampling
  // You may also need to modify other functions marked with "Task 2".

  // To do a unit area box filter we just average over the appropriate squares 

  // Let's use a mask vector and do a dot product, will reduce code clutter (hopefully)


  


  // CHECK! Come back and change filter and supersample types to be unsigned char to optimize
  std::vector<unsigned char> filter_mask;
  printf("width is %d\n", (int)supersample_width);
  filter_mask.resize( 4*sample_rate*(supersample_width + 1) ); // simplified from sample_rate*supersample_width + sample_rate*4

  // std::fill(filter_mask.begin(), filter_mask.end(), 0);
  // filter_mask = {0}; // Initialize all elemenets to 0

  // Now create the filter mask
  // If sample rate == 2, then we have an adjacent set of pixels, plus 2 at 4*supersample_width away from each one
  // If sample rate == 2, have 2 adjacent pixels, plus 3 at 4*supersample_width, plus another 3 at 4*supersample_width*2
  // etc

  // This will loop over sample_rate^2 points, since thats how many extra samples per pixels we have 
  for (int i = 0; i < sample_rate; i++) {
    for (int j = 0; j < sample_rate; j++) {
      filter_mask[4*(j + supersample_width*i)] = 1;
      printf("curr i is %d\n", i);
      printf("curr j is %d\n", j);
      printf("curr pos is %d\n", 4*j + supersample_width*i);
    }
  }

  // Now calculate inner product over the supersample_render_target by shifting the start point of the dot product
  int render_target_length = 4 * target_w * target_h;
  int filter_len = filter_mask.size();

  printf("filter length: %d\n", filter_len);


  printf("\n");

  // for (int letsloop = 0; letsloop < 12; letsloop ++) {
  //     supersample_render_target[letsloop] = 200;
  // }

  // int tempSum = 0;
  int indx_target = 0;
  int indx_supsamp_target = 0;
  int sample_rate_sq = sample_rate*sample_rate;
  for(int i = 0; i < target_h; i++) {
    for(int j = 0; j < target_w; j++) {





      indx_target = 4*(i*target_w + j);
      indx_supsamp_target = 4*sample_rate*(i*supersample_width + j);

    
      render_target[indx_target] = std::inner_product(supersample_render_target.begin()+indx_supsamp_target, supersample_render_target.begin()+indx_supsamp_target+filter_len, filter_mask.begin(), 0) / sample_rate_sq;
      render_target[indx_target+1] = std::inner_product(supersample_render_target.begin()+indx_supsamp_target+1, supersample_render_target.begin()+indx_supsamp_target+1+filter_len, filter_mask.begin(), 0) / sample_rate_sq;
      render_target[indx_target+2] = std::inner_product(supersample_render_target.begin()+indx_supsamp_target+2, supersample_render_target.begin()+indx_supsamp_target+2+filter_len, filter_mask.begin(), 0) / sample_rate_sq;
      render_target[indx_target+3] = std::inner_product(supersample_render_target.begin()+indx_supsamp_target+3, supersample_render_target.begin()+indx_supsamp_target+3+filter_len, filter_mask.begin(), 0) / sample_rate_sq;


      // render_target[i] = supersample_render_target[4*sample_rate_sq*i];
      // render_target[i+1] = supersample_render_target[4*sample_rate_sq*i + 1];
      // render_target[i] = supersample_render_target[4*sample_rate_sq*i + 2];
      // render_target[i] = supersample_render_target[4*sample_rate_sq*i + 3];

    }
  }

  std::fill(supersample_render_target.begin(), supersample_render_target.end(), 0); // Fill super_render_target with 0s so we can reset it!

    
  return;

}

Color SoftwareRendererImp::alpha_blending(Color pixel_color, Color color)
{
  // Task 5
  // Implement alpha compositing
  return pixel_color;
}


} // namespace CS248

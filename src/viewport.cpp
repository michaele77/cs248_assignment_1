#include "viewport.h"

#include "CS248.h"

namespace CS248 {

void ViewportImp::set_viewbox( float x, float y, float span ) {

  // Task 3 (part 2): 
  // Set svg to normalized device coordinate transformation. Your input
  // arguments are defined as SVG canvans coordinates.


  // Translation matrix: identity matrix + bx @ [2][0] + by @ [2][1]
  Matrix3x3 transmat_translation = Matrix3x3::identity(); // Get identity matrix
  transmat_translation[2][0] = -1*(x - span)/(2*span); // right shift by x - span
  transmat_translation[2][1] = -1*(y - span)/(2*span); // downshift by y - span

  // Scale matrix: identity matrix + sx @ [0][0] + sy @ [1][1]
  Matrix3x3 transmat_scale = Matrix3x3::identity(); // Get identity matrix
  transmat_scale[0][0] = 1/(2*span);
  transmat_scale[1][1] = 1/(2*span);


  // Matrix3x3 transmat = transmat_reflection * transmat_translation * transmat_scale;
  Matrix3x3 transmat = transmat_translation*transmat_scale;

  this->svg_2_norm = transmat;


  this->x = x;
  this->y = y;
  this->span = span; 

}

void ViewportImp::update_viewbox( float dx, float dy, float scale ) { 
  
  this->x -= dx;
  this->y -= dy;
  this->span *= scale;
  set_viewbox( x, y, span );
}

} // namespace CS248

#include "viewport.h"

#include "CS248.h"

namespace CS248 {

void ViewportImp::set_viewbox( float x, float y, float span ) {

  // Task 3 (part 2): 
  // Set svg to normalized device coordinate transformation. Your input
  // arguments are defined as SVG canvans coordinates.

  // CHECK! How did this work before we even did any coordinate transformations??

  // Reflect about x matrix: identity matrix, but y is negative
  Matrix3x3 transmat_reflection = Matrix3x3::identity(); // Get identity matrix
  transmat_reflection[1][1] = -1; // x axis reflection

  // Translation matrix: identity matrix + bx @ [0][2] + by @ [1][2]
  Matrix3x3 transmat_translation = Matrix3x3::identity(); // Get identity matrix
  transmat_translation[0][2] = x - span; // right shift by x - span
  transmat_translation[1][2] = -y + span; // downshift by y - span

  // Scale matrix: identity matrix + sx @ [0][0] + sy @ [1][1]
  Matrix3x3 transmat_scale = Matrix3x3::identity(); // Get identity matrix
  transmat_scale[0][0] = 2*span; // scale by 2*span in x
  transmat_scale[1][1] = 2*span; // scale by 2*span in y


  Matrix3x3 transmat = transmat_reflection * transmat_translation * transmat_scale;

  Vector3D u( x, y, 1.0 );

  // apply projective space transformation
  u = transmat * u; 

  // project back to 2D Euclidean plane
  Vector2D u_out(u.x / u.z, u.y / u.z);



  this->x = u_out.x;
  this->y = u_out.y;
  this->span = span; 

}

void ViewportImp::update_viewbox( float dx, float dy, float scale ) { 
  
  this->x -= dx;
  this->y -= dy;
  this->span *= scale;
  set_viewbox( x, y, span );
}

} // namespace CS248

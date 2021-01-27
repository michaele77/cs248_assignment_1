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
  transmat_translation[0][2] = -1*(x - span); // right shift by x - span
  transmat_translation[1][2] = -1*(y - span); // downshift by y - span

  // Scale matrix: identity matrix + sx @ [0][0] + sy @ [1][1]
  Matrix3x3 transmat_scale = Matrix3x3::identity(); // Get identity matrix
  // transmat_scale[0][0] = 2*span; // scale by 2*span in x
  // transmat_scale[1][1] = 2*span; // scale by 2*span in y
  transmat_scale[0][0] = 1/(2*span);
  transmat_scale[1][1] = 1/(2*span);
  // transmat_scale[0][2] /= (2*span);
  // transmat_scale[1][2] /= (2*span);

  // Matrix3x3 transmat = transmat_reflection * transmat_translation * transmat_scale;
  // Matrix3x3 transmat = transmat_translation*transmat_scale;

  Matrix3x3 transmat = Matrix3x3::identity();
  transmat[0][2] = -1*(x - span); 
  transmat[1][2] = -1*(y - span);


  transmat[0][0] = 1/(2*span);
  transmat[1][1] = 1/(2*span);
  transmat[0][2] /= 2*span;
  transmat[1][2] /= 2*span;


  Matrix3x3 printg_mat  = this->svg_2_norm;

  printf("Checking svg_2_norm mat:\n");
  printf("%f, %f, %f\n", printg_mat[0][0], printg_mat[0][1], printg_mat[0][2]);
  printf("%f, %f, %f\n", printg_mat[1][0], printg_mat[1][1], printg_mat[1][2]);
  printf("%f, %f, %f\n", printg_mat[2][0], printg_mat[2][1], printg_mat[2][2]);

  printg_mat  = transmat_reflection;
  printf("Checking reflection mat:\n");
  printf("%f, %f, %f\n", printg_mat[0][0], printg_mat[0][1], printg_mat[0][2]);
  printf("%f, %f, %f\n", printg_mat[1][0], printg_mat[1][1], printg_mat[1][2]);
  printf("%f, %f, %f\n", printg_mat[2][0], printg_mat[2][1], printg_mat[2][2]);

  printg_mat  = transmat_translation;
  printf("Checking translation mat:\n");
  printf("%f, %f, %f\n", printg_mat[0][0], printg_mat[0][1], printg_mat[0][2]);
  printf("%f, %f, %f\n", printg_mat[1][0], printg_mat[1][1], printg_mat[1][2]);
  printf("%f, %f, %f\n", printg_mat[2][0], printg_mat[2][1], printg_mat[2][2]);

  printg_mat  = transmat_scale;
  printf("Checking scale mat:\n");
  printf("%f, %f, %f\n", printg_mat[0][0], printg_mat[0][1], printg_mat[0][2]);
  printf("%f, %f, %f\n", printg_mat[1][0], printg_mat[1][1], printg_mat[1][2]);
  printf("%f, %f, %f\n", printg_mat[2][0], printg_mat[2][1], printg_mat[2][2]);

  printg_mat  = transmat;
  printf("Checking fullscale mat:\n");
  printf("%f, %f, %f\n", printg_mat[0][0], printg_mat[0][1], printg_mat[0][2]);
  printf("%f, %f, %f\n", printg_mat[1][0], printg_mat[1][1], printg_mat[1][2]);
  printf("%f, %f, %f\n", printg_mat[2][0], printg_mat[2][1], printg_mat[2][2]);

  printf("~~~~~~~~~~~~~");

  
  // this->svg_2_norm = transmat_reflection;
  // this->svg_2_norm = transmat_scale*this->svg_2_norm;
  // this->svg_2_norm = transmat_translation*this->svg_2_norm;


  printg_mat  = this->svg_2_norm;
  printf("Checking copy mat:\n");
  printf("%f, %f, %f\n", printg_mat[0][0], printg_mat[0][1], printg_mat[0][2]);
  printf("%f, %f, %f\n", printg_mat[1][0], printg_mat[1][1], printg_mat[1][2]);
  printf("%f, %f, %f\n", printg_mat[2][0], printg_mat[2][1], printg_mat[2][2]);


  // CHECK!
  // issue: MODIFY MATRIX svg_2_norm NOT X,Y YOU DUMMY
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

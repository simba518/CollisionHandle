#ifndef _TETRAHEDRONVERTEXCON_H_
#define _TETRAHEDRONVERTEXCON_H_

class TetrahedronVertexCon{
  
public:
  static double functionValue(const double *dx, const double *vx);
  static void gradient(const double *dx, const double *vx, double *grad);
  
};

#endif /* _TETRAHEDRONVERTEXCON_H_ */

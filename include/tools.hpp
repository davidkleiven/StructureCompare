#ifndef TOOLS_H
#define TOOLS_H
#include <cmath>
#include "linalg.hpp"

namespace tools
{
  double dot( const Vector &v1, const Vector &v2 );

  void dot( const Matrix &m1,  const Matrix &m2, Matrix &out );

  void dot( const Matrix &mat, const Vector &vec, Vector &out );

  double length( const Vector &vec );

  double angle( const Vector &vec1, const Vector &vec2 );

  double det2x2( double a11, double a12, double a21, double a22 );

  double det3x3( const Matrix &mat );

  void inv3x3( const Matrix &mat, Matrix &out );
}

#endif

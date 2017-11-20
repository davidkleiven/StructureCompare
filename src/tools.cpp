#include "tools.hpp"
#include <cassert>

double tools::dot( const Vector &v1, const Vector &v2 )
{
  double dotprod = 0.0;
  for ( unsigned int i=0;i<v1.size();i++ )
  {
    dotprod += v1(i)*v2(i);
  }
  return dotprod;
};

void tools::dot( const Matrix &m1,  const Matrix &m2, Matrix &out )
{
  assert( m1.get_ncols() == m2.get_nrows() );
  assert( out.get_nrows() == m1.get_nrows() );
  assert( out.get_ncols() == m2.get_ncols() );
  
  for ( unsigned int i=0;i<out.get_nrows();i++ )
  {
    for ( unsigned int j=0;j<out.get_ncols();j++ )
    {
      out(i,j) = 0.0;
      for ( unsigned int k=0;k<m1.get_ncols();k++ )
      {
        out(i,j) += m1(i,k)*m2(k,j);
      }
    }
  }
};

void tools::dot( const Matrix &mat, const Vector &vec, Vector &out )
{
  for ( unsigned int i=0;i<out.size();i++ )
  {
    out(i) = 0.0;
    for ( unsigned int j=0;j<mat.get_ncols();j++ )
    {
      out(i) += mat(i,j)*vec(j);
    }
  }
};

double tools::length( const Vector &vec )
{
  return sqrt( dot(vec,vec) );
};

double tools::angle( const Vector &vec1, const Vector &vec2 )
{
  double dotprod = dot(vec1,vec2);
  double l1 = length(vec1);
  double l2 = length(vec2);
  return acos( dotprod/(l1*l2) );
};

double tools::det2x2( double a11, double a12, double a21, double a22 )
{
  return a11*a22 - a12*a21;
};

double tools::det3x3( const Matrix &mat )
{
  double V = 0.0;
  V += mat(0,0)*( mat(1,1)*mat(2,2) - mat(1,2)*mat(2,1) );
  V -= mat(0,1)*( mat(1,0)*mat(2,2) - mat(2,0)*mat(1,2) );
  V += mat(0,2)*( mat(1,0)*mat(2,1) - mat(2,0)*mat(1,1) );
  return V;
};

void tools::inv3x3( const Matrix &mat, Matrix &out )
{
  double det = det3x3( mat );
  out(0,0) = det2x2( mat(1,1), mat(1,2), mat(2,1), mat(2,2) )/det;
  out(0,1) = det2x2( mat(0,2), mat(0,1), mat(2,2), mat(2,1) )/det;
  out(0,2) = det2x2( mat(0,1), mat(0,2), mat(1,1), mat(1,2) )/det;

  out(1,0) = det2x2( mat(1,2), mat(1,0), mat(2,2), mat(2,0) )/det;
  out(1,1) = det2x2( mat(0,0), mat(0,2), mat(2,0), mat(2,2) )/det;
  out(1,2) = det2x2( mat(0,2), mat(0,0), mat(1,2), mat(1,0) )/det;

  out(2,0) = det2x2( mat(1,0), mat(1,1), mat(2,0), mat(2,1) )/det;
  out(2,1) = det2x2( mat(0,1), mat(0,0), mat(2,1), mat(2,0) )/det;
  out(2,2) = det2x2( mat(0,0), mat(0,1), mat(1,0), mat(1,1) )/det;
};

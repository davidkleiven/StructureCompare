#ifndef TOOLS_H
#define TOOLS_H
#include <cmath>

namespace tools
{
  template<unsigned int N>
  double dot( const std::array<double,N> v1, const std::array<double,N> &v2 )
  {
    double dot = 0.0;
    for ( unsigned int i=0;i<N;i++ )
    {
      dot += v1[i]*v2[i];
    }
  }

  template<unsigned int N>
  double dot( const std::array< std::array<double,N>, N> &m1,  const std::array< std::array<double,N>, N> &m2, std::array< std::array<double,N>, N> &out )
  {
    for ( unsigned int i=0;i<N;i++ )
    {
      for ( unsigned int j=0;j<N;j++ )
      {
        out[i][j] = 0.0;
        for ( unsigned int k=0;k<N;k++ )
        {
          out[i][j] += m1[i][k]*m2[k][j];
        }
      }
    }
  }

  template<unsigned int N>
  double length( const std::array<double,N> &vec )
  {
    return sqrt( dot<N>(vec,vec) );
  }

  template<unsigned int N>
  double angle( const std::array<double,N> &vec1, const std::array<double,N> &vec2 )
  {
    double dotprod = dot<N>(vec1,vec2);
    double l1 = length<N>(vec1);
    double l2 = length<N>(vec2);
    return acos( dotprod/(l1*l2) );
  }

  double det2x2( double a11, double a12, double a21, double a22 )
  {
    return a11*a22 - a12*a21;
  }

  double det3x3( const std::array< std::array<double,3>, 3> &mat )
  {
    double V = 0.0;
    V += mat[0][0]*( mat[1][1]*mat[2][2] - mat[1][2]*mat[2][1] );
    V -= mat[0][1]*( mat[1][0]*mat[2][2] - mat[2][0]*mat[1][2] );
    V += mat[0][2]*( mat[1][0]*mat[2][1] - mat[2][0]*mat[1][1] );
    return V;
  }

  double inv3x3( const std::array< std::array<double,3>, 3> &mat, std::array< std::array<double,3>, 3> &out )
  {
    double det = det3x3( mat );
    out[0][0] = det2x2( mat[1][1], mat[1][2],mat[2][1], mat[2][2] )/det;
    out[0][1] = det2x2( mat[0][2], mat[0][1],mat[2][2], mat[2][1] )/det;
    out[0][2] = det2x2( mat[0][1], mat[0][2],mat[1][1], mat[1][2] )/det;

    out[1][0] = det2x2( mat[1][2], mat[1][0],mat[2][2], mat[2][0] )/det;
    out[1][1] = det2x2( mat[0][0], mat[0][2],mat[2][0], mat[2][2] )/det;
    out[1][2] = det2x2( mat[0][2], mat[0][0],mat[1][2], mat[1][0] )/det;

    out[2][0] = det2x2( mat[1][0], mat[1][1],mat[2][0], mat[2][1] )/det;
    out[2][1] = det2x2( mat[0][1], mat[0][0],mat[2][1], mat[2][0] )/det;
    out[2][2] = det2x2( mat[0][0], mat[0][1],mat[1][0], mat[1][1] )/det;
  }
};

#endif

#ifndef LINALG_H
#define LINALG_H

#include <iostream>

class Vector;
class Matrix
{
public:
  Matrix( unsigned int nrows, unsigned int ncols );
  Matrix(){};
  Matrix(const Matrix &other);
  ~Matrix();

  void set_size( unsigned int nrows, unsigned int ncols );

  const double& operator()(unsigned int n, unsigned int m ) const;
  double& operator()(unsigned int n, unsigned int m );

  Matrix& operator +=( double num );
  Matrix& operator -=(double num);
  Matrix& operator=( const Matrix &other );

  Matrix T();

  Vector row( unsigned int row ) const;
  Vector col( unsigned int col ) const;

  void set_row( const Vector &vec, unsigned int row );
  void set_col( const Vector &vec, unsigned int col );

  unsigned int get_nrows() const { return nrows; };
  unsigned int get_ncols() const { return ncols; };
  friend std::ostream& operator <<( std::ostream& stream, const Matrix &mat );
protected:
  unsigned int nrows{0};
  unsigned int ncols{0};
  double *data{nullptr};
};

class Vector: public Matrix
{
public:
  Vector( unsigned int N ): Matrix(N,1){};

  const double& operator()(unsigned int n) const;
  double& operator()( unsigned int n );
  double& operator[](unsigned int n){ return (*this)(n);};
  const double& operator[](unsigned int n ) const { return (*this)(n);};

  unsigned int size() const { return nrows; }
};


#endif

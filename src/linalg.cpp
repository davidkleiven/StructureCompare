#include "linalg.hpp"

using namespace std;

Matrix::Matrix( unsigned int nrows, unsigned int ncols ):nrows(nrows),ncols(ncols)
{
  data = new double[nrows*ncols];
}

Matrix::~Matrix()
{
  delete [] data;
}

void Matrix::set_size( unsigned int new_nrows, unsigned int new_ncols )
{
  nrows = new_nrows;
  ncols = new_ncols;
  delete [] data;
  data = new double[nrows*ncols];

}

Matrix::Matrix( const Matrix& other )
{
  this->set_size( other.nrows, other.ncols );
  for ( unsigned int i=0;i<nrows*ncols;i++ )
  {
    data[i] = other.data[i];
  }
}

Matrix& Matrix::operator=( const Matrix &other )
{
  set_size( other.nrows, other.ncols );
  for ( unsigned int i=0;i<nrows*ncols;i++ )
  {
    data[i] = other.data[i];
  }
  return *this;
}

Matrix Matrix::T()
{
  Matrix newmat(*this);
  for ( unsigned int i=0;i<ncols;i++ )
  {
    for ( unsigned int j=0;j<nrows;j++ )
    {
      newmat(j,i) = (*this)(i,j);
    }
  }
  return newmat;
}

const double& Matrix::operator()( unsigned int n, unsigned int m ) const
{
  return data[nrows*m+n];
}

double& Matrix::operator()( unsigned int n, unsigned int m )
{
  return data[nrows*m+n];
}

Vector Matrix::row( unsigned int row ) const
{
  Vector newrow(nrows);
  for ( unsigned int i=0;i<ncols;i++ )
  {
    newrow(i) = (*this)(row,i);
  }
  return newrow;
}

void Matrix::set_row( const Vector &vec, unsigned int row )
{
  for ( unsigned int i=0;i<ncols;i++ )
  {
    (*this)(row,i) = vec(i);
  }
}

Matrix& Matrix::operator+=(double num)
{
  for ( unsigned int i=0;i<ncols;i++)
  {
    for ( unsigned int j=0;j<nrows;j++ )
    {
      (*this)(j,i) += num;
    }
  }
  return *this;
}

Matrix& Matrix::operator-=(double num)
{
  return *this += (-num);
}

ostream& operator << (ostream& out, const Matrix &mat )
{
  for ( unsigned int i=0;i<mat.nrows;i++ )
  {
    out << "[";
    for ( unsigned int j=0;j<mat.ncols;j++ )
    {
      out << mat(i,j) << " ";
    }
    out <<"]\n";
  }
  return out;
}

void Matrix::set_col( const Vector &vec, unsigned int col )
{
  for ( unsigned int i=0;i<nrows;i++ )
  {
    (*this)(i,col) = vec(i);
  }
}

Vector Matrix::col( unsigned int column ) const
{
  Vector vec(nrows);
  for ( unsigned int i=0;i<nrows;i++ )
  {
    vec(i) = (*this)(i,column);
  }
  return vec;
}

const double& Vector::operator()( unsigned int n ) const
{
  return data[n];
}

double& Vector::operator()( unsigned int n )
{
  return data[n];
}

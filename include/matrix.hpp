#ifndef INCLUDE_MATRIX_HPP_
#define INCLUDE_MATRIX_HPP_
#pragma once
#include<math.h>

template <class T>
class Matrix
{
private:
    T **matr;
    int rows, cols;
public:
    ~Matrix();
    int get_rows() const;
    int get_cols() const;
    Matrix(int rows_m, int cols_m);
    Matrix(const Matrix& copy);
    Matrix& operator=(const Matrix<T>& rhs);
    Matrix<T> operator+(Matrix<T>& m2);
    Matrix operator -(Matrix &m2);
    Matrix operator *(Matrix &m2);
    double det(Matrix m);
    Matrix deleterows_cols(Matrix m, int nrows, int ncols);
    Matrix Inversion();
    template<class Vec1>
    friend bool operator ==(const Matrix<Vec1> &m1, const Matrix<Vec1> &m2);
    template<class Vec2>
    friend bool operator !=(const Matrix<Vec2> &m1, const Matrix<Vec2> &m2);
    T* operator [](size_t i) const;
};

template<class T>
Matrix<T>::~Matrix(){
    for (int i = 0; i < rows; i++){
        delete[] matr[i];
    }
    delete[] matr;
}

template<class T>
int Matrix<T>::get_rows() const {
 return rows;
}

template<class T>
int Matrix<T>::get_cols() const {
 return cols;
}

template<class T>
Matrix<T>::Matrix(int rows_m, int cols_m)
    {
        
        rows = rows_m;
        cols = cols_m;
        matr = new T*[rows_m];
        for (int i = 0; i < rows_m; i++)
            matr[i] = new T[cols_m];
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols_m; j++)
                matr[i][j] = 0;
    }

template<class T>
    Matrix<T>::Matrix(const Matrix& cp)
        {
            rows = cp.rows;
            cols = cp.cols;
            matr = new T*[rows];
            for (int i = 0; i < rows; i++)
                matr[i] = new T[cols];
            for (int i = 0; i < rows; i++)
                for (int j = 0; j < cols; j++)
                    matr[i][j] = cp.matr[i][j];
}

template<class T>
T* Matrix<T>::operator [](size_t i) const{
 return matr[i];
}

template<class T>
Matrix<T>& Matrix<T>::operator=(const Matrix<T>& cp)
    {
        if (&cp != this)
            {
                rows = cp.rows;
                cols = cp.cols;
                for (int i = 0; i < rows; i++)
                    for (int j = 0; j < cols; j++)
                        matr[i][j] = cp[i][j];
            }
        return *this;
    }


template<class T>
Matrix<T> Matrix<T>::operator+(Matrix<T>& m2)
    {
        Matrix<T> m(0, 0);
        if (rows == m2.rows && cols == m2.cols)
            {
                Matrix<T> m(rows, cols);
                for (int i = 0; i < rows; i++)
                    for (int j = 0; j < cols; j++)
                        m[i][j] = matr[i][j] + m2[i][j];
                return m;
            }
        return m;
    }


template <class T>
Matrix<T> Matrix<T>::operator -(Matrix& m2){
    if (this->get_cols() != m2.get_cols() || this->get_rows() != m2.get_rows()){
  Matrix<T> res(0, 0);
  return res;
 }
 Matrix<T> res(rows, cols);
 for (int i = 0; i < rows ; i++){
  for (int j = 0 ; j < cols; j++){
   res[i][j] = (*this)[i][j] - m2[i][j];
  }
 }
 return res;
}

template <class T>
Matrix<T> Matrix<T>::operator *(Matrix& m2){
 if (this->get_cols() != m2.get_rows()){
  Matrix<T> res(0, 0);
  return res;
 }
 Matrix<T> res(rows, m2.get_cols());
 for (int i = 0 ; i < res.get_rows(); i++){
  for (int j = 0 ; j < res.get_cols(); j++){
   for (int k = 0 ; k < cols ; k++){
    res[i][j] += (*this)[i][k] * m2[k][j];
   }
  }
 }
 return res;
}

template<class T>
bool operator==(const Matrix<T> &m1, const Matrix<T> &m2){
 if (std::is_floating_point<T>::value){
  for (int i = 0; i < m1.get_rows(); i++){
   for (int j = 0; j < m1.get_cols(); j++){
    if (abs(m1[i][j] - m2[i][j]) > std::numeric_limits<double>::epsilon()){
     return false;
    }
   }
  }
  return true;
 }else {
  for (int i = 0; i < m1.get_rows(); i++){
   for (int j = 0; j < m1.get_cols(); j++){
    if (m1[i][j] != m2[i][j]){
     return false;
    }
   }
  }
  return true;
 }
}

template<>
bool operator==(const Matrix<float> &m1, const Matrix<float> &m2){
 for (int i = 0; i < m1.get_rows(); i++){
  for (int j = 0; j < m1.get_cols(); j++){
   if (abs(m1[i][j] - m2[i][j]) > std::numeric_limits<float>::epsilon()){
    return false;
   }
  }
 }
 return true;
}

template<class T>
bool operator!=(const Matrix<T> &m1, const Matrix<T> &m2){
 return !(m1 == m2);
}

template<class T>
Matrix<T> Matrix<T>::deleterows_cols(Matrix<T> mat, int nrows, int ncols){
 Matrix<T> res(mat.get_rows() - 1, mat.get_cols() - 1);
 int numrow = 0;
 int numcol = 0;
 for (int i = 0 ; i < mat.get_rows() ; i ++){
  if (i != nrows){
   for (int j = 0 ; j < mat.get_cols() ; j++){
    if (j != ncols){
     res[numrow][numcol] = mat[i][j];
     numcol += 1;
    }else{
     continue;
    }
   }
   numrow += 1;
   numcol = 0;
  }else{
   continue;
  }
 }
 return res;
}

template<class T>
double Matrix<T>::det(Matrix<T> matr){
 double db = 0;
 if (matr.get_rows() > 2){
  for (int i = 0 ; i < matr.get_rows() ; i++){
   db += pow(-1, i) * matr[0][i] * det(deleterows_cols(matr, 0, i));
  }
 }else{
  db = matr[0][0] * matr[1][1] - matr[0][1]*matr[1][0];
 }
 return db;
}


template<class T>
Matrix<T> Matrix<T>::Inversion(){
 Matrix<T> inv(this->rows, this->cols);
 double det_m = det(*this);
 for (int i = 0; i < inv.get_rows() ; i++){
  for (int j = 0; j < inv.get_cols() ; j++){
   inv[i][j] = pow(-1, i+j)*det(deleterows_cols(*this, i, j));
  }
 }
 Matrix<T> invT(inv.get_rows(), inv.get_cols());
 for (int i = 0 ; i < invT.get_rows(); i++){
  for (int j = 0 ; j < invT.get_cols(); j++){
   T temp = inv[j][i];
   T det_temp = temp / det_m;
   invT[i][j] = det_temp;
  }
 }
 return invT;
}

#endif // INCLUDE_MATRIX_HPP_

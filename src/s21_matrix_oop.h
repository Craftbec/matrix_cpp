#ifndef SRC_S21_MATRIX_H_
#define SRC_S21_MATRIX_H_

#define EPS 1e-7
#include <cmath>
#include <iostream>

class S21Matrix {
 private:
  int rows_, cols_;
  double** matrix_;
  S21Matrix S21_decrease_matrix(int irows, int jcolumns) const;

 public:
  S21Matrix();
  S21Matrix(int rows, int cols);
  S21Matrix(const S21Matrix& other);
  S21Matrix(S21Matrix&& other);
  ~S21Matrix();
  S21Matrix& operator=(const S21Matrix& other);
  S21Matrix& operator=(S21Matrix&& other) noexcept;

  int getRows() const;
  int getCols() const;
  void setRows(const int rows);
  void setCols(const int cols);

  bool EqMatrix(const S21Matrix& other) const;
  void SumMatrix(const S21Matrix& other);
  void SubMatrix(const S21Matrix& other);
  void MulNumber(const double num);
  void MulMatrix(const S21Matrix& other);
  S21Matrix Transpose() const;
  S21Matrix CalcComplements() const;
  double Determinant() const;
  S21Matrix InverseMatrix() const;

  void operator+=(const S21Matrix& other);
  void operator-=(const S21Matrix& other);
  void operator*=(const S21Matrix& other);
  void operator*=(double num);

  bool operator==(const S21Matrix& other) const;
  double& operator()(int i, int j);

  S21Matrix operator+(const S21Matrix& other) const;
  S21Matrix operator-(const S21Matrix& other) const;
  S21Matrix operator*(const S21Matrix& other) const;
  S21Matrix operator*(double num) const;

 private:
  void allocate_memory();
  void clean_memory();
  void copy_matrix(const S21Matrix& other);
};
S21Matrix operator*(double num, const S21Matrix& other);

#endif  // SRC_S21_MATRIX_H_
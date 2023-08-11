#include "s21_matrix_oop.h"

S21Matrix::S21Matrix() : rows_(5), cols_(5) { allocate_memory(); }

S21Matrix::S21Matrix(int rows, int cols) : rows_(rows), cols_(cols) {
  if (rows <= 0 || cols <= 0) {
    throw std::out_of_range("incorrect matrix");
  }
  allocate_memory();
}

S21Matrix::S21Matrix(const S21Matrix &other)
    : rows_(other.rows_), cols_(other.cols_) {
  allocate_memory();
  copy_matrix(other);
}

S21Matrix::S21Matrix(S21Matrix &&other)
    : rows_(other.rows_), cols_(other.cols_), matrix_(other.matrix_) {
  other.cols_ = other.rows_ = 0;
  other.matrix_ = nullptr;
}

S21Matrix::~S21Matrix() {
  clean_memory();
  cols_ = rows_ = 0;
}
int S21Matrix::getRows() const { return rows_; }

int S21Matrix::getCols() const { return cols_; }

void S21Matrix::setRows(const int rows) {
  if (rows <= 0) {
    throw std::out_of_range("Invalid value");
  }
  S21Matrix tmp(rows, cols_);
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols_; j++) {
      if (i >= rows_) {
        tmp.matrix_[i][j] = 0.0;
      } else {
        tmp.matrix_[i][j] = matrix_[i][j];
      }
    }
  }
  std::swap(*this, tmp);
}

void S21Matrix::setCols(const int cols) {
  if (cols <= 0) {
    throw std::out_of_range("Invalid value");
  }
  S21Matrix tmp(rows_, cols);
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols; j++) {
      if (j >= cols_) {
        tmp.matrix_[i][j] = 0;
      } else {
        tmp.matrix_[i][j] = matrix_[i][j];
      }
    }
  }
  std::swap(*this, tmp);
}

bool S21Matrix::EqMatrix(const S21Matrix &other) const {
  bool res = false;
  if (rows_ == other.rows_ && cols_ == other.cols_) {
    res = true;
    for (int i = 0; i < rows_; i++) {
      for (int j = 0; j < cols_; j++) {
        if (fabs(matrix_[i][j] - other.matrix_[i][j]) > EPS) {
          res = false;
          break;
        }
      }
    }
  }
  return res;
}

void S21Matrix::SumMatrix(const S21Matrix &other) {
  if (rows_ != other.rows_ || cols_ != other.cols_) {
    throw std::out_of_range("Different dimensions of Matrix");
  }
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] = matrix_[i][j] + other.matrix_[i][j];
    }
  }
}

void S21Matrix::SubMatrix(const S21Matrix &other) {
  if (rows_ != other.rows_ || cols_ != other.cols_) {
    throw std::out_of_range("Different dimensions of Matrix");
  }
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] = matrix_[i][j] - other.matrix_[i][j];
    }
  }
}

void S21Matrix::MulNumber(const double num) {
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] = matrix_[i][j] * num;
    }
  }
}

void S21Matrix::MulMatrix(const S21Matrix &other) {
  if (cols_ != other.rows_) {
    throw std::out_of_range(
        "The number of columns of the first matrix is not equal to the number "
        "of rows of the second matrix");
  }
  S21Matrix result(rows_, other.cols_);
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < other.cols_; j++) {
      result.matrix_[i][j] = 0;
      for (int k = 0; k < cols_; k++) {
        result.matrix_[i][j] += matrix_[i][k] * other.matrix_[k][j];
      }
    }
  }
  std::swap(*this, result);
}

S21Matrix S21Matrix::Transpose() const {
  S21Matrix result(cols_, rows_);
  for (int i = 0; i < result.rows_; i++) {
    for (int j = 0; j < result.cols_; j++) {
      result.matrix_[i][j] = matrix_[j][i];
    }
  }
  return result;
}

S21Matrix S21Matrix::S21_decrease_matrix(int irows, int jcolumns) const {
  S21Matrix tmp(rows_ - 1, cols_ - 1);
  for (int i = 0, a = 0; i < rows_; i++) {
    if (i != irows) {
      for (int j = 0, b = 0; j < cols_; j++) {
        if (j != jcolumns) {
          tmp.matrix_[a][b] = matrix_[i][j];
          b++;
        }
      }
      a++;
    }
  }
  return tmp;
}

S21Matrix S21Matrix::CalcComplements() const {
  if (rows_ != cols_) {
    throw std::out_of_range("Matrix isn ot square");
  }
  double s;
  S21Matrix result(rows_, cols_);
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      if (rows_ > 2) {
        S21Matrix tmp = S21_decrease_matrix(i, j);
        s = tmp.Determinant();
      } else {
        s = Determinant();
      }
      result(i, j) = pow(-1, i + j) * s;
    }
  }
  return result;
}

double S21Matrix::Determinant() const {
  if (rows_ != cols_) {
    throw std::out_of_range("Matrix isn ot square");
  }
  double result = 0;
  if (rows_ == 1) {
    result = matrix_[0][0];
  } else if (rows_ == 2) {
    result = matrix_[0][0] * matrix_[1][1] - matrix_[0][1] * matrix_[1][0];
  } else {
    for (int i = 0; i < rows_; i++) {
      S21Matrix tmp = S21_decrease_matrix(i, 0);
      double s = tmp.Determinant();
      result += matrix_[i][0] * s * pow(-1, i);
    }
  }
  return result;
}

S21Matrix S21Matrix::InverseMatrix() const {
  if (rows_ != cols_) {
    throw std::out_of_range("Matrix isn ot square");
  }
  double det = Determinant();
  if (fabs(det) <= EPS) {
    throw std::out_of_range("Matrix determinant is 0");
  }
  S21Matrix complement(rows_, cols_);
  complement = CalcComplements();
  S21Matrix result(cols_, rows_);
  result = complement.Transpose();
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      result.matrix_[i][j] *= (1 / det);
    }
  }
  return result;
}

bool S21Matrix::operator==(const S21Matrix &other) const {
  return EqMatrix(other);
}

void S21Matrix::operator+=(const S21Matrix &other) { SumMatrix(other); }

void S21Matrix::operator-=(const S21Matrix &other) { SubMatrix(other); }

void S21Matrix::operator*=(const S21Matrix &other) { MulMatrix(other); }

void S21Matrix::operator*=(double num) { MulNumber(num); }

S21Matrix &S21Matrix::operator=(const S21Matrix &other) {
  if (this == &other) return *this;
  clean_memory();
  rows_ = other.rows_;
  cols_ = other.cols_;
  allocate_memory();
  copy_matrix(other);
  return *this;
}

S21Matrix &S21Matrix::operator=(S21Matrix &&other) noexcept {
  if (this == &other) return *this;
  clean_memory();
  rows_ = other.rows_;
  cols_ = other.cols_;
  matrix_ = other.matrix_;
  other.rows_ = other.cols_ = 0;
  other.matrix_ = nullptr;
  return *this;
}

double &S21Matrix::operator()(int i, int j) {
  if (i >= rows_ || j >= cols_) {
    throw std::out_of_range("Index outside matrix");
  }

  return matrix_[i][j];
}

S21Matrix S21Matrix::operator+(const S21Matrix &other) const {
  S21Matrix result(*this);
  result.SumMatrix(other);
  return result;
}

S21Matrix S21Matrix::operator-(const S21Matrix &other) const {
  S21Matrix result(*this);
  result.SubMatrix(other);
  return result;
}

S21Matrix S21Matrix::operator*(const S21Matrix &other) const {
  S21Matrix result(*this);
  result.MulMatrix(other);
  return result;
}

S21Matrix S21Matrix::operator*(double num) const {
  S21Matrix result(*this);
  result.MulNumber(num);
  return result;
}

S21Matrix operator*(double num, const S21Matrix &other) {
  S21Matrix result(other);
  result.MulNumber(num);
  return result;
}

void S21Matrix::allocate_memory() {
  matrix_ = new double *[rows_];
  for (int i = 0; i < rows_; i++) {
    matrix_[i] = new double[cols_];
  }
}

void S21Matrix::clean_memory() {
  if (matrix_) {
    for (int i = 0; i < rows_; i++) {
      delete[] matrix_[i];
    }
    delete[] matrix_;
    matrix_ = nullptr;
  }
}

void S21Matrix::copy_matrix(const S21Matrix &other) {
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] = other.matrix_[i][j];
    }
  }
}
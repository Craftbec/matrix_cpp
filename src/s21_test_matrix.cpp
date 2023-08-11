#include <gtest/gtest.h>

#include "s21_matrix_oop.h"

#define SUCCESS 1
#define FAILURE 0

unsigned int seed = 0;
using namespace std;

TEST(Constructor_matrix, constructor) {
  const int rows = rand_r(&seed) % 100 + 1;
  const int cols = rand_r(&seed) % 100 + 1;
  S21Matrix A(rows, cols);

  double k = 0.1;
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      double rand_val = rand_r(&seed) % 2001 - 1000 + k;
      A(i, j) = rand_val;
      k += 0.000001;
    }
  }
  S21Matrix Copy_A(A);
  S21Matrix B = std::move(A);
  ASSERT_TRUE(Copy_A.EqMatrix(B) == SUCCESS);
}

TEST(Set_matrix, setrows_more) {
  S21Matrix A(2, 2);
  A(0, 0) = 1;
  A(0, 1) = 2;
  A(1, 0) = 3;
  A(1, 1) = 4;
  A.setRows(3);
  S21Matrix B(3, 2);
  B(0, 0) = 1;
  B(0, 1) = 2;
  B(1, 0) = 3;
  B(1, 1) = 4;
  B(2, 0) = 0;
  B(2, 1) = 0;
  ASSERT_TRUE(A.EqMatrix(B) == SUCCESS);
}

TEST(Set_matrix, setrows_less) {
  S21Matrix A(2, 2);
  A(0, 0) = 1;
  A(0, 1) = 2;
  A(1, 0) = 3;
  A(1, 1) = 4;
  A.setRows(1);
  S21Matrix B(1, 2);
  B(0, 0) = 1;
  B(0, 1) = 2;
  ASSERT_TRUE(A.EqMatrix(B) == SUCCESS);
}

TEST(Set_matrix, setcols_more) {
  S21Matrix A(2, 2);
  A(0, 0) = 1;
  A(0, 1) = 2;
  A(1, 0) = 3;
  A(1, 1) = 4;
  A.setCols(3);
  S21Matrix B(2, 3);
  B(0, 0) = 1;
  B(0, 1) = 2;
  B(0, 2) = 0;
  B(1, 0) = 3;
  B(1, 1) = 4;
  B(1, 2) = 0;

  ASSERT_TRUE(A.EqMatrix(B) == SUCCESS);
}

TEST(Set_matrix, setcols_less) {
  S21Matrix A(2, 2);
  A(0, 0) = 1;
  A(0, 1) = 2;
  A(1, 0) = 3;
  A(1, 1) = 4;
  A.setCols(1);
  S21Matrix B(2, 1);
  B(0, 0) = 1;
  B(1, 0) = 3;
  ASSERT_TRUE(A.EqMatrix(B) == SUCCESS);
}

TEST(Operator_matrix, eq_matrix) {
  const int rows = rand_r(&seed) % 100 + 1;
  const int cols = rand_r(&seed) % 100 + 1;

  S21Matrix A;
  A.setRows(rows);
  A.setCols(cols);
  S21Matrix B(rows, cols);
  double k = 0.1;
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      double rand_val = rand_r(&seed) % 2001 - 1000 + k;
      A(i, j) = rand_val;
      B(i, j) = rand_val;
      k += 0.000001;
    }
  }
  ASSERT_TRUE(A.EqMatrix(B) == SUCCESS);
}

TEST(Operator_matrix, not_eq_matrix) {
  const int rows = rand_r(&seed) % 100 + 1;
  const int cols = rand_r(&seed) % 100 + 1;

  S21Matrix A(rows, cols);
  S21Matrix B(rows, cols);
  double k = 0.1;
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      double rand_val = rand_r(&seed) % 2001 - 1000 + k;
      A(i, j) = rand_val;
      B(i, j) = 2 * rand_val;
      k += 0.000001;
    }
  }
  ASSERT_TRUE(A.EqMatrix(B) == FAILURE);
}

TEST(Operator_matrix, not_rows_eq_matrix) {
  const int rows = rand_r(&seed) % 100 + 1;
  const int cols = rand_r(&seed) % 100 + 1;

  S21Matrix A(rows, cols);
  S21Matrix B((rows + 10), cols);

  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      double rand_val = i * 10.54;
      A(i, j) = rand_val;
    }
  }
  for (int i = 0; i < rows + 10; i++) {
    for (int j = 0; j < cols; j++) {
      double rand_val = i * 10.54;
      B(i, j) = rand_val;
    }
  }
  ASSERT_TRUE(A.EqMatrix(B) == FAILURE);
}

TEST(Operator_matrix, sum_matrix) {
  const int rows = rand_r(&seed) % 100 + 1;
  const int cols = rand_r(&seed) % 100 + 1;

  S21Matrix A(rows, cols);
  S21Matrix B(rows, cols);
  S21Matrix Check_m(rows, cols);
  double k = 0.1;
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      double rand_val = rand_r(&seed) % 2001 - 1000 + k;
      A(i, j) = rand_val;
      B(i, j) = 2 * rand_val + 0.015;
      Check_m(i, j) = 3 * rand_val + 0.015;
      k += 0.000001;
    }
  }
  A.SumMatrix(B);
  ASSERT_TRUE(A.EqMatrix(Check_m) == SUCCESS);
}

TEST(incorrect, calc_incorr_matrix) {
  S21Matrix A(1, 3);
  A(0, 0) = 1;
  try {
    A.CalcComplements();
  } catch (const std::exception& e) {
    A(0, 0) = 100;
  }

  GTEST_ASSERT_EQ(A(0, 0), 100);
}

TEST(incorrect, calc_one_to_one) {
  S21Matrix A(1, 1);
  S21Matrix B(1, 1);
  A(0, 0) = 4;
  B = A.CalcComplements();

  GTEST_ASSERT_EQ(B(0, 0), 4);
}

TEST(incorrect, determ_not_square) {
  S21Matrix A(1, 2);
  A(0, 0) = 1;
  A(0, 1) = 1;
  try {
    A.Determinant();
  } catch (const std::exception& e) {
    A(0, 0) = 100;
  }
  GTEST_ASSERT_EQ(A(0, 0), 100);
}

TEST(incorrect, determ_one_x_one) {
  S21Matrix A(1, 1);
  A(0, 0) = 1;
  int res = A.Determinant();
  GTEST_ASSERT_EQ(res, 1);
}

TEST(incorrect, determ_one_x_two) {
  S21Matrix A(2, 2);
  A(0, 0) = 1;
  A(0, 1) = 2;
  A(1, 0) = 3;
  A(1, 1) = 4;
  int res = A.Determinant();
  GTEST_ASSERT_EQ(res, -2);
}

TEST(incorrect, incorrect_invers) {
  S21Matrix A(1, 1);
  A(0, 0) = 0;
  try {
    A.InverseMatrix();
  } catch (const std::exception& e) {
    A(0, 0) = 100;
  }
  GTEST_ASSERT_EQ(A(0, 0), 100);
}

TEST(incorrect, incorrect_iterator) {
  S21Matrix A(1, 1);
  A(0, 0) = 18;
  try {
    A(0, 9) = 0;
  } catch (const std::exception& e) {
    A(0, 0) = 100;
  }
  GTEST_ASSERT_EQ(A(0, 0), 100);
}

TEST(incorrect, sum_incorr_matrix) {
  const int rows = 1, cols = 1;
  S21Matrix A(rows, cols);
  S21Matrix B(rows + 1, cols);
  A(0, 0) = 1;
  B(0, 0) = 1;
  try {
    A.SumMatrix(B);
  } catch (const std::exception& e) {
  }

  GTEST_ASSERT_EQ(A(0, 0), 1);
}

TEST(incorrect, sub_incorr_matrix) {
  const int rows = 1, cols = 1;
  S21Matrix A(rows, cols);
  S21Matrix B(rows + 1, cols);
  A(0, 0) = 1;
  B(0, 0) = 1;
  try {
    A.SubMatrix(B);
  } catch (const std::exception& e) {
  }

  GTEST_ASSERT_EQ(A(0, 0), 1);
}

TEST(incorrect, mul_incorr_matrix) {
  const int rows = 1, cols = 1;
  S21Matrix A(rows, cols);
  S21Matrix B(rows + 1, cols);
  A(0, 0) = 1;
  B(0, 0) = 1;
  try {
    A.MulMatrix(B);
  } catch (const std::exception& e) {
  }

  GTEST_ASSERT_EQ(A(0, 0), 1);
}

TEST(Operator_matrix, sub_matrix) {
  const int rows = rand_r(&seed) % 100 + 1;
  const int cols = rand_r(&seed) % 100 + 1;

  S21Matrix A(rows, cols);
  S21Matrix B(rows, cols);
  S21Matrix Check_m(rows, cols);
  double k = 0.1;
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      double rand_val = rand_r(&seed) % 2001 - 1000 + k;
      A(i, j) = rand_val;
      B(i, j) = 2 * rand_val + 0.015;
      Check_m(i, j) = -rand_val - 0.015;
      k += 0.000001;
    }
  }
  A.SubMatrix(B);
  ASSERT_TRUE(A.EqMatrix(Check_m) == SUCCESS);
}

TEST(Operator_matrix, mult_number) {
  const int rows = rand_r(&seed) % 100 + 1;
  const int cols = rand_r(&seed) % 100 + 1;
  S21Matrix A(rows, cols);
  S21Matrix Check_m(rows, cols);

  double k = 0.1;
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      double rand_val = rand_r(&seed) % 2001 - 1000 + k;
      A(i, j) = rand_val;
      Check_m(i, j) = rand_val * 0.345;
      k += 0.000001;
    }
  }
  A.MulNumber(0.345);
  ASSERT_TRUE(A.EqMatrix(Check_m) == SUCCESS);
}

TEST(Operator_matrix, mult_matrix) {
  const int rows = rand_r(&seed) % 100 + 1;
  const int cols = rand_r(&seed) % 100 + 1;
  S21Matrix A(rows, cols);
  S21Matrix B(cols, rows + 10);
  S21Matrix Check_m(rows, rows + 10);

  double k = 0.1;
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      double rand_val = rand_r(&seed) % 2001 - 1000 + k;
      A(i, j) = rand_val;
      k += 0.000001;
    }
  }
  for (int i = 0; i < cols; i++) {
    for (int j = 0; j < rows + 10; j++) {
      double rand_val = rand_r(&seed) % 2001 - 1000 + k;
      B(i, j) = rand_val + k + 0.32;
      k += 0.000001;
    }
  }
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < rows + 10; j++) {
      Check_m(i, j) = 0;
      for (int s = 0; s < cols; s++) {
        Check_m(i, j) += A(i, s) * B(s, j);
      }
    }
  }

  A.MulMatrix(B);
  ASSERT_TRUE(A.EqMatrix(Check_m) == SUCCESS);
}

TEST(Operator_matrix, transpose_matrix) {
  const int rows = rand_r(&seed) % 100 + 1;
  const int cols = rand_r(&seed) % 100 + 1;
  S21Matrix A(rows, cols);
  S21Matrix A_trans(cols, rows);
  S21Matrix Check_m(cols, rows);

  double k = 0.1;
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      double rand_val = rand_r(&seed) % 2001 - 1000 + k;
      A(i, j) = rand_val;
      Check_m(j, i) = rand_val;
      k += 0.000001;
    }
  }
  A_trans = A.Transpose();
  ASSERT_TRUE(A_trans.EqMatrix(Check_m) == SUCCESS);
}

TEST(Operator_matrix, calc_complements_mat) {
  S21Matrix A(3, 3);
  S21Matrix A_calc(3, 3);
  S21Matrix check(3, 3);

  A(0, 0) = 1;
  A(0, 1) = 2;
  A(0, 2) = 3;
  A(1, 0) = 0;
  A(1, 1) = 4;
  A(1, 2) = 2;
  A(2, 0) = 5;
  A(2, 1) = 2;
  A(2, 2) = 1;

  check(0, 0) = 0;
  check(0, 1) = 10;
  check(0, 2) = -20;
  check(1, 0) = 4;
  check(1, 1) = -14;
  check(1, 2) = 8;
  check(2, 0) = -8;
  check(2, 1) = -2;
  check(2, 2) = 4;

  A_calc = A.CalcComplements();
  ASSERT_TRUE(A_calc.EqMatrix(check) == SUCCESS);
}

TEST(Operator_matrix, determinant_mat) {
  S21Matrix A(3, 3);
  A(0, 0) = 1;
  A(0, 1) = 2;
  A(0, 2) = 3;
  A(1, 0) = 0;
  A(1, 1) = 4;
  A(1, 2) = 2;
  A(2, 0) = 5;
  A(2, 1) = 2;
  A(2, 2) = 1;

  double result = A.Determinant();
  ASSERT_TRUE(result == -40);
}

TEST(Operator_matrix, inverse_mat) {
  S21Matrix A(3, 3);
  S21Matrix A_inverse(3, 3);
  S21Matrix check(3, 3);

  A(0, 0) = 1;
  A(0, 1) = 2;
  A(0, 2) = 3;
  A(1, 0) = 0;
  A(1, 1) = 4;
  A(1, 2) = 2;
  A(2, 0) = 5;
  A(2, 1) = 2;
  A(2, 2) = 1;

  check(0, 0) = 0;
  check(0, 1) = -0.1;
  check(0, 2) = 0.2;
  check(1, 0) = -0.25;
  check(1, 1) = 0.35;
  check(1, 2) = 0.05;
  check(2, 0) = 0.5;
  check(2, 1) = -0.2;
  check(2, 2) = -0.1;
  A_inverse = A.InverseMatrix();
  ASSERT_TRUE(A_inverse.EqMatrix(check) == SUCCESS);
}

TEST(Overload_matrix, operator_plus) {
  const int rows = rand_r(&seed) % 100 + 1;
  const int cols = rand_r(&seed) % 100 + 1;

  S21Matrix A(rows, cols);
  S21Matrix B(rows, cols);
  S21Matrix Check_m(rows, cols);
  double k = 0.1;
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      double rand_val = rand_r(&seed) % 2001 - 1000 + k;
      A(i, j) = rand_val;
      B(i, j) = 2 * rand_val + 0.015;
      k += 0.000001;
    }
  }
  S21Matrix Copy_A(A);
  Check_m = Copy_A + B;
  A.SumMatrix(B);
  ASSERT_TRUE(A.EqMatrix(Check_m) == SUCCESS);
}

TEST(Overload_matrix, operator_minus) {
  const int rows = rand_r(&seed) % 100 + 1;
  const int cols = rand_r(&seed) % 100 + 1;

  S21Matrix A(rows, cols);
  S21Matrix B(rows, cols);
  S21Matrix Check_m(rows, cols);
  double k = 0.1;
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      double rand_val = rand_r(&seed) % 2001 - 1000 + k;
      A(i, j) = rand_val;
      B(i, j) = 2 * rand_val + 0.015;
      k += 0.000001;
    }
  }
  S21Matrix Copy_A(A);
  Check_m = Copy_A - B;
  A.SubMatrix(B);
  ASSERT_TRUE(A.EqMatrix(Check_m) == SUCCESS);
}

TEST(Overload_matrix, operator_mult_matrix) {
  const int rows = rand_r(&seed) % 100 + 1;
  const int cols = rand_r(&seed) % 100 + 1;
  S21Matrix A(rows, cols);
  S21Matrix B(cols, rows + 10);
  S21Matrix Check_m(rows, rows + 10);

  double k = 0.1;
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      double rand_val = rand_r(&seed) % 2001 - 1000 + k;
      A(i, j) = rand_val;
      k += 0.000001;
    }
  }
  for (int i = 0; i < cols; i++) {
    for (int j = 0; j < rows + 10; j++) {
      double rand_val = rand_r(&seed) % 2001 - 1000 + k;
      B(i, j) = rand_val + k + 0.32;
      k += 0.000001;
    }
  }

  S21Matrix Copy_A(A);
  Check_m = Copy_A * B;
  A.MulMatrix(B);
  ASSERT_TRUE(A.EqMatrix(Check_m) == SUCCESS);
}

TEST(Overload_matrix, operator_mult_number_1) {
  const int rows = rand_r(&seed) % 100 + 1;
  const int cols = rand_r(&seed) % 100 + 1;
  S21Matrix A(rows, cols);
  S21Matrix Check_m(rows, cols);

  double k = 0.1;
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      double rand_val = rand_r(&seed) % 2001 - 1000 + k;
      A(i, j) = rand_val;
      k += 0.000001;
    }
  }
  Check_m = A * 0.345;
  A.MulNumber(0.345);
  ASSERT_TRUE(A.EqMatrix(Check_m) == SUCCESS);
}

TEST(Overload_matrix, operator_eq_1) {
  const int rows = rand_r(&seed) % 100 + 1;
  const int cols = rand_r(&seed) % 100 + 1;

  S21Matrix A(rows, cols);
  S21Matrix B(rows, cols);
  double k = 0.1;
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      double rand_val = rand_r(&seed) % 2001 - 1000 + k;
      A(i, j) = rand_val;
      B(i, j) = rand_val;
      k += 0.000001;
    }
  }
  ASSERT_TRUE(A.EqMatrix(B) == (A == B));
}

TEST(Overload_matrix, operator_eq_2) {
  const int rows = rand_r(&seed) % 100 + 1;
  const int cols = rand_r(&seed) % 100 + 1;

  S21Matrix A(rows, cols);
  S21Matrix B(rows, cols);
  double k = 0.1;
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      double rand_val = rand_r(&seed) % 2001 - 1000 + k;
      A(i, j) = rand_val;
      B(i, j) = 2 * rand_val;
      k += 0.000001;
    }
  }
  ASSERT_TRUE(A.EqMatrix(B) == (A == B));
}

TEST(Overload_matrix, operator_plus_eq) {
  const int rows = rand_r(&seed) % 100 + 1;
  const int cols = rand_r(&seed) % 100 + 1;

  S21Matrix A(rows, cols);
  S21Matrix B(rows, cols);
  S21Matrix Check_m(rows, cols);
  double k = 0.1;
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      double rand_val = rand_r(&seed) % 2001 - 1000 + k;
      A(i, j) = rand_val;
      B(i, j) = 2 * rand_val + 0.015;
      Check_m(i, j) = 3 * rand_val + 0.015;
      k += 0.000001;
    }
  }
  S21Matrix Copy_A(A);
  Copy_A += B;
  A.SumMatrix(B);
  ASSERT_TRUE((Copy_A == A) == SUCCESS);
}

TEST(Overload_matrix, operator_minus_eq) {
  const int rows = rand_r(&seed) % 100 + 1;
  const int cols = rand_r(&seed) % 100 + 1;

  S21Matrix A(rows, cols);
  S21Matrix B(rows, cols);
  S21Matrix Check_m(rows, cols);
  double k = 0.1;
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      double rand_val = rand_r(&seed) % 2001 - 1000 + k;
      A(i, j) = rand_val;
      B(i, j) = 2 * rand_val + 0.015;
      Check_m(i, j) = 3 * rand_val + 0.015;
      k += 0.000001;
    }
  }
  S21Matrix Copy_A(A);
  Copy_A -= B;
  A.SubMatrix(B);
  ASSERT_TRUE((Copy_A == A) == SUCCESS);
}

TEST(Overload_matrix, operator_mult_number_eq) {
  const int rows = rand_r(&seed) % 100 + 1;
  const int cols = rand_r(&seed) % 100 + 1;
  S21Matrix A(rows, cols);
  S21Matrix Check_m(rows, cols);

  double k = 0.1;
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      double rand_val = rand_r(&seed) % 2001 - 1000 + k;
      A(i, j) = rand_val;
      Check_m(i, j) = rand_val * 0.345;
      k += 0.000001;
    }
  }
  S21Matrix Copy_A(A);
  Copy_A *= 0.345;
  A.MulNumber(0.345);
  ASSERT_TRUE((A == Copy_A) == SUCCESS);
}

TEST(Overload_matrix, operator_mult_matrix_eq) {
  const int rows = rand_r(&seed) % 100 + 1;
  const int cols = rand_r(&seed) % 100 + 1;
  S21Matrix A(rows, cols);
  S21Matrix B(cols, rows + 10);
  S21Matrix Check_m(rows, rows + 10);

  double k = 0.1;
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      double rand_val = rand_r(&seed) % 2001 - 1000 + k;
      A(i, j) = rand_val;
      k += 0.000001;
    }
  }
  for (int i = 0; i < cols; i++) {
    for (int j = 0; j < rows + 10; j++) {
      double rand_val = rand_r(&seed) % 2001 - 1000 + k;
      B(i, j) = rand_val + k + 0.32;
      k += 0.000001;
    }
  }

  S21Matrix Copy_A(A);
  Copy_A *= B;
  A.MulMatrix(B);
  ASSERT_TRUE((A == Copy_A) == SUCCESS);
}

TEST(Get_matrix, get_rows) {
  const int rows = rand_r(&seed) % 100 + 1;
  const int cols = rand_r(&seed) % 100 + 1;
  S21Matrix A(rows, cols);

  ASSERT_TRUE(rows == A.getRows());
}

TEST(Get_matrix, get_cols) {
  const int rows = rand_r(&seed) % 100 + 1;
  const int cols = rand_r(&seed) % 100 + 1;
  S21Matrix A(rows, cols);

  ASSERT_TRUE(cols == A.getCols());
}

TEST(Set_matrix, set_rows) {
  const int rows = rand_r(&seed) % 100 + 1;
  const int cols = rand_r(&seed) % 100 + 1;
  S21Matrix A(rows, cols);
  A.setRows(rows + 10);

  ASSERT_TRUE(rows + 10 == A.getRows());
}

TEST(Set_matrix, set_cols) {
  const int rows = rand_r(&seed) % 100 + 1;
  const int cols = rand_r(&seed) % 100 + 1;
  S21Matrix A(rows, cols);
  A.setCols(cols + 10);

  ASSERT_TRUE(cols + 10 == A.getCols());
}

TEST(incorrect, set_incorrect_rows) {
  S21Matrix test(3, 3);
  try {
    test.setRows(-1);
  } catch (const std::exception& e) {
  }
  GTEST_ASSERT_EQ(test.getRows(), 3);
}

TEST(incorrect, set_incorrect_cols) {
  S21Matrix test(3, 3);
  try {
    test.setCols(-1);
  } catch (const std::exception& e) {
  }
  GTEST_ASSERT_EQ(test.getRows(), 3);
}

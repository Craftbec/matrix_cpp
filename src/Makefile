CC			=	g++ -Wall -Wextra -Werror -std=c++17
C++			=   -lstdc++ -lm
TEST_FLAGS 	= 	-lcheck -lgtest -lgtest_main -lpthread 
SRC			=	s21_matrix_oop.cpp
LIB			=	s21_matrix_oop.a
OBJ			= 	$(SRC:%.cpp=%.o)

all:	clean s21_matrix_oop.a  test

%.o:	%.cpp
		$(CC) -c $<

s21_matrix_oop.a:	$(OBJ) s21_matrix_oop.h
		ar rvs $(LIB) $^

clean:
		rm -rf $(LIB)
		rm -rf ./test
		rm -rf test.info test_output/
		rm -rf *.o  *.gcno *.gcda
		rm -rf RESULT_VALGRIND.txt test.dSYM

test:
		$(CC) -lgcov -coverage s21_test_matrix.cpp $(SRC) $(TEST_FLAGS) $(C++) -o test
		./test

gcov_report: test
		lcov -t "test" -o test.info -c -d .
		genhtml -o test_output test.info
		@cd ./test_output ; open index.html

check:
	cppcheck --enable=all --language=c++ --suppress=missingIncludeSystem . *.cpp *.h
	cp ../materials/linters/.clang-format .
	clang-format -n *.cpp *.h

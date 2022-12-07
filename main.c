#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include "tinyexpr.c"
#include "tinyexpr.h"

#define MAX_SIZE 100 //Максимальный размер вводимой строки

const double INIT_APP = 0.5; //Начальные значения иксов
const double EPS = 0.0000001;
const int CNT = 1000; //Количество итераций
char * X_N[MAX_SIZE];

typedef double ** num_type;
typedef char ** str_type;

void allocate_memory(void *** matrix, int row, int col, int size_type, int size_ptr) {
    (*matrix) = malloc(size_ptr * (row));
    for (int i = 0; i < row; i++) {
        (*matrix)[i] = malloc(col * size_type);
    }
}

num_type mat_mat_mul(num_type matrix_a, int row_a, int col_a,
                     num_type matrix_b, int row_b, int col_b) {
    num_type ans;
    allocate_memory(&ans, row_a, col_b, sizeof(double), sizeof(double *));
    for (int i = 0; i < row_a; ++i) {
        for (int j = 0; j < col_b; ++j) {
            ans[i][j] = 0;
            for (int k = 0; k < col_a; ++k) {
                ans[i][j] += matrix_a[i][k] * matrix_b[k][j];
            }
        }
    }
    return ans;
}

num_type mat_num_mul(num_type matrix, int row, int col, double number) {
    num_type ans;
    allocate_memory(&ans, row, col, sizeof(double), sizeof(double *));
    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < col; ++j) {
            ans[i][j] = matrix[i][j] * number;
        }
    }
    return ans;
}

num_type mat_mat_add(num_type matrix_a, num_type matrix_b, int row, int col) {
    num_type ans;
    allocate_memory(&ans, row, col, sizeof(double), sizeof(double *));
    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < col; ++j) {
            ans[i][j] = matrix_a[i][j] + matrix_b[i][j];
        }
    }
    return ans;
}

num_type get_transposed_matrix(num_type matrix, int row, int col) {
    num_type ans;
    allocate_memory(&ans, col, row, sizeof(double), sizeof(double *));
    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < col; ++j) {
            ans[j][i] = matrix[i][j];
        }
    }
    return ans;
}

num_type get_cofactor(num_type matrix, int order, int row, int col) {
    num_type ans;
    allocate_memory(&ans, order - 1, order - 1, sizeof(double), sizeof(double *));
    for (int i = 0; i < order - 1; ++i) {
        for (int j = 0; j < order - 1; ++j) {
            int x, y;
            x = i + (i >= row);
            y = j + (j >= col);
            ans[i][j] = matrix[x][y];
        }
    }
    return ans;
}

double get_determinant(num_type matrix, int order) {
    if (order == 1) {
        return matrix[0][0];
    }
    double ans = 0;
    for (int i = 0; i < order; ++i) {
        int sign = ((i & 1) ? -1 : 1);
        ans += (double)sign * matrix[0][i] * get_determinant(get_cofactor(matrix, order, 0, i), order - 1);
    }
    return ans;
}

num_type get_union_matrix(num_type matrix, int order) {
    num_type ans;
    allocate_memory(&ans, order, order, sizeof(double), sizeof(double *));
    if (order == 1) {
        ans[0][0] = matrix[0][0];
        return ans;
    }
    for (int i = 0; i < order; ++i) {
        for (int j = 0; j < order; ++j) {
            int sign = ((i + j & 1) ? -1 : 1);
            ans[i][j] = (double)sign * get_determinant(get_cofactor(matrix, order, i, j), order - 1);
        }
    }
    ans = get_transposed_matrix(ans, order, order);
    return ans;
}

num_type get_inversed_matrix(num_type matrix, int order) {
    num_type ans;
    allocate_memory(&ans, order, order, sizeof(double), sizeof(double *));
    ans = mat_num_mul(get_union_matrix(matrix, order), order, order, 1.0 / get_determinant(matrix, order));
    return ans;
}

num_type get_func_value(num_type vector_x, str_type vector_func, int n) {
    num_type ans;
    allocate_memory(&ans, n, 1, sizeof(double), sizeof(double *));
    te_variable vars[n];
    for (int i = 0; i < n; ++i) {
        vars[i] = (te_variable){X_N[i], &(vector_x[i][0])};
    }
    for (int i = 0; i < n; ++i) {
        int err = 0;
        te_expr *expr = te_compile(vector_func[i], vars, n, &err);
        ans[i][0] = te_eval(expr);
    }
    return ans;
}

double get_partial_derivative(num_type vector_x, str_type vector_func, int row, int col, int n) {
    te_variable vars[n];
    for (int i = 0; i < n; ++i) {
        vars[i] = (te_variable){X_N[i], &(vector_x[i][0])};
    }
    int err;
    te_expr *expr = te_compile(vector_func[row], vars, n, &err);
    const double upper = te_eval(expr);
    vector_x[col][0] -= EPS;
    const double lower = te_eval(expr);
    te_free(expr);
    return (upper - lower) / EPS;
}

num_type get_mat_yakobi(int n, num_type vector_x, str_type vector_func) {
    num_type ans;
    allocate_memory(&ans, n, n, sizeof(double), sizeof(double *));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            ans[i][j] = get_partial_derivative(vector_x, vector_func, i, j, n);
        }
    }
    return ans;
}

char * get_str_x(int id) {
    char * ans;
    int len = 1, tmp = id;
    while (tmp /= 10) len++;
    ans = (char *)malloc(sizeof(char) * (len + 2));
    ans[0] = 'x';
    ans[len + 1] = '\0';
    while (id) {
        ans[len--] = id % 10 + '0';
        id /= 10;
    }
    return ans;
}

void make_xes() {
    for (int i = 0; i < MAX_SIZE; ++i) {
        X_N[i] = get_str_x(i + 1);
    }
}

void init(int * n, num_type * matrix_numb, str_type * matrix_str) {
    make_xes();
    printf("Input n:\n");
    scanf("%d", n);
    while (!getchar());
    allocate_memory(matrix_numb, *n, 1, sizeof(double), sizeof(double *));
    allocate_memory(matrix_str, *n, MAX_SIZE, sizeof(char), sizeof(char *));
    for (int i = 0; i < *n; ++i) {
        (*matrix_numb)[i][0] = INIT_APP;
    }
    printf("Input %d functions separated by enters\n(each must contain arguments of the form xN, N belongs to the set of natural numbers):\n", *n);
    for (int i = 0; i < *n; ++i) {
        int id = 0;
        printf("f%d(x1..x%d)=", i + 1, *n);
        do {
            (*matrix_str)[i][id] = getchar();
        } while ((*matrix_str)[i][id++] != '\n');
        (*matrix_str)[i][id - 1] = '\0';
    }
}

void output_matrix(num_type matrix, int row, int col) {
    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < col; ++j) {
            if (j == 0) {
                if (row == 1)
                    printf("( ");
                else if (i == 0)
                    printf("/ ");
                else if (i == row - 1)
                    printf("\\ ");
                else
                    printf("| ");
            }
            printf("%.6lf", matrix[i][j]);
            if (i == row - 1 && j == col - 1)
                printf("  ");
            else
                printf(", ");
            if (j == col - 1) {
                if (row == 1)
                    printf(")");
                else if (i == 0)
                    printf("\\");
                else if (i == row - 1)
                    printf("/");
                else
                    printf("|");
            }
        }
        puts("");
    }
}

num_type iteration(int n, num_type vector_x, str_type vector_func) {
    num_type ans;
    allocate_memory(&ans, n, 1, sizeof(double), sizeof(double *));
    num_type m_yakobi = get_mat_yakobi(n, vector_x, vector_func);
    ans = mat_mat_add(vector_x, mat_num_mul(mat_mat_mul(get_inversed_matrix(m_yakobi, n), n, n, get_func_value(vector_x, vector_func, n), n, 1), n, 1, -1), n, 1);
    return ans;
}

int main() {
    num_type x_k;
    str_type functions_str;
    int n;
    init(&n, &x_k, &functions_str);
    for (int i = 0; i < CNT; ++i) {
        x_k = iteration(n, x_k, functions_str);
    }
    puts("The answer is:\n X = ");
    output_matrix(x_k, n, 1);
}

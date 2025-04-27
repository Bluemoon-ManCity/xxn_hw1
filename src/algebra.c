#include "algebra.h"
#include <stdio.h>
#include <math.h>

Matrix create_matrix(int row, int col)
{
    Matrix m;
    m.rows = row;
    m.cols = col;
    return m;
}

Matrix add_matrix(Matrix a, Matrix b)
{
    int i, j;
    if (a.rows != b.rows || a.cols != b.cols)
    {
        printf("Error: Matrix a and b must have the same rows and cols.\n");
        return create_matrix(0, 0);
    }
    Matrix c = create_matrix(a.rows, a.cols);
    for (i = 0; i < a.rows; i++)
    {
        for (j = 0; j < a.cols; j++)
        {
            c.data[i][j] = a.data[i][j] + b.data[i][j];
        }
    }
    return c;
}

Matrix sub_matrix(Matrix a, Matrix b)
{
    int i, j;
    if (a.rows != b.rows || a.cols != b.cols)
    {
        printf("Error: Matrix a and b must have the same rows and cols.\n");
        return create_matrix(0, 0);
    }
    Matrix c = create_matrix(a.rows, a.cols);
    for (i = 0; i < a.rows; i++)
    {
        for (j = 0; j < a.cols; j++)
        {
            c.data[i][j] = a.data[i][j] - b.data[i][j];
        }
    }
    return c;
}

Matrix mul_matrix(Matrix a, Matrix b)
{
    int i, j;
    if (a.cols != b.rows)
    {
        printf("Error: The number of cols of matrix a must be equal to the number of rows of matrix b.\n");
        return create_matrix(0, 0);
    }
    Matrix c = create_matrix(a.rows, b.cols);
    for (i = 0; i < a.rows; i++)
    {
        for (j = 0; j < b.cols; j++)
        {
            c.data[i][j] = 0;
            for (int k = 0; k < a.cols; k++)
            {
                c.data[i][j] += a.data[i][k] * b.data[k][j];
            }
        }
    }
    return c;
}

Matrix scale_matrix(Matrix a, double k)
{
    int i, j;
    Matrix c = create_matrix(a.rows, a.cols);
    for (i = 0; i < a.rows; i++)
    {
        for (j = 0; j < a.cols; j++)
        {
            c.data[i][j] = a.data[i][j] * k;
        }
    }
    return c;
}

Matrix transpose_matrix(Matrix a)
{
    int i, j;
    Matrix c = create_matrix(a.cols, a.rows);
    for (i = 0; i < a.rows; i++)
    {
        for (j = 0; j < a.cols; j++)
        {
            c.data[j][i] = a.data[i][j];
        }
    }
    return c;
}

double det_matrix(Matrix a)
{
    int i, j, k;
    if (a.rows != a.cols)
    {
        printf("Error: The matrix must be a square matrix.\n");
        return 0;
    }
    if (a.rows == 1)
    {
        return a.data[0][0];
    }
    if (a.rows == 2)
    {
        return a.data[0][0] * a.data[1][1] - a.data[0][1] * a.data[1][0];
    }
    double det = 0;
    for (i = 0; i < a.cols; i++)
    {
        Matrix sub_matrix = create_matrix(a.rows - 1, a.cols - 1);
        for (j = 1; j < a.rows; j++)
        {
            for (k = 0; k < a.cols; k++)
            {
                if (k < i)
                {
                    sub_matrix.data[j - 1][k] = a.data[j][k];
                }
                else if (k > i)
                {
                    sub_matrix.data[j - 1][k - 1] = a.data[j][k];
                }
            }
        }
        det += pow(-1, i) * a.data[0][i] * det_matrix(sub_matrix);
    }
    return det; 
}

Matrix inv_matrix(Matrix a)
{
    int i, j, k, l;
    if (a.rows != a.cols)
    {
        printf("Error: The matrix must be a square matrix.\n");
        return create_matrix(0, 0);
    }
    double det = det_matrix(a);
    if (det == 0)
    {
        printf("Error: The matrix is singular.\n");
        return create_matrix(0, 0);
    }
    Matrix c = create_matrix(a.rows, a.cols);
    for (i = 0; i < a.rows; i++)
    {
        for (j = 0; j < a.cols; j++)
        {
            Matrix sub_matrix = create_matrix(a.rows - 1, a.cols - 1);
            for (k = 0; k < a.rows; k++)
            {
                for (l = 0; l < a.cols; l++)
                {
                    if (k != i && l != j)
                    {
                        sub_matrix.data[k - (k > i)][l - (l > j)] = a.data[k][l];
                    }
                }
            }
            c.data[j][i] = pow(-1, i + j) * det_matrix(sub_matrix) / det;
        }
    }
    return c;
}

int rank_matrix(Matrix a)
{
    int i, j, k;
    double ratio;
    if (a.rows != a.cols)
    {
        printf("Error: The matrix must be a square matrix.\n");
        return 0;
    }
    int rank = a.rows;
    for (i = 0; i < rank; i++)
    {
        if (a.data[i][i] != 0)
        {
            for (j = i + 1; j < rank; j++)
            {
                ratio = a.data[j][i] / a.data[i][i];
                for (k = i; k < rank; k++)
                {
                    a.data[j][k] -= ratio * a.data[i][k];
                }
            }
        }
        else
        {
            rank--;
            for (j = i; j < rank; j++)
            {
                for (k = 0; k < a.cols; k++)
                {
                    a.data[j][k] = a.data[j + 1][k];
                }
            }
            i--;
        }
    }
    return rank;
}

double trace_matrix(Matrix a)
{
    int i;
    if (a.rows != a.cols)
    {
        printf("Error: The matrix must be a square matrix.\n");
        return 0;
    }
    double trace = 0;
    for (i = 0; i < a.rows; i++)
    {
        trace += a.data[i][i];
    }
    return trace;
}

void print_matrix(Matrix a)
{
    int i, j;
    for (i = 0; i < a.rows; i++)
    {
        for (j = 0; j < a.cols; j++)
        {
            // 按行打印，每个元素占8个字符的宽度，小数点后保留2位，左对齐
            printf("%-8.2f", a.data[i][j]);
        }
        printf("\n");
    }
}
#ifndef KABSCH_RMSD
#define KABSCH_RMSD

#include <valarray>

#include <cblas.h>
#include <lapacke.h>

namespace kabsch
{


template <class M>
double rmsd(
        const M P,
        const M Q)
{
    double rmsd {0.0};
    const unsigned int D {3};
    const unsigned int size = P.size();
    const unsigned int N = size / D;

    for(unsigned int i = 0; i < size; ++i) {
        rmsd += (P[i] - Q[i])*(P[i] - Q[i]);
    }

    return sqrt(rmsd/N);
}


template <class M>
M centroid(M coordinates)
{

    double x {0};
    double y {0};
    double z {0};
    unsigned int size = coordinates.size();
    unsigned int n_atoms = size / 3;

    unsigned int i = 0;
    while(i<size)
    {
        x += coordinates[i++];
        y += coordinates[i++];
        z += coordinates[i++];
    }

    x /= n_atoms;
    y /= n_atoms;
    z /= n_atoms;

    return M {x, y, z};
}


template <class Matrix>
Matrix multiply(Matrix A, Matrix B,
    const int M,
    const int N,
    const int K)
{
    double one = 1.0;

    Matrix C(M*N);

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
        M, N, K, one,
        begin(A), K,
        begin(B), N, 1.0,
        begin(C), N);

    return C;
}


template <class Matrix>
Matrix transpose_multiply(Matrix A, Matrix B,
    const int M,
    const int N,
    const int K)
{
    double one = 1.0;

    Matrix C(M*N);

    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans,
        M, N, K, one,
        begin(A), M,
        begin(B), N, one,
        begin(C), N);

    return C;
}


template <class Matrix>
std::tuple<Matrix, Matrix, Matrix> matrix_svd(Matrix A, int rows, int cols)
{
    // lapack_int LAPACKE_dgesvd( int matrix_layout, char jobu, char jobvt,
    //     lapack_int m, lapack_int n,
    //     double* a, lapack_int lda,
    //     double* s, double* u, lapack_int ldu,
    //     double* vt, lapack_int ldvt,
    //     double* superb );

    Matrix U(cols*rows);
    Matrix S(rows);
    Matrix VT(cols*rows);
    Matrix superb(cols*rows);
    int info;

    info = LAPACKE_dgesvd(LAPACK_ROW_MAJOR, 'A', 'A',
        rows, cols,
        begin(A), rows,
        begin(S),
        begin(U), rows,
        begin(VT), rows,
        begin(superb));

    // TODO check if SVD is unconverged
    if(info > 0)
    {
        // Failed
    }

    return make_tuple(U, S, VT);
}


template <class M>
double determinant3x3(M A)
{
    // determinant of a square 3x3 matrix
    double det = A[0]*A[4]*A[8]
        +A[1]*A[5]*A[6]
        +A[2]*A[3]*A[7]
        -A[2]*A[4]*A[6]
        -A[1]*A[3]*A[8]
        -A[0]*A[5]*A[7];

    return det;
}


template <class M, class T>
M kabsch(M P, M Q, const T n_atoms)
{
    // const unsigned int L = P.size();
    // const unsigned int D = 3;

    M U;
    M V;
    M S;

    M C = transpose_multiply(P, Q, 3, 3, n_atoms);

    tie(U, S, V) = matrix_svd(C, 3, 3);

    // Getting the sign of the det(U)*(V) to decide whether we need to correct
    // our rotation matrix to ensure a right-handed coordinate system.
    if(determinant3x3(U)*determinant3x3(V) < 0.0)
    {
        // TODO More numpy'ish way to do this?
        U[std::slice( 2, 3, 3 )] = M({-U[3*0+2], -U[3*1+2], -U[3*2+2]});
    }

    M rotation = multiply(U, V, 3, 3, 3);

    return rotation;
}


template <class M, class T>
M kabsch_rotate(
        const M P,
        const M Q,
        const T n_atoms)
{
    M U = kabsch(P, Q, n_atoms);
    M product = multiply(P, U, n_atoms, 3, 3);
    return product;
}


template <class M, class T>
double kabsch_rmsd(
        const M P,
        const M Q,
        const T n_atoms)
{
    M P_rotated = kabsch_rotate(P, Q, n_atoms);
    return rmsd(P_rotated, Q);
}


} // namespace rmsd

#endif

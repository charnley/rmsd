
#include <iostream>
#include <sstream>
#include <fstream>
#include <iterator>
#include <vector>
#include <string.h>

#include <tuple>

#include <cblas.h>
#include <lapacke.h>

#include <valarray>

#include "matrix.h"

typedef std::valarray<double> matrix;

template <class M>
double matrix_det3x3(M A)
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

matrix multiply(matrix A, matrix B,
    const int M,
    const int N,
    const int K)
{
    double one = 1.0;
    char no = 'N';

    matrix C(M*N);

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
        M, N, K, one,
        begin(A), K,
        begin(B), N, 1.0,
        begin(C), N);

    return C;
}

// Should be
// template <class M>
// M transpose_multiply(M A, M B)

matrix transpose_multiply(matrix A, matrix B,
    const int M,
    const int N,
    const int K)
{
    double one = 1.0;
    char no = 'N';
    char yes = 'T';

    std::valarray<double> C(M*N);

    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans,
        M, N, K, one,
        begin(A), M,
        begin(B), N, one,
        begin(C), N);

    return C;
}


std::tuple<matrix, std::valarray<double>, matrix> matrix_svd(matrix A, int rows, int cols)
{
    // lapack_int LAPACKE_dgesvd( int matrix_layout, char jobu, char jobvt,
    //     lapack_int m, lapack_int n,
    //     double* a, lapack_int lda,
    //     double* s, double* u, lapack_int ldu,
    //     double* vt, lapack_int ldvt,
    //     double* superb );

    matrix U(cols*rows);
    std::valarray<double> S(rows);
    matrix VT(cols*rows);
    matrix superb(cols*rows);
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


template<typename Out>
void split(const std::string &s, char delim, Out result)
{
    std::stringstream ss;
    ss.str(s);
    std::string item;
    while(getline(ss, item, delim))
    {
        // Ignore multiple delimeters in a row
        if (item.empty())
            continue;

        *(result++) = item;
    }
}


std::vector<std::string> split(const std::string &s, char delim)
{
    std::vector<std::string> elems;
    split(s, delim, back_inserter(elems));
    return elems;
}


void print_coordinates(
    const matrix coord,
    const std::valarray<std::string> atoms,
    const unsigned int n_atoms)
{
    for(unsigned int i=0; i<n_atoms; i++)
    {
        std::cout << atoms[i] << " ";
        for(unsigned int j=0; j<3; j++)
        {
            std::cout << coord[3*i + j] << " ";
        }
        std::cout << std::endl;
    }
}


void print_matrix(const matrix D, unsigned int rows, unsigned int cols)
{
    for(unsigned int i=0; i<rows; i++)
    {
        for(unsigned int j=0; j<cols; j++)
        {
            std::cout << D[cols*i + j] << " ";
        }
        std::cout << std::endl;
    }
}


std::tuple<std::valarray<std::string>, matrix, int> get_coordinates(const std::string filename)
{
    std::string line;
    std::vector<std::string> coordinate_line;
    std::ifstream xyzfile (filename);

    unsigned int n_atoms {0};
    std::string atom;
    std::vector<double> coordinate {};
    double x, y, z;

    if(xyzfile.is_open())
    {
        // Get number of atoms
        getline(xyzfile, line);
        n_atoms = stoi(line);

        // Allocate atoms and coordinates,
        // now that we know how many atoms are included
        std::valarray<std::string> atoms(n_atoms);
        matrix coordinates(3 * n_atoms);

        // Empty line
        getline(xyzfile, line);

        // Read coordinates
        int i = 0;
        while(getline(xyzfile, line))
        {
            coordinate_line = split(line, ' ');

            atom = coordinate_line[0];
            x = std::stod(coordinate_line[1]);
            y = std::stod(coordinate_line[2]);
            z = std::stod(coordinate_line[3]);

            // n_cols * row + col
            coordinates[3*i + 0] = x;
            coordinates[3*i + 1] = y;
            coordinates[3*i + 2] = z;

            atoms[i] = atom;
            i++;
        }

        xyzfile.close();

        return make_tuple(atoms, coordinates, n_atoms);
    }
    else
    {
        std::cout << "Unable to open " << filename << std::endl;
        exit(0);
    }
}


double rmsd(const matrix P, const matrix Q)
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


matrix kabsch(matrix P, matrix Q, const unsigned int n_atoms)
{
    const unsigned int L = P.size();
    const unsigned int D = 3;

    matrix U;
    matrix V;
    std::valarray<double> S;

    matrix C = transpose_multiply(P, Q, 3, 3, n_atoms);

    tie(U, S, V) = matrix_svd(C, 3, 3);

    // Getting the sign of the det(U)*(V) to decide whether we need to correct
    // our rotation matrix to ensure a right-handed coordinate system.
    if(matrix_det3x3(U)*matrix_det3x3(V) < 0.0)
    {
        // TODO More numpy'ish way to do this?
        U[std::slice( 2, 3, 3 )] = {-U[3*0+2], -U[3*1+2], -U[3*2+2]};
    }

    matrix rotation = multiply(U, V, 3, 3, 3);

    return rotation;
}


matrix kabsch_rotate(const matrix P, const matrix Q, const unsigned int n_atoms)
{
    matrix U = kabsch(P, Q, n_atoms);
    matrix product = multiply(P, U, n_atoms, 3, 3);
    return product;
}


double kabsch_rmsd(const matrix P, const matrix Q, const unsigned int n_atoms)
{
    matrix P_rotated = kabsch_rotate(P, Q, n_atoms);
    return rmsd(P_rotated, Q);
}

std::valarray<double> centroid(const matrix coordinates)
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

    return std::valarray<double> {x, y, z};
}


int main(int argc, char* argv[])
{

    if(argc == 1)
    {
        printf("Usage: rmsd structure_a.xyz structure_b.xyz\n");
        exit(0);
    }

    if(argc < 3)
    {
        printf ("Not enough arguments\n");
        exit(EXIT_FAILURE);
    }

    bool output = false;
    int a = 0;
    int b = 0;

    if(argc > 3)
    {
        for(int i=1; i<argc; i++)
        {
            if(strcmp(argv[i], "-o") == 0)
            {
                output = true;
            }
            else if (a == 0)
            {
                a = i;
            }
            else
            {
                b = i;
            }

        }
    }
    else
    {
        a = 1;
        b = 2;
    }

    const std::string p_file = argv[a];
    const std::string q_file = argv[b];

    unsigned int n_atoms {0};
    unsigned int p_no_atoms;
    unsigned int q_no_atoms;

    std::valarray<std::string> p_atoms;
    std::valarray<std::string> q_atoms;

    matrix q_coord;
    matrix p_coord;

    std::valarray<double> p_center;
    std::valarray<double> q_center;

    tie(p_atoms, p_coord, p_no_atoms) = get_coordinates(p_file);
    tie(q_atoms, q_coord, q_no_atoms) = get_coordinates(q_file);

    n_atoms = p_no_atoms;

    // Calculate the center of the molecules
    p_center = centroid(p_coord);
    q_center = centroid(q_coord);
    print_coordinates(p_coord, p_atoms, n_atoms);
    printf("\n");

    // Recenter molecules in origin
    for(unsigned int i = 0; i < n_atoms; i++) {
        for(unsigned int d = 0; d < 3; d++) {
            p_coord[3*i + d] -= p_center[d];
            q_coord[3*i + d] -= q_center[d];
        }
    }

    if (!output)
    {
        double krmsd = kabsch_rmsd(p_coord, q_coord, n_atoms);
        std::cout << krmsd << std::endl;
    }
    else
    {

        printf("p_coord\n");
        print_matrix(p_coord, n_atoms, 3);
        printf("\n");

        matrix U = kabsch(p_coord, q_coord, n_atoms);
        matrix product = multiply(p_coord, U, n_atoms, 3, 3);

        printf("\n");

        print_matrix(U, 3, 3);

        printf("\n");

        for(unsigned int i = 0; i < n_atoms; i++) {
            for(unsigned int d = 0; d < 3; d++) {
                product[3*i + d] += q_center[d];
            }
        }

        print_coordinates(product, p_atoms, n_atoms);
    }
}


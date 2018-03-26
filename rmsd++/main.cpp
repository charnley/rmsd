#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h>
#include <sys/stat.h>
#include <tuple>
#include <vector>

#include <valarray>

#include "kabsch.h"

typedef std::valarray<double> matrix;

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


inline bool exists (const std::string& name)
{
    struct stat buffer;
    return (stat (name.c_str(), &buffer) == 0);
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
        std::cerr << "Unable to open " << filename << std::endl;
        exit(0);
    }
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
            printf("%10.6lf ", coord[3*i+j]);
        }
        std::cout << std::endl;
    }
}


int main(int argc, char* argv[])
{

    if(argc == 1)
    {
        printf("Usage: calculate_rmsd [options] alpha.xyz beta.xyz\n");
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

    for(int i=1; i < argc; i++)
    {

        if(strcmp(argv[i], "-o") == 0 || strcmp(argv[i], "--output") == 0)
        {
            output = true;
            continue;
        }

        if(a==0)
        {
            a = i;
        }
        else
        {
            b = i;
        }

    }

    const std::string p_file = argv[a];
    const std::string q_file = argv[b];

    unsigned int n_atoms {0};

    std::valarray<std::string> p_atoms;
    std::valarray<std::string> q_atoms;

    matrix p_coord;
    matrix q_coord;

    std::valarray<double> p_center;
    std::valarray<double> q_center;

    std::tie(p_atoms, p_coord, n_atoms) = get_coordinates(p_file);
    std::tie(q_atoms, q_coord, n_atoms) = get_coordinates(q_file);

    // find the center of the structures
    p_center = kabsch::centroid(p_coord);
    q_center = kabsch::centroid(q_coord);

    // Recenter molecules in origin
    for(unsigned int i = 0; i < n_atoms; i++) {
        for(unsigned int d = 0; d < 3; d++) {
            p_coord[3*i + d] -= p_center[d];
            q_coord[3*i + d] -= q_center[d];
        }
    }

    if(output)
    {
        matrix p_rotated = kabsch::kabsch_rotate(p_coord, q_coord, n_atoms);

        // Center coordinates P on Q
        for(unsigned int i = 0; i < n_atoms; i++) {
            for(unsigned int d = 0; d < 3; d++) {
                p_rotated[3*i + d] += q_center[d];
            }
        }

        // print the structure in XYZ format
        std::cout << n_atoms << std::endl << std::endl;
        print_coordinates(p_rotated, p_atoms, n_atoms);

        return 0;
    }

    // rotate and calculate rmsd
    double krmsd = kabsch::kabsch_rmsd(p_coord, q_coord, n_atoms);
    printf( "%.10lf \n", krmsd);

} // main


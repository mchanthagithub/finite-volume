//
// Created by maytee on 4/26/18.
//
// Structured Carteresian Grid
//     _______F11______________F12_______
//    |                 |                |
//  F8|        2      F9|       3     F10|
//    |                 |                |
//    |_________________|________________|
//    |        F6       |       F7       |
//  F3|               F4|              F5|
//    |        0        |       1        |
//    |_________________|________________|
//           F1                F2
//
//    3______N=2________2
//    |                 |
// W=3|                 |E=1
//    |_________________|
//    0       S=0        1
// -Local node numbering is assumed to go from to 0 to 3 counterclockwise, starting from bottom left
// -Cells are numbered from right to left, bottom up, so bottom-left most node is 0, bottom-right most node is Nx-1
//      and top-rightmost is numCells-1
// -Global nodes are numbered similarly to cells
// -Local face numbering goes 0 = S, 1 = E, 2 = N and 3 = W
// -Global face numbering starts from bottom on horizontal faces, then goes vertical faces, then goes next horizontal, etc.
//


#ifndef PROJECT_STRUCTUREDGRID_H
#define PROJECT_STRUCTUREDGRID_H
#include "Grid.h"

class CartesianGrid : public Grid
{
public:
    // Takes number of elements in X and Y and length of X and Y sides
    CartesianGrid();
    CartesianGrid(int Nx, int Ny, double Lx, double Ly, double smallX, double smallY);
    void setVariableSizes(int nCells, int nFaces, int nCorners);
    void generateMesh();

    void setBCs() override;
    void setInitialValues() override;

    int Nx;
    int Ny;
    double Lx;
    double Ly;
    double delX;
    double delY;
    double minX;
    double minY;
};


#endif //PROJECT_STRUCTUREDGRID_H

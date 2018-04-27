//
// Created by maytee on 4/26/18.
//
// Structured Carteresian Grid

#ifndef PROJECT_STRUCTUREDGRID_H
#define PROJECT_STRUCTUREDGRID_H
#include "Grid.h"

class CartesianGrid : public Grid
{
public:
    // Takes number of elements in X and Y and length of X and Y sides
    CartesianGrid(int Nx, int Ny, double Lx, double Ly);
    void setVariableSizes(int nCells, int nFaces, int nCorners);
    void generateMesh();

    void setBCs();
    void setInitialValues();

    int Nx;
    int Ny;
    double Lx;
    double Ly;
    double delX;
    double delY;
};


#endif //PROJECT_STRUCTUREDGRID_H

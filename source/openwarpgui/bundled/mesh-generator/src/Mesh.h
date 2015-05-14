/*

    Copyright (C) 2014 TopCoder Inc., All Rights Reserved.

    @author TCSASSEMBLER, caoweiquan322
    @version 1.1
  */
#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <cmath>
#include <list>
#include <map>
#include <algorithm>
#include <functional>
#include <cctype>
#include <locale>
#include <set>

#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkPolyData.h>
#include <vtkPolyDataWriter.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkMeshQuality.h>
#include <vtkDoubleArray.h>
#include <vtkCellData.h>

#define WNT
    #include "TopTools_IndexedMapOfShape.hxx"
namespace nglib{
    #define OCCGEOMETRY
    #include "nglib.h"
}
using namespace std;
using namespace nglib;
#include "Vector.h"
#include "MeshStats.h"



#define PI 3.14159265359


class Mesh{
    /*
    This class heavily uses Ng and vtk library.
    A program must call Ng_Init() once before calling any Ng function.
    It must call Ng_Exit() once after all Ng function call and before program terminates.
    */
public:
    // name of mesh
    string name;

    // should contain unique point in space
    // the only exception is after loadSTLFile
    vector<Vector> point;

    // array of array of index_of_point
    // each array of index_of_point specify one surface
    // index_of_point identify index in "vector<Vector> point"
    vector<vector<int>> surface;

    // feed Ng_Mesh to this object
    // param mesh: a valid Ng_Mesh
    void feedNgMesh(Ng_Mesh* mesh);

    // feed Ng_STL_Geometry to this object
    // param geo: a valid Ng_STL_Geometry object
    void feedNgSTLGeo(Ng_STL_Geometry* geo);

    // save this object as STL file
    // param path: path to file
    void saveAsSTL(const string& path);//STL(STereoLithography) ascii format

    //save this object as GDF file,called after forceTriToQuad()
    //param path: path to file
    void saveAsGDF(const string& path);

    // save this object as VTK file
    // param path: path to file
    void saveAsVTK(const string& path);//vtk ascii format

    // save this object as VTP file
    // param path: path to file
    void saveAsVTP(const string& path);//vtk xml format

    // save this object as DAT file
    // param path: path to file
    void saveAsDAT(const string& path);//Nemoh ascii format

    // load ascii STL file
    // param path: path to file
    // postcondition:
    // the resulting Mesh does not contain unique points.
    // one point in space may be duplicated in "vector<Vector> point"
    void loadSTLFile(const string& path);// load ascii STL format

    // print vector textually to out
    // param out: ostream object as output
    // param v: an object Vector
    void print(ostream& out, const Vector& v);

    // calculate total surface area
    // return: total surface area
    double getSurfaceArea();

    // calculate volume based on ref as point of reference
    // param ref: point as reference to calculate volume.
    // because of the nature of the algorithm,
    //      for a closed surface the return value will be independent
    //      of ref. But the value will depend on ref if
    //      surface is not closed
    double getVolume(const Vector& ref = Vector());

    // get edge count
    // return: number of edge
    // precondition: point in mesh should be unique
    int getEdgeCount();


    // check if this mesh is closed
    // return : return true if this mesh is closed, false otherwise
    // a surface is closed (contains volume) if volume calculated by getVolume(ref) does not depends on
    // point of reference. Otherwise it is an open surface (shell elements).
    bool isClosed(){
        return abs(getVolume() - getVolume(Vector(10,10,10)))<0.01;
    }

    // convert the mesh to VtkPolyData object
    // return: the vtkPolyData object
    // user must destroy the return value manually, using value->Delete();
    vtkPolyData* toVtkPolyData();

    // convert the mesh to Ng_Mesh object.
    // return the Ng_Mesh object.
    // user must destroy the return value manually, using nglib::Ng_DeleteMesh(value);
    Ng_Mesh* toNgMesh();

    // convert the mesh to Ng_STL_Geometry object.
    // return: the Ng_STL_Geometry object.
    // no need to delete, since there is no Ng_DeleteGeometry
    Ng_STL_Geometry* toSTLGeometry();

    // force the remaining tri into quad
    // split tri into 3 quad, propagate new points to split the remaining quad/tri
    // this will half the dimension of each tri/quad.
    // precondition:
    //      -. surface must only contain tri or quad.
    //      -. point with the same coordinate already merged
    // post condition: all surfaces are quad
    // return: is any tri converted?
    bool forceTriToQuad();

    //calculate statistics
    // param: object that will contain statistics
    void stats(MeshStats*);

private:
    // Serialize one surface to Nemoh dat file.
    // Note that we assume it MUST be a QUADRANGLE surface.
    // param out: ostream object as output
    // param index: array of index of vertices specifying a polygon
    //      index of vertices as in vector<Vector> point
    // param points: the points
    void serializeDatSurface(ostream& out, const vector<int>& index,
        const vector<Vector>& points);
    // surface area vector of a polygon
    // to get value of area use .length()
    // param out: surface area Vector
    // param index: array of index of vertices specifying a polygon
    //      index of vertices as in vector<Vector> point
    void surfaceAreaVector(Vector* out, const vector<int>& index);

    // center of a polygon
    // param out: Vector of center of the polygon
    // param index: array of index of vertices specifying a polygon
    //      index of vertices as in vector<Vector> point
    void center(Vector* out, const vector<int>& index);

    // equiangular value based on
    // http://en.wikipedia.org/wiki/Types_of_mesh#Equiangular_skew
    // param index: array of index of vertices specifying a polygon
    //      index of vertices as in vector<Vector> point
    // return: equiangular value of the polygon
    double equiangularSkew(const vector<int>& index);

    // aspect ratio based on
    // http://en.wikipedia.org/wiki/Types_of_mesh#Aspect_ratio
    // param index: array of index of vertices specifying a polygon
    //      index of vertices as in vector<Vector> point
    // return: aspect ratio of input polygon
    double aspectRatio(const vector<int>& index);

    // reverse map,
    // collect for each point which surface do they belong
    // param out: reverse mapping, that is from point to surface
    //      given index of a point, this return value can give
    //      list of surface that this point belong to.
    void collectPointToSurfaceMapping(vector<vector<int>>* out);

    // param out: array of array of neighbour id
    //      each correspond to an edge
    //      if that edge does not have neighbor, then its value is -1
    // param pts: point to surface mapping
    void collectNeighborMapping(vector<vector<int>>* out, const vector<vector<int>>& pts);

    // divide one polygon, and propagate the new point created to the remaining
    // param startPolyId: valid id of either tri or quad
    // param neigh: output from collectNeighborMapping().
    void devideToQuadAndPropagate(int startPolyId, const vector<vector<int>>& nigh);
};

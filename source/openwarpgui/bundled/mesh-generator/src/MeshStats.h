/*

    Copyright (C) 2014 TopCoder Inc., All Rights Reserved.

    @author TCSASSEMBLER
    @version 1.0
  */
#pragma once

#include <vector>
#include <string>
#include <sstream>
#include <algorithm>

using namespace std;

/*
Utility class to simplify stats calculation
*/
class MeshStats{
public:
    // minimal edge length
    double minEdge;

    // maximal edge length
    double maxEdge;

    // minimal skewness of polygon
    double minSkew;

    // maximal skewness of polygon
    double maxSkew;

    // minimal aspect ratio of polygon
    double minRatio;

    // maximal aspect ratio of polygon
    double maxRatio;

    // constructor
    MeshStats(){
        minEdge = minSkew = minRatio = 1e15;
        maxEdge = maxSkew = maxRatio = -1e15;
    }
    // feed this object with length of an edge
    // param in: edge length
    void feedEdge(double in){
        minEdge = min(minEdge,in);
        maxEdge = max(maxEdge,in);
    }
    // feed this object with skewness of a polygon
    // param in: skewness
    void feedSkew(double in){
        minSkew = min(minSkew,in);
        maxSkew = max(maxSkew,in);
    }
    // feed this object with aspect ratio of a polygon
    // param in: aspect ratio
    void feedRatio(double in){
        minRatio = min(minRatio,in);
        maxRatio = max(maxRatio,in);
    }
    // convert to readable string
    // return: readable string representing this object
    string tostring(){
        ostringstream out;
        out << "Edge(resolution) (min -> max): " << minEdge << " -> " << maxEdge << endl;
        out << "Skew (max -> min): " << maxSkew << " -> " << minSkew << endl;
        out << "AspectRatio (max -> min): " << maxRatio << " -> " << minRatio << endl;
        return out.str();
    }
};


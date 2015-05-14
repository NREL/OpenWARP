/*

    Copyright (C) 2014 TopCoder Inc., All Rights Reserved.

    @author TCSASSEMBLER, caoweiquan322
    @version 1.1
  */
#pragma once

#include <vector>
#include <string>
#include <sstream>
#include <cmath>

using namespace std;

/*
3 dimensional vector
each method is self describing
*/
class Vector{
public:
    double x, y, z;

    // constructor
    // param x_: the x component
    // param y_: the y component
    // param z_: the z component
    Vector(double x_=0, double y_=0, double z_=0):x(x_),y(y_),z(z_){
    }
    // commonly understood: this -= b
    // param b: input Vector
    // return: this
    Vector& operator-=(const Vector& b){
        x-=b.x;
        y-=b.y;
        z-=b.z;
        return *this;
    }
    // commonly understood: mid point of between this and b
    // param b: input Vector
    // return: mid vector between this and b
    Vector mid(const Vector& b){
        return Vector((x+b.x)/2,(y+b.y)/2,(z+b.z)/2);
    }
    // return: length of this vector
    double length() const{
        return sqrt(x*x+y*y+z*z);
    }
    // commonly understood: this += b
    // param b: input Vector
    // return: this
    Vector& operator+=(const Vector& in){
        x += in.x;
        y += in.y;
        z += in.z;
        return *this;
    }
    // commonly understood: this *= b
    // param b: input Vector
    // return: this
    Vector& operator*=(double in){
        x*=in;
        y*=in;
        z*=in;
        return *this;
    }
    // return: string representation of this vector
    string tostring(){
        ostringstream out;
        out << "<" << x << "," << y << "," << z << ">";
        return out.str();
    }

    // Calculate cross production between this vector and another
    // param in: the vector to do cross production with
    // return: cross product of this vector with in
    Vector cross(const Vector& in){
        return Vector(y*in.z-z*in.y,z*in.x-x*in.z,x*in.y-y*in.x);
    }

    // Calculate dot production between this vector and another
    // param in: the vector to do dot production with
    // return: dot product of this vector with in
    double dot(const Vector& in){
        return x*in.x+y*in.y+z*in.z;
    }
    // get consine of angle wrt in
    // param in: vector to which angle will be measured
    double cosine(const Vector& in){
        return dot(in)/length()/in.length();
    }
};

// Calculate summation of two vectors
// param lhs: the input vector
// param rhs: the vector to be added
// return: lhs + rhs
inline Vector operator+(Vector lhs, const Vector& rhs){
    lhs += rhs;
    return lhs;
}

// Calculate subtraction of two vectors
// param lhs: the input vector
// param rhs: the vector to be subtracted
// return: lhs - rhs
inline Vector operator-(Vector lhs, const Vector& rhs){
    lhs -= rhs;
    return lhs;
}

// Calculate product of one vector and one double value
// param lhs: the vector
// param rhs: the double value
// return: lhs * rhs
inline Vector operator*(Vector lhs, double rhs){
    lhs *= rhs;
    return lhs;
}


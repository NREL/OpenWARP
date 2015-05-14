/*

    Copyright (C) 2014 TopCoder Inc., All Rights Reserved.

    @author TCSASSEMBLER, caoweiquan322
    @version 1.1
  */

#include "stdafx.h"
#include <exception>
#include "Mesh.h"

void Mesh::serializeDatSurface(ostream& out, const vector<int>& index,
                               const vector<Vector>& points) {
    Vector a = points[index[2]] - points[index[0]];
    Vector b = points[index[3]] - points[index[1]];
    Vector sum;
    for (int x : index)
    {
        sum += points[x];
    }
    // Nemoh .dat file requires the normal vector NOT point to x axis.
    if ((a.cross(b).y) * (sum.y) < 0.0)
    {
        out << (index[3] + 1) << " " << (index[2] + 1) << " "
            << (index[1] + 1) << " " << (index[0] + 1) << endl;
    }
    else
    {
        out << (index[0] + 1) << " " << (index[1] + 1) << " "
            << (index[2] + 1) << " " << (index[3] + 1) << endl;
    }
}

void Mesh::surfaceAreaVector(Vector *out, const vector<int>& index){
    Vector sum;
    for(int x = 0; x < index.size(); x++){
        sum += point[index[x]].cross(point[index[(x+1)%index.size()]]);
    }
    (*out) = sum * 0.5;
}
void Mesh::center(Vector* out, const vector<int>& index){
    Vector sum;
    for(int x:index){
        sum += point[x];
    }
    (*out) = sum * (1.0/index.size());
}
double Mesh::getSurfaceArea(){
    double all = 0;
    Vector out;
    for(auto index:surface){
        surfaceAreaVector(&out, index);
        all += out.length();
    }
    return all;
}
double Mesh::getVolume(const Vector& ref){
    double all = 0;
    Vector area, out;
    for(auto index:surface){
        surfaceAreaVector(&area, index);
        center(&out, index);
        all += area.dot(out-ref);
    }
    return all*(1.0/3);
}
int Mesh::getEdgeCount(){
    // this algorithm applicable to both closed & open surface (shell elements)
    vector<set<int>> edge(point.size());
    for(auto vertexes:surface){
        for(int i=0;i<vertexes.size();i++){
            int a = vertexes[i];
            int b = vertexes[(i+1)%vertexes.size()];
            if(a>b){
                swap(a,b);
            }
            edge[a].insert(b);
        }
    }
    int out = 0;
    for(auto e:edge){
        out += (int) e.size();
    }
    return out;
}
void Mesh::print(ostream& out, const Vector& v){
    out << "vertex " << v.x << " " << v.y << " " << v.z << endl;
}

void Mesh::saveAsGDF(const string& path) {
    ofstream out(path);
    out << "WAMIT file export (mesh)" << endl;
    out << "1 9.80665   ULEN GRAV" << endl;
    out << "0  0    ISX  ISY" << endl;
    out << surface.size() << endl;
    for(vector<int> index:surface)
        if(index.size() == 4 || index.size() == 3)
        {
            for(int i = 0 ; i < index.size(); ++i)
                out << point[index[i]].x<<" "<< point[index[i]].y <<" " << point[index[i]].z << endl;
            if(index.size() == 3)
                out << point[index[2]].x<<" "<< point[index[2]].y <<" " << point[index[2]].z << endl;
        }
    out.flush();
    out.close();
}

void Mesh::saveAsSTL(const string& path){
    ofstream out(path);
    if(name.length()>0){
        out << "solid " << name << endl;
    }else{
        out << "solid " << "Dummy" << endl;
    }
    for(vector<int> index:surface){
        out << "facet normal 1.0 1.0 1.0" << endl;
        out << "outer loop" << endl;
        print(out,point[index[0]]);
        print(out,point[index[1]]);
        print(out,point[index[2]]);
        out << "endloop" << endl;
        out << "endfacet" << endl;
        if(index.size()==4){
            out << "facet normal 1.0 1.0 1.0" << endl;
            out << "outer loop" << endl;
            print(out,point[index[0]]);
            print(out,point[index[2]]);
            print(out,point[index[3]]);
            out << "endloop" << endl;
            out << "endfacet" << endl;
        }
    }
    out << "endsolid" << endl;
    out.flush();
    out.close();
}
vtkPolyData* Mesh::toVtkPolyData(){
    vtkPoints* points = vtkPoints::New();
    for(int i=0;i<point.size();i++){
        auto p = point[i];
        points->InsertPoint(i,p.x,p.y,p.z);
    }
    vtkPolyData* solid = vtkPolyData::New();
    solid->SetPoints(points);
    points->Delete();

    vtkCellArray* polys = vtkCellArray::New();
    vtkDoubleArray* skew = vtkDoubleArray::New();
    solid->GetCellData()->SetScalars(skew);//put skewness value
    for(auto index:surface){
        vector<vtkIdType> buf;
        for(auto x:index){
            buf.push_back(x);
        }
        vtkIdType cell = polys->InsertNextCell(buf.size(),&buf[0]);
        skew->InsertTuple1(cell,equiangularSkew(index));
    }
    solid->SetPolys(polys);
    polys->Delete();
    skew->Delete();
    return solid;
}
Ng_Mesh* Mesh::toNgMesh(){
    Ng_Mesh* out = nglib::Ng_NewMesh();
    for(auto v:point){
        double p[] = {v.x,v.y,v.z};
        nglib::Ng_AddPoint(out,p);
    }
    for(auto x:surface){
        if(x.size()==3){
            int index[] = {x[0]+1,x[1]+1,x[2]+1};//+1
            nglib::Ng_AddSurfaceElement(out,nglib::Ng_Surface_Element_Type::NG_TRIG,index);
        }else if(x.size()==4){
            int index[] = {x[0]+1,x[1]+1,x[2]+1,x[3]+1};//+1
            nglib::Ng_AddSurfaceElement(out,nglib::Ng_Surface_Element_Type::NG_QUAD,index);
        }
    }
    return out;
}
Ng_STL_Geometry* Mesh::toSTLGeometry(){
    // convert this mesh to Ng_STL_Geometry
    // as preparation to-mesh the surface
    Ng_STL_Geometry* geo = nglib::Ng_STL_NewGeometry();
    for(auto x:surface){
        double p1[] = {point[x[0]].x,point[x[0]].y,point[x[0]].z};
        double p2[] = {point[x[1]].x,point[x[1]].y,point[x[1]].z};
        double p3[] = {point[x[2]].x,point[x[2]].y,point[x[2]].z};
        nglib::Ng_STL_AddTriangle(geo,p1,p2,p3);
        nglib::Ng_STL_AddEdge(geo,p1,p2);
        nglib::Ng_STL_AddEdge(geo,p2,p3);
        nglib::Ng_STL_AddEdge(geo,p3,p1);
        if(x.size()==4){
            double p1[] = {point[x[2]].x,point[x[2]].y,point[x[2]].z};
            double p2[] = {point[x[3]].x,point[x[3]].y,point[x[3]].z};
            double p3[] = {point[x[0]].x,point[x[0]].y,point[x[0]].z};
            nglib::Ng_STL_AddTriangle(geo,p1,p2,p3);
            nglib::Ng_STL_AddEdge(geo,p1,p2);
            nglib::Ng_STL_AddEdge(geo,p2,p3);
            nglib::Ng_STL_AddEdge(geo,p3,p1);
        }
    }
    // initialize after adding triangles
    nglib::Ng_STL_InitSTLGeometry(geo);
    return geo;
}
void Mesh::saveAsVTK(const string& path){
    vtkPolyData* solid = toVtkPolyData();
    vtkPolyDataWriter* writer = vtkPolyDataWriter::New();
    writer->SetFileName(path.c_str());
    writer->SetInput(solid);
    writer->Write();
    solid->Delete();
    writer->Delete();
}
void Mesh::saveAsVTP(const string& path){
    vtkPolyData* solid = toVtkPolyData();
    vtkXMLPolyDataWriter* writer = vtkXMLPolyDataWriter::New();
    writer->SetDataModeToAscii();
    writer->SetFileName(path.c_str());
    writer->SetInput(solid);
    writer->Write();
    solid->Delete();
    writer->Delete();
}
void Mesh::saveAsDAT(const string& path) {
    // Check if all surfaces are quadrangle surfaces
    for (vector<int> index : surface) {
        // Nemoh is only able to process quadrangle surface
        if (index.size() != 4) {
            cout<< "All the surfaces should be quadrangle as input of Nemoh .dat files." << endl;
            cout<< "This process would be canceled before finished." << endl;
            return;
        }
    }

    ofstream out(path);
    out << "2 0" << endl;
    // Output points
    int itr = 0;
    for (auto v : point)
    {
        out << (++itr) << " " << v.x << " " << v.y << " " << v.z << endl;
    }
    out << "0 0.00 0.00 0.00" << endl;

    // Output surfaces
    for (vector<int> index : surface) {
        // Nemoh is only able to process quadrangle surface
        if (index.size() == 4) {
            serializeDatSurface(out, index, point);
        }
    }
    out << "0 0 0 0" << endl;
    out.flush();
    out.close();
}
void Mesh::feedNgMesh(Ng_Mesh* mesh){
    // feed a valid Ng_Mesh
    point.clear();
    surface.clear();
    int pcount = nglib::Ng_GetNP(mesh);
    double* x = new double[3];
    for(int i=0;i<pcount;i++){
        nglib::Ng_GetPoint(mesh,i+1,x);//+1
        point.push_back(Vector(x[0],x[1],x[2]));
    }
    delete[] x;
    int scount = nglib::Ng_GetNSE(mesh);
    int* ip = new int[4];
    for(int i=0;i<scount;i++){
        auto tipe = nglib::Ng_GetSurfaceElement(mesh,i+1,ip);//+1
        vector<int> index;
        bool indexOverflow = false;
        if(tipe == nglib::Ng_Surface_Element_Type::NG_TRIG ){
            for(int j=0;j<3;j++){
                index.push_back(ip[j]-1);//-1
                if(ip[j]>point.size()){
                    indexOverflow = true;
                }
            }
        }else if(tipe == nglib::Ng_Surface_Element_Type::NG_QUAD){
            for(int j=0;j<4;j++){
                index.push_back(ip[j]-1);//-1
                if(ip[j]>point.size()){
                    indexOverflow = true;
                }
            }
        }
        if(!indexOverflow){
            surface.push_back(index);
        }else{
            cout << "surface " << i << " index overflow" << endl;
            for(auto x:index){
                cout << "  " << x;
            }
            cout << endl;
        }
    }
    delete[] ip;
}
// utility used to loadSTLFile
string getToken(const string& text, int& index){

    while(index < text.length() && isspace(text[index])){
        index ++;
    }
    int start = index;
    while(index < text.length() && !isspace(text[index])){
        index ++;
    }
    return text.substr(start,index-start);
}
void Mesh::loadSTLFile(const string& path){
    // parse ASCII STL file directly
    // this is needed because the alternative:
    // 1. use nglib to load STL results in Ng_STL_Geometry
    //      and nglib does not provide support to read each data in Ng_STL_Geometry
    //      only after we converted to Ng_Mesh that we can read it.
    //      Ng_Mesh is an already meshed version.
    //      Ng_STL_Geometry need to be closed-surface in order to validly converted
    //      to Ng_Mesh
    // 2. use vtk lib.
    //      some version of vtk lib error in reading ASCII STL
    //      the version that checked to be error: vtk5.10.1, vtk6.0.0
    //      the version that checked to be valid: either 5.6 or 5.8 (i forget)
    //      but then those versions have some other problems for not able generate all library

    ifstream fin(path);
    string text;
    getline(fin,text,char(-1));
    fin.close();
    const int STACKSIZE = 10;
    string stack[STACKSIZE];
    int stackpos = 0;
    int index = 0;
    point.clear();
    surface.clear();
    vector<int> pointIndex;
    while(index < text.length() && stackpos<STACKSIZE-1){
        string token = getToken(text,index);
        if(token.length() == 0){
            break;
        }
        if(stackpos == 0 && token.compare("solid") == 0){
            token = getToken(text,index);
            name = token;
            stackpos += 1;
        }else if(stackpos == 1 && token.compare("endsolid") == 0){
            stackpos -= 1;
        }else if(stackpos == 1 && token.compare("facet") == 0){
            token = getToken(text,index); // normal
            if(token.compare("normal")){
                throw new runtime_error("required text 'normal' after facet");
            }
            double x = atof(getToken(text,index).c_str());
            double y = atof(getToken(text,index).c_str());
            double z = atof(getToken(text,index).c_str());
            Vector normal = Vector(x,y,z);
            // normal not saved
            stackpos += 1;
        }else if(stackpos == 2 && token.compare("endfacet") == 0){
            surface.push_back(pointIndex);
            pointIndex.clear();
            stackpos -= 1;
        }else if(stackpos == 2 && token.compare("outer") == 0){
            token = getToken(text,index); // loop
            if(token.compare("loop")){
                throw new runtime_error("required text 'loop' after outer");
            }
            stackpos += 1;
        }else if(stackpos == 3 && token.compare("endloop") == 0){
            stackpos -= 1;
        }else if(stackpos == 3 && token.compare("vertex") == 0){
            double x = atof(getToken(text,index).c_str());
            double y = atof(getToken(text,index).c_str());
            double z = atof(getToken(text,index).c_str());
            point.push_back(Vector(x,y,z));
            pointIndex.push_back((int) point.size()-1);
        }else if(stackpos < 0){
            throw new runtime_error("invalid file format: stackpos < 0");
        }else{
            throw new runtime_error("invalid file format: unknown");
        }
    }
    if(stackpos != 0){
        throw new runtime_error("invalid file format: stackpos != 0");
    }
}
void Mesh::stats(MeshStats* out){
    // statistics counter
    for(auto s:surface){
        for(int i=0;i<s.size();i++){
            out->feedEdge((point[s[i]]-point[s[(i+1)%s.size()]]).length());
        }
        out->feedSkew(equiangularSkew(s));
        out->feedRatio(aspectRatio(s));
    }
}
double Mesh::equiangularSkew(const vector<int>& index){
    // http://en.wikipedia.org/wiki/Types_of_mesh#Equiangular_skew
    // angular skew is max((theta_max-theta_e)/(pi-theta_e),(theta_e-theta_min)/theta_e)
    // theta_e = PI/3 for tri, = PI/2 for quad
    double thetaMin = 10000;
    double thetaMax = -10000;
    double thetaE = 90;
    if(index.size()==3){
        thetaE=PI/3;
    }else if(index.size() ==4){
        thetaE = PI/2;
    }
    for(int i=0;i<index.size();i++){
        Vector v1 = point[index[(i+1)%index.size()]]-point[index[i]];
        Vector v2 = point[index[(i+2)%index.size()]]-point[index[(i+1)%index.size()]];
        double cosValue = - v1.dot(v2)/(v1.length()*v2.length());
        double theta = acos(cosValue);
        thetaMin = min(theta,thetaMin);
        thetaMax = max(theta,thetaMax);
    }
    return max((thetaMax-thetaE)/(PI-thetaE),(thetaE-thetaMin)/thetaE);
}
double Mesh::aspectRatio(const vector<int>& index){
    // http://en.wikipedia.org/wiki/Types_of_mesh#Aspect_ratio
    // aspect ratio of a polygon is max(edgeLength)/min(edgeLength)
    double low = 1000000;
    double top = -1000000;
    for(int i=0;i<index.size();i++){
        Vector v1 = point[index[i]];
        Vector v2 = point[index[(i+1)%index.size()]];
        double len = (v1-v2).length();
        low = min(low,len);
        top = max(top,len);
    }
    return top/low;
}

/*
convert all tri to quad, adjust the remaining quad.

algo:
    if exist a tri:
        1. divide the tri to quad
            this create 3 new points
        2. propagate the three new point to divide the neighbor quad/tri
*/
bool Mesh::forceTriToQuad(){
    // does a tri exist?
    bool hasTri = false;
    int tri;
    for(int i=0;i<surface.size();i++){
        if(surface[i].size()==3){
            hasTri = true;
            tri = i;
            break;
        }
    }
    if(!hasTri){
        return false;
    }
    // create point to poly mapping
    vector<vector<int>> pointToSurface;
    collectPointToSurfaceMapping(&pointToSurface);

    // create poly to neighbor mapping
    vector<vector<int>> neighbor;
    collectNeighborMapping(&neighbor, pointToSurface);

    devideToQuadAndPropagate(tri,neighbor);
    return true;
}


// utility class used only in converting tri to quad
class Poly{
public:
    int id;
    vector<int> neighborId;
    vector<int> vertexId;
    vector<int> childId;
    vector<int> newPointId;
    bool queued;
    // utility function used to find
    // index of point that equals to input p
    // param p: point to search
    // param point: array of Vector, index of array is id of that vector
    // return : index of that point if exist, -1 if not
    int getNewPointIdEqualTo(const Vector& p,const vector<Vector>& point){
        double eps = 1e-10;
        for(auto x:newPointId){
            if((p-point[x]).length() < eps){
                return x;
            }
        }
        return -1;
    }
};


void Mesh::devideToQuadAndPropagate(int startPolyId, const vector<vector<int>>& neigh){

    vector<Poly> polys;
    vector<vector<int>> newSurface;
    for(int i=0;i<surface.size();i++){
        Poly p;
        p.id = i;
        p.neighborId = neigh[i];
        p.vertexId = surface[i];
        p.queued = false;
        polys.push_back(p);
    }
    list<int> idQueue;
    idQueue.push_back(startPolyId);
    polys[startPolyId].queued = true;
    while(idQueue.size()>0){
        Vector center;
        vector<int> newPointId;
        int cId = idQueue.front();
        Poly& cPoly = polys[cId];
        idQueue.pop_front();
        for(int i=0;i<cPoly.vertexId.size();i++){
            Vector a = point[cPoly.vertexId[i]];
            Vector b = point[cPoly.vertexId[(i+1)%cPoly.vertexId.size()]];
            center += a;
            if(cPoly.neighborId[i]==-1){
                newPointId.push_back((int) point.size());
                point.push_back((a+b)*0.5);
            }else{
                Poly np = polys[cPoly.neighborId[i]];
                int pointId = np.getNewPointIdEqualTo((a+b)*0.5,point);
                if(pointId>=0){
                    newPointId.push_back(pointId);
                }else{
                    newPointId.push_back((int) point.size());
                    point.push_back((a+b)*0.5);
                }
            }
        }
        center *= 1.0/cPoly.vertexId.size();
        int centerId = (int) point.size();
        newPointId.push_back(centerId);
        point.push_back(center);
        cPoly.newPointId = newPointId;
        for(int i=0;i<cPoly.vertexId.size();i++){
            int iPlus1 = (i+1)%cPoly.vertexId.size();
            vector<int> cSurface;
            cSurface.push_back(newPointId[i]);
            cSurface.push_back(cPoly.vertexId[iPlus1]);
            cSurface.push_back(newPointId[iPlus1]);
            cSurface.push_back(centerId);
            newSurface.push_back(cSurface);
        }
        for(int i=0;i<cPoly.vertexId.size();i++){
            if(cPoly.neighborId[i] != -1 &&
                polys[cPoly.neighborId[i]].newPointId.size()==0 &&
                polys[cPoly.neighborId[i]].queued == false){
                idQueue.push_back(cPoly.neighborId[i]);
                polys[cPoly.neighborId[i]].queued = true;
            }
        }
    }
    surface = newSurface;
}
void Mesh::collectPointToSurfaceMapping(vector<vector<int>>* out){
    out->resize(point.size());
    for(int s=0;s<surface.size();s++){
        for(int v:surface[s]){
            (*out)[v].push_back(s);
        }
    }
}

void Mesh::collectNeighborMapping(vector<vector<int>>* out, const vector<vector<int>>& pts){
    for(int sid=0;sid<surface.size();sid++){
        auto s = surface[sid];
        // find all side neighbors
        vector<int> neigh;
        for(int i=0;i<s.size();i++){
            // it is a neighbor side AB if it also has point AB
            int a = s[i];
            int b = s[(i+1)%s.size()];
            bool found = false;
            for(auto cA:pts[a]){
                if(cA == sid){
                    continue;
                }
                // does it have point B as well?
                for(auto cB:pts[b]){
                    if(cA==cB){
                        neigh.push_back(cA);
                        found = true;
                        break;
                    }
                }
                if(found){
                    break;
                }
            }
            if(!found){
                neigh.push_back(-1);
            }
        }
        out->push_back(neigh);
    }
}

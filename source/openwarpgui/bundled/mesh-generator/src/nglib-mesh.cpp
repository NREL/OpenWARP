/*

    Copyright (C) 2014 TopCoder Inc., All Rights Reserved.

    @author TCSASSEMBLER, caoweiquan322
    @version 1.1
  */

// nglib-mesh.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

#include "Mesh.h"


#include "TopoDS.hxx"
#include "TopoDS_Face.hxx"
#include "TopoDS_Shape.hxx"
#include "GProp_GProps.hxx"
#include "BRepGProp.hxx"
// param mesh: input Mesh
void printStatistics(Mesh& mesh){
    // print some statistics
    cout << "statistics:" << endl;
    cout << "total area:" << mesh.getSurfaceArea() << endl;
    cout << "total volume:" << mesh.getVolume() << endl;
    cout << "points count: " << mesh.point.size() << endl;
    cout << "edges count:" << mesh.getEdgeCount() << endl;
    cout << "triangle count:" << count_if(mesh.surface.begin(),mesh.surface.end(),
        [](const vector<int>& v){return v.size()==3;}) << endl;
    cout << "quad count:" << count_if(mesh.surface.begin(),mesh.surface.end(),
        [](const vector<int>& v){return v.size()==4;}) << endl;
    cout << "surface (elements) count: " << mesh.surface.size() << endl;
    cout << "is closed surface: " << mesh.isClosed() << endl;
    MeshStats out;
    mesh.stats(&out);
    cout << out.tostring() << endl;

}
// this function is a meshing demo of step format and statistics print
// param: argument to config
void polyDataAndStatisticsSample_step(map<string,string>& config) {

    nglib::Ng_Init();
    nglib::Ng_Mesh * occ_mesh;
    nglib::Ng_Meshing_Parameters mp=nglib::Ng_Meshing_Parameters();

    auto geo = nglib::Ng_OCC_Load_STEP(config["infile"].c_str());
    occ_mesh = nglib::Ng_NewMesh();


    //mp.quad_dominated = 1;
    mp.maxh = atof(config["maxh"].c_str());
    mp.minh = atof(config["minh"].c_str());
    mp.fineness = atof(config["fineness"].c_str());
    mp.grading = atof(config["grading"].c_str());

    nglib::Ng_OCC_SetLocalMeshSize(geo,occ_mesh,&mp);
    cout<<"start edgemesh"<<endl;
    nglib::Ng_OCC_GenerateEdgeMesh(geo,occ_mesh,&mp);

    cout<<"edgemesh ok!!!"<<endl;
    cout<<"start surfacemesh"<<endl;
    nglib::Ng_OCC_GenerateSurfaceMesh(geo,occ_mesh, &mp);
    cout<<"surfacemesh ok!!!"<<endl;

    Mesh m;
    m.feedNgMesh(occ_mesh);
    nglib::Ng_DeleteMesh(occ_mesh);
    m.saveAsSTL(config["outfile"]+".stl");
    m.saveAsVTK(config["outfile"]+".vtk");
    m.saveAsVTP(config["outfile"]+".vtp");
    m.saveAsGDF(config["outfile"]+".gdf");
    // statistics
    printStatistics(m);

    cout << "after call to .forceTriToQuad()" << endl;
    m.forceTriToQuad();
    m.saveAsSTL(config["outfile"]+"-quad.stl");
    m.saveAsVTK(config["outfile"]+"-quad.vtk");
    m.saveAsVTP(config["outfile"]+"-quad.vtp");
    m.saveAsGDF(config["outfile"]+"-quad.gdf");
    m.saveAsDAT(config["outfile"]+"-quad.dat");
    // statistics
    printStatistics(m);

    cout<<"print statistics OK!!!"<<endl;
}

// this function is a meshing demo of igs format and statistics print
// param: argument to config
void polyDataAndStatisticsSample_iges(map<string,string>& config) {

    nglib::Ng_Init();
    nglib::Ng_Mesh * occ_mesh;
    nglib::Ng_Meshing_Parameters mp=nglib::Ng_Meshing_Parameters();

    auto geo = nglib::Ng_OCC_Load_IGES(config["infile"].c_str());
    occ_mesh = nglib::Ng_NewMesh();

    //mp.quad_dominated = 1;
    mp.maxh = atof(config["maxh"].c_str());
    mp.minh = atof(config["minh"].c_str());
    mp.fineness = atof(config["fineness"].c_str());
    mp.grading = atof(config["grading"].c_str());

    nglib::Ng_OCC_SetLocalMeshSize(geo,occ_mesh,&mp);
    cout<<"start edgemesh"<<endl;
    nglib::Ng_OCC_GenerateEdgeMesh(geo,occ_mesh,&mp);

    cout<<"edgemesh ok!!!"<<endl;
    cout<<"start surfacemesh"<<endl;
    nglib::Ng_OCC_GenerateSurfaceMesh(geo,occ_mesh, &mp);
    cout<<"surfacemesh ok!!!"<<endl;

    Mesh m;
    m.feedNgMesh(occ_mesh);
    nglib::Ng_DeleteMesh(occ_mesh);
    m.saveAsSTL(config["outfile"]+".stl");
    m.saveAsVTK(config["outfile"]+".vtk");
    m.saveAsVTP(config["outfile"]+".vtp");
    m.saveAsGDF(config["outfile"]+".gdf");
    // statistics
    printStatistics(m);

    cout << "after call to .forceTriToQuad()" << endl;
    m.forceTriToQuad();
    m.saveAsSTL(config["outfile"]+"-quad.stl");
    m.saveAsVTK(config["outfile"]+"-quad.vtk");
    m.saveAsVTP(config["outfile"]+"-quad.vtp");
    m.saveAsGDF(config["outfile"]+"-quad.gdf");
    m.saveAsDAT(config["outfile"]+"-quad.dat");
    // statistics
    printStatistics(m);

    cout<<"print statistics OK!!!"<<endl;
}

// this function is a meshing demo of stl format and statistics print
// param: argument to config
void polyDataAndStatisticsSample(map<string,string>& config){
    // meshing demo and statistics
    // this also include forceTriToQuad()
    nglib::Ng_Init();
    auto geo = nglib::Ng_STL_LoadGeometry(config["infile"].c_str());
    nglib::Ng_Mesh* mesh = nglib::Ng_NewMesh();
    nglib::Ng_Meshing_Parameters mp = nglib::Ng_Meshing_Parameters();
    //mp.quad_dominated = 1;
    mp.maxh = atof(config["maxh"].c_str());
    mp.minh = atof(config["minh"].c_str());
    mp.fineness = atof(config["fineness"].c_str());
    mp.grading = atof(config["grading"].c_str());


    nglib::Ng_STL_InitSTLGeometry(geo);
    nglib::Ng_STL_MakeEdges(geo,mesh,&mp);

    cout << "Start Surface Meshing...." << endl;
    nglib::Ng_Result result = nglib::Ng_STL_GenerateSurfaceMesh(geo,mesh,&mp);
    if(result != nglib::NG_OK){
        cout << "Error in Surface Meshing....Aborting!!" << endl;
    }

    Mesh m;
    m.feedNgMesh(mesh);
    nglib::Ng_DeleteMesh(mesh);
    m.saveAsSTL(config["outfile"]+".stl");
    m.saveAsVTK(config["outfile"]+".vtk");
    m.saveAsVTP(config["outfile"]+".vtp");

    // statistics
    printStatistics(m);

    cout << "after call to .forceTriToQuad()" << endl;
    m.forceTriToQuad();
    m.saveAsSTL(config["outfile"]+"-quad.stl");
    m.saveAsVTK(config["outfile"]+"-quad.vtk");
    m.saveAsVTP(config["outfile"]+"-quad.vtp");
    m.saveAsGDF(config["outfile"]+"-quad.gdf");
    m.saveAsDAT(config["outfile"]+"-quad.dat");
    // statistics
    printStatistics(m);

    cout<<"print statistics OK!!!"<<endl;
    Mesh original;
    original.loadSTLFile(config["infile"]);
    cout << "original volume: " << original.getVolume() << endl;
    double oriVolume = original.getVolume();
    if(config["usetolerance"].compare("1") == 0){
        double tolerance = 0.01*atof(config["tolerance"].c_str());
        for(int i=0;i<10;i++){
            double cVolume = m.getVolume();
            if(abs((oriVolume-cVolume)/oriVolume) > tolerance){
                mp.minh /= 2;
                mp.maxh /= 2;
                auto geo = nglib::Ng_STL_LoadGeometry(config["infile"].c_str());
                nglib::Ng_Mesh* mesh = nglib::Ng_NewMesh();
                nglib::Ng_STL_InitSTLGeometry(geo);
                nglib::Ng_STL_MakeEdges(geo,mesh,&mp);
                nglib::Ng_Result result = nglib::Ng_STL_GenerateSurfaceMesh(geo,mesh,&mp);
                if(!result){
                    cout << "minh:" << mp.minh << ", maxh:" << mp.maxh <<
                        "result OK." << endl;
                }
                m.feedNgMesh(mesh);
                nglib::Ng_DeleteMesh(mesh);
                m.forceTriToQuad();
            }else{
                printStatistics(m);
                cout << "tolerance:" << abs((oriVolume-cVolume)/oriVolume) << endl;
                cout << "ori volume:" << oriVolume << endl;
                break;
            }

        }
    }

    nglib::Ng_Exit();
}



// trimming utilities
//http://stackoverflow.com/questions/216823/whats-the-best-way-to-trim-stdstring
// trim from start
static inline std::string &ltrim(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
    return s;
}
// trim from end
static inline std::string &rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
    return s;
}
// trim from both ends
// method declaration was modified from:
// static inline std::string &trim(std::string &s)
// to:
// static inline std::string &trim(std::string &s)
// because the original signature with non-const references cannot be used with r-values
// and caused compile error in mac
static inline std::string trim(std::string s) {
    return ltrim(rtrim(s));
}

// load configuration, assume config well written
void loadConfig(const string& fn, map<string,string>& config){
    ifstream fin(fn);
    string line;
    while(getline(fin,line)){
        line = trim(line);
        if(line.length()>0 && line[0] != '#'){
            int index = (int) line.find(':');
            if(index>0){
                string name = line.substr(0,index);
                string value = trim(line.substr(index+1));
                config[name] = value;
            }
        }
    }

    string infile = config["infile"];
    string infiletype = trim(infile.substr(infile.find('.')+1));
    config["infiletype"] = infiletype;
    fin.close();
}

//checks if a string is number or not
//returns true if string is a number otherwise false.
bool isNumber( string myString )
{
    std::istringstream iss(myString);
    float f;
    iss >>f; // noskipws considers leading whitespace invalid
    // Check the entire string was consumed and if either failbit or badbit is set
    return iss.eof() && !iss.fail();
}


//This function checks whether given file exists or not
// return True if file exists otherwise False;
bool checkFileExists(string const &fileName){
    ifstream fileStream(fileName.c_str());
    if(fileStream.good()){
        fileStream.close();
        return true;
    }
    else{
        fileStream.close();
        return false;
    }
}

//This function checks all the parameters of configuration file.
// If there are invalid parameters, it overwrites them with default values and report errors.
// if input file does not exist or filetype is wrong - returns false
bool handleConfigParams(map<string,string> &config) {
    if(!checkFileExists(config["infile"])){
        cout<<"Config parameter: infile "<<config["infile"]<<" does not exist"<<endl;
        //config["infile"]="FlapSingle.stl"; //set the value to the default value
        return false;
    }
	if(config["infiletype"] != "stl" && config["infiletype"] != "step" && config["infiletype"] != "igs")
	{
        cout<<"config parameter :infiletype is not one of stl, step, igs"<<endl;
		//config["infiletype"] = "stl";
        return false;
	}
    if(!isNumber(config["maxh"]))
    {
        if(config["maxh"].empty())
            cout<<"config parameter :maxh is empty"<<endl;
        else
            cout<<"config parameter :maxh "<<config["maxh"]<<" is not a valid number "<<endl;
        config["maxh"]="2.0";
    }
    if(!isNumber(config["minh"]))
    {
        if(config["minh"].empty())
            cout<<"config parameter :minh is empty"<<endl;
        else
            cout<<"config parameter :minh "<<config["minh"]<<" is not a valid number "<<endl;
        config["minh"]="0.5";
    }
    if(!isNumber(config["fineness"]))
    {
        if(config["fineness"].empty())
            cout<<"config parameter :fineness is empty"<<endl;
        else
            cout<<"config parameter :fineness "<<config["fineness"]<<" is not a valid number "<<endl;
        config["fineness"]="0.5";
    }
    if(!isNumber(config["grading"]))
    {
        if(config["grading"].empty())
            cout<<"config parameter :grading is empty"<<endl;
        else
            cout<<"config parameter :grading "<<config["grading"]<<" is not a valid number "<<endl;
        config["grading"]="0.5";
    }
    if(!isNumber(config["usetolerance"]))
    {
        if(config["usetolerance"].empty())
            cout<<"config parameter :usetolerance is empty"<<endl;
        else
            cout<<"config parameter :usetolerance "<<config["usetolerance"]<<" is not a valid number "<<endl;
        config["usetolerance"]="0";
    }
    if(!isNumber(config["tolerance"]))
    {
        if(config["tolerance"].empty())
            cout<<"config parameter :tolerance is empty"<<endl;
        else
            cout<<"config parameter :tolerance "<<config["tolerance"]<<" is not a valid number "<<endl;
        config["tolerance"]="1";
    }
    
    return true;
}

int main(int argc, char* argv[])
{

    map<string,string> config;
    config["infiletype"] = "stl";
    config["infile"] = "FlapSingle.stl";
    config["outfile"] = "flap-meshed";
    config["maxh"] = "2.0";
    config["minh"] = "0.5";
    config["fineness"] = "0.5";
    config["grading"] = "0.5";
    config["usetolerance"] = "0";
    config["tolerance"] = "1";
    
    if (checkFileExists("config.txt")) {
        loadConfig("config.txt",config);
    } else {
        cout << "unable to find config.txt, using default values" << endl;
    }
    
    if (!handleConfigParams(config))
        return 1;

    if(config["infiletype"] == "stl")
        polyDataAndStatisticsSample(config);
    else if(config["infiletype"] == "step" )
        polyDataAndStatisticsSample_step(config);
    else if(config["infiletype"] == "igs")
        polyDataAndStatisticsSample_iges(config);

    return 0;
}


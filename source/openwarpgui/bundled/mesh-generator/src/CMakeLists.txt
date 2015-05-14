cmake_minimum_required(VERSION 2.8)

PROJECT(nglib_mesh)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wno-deprecated -fpermissive") 
find_package(VTK REQUIRED)
find_package(OCE REQUIRED COMPONENTS TKernel TKMath TKAdvTools TKG2d TKG3d TKGeomBase TKBRep TKGeomAlgo TKTopAlgo TKPrim TKBO TKHLR TKMesh TKShHealing TKXMesh TKBool TKFillet TKFeat TKOffset TKService TKV3d TKOpenGl TKMeshVS TKNIS TKVoxel TKCDF PTKernel TKLCAF FWOSPlugin TKPShape TKBinL TKXmlL TKPLCAF TKTObj TKShapeSchema TKStdLSchema TKCAF TKBin TKXml TKPCAF TKBinTObj TKXmlTObj TKStdSchema TKSTL TKXSBase TKSTEPBase TKIGES TKSTEPAttr TKSTEP209 TKSTEP TKVRML TKXCAF TKXCAFSchema TKXmlXCAF TKBinXCAF TKXDEIGES TKXDESTEP TKTopTest TKViewerTest  )

#this section can be moved to FindNETGEN.cmake
SET(NETGEN_INCLUDE_DIRS "/usr/include")
FIND_LIBRARY(NGLIB_LIBRARIES nglib /usr/lib)

include(${VTK_USE_FILE})
include_directories(${OCE_INCLUDE_DIRS})
include_directories(${NETGEN_INCLUDE_DIRS})


add_executable(nglib_mesh nglib-mesh.cpp MeshStats.cpp Vector.cpp stdafx.cpp Mesh.cpp)
add_definitions( -DHAVE_NO_DLL=1 )

target_link_libraries(nglib_mesh ${VTK_LIBRARIES} TKernel TKMath TKAdvTools TKG2d TKG3d TKGeomBase TKBRep TKGeomAlgo TKTopAlgo TKPrim TKBO TKHLR TKMesh TKShHealing TKXMesh TKBool TKFillet TKFeat TKOffset TKService TKV3d TKOpenGl TKMeshVS TKNIS TKVoxel TKCDF PTKernel TKLCAF FWOSPlugin TKPShape TKBinL TKXmlL TKPLCAF TKTObj TKShapeSchema TKStdLSchema TKCAF TKBin TKXml TKPCAF TKBinTObj TKXmlTObj TKStdSchema TKSTL TKXSBase TKSTEPBase TKIGES TKSTEPAttr TKSTEP209 TKSTEP TKVRML TKXCAF TKXCAFSchema TKXmlXCAF TKBinXCAF TKXDEIGES TKXDESTEP TKTopTest TKViewerTest ${NGLIB_LIBRARIES})

message(WARNING "${VTK_LIBRARIES} ${OCE_LIBRARIES} ${NGLIB_LIBRARIES}")

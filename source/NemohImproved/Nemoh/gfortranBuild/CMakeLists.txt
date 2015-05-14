# CMake project file for NEMOH

cmake_minimum_required (VERSION 2.8)
project (NEMOH Fortran)
enable_language (Fortran)

if(NOT CMAKE_BUILD_TYPE)
	SET(CMAKE_BUILD_TYPE Release CACHE STRING "default to debug" FORCE)
endif(NOT CMAKE_BUILD_TYPE)

if(WIN32)

	if ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel")
		set(CMAKE_Fortran_FLAGS "/stand /Qopenmp /Qmkl ")
		set(CMAKE_Fortran_FLAGS_DEBUG "/stand /Qopenmp /Qmkl /fp:source /Od /traceback /check:all /warn:all /debug:all ")
		set(CMAKE_Fortran_FLAGS_RELEASE "/stand /Qopenmp /Qmkl ")
	
		# You can uncomment the following lines for further optimisation but doing so results will not be completely exact with original version.
		#set(CMAKE_Fortran_FLAGS "/stand /Qopenmp /Qmkl /fast")
		
	elseif ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
		set(CMAKE_Fortran_FLAGS "-W -Wall -std=f2008 -fopenmp -fmax-errors=0  -O3 -w")
		set(CMAKE_Fortran_FLAGS_DEBUG "-W -Wall -std=f2008 -fopenmp   -fmax-errors=0 -O3 -w")
		set(CMAKE_Fortran_FLAGS_RELEASE "-W -Wall -std=f2008 -fopenmp   -fmax-errors=0 -w")
	else()
	
		message(SEND_ERROR "Only gfortran and Intel ifort compiler are currently supported.")
	
	endif()
endif (WIN32)

if (UNIX)
	if ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel")
		set(CMAKE_Fortran_FLAGS "-stand -openmp -mkl -heap-arrays 1024")
		set(CMAKE_Fortran_FLAGS_DEBUG "-stand -openmp -mkl -traceback -g -fp-model source -check all -warn all -debug all -heap-arrays 1024")
		set(CMAKE_Fortran_FLAGS_RELEASE "-stand -openmp -mkl -heap-arrays 1024")
	
		# You can uncomment the following lines for further optimisation but doing so results will not be completely exact with original version.
		#set(CMAKE_Fortran_FLAGS "-stand -openmp -mkl -fast")
		
	elseif ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
		set(CMAKE_Fortran_FLAGS "-W -Wall -std=f2008 -fopenmp   -fmax-errors=0  -O3 -w ")
		set(CMAKE_Fortran_FLAGS_DEBUG "-W -Wall -std=f2008 -fopenmp   -fmax-errors=0  -O3 -w ")
		set(CMAKE_Fortran_FLAGS_RELEASE "-W -Wall -std=f2008 -fopenmp   -fmax-errors=0 -O3 -w ")
	
	else()
	
		message(SEND_ERROR "Only gfortran and Intel ifort compiler are currently supported.")
	
	endif()

	find_package(BLAS REQUIRED)
	find_package(LAPACK REQUIRED)
endif (UNIX)





set(CMAKE_EXE_LINKER_FLAGS "")

add_library(nemoh SHARED "../Solver/NEMOH.f90" "../Solver/Core/ALLOCATE_DATA.f90" "../Solver/Core/Bodyconditions.f90" "../Solver/Core/COMPUTE_GREEN_FD.f90" "../Solver/Core/COMPUTE_GREEN_FREESURFACE.f90" "../Solver/Core/COMPUTE_GREEN_INFD.f90" "../Solver/Core/COMPUTE_KOCHIN.f90" "../Solver/Core/COMPUTE_POTENTIAL_DOMAIN.f90" "../Solver/Core/COM_VAR.f90" "../Solver/Core/DEALLOCATE_DATA.f90" "../Solver/Core/ELEMENTARY_FNS.f90" "../Solver/Core/INITIALIZATION.f90" "../Solver/Core/M_SOLVER.f90" "../Solver/Core/OUTPUT.f90" "../Solver/Core/PREPARE_MESH.f90" "../Solver/Core/SOLVE_BEM.f90" "../Solver/Core/SOLVE_BEM_FD_DIRECT.f90"  "../Solver/Core/SOLVE_BEM_INFD_DIRECT.f90"  "../Solver/Core/GMRES.f90" "../Solver/Core/SOLVE_BEM_INFD_ITERATIVE.f90"  "../Solver/Core/SOLVE_BEM_FD_ITERATIVE.f90" "../commonFiles/Environment.f90" "../commonFiles/Identification.f90" "../commonFiles/Results.f90" "../commonFiles/Mesh.f90" "../Solver/Core/UTILITY.f90" "../Solver/Core/COMPUTE_MEAN_DRIFT_FORCES_FIN.f90" "../Solver/Core/COMPUTE_MEAN_DRIFT_FORCES_INF.f90"
"../Solver/Core/COMPUTE_YAW_MOMENT_FIN.f90" "../Solver/Core/COMPUTE_YAW_MOMENT_INF.f90")
#add_executable(preProcessor "preProcessor/Main.f90" "commonFiles/Environment.f90" "commonFiles/Identification.f90" "preProcessor/Mesh.f90" "preProcessor/BodyConditions.f90" "preProcessor/Integration.f90" )

#add_executable(postProcessor "postProcessor/Main.f90" "commonFiles/Identification.f90" "commonFiles/Environment.f90" "commonFiles/Results.f90" "postProcessor/IRF.f90" "postProcessor/Plot_WaveElevation.f90" "postProcessor/Compute_RAOs.f90" )

#add_executable(Mesh "commonFiles/Identification.f90" "Mesh/mesh.f90" "Mesh/Netting.f90"  "Mesh/hydre.f90")

#add_executable(Solver "Solver/NEMOH.f90" "Solver/Core/ALLOCATE_DATA.f90" "Solver/Core/Bodyconditions.f90" "Solver/Core/COMPUTE_GREEN_FD.f90" "Solver/Core/COMPUTE_GREEN_FREESURFACE.f90" "Solver/Core/COMPUTE_GREEN_INFD.f90" "Solver/Core/COMPUTE_KOCHIN.f90" "Solver/Core/COMPUTE_POTENTIAL_DOMAIN.f90" "Solver/Core/COM_VAR.f90" "Solver/Core/DEALLOCATE_DATA.f90" "Solver/Core/ELEMENTARY_FNS.f90" "Solver/Core/INITIALIZATION.f90" "Solver/Core/M_SOLVER.f90" "Solver/Core/OUTPUT.f90" "Solver/Core/PREPARE_MESH.f90" "Solver/Core/SOLVE_BEM.f90" "Solver/Core/SOLVE_BEM_FD_DIRECT.f90"  "Solver/Core/SOLVE_BEM_INFD_DIRECT.f90"  "Solver/Core/GMRES.f90" "Solver/Core/SOLVE_BEM_INFD_ITERATIVE.f90"  "Solver/Core/SOLVE_BEM_FD_ITERATIVE.f90" "commonFiles/Environment.f90" "commonFiles/Identification.f90" "commonFiles/Results.f90" "commonFiles/Mesh.f90")


if ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")

	SET_TARGET_PROPERTIES(nemoh PROPERTIES COMPILE_FLAGS "-fopenmp" LINK_FLAGS "-fopenmp")
	SET_TARGET_PROPERTIES(nemoh PROPERTIES COMPILE_FLAGS "-fopenmp -lblas -llapack" LINK_FLAGS "-fopenmp -lblas -llapack")
	target_link_libraries(nemoh ${BLAS_LIBRARIES})
	target_link_libraries(nemoh ${LAPACK_LIBRARIES})
	target_link_libraries(nemoh ${OpenMP_LIBRARIES})
endif()
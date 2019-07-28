cmake_minimum_required(VERSION 3.0.2)


SET (PROJECTNAME "esther")

PROJECT (${PROJECTNAME} Fortran)

ENABLE_LANGUAGE(C) 

#FIND_PACKAGE(LAPACK REQUIRED)

#FIND_PACKAGE(HDF5 COMPONENTS Fortran REQUIRED)
#INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIR})


#SET (CMAKE_Fortran_FLAGS "-O3 -cpp -mtune=native -fopenmp")
set(CMAKE_BUILD_TYPE Debug
        CACHE STRING "Build type (Debug, Release)" FORCE)
#set(warn "-Wall -Wextra -Wimplicit-interface -fPIC")# -fmax-errors=1")
set(CMAKE_Fortran_FLAGS_DEBUG   "${common} ${warn} -g -fbounds-check -fcheck-array-temporaries -fbacktrace")


ADD_EXECUTABLE(${PROJECTNAME}
        blf2in3.f  caxsrefl3d.f  gen2in3vec.f  pccslabc.f  secularc.f  stack.f
        blmy.f     dlmsf2in3.f   order2dl.f    pccslab.f   secular.f   tmtaxsp.f

        "../chebev.f"
        "../beschb.f"
        "../../../3d/biga.f"
        "../bessik.f"
        "../bessjy.f"
        "../../../3d/gnricbessh.f"
        "../gnzbess.f"
        "../zbessf.f"

        )


target_link_libraries(${PROJECTNAME} #${BLAS_LIBRARIES}
        #${LAPACK_LIBRARIES} ${HDF5_Fortran_LIBRARIES}
        )



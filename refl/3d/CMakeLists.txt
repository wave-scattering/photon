cmake_minimum_required(VERSION 3.0.2)


SET (PROJECTNAME "c1axsrefl3d")

PROJECT (${PROJECTNAME} Fortran)

ENABLE_LANGUAGE(C) 

#FIND_PACKAGE(LAPACK REQUIRED)

#FIND_PACKAGE(HDF5 COMPONENTS Fortran REQUIRED)
#INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIR})


#SET (CMAKE_Fortran_FLAGS "-O3 -cpp -mtune=native -fopenmp")


ADD_EXECUTABLE(${PROJECTNAME} "c1axsrefl3d.f"
#        "bess.f"
        "chebev.f"
        "beschb.f"
        "../../3d/biga.f"
        "bessik.f"
        "bessjy.f"
        "../../3d/gnricbessh.f"
        "gnzbess.f"
        "zbessf.f"
#        "ampld.par.f"
        "blf2in3.f"
        "blmy.f"
        "dlmsf2in3.f"
        "gen2in3vec.f"
        "latgen2d.f"
        "../../3d/mediumn.f"
        "order2dl.f"
        "pccslab.f"
        "pccslabc.f"
        "../../3d/readmat.f"
        "secular.f"
        "secularc.f"
        "sordalc.f"
        "stack.f"
        "tmtaxsp1.f"
        "../../zartan.f"
        "zcopya.f"
        "zcopyv.f"
        "znsrefind.f")


target_link_libraries(${PROJECTNAME} #${BLAS_LIBRARIES}
        #${LAPACK_LIBRARIES} ${HDF5_Fortran_LIBRARIES}
        )



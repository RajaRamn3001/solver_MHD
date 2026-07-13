add_library(read_input   read_input.f90)
add_library(grid_generate   grid_generate.f90)
add_library(array_processes   array_processes.f90)
add_library(build_linear_system   build_linear_system.f90)
add_library(matrix_vector_op   matrix_vector_op.f90)
add_library(linear_solvers   linear_solvers.f90)
add_library(green_function_sub   green_function_sub.f90)
add_library(symbolic_LU_op   symbolic_LU_op.f90)


# Express library-on-library deps
target_link_libraries(symbolic_LU_op 
    PUBLIC 
        array_processes        
        build_linear_system   
)

target_link_libraries(linear_solvers 
    PUBLIC 
        symbolic_LU_op
        array_processes
        matrix_vector_op        
        build_linear_system   
)



add_executable(myexe main_DH_BEM.f90)
#add_executable(myexe main_DH_pseudo_vac.f90)
#add_executable(myexe main_test.f90)

# Usually enough to link only what you use directly
target_link_libraries(myexe 
             PRIVATE 
                    read_input
                    grid_generate 
                    build_linear_system 
                    matrix_vector_op
                    linear_solvers
                    green_function_sub
                   
 )

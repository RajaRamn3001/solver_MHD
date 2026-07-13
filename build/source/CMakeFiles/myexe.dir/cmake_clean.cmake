file(REMOVE_RECURSE
  "../run/myexe"
  "../run/myexe.pdb"
  "CMakeFiles/myexe.dir/main_test.f90.o"
)

# Per-language clean rules from dependency scanning.
foreach(lang Fortran)
  include(CMakeFiles/myexe.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()

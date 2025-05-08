add_library(eigen_types STATIC ${PROJECT_SOURCE_DIR}/src/math/eigen_types.h)
target_link_libraries(eigen_types Eigen3:Eigen)
set_target_properties(eigen_types PROPERTIES LINKER_LANGUAGE CXX)

add_subdirectory(integrator)
add_subdirectory(unscented_transform)

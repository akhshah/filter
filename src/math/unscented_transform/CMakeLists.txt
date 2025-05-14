add_library(unscented_transform SHARED
    ${PROEJCT_SOURCE_DIR}/src/math/unscented_transform/unscented_transform.h)
target_link_libraries(unscented_transform eigen_types Eigen3:Eigen)
set_target_properties(uscented_transform PROPERTIES LINKER_LANGUAGE CXX)

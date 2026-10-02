cmake_minimum_required(VERSION 3.18)
project(qsiNative LANGUAGES C CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
find_package(ITK 5 REQUIRED COMPONENTS
  ITKCommon ITKTransform ITKImageGrid ITKImageIntensity
  ITKIONIFTI ITKIOTransformHDF5 ITKIOTransformMatlab ITKIOTransformInsightLegacy)
include(${ITK_USE_FILE})
add_library(qsiNative SHARED "${QSI_NATIVE_SOURCE}" "${RCPP_EXPORTS}")
target_include_directories(qsiNative PRIVATE "${R_INCLUDE_DIR}" "${RCPP_INCLUDE_DIR}")
target_link_libraries(qsiNative PRIVATE ${ITK_LIBRARIES} "${R_LIBRARY}")
set_target_properties(qsiNative PROPERTIES PREFIX "" SUFFIX "${R_DLL_SUFFIX}"
  BUILD_WITH_INSTALL_RPATH TRUE INSTALL_RPATH_USE_LINK_PATH TRUE)
install(TARGETS qsiNative LIBRARY DESTINATION . RUNTIME DESTINATION .)

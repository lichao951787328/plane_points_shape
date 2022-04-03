cmake_minimum_required(VERSION 3.0.2)
project(use_py_to_cpp)

## Compile as C++11, supported in ROS Kinetic and newer
add_compile_options(-std=c++17)
find_package(Eigen3 REQUIRED)
find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
find_package(Python3 COMPONENTS NumPy)
find_package(PythonLibs REQUIRED)
# find_package(Python3 COMPONENTS Development NumPy)
# target_include_directories(myproject PRIVATE ${Python2_INCLUDE_DIRS} ${Python2_NumPy_INCLUDE_DIRS})
# target_link_libraries(myproject Python2::Python Python2::NumPy)
find_package(PCL REQUIRED)
# find_package(PythonLibs REQUIRED)



include_directories(
  ${EIGEN3_INCLUDE_DIRS}
  ${PYTHON_INCLUDE_DIRS}
#   ${PYTHON_INCLUDE_DIRS}
  ${PROJECT_SOURCE_DIR}/include
)
# message("------------${PROJECT_SOURCE_DIR}")
# file(GLOB USER_LIBS_PATH ${PROJECT_SOURCE_DIR}/dt/*.cpp)
# aux_source_directory(${PROJECT_SOURCE_DIR}/dt DIR_SRCS)

add_executable(main main.cpp dt/delaunay.cpp dt/edge.cpp dt/triangle.cpp dt/vector2.cpp)
add_executable(basic basic.cpp)

target_link_libraries(main
# ${PYTHON_LIBRARIES}
${EIGEN3_LIBRARIES}
${PYTHON_LIBRARIES}
)

target_link_libraries(basic
# ${PYTHON_LIBRARIES}
${EIGEN3_LIBRARIES}
${PYTHON_LIBRARIES}
)
# target_link_libraries(test_new
# # ${PYTHON_LIBRARIES}
# # ${PCL_LIBRARIES}
# )
cmake_minimum_required(VERSION 3.16)
project(example)

list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Libigl
option(LIBIGL_WITH_OPENGL            "Use OpenGL"         OFF)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           OFF)
option(LIBIGL_WITH_CGAL       "Use CGAL"           ON)
include(libigl)

# Add your project files
file(GLOB SRC_FILES *.cpp)
add_executable(${PROJECT_NAME} ${SRC_FILES})
target_link_libraries(${PROJECT_NAME} PUBLIC igl::core igl::cgal)

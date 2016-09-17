
cmake_minimum_required(VERSION 2.8)

PROJECT(AnglePicture)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(AnglePicture MACOSX_BUNDLE AnglePicture)

if(VTK_LIBRARIES)
  target_link_libraries(AnglePicture ${VTK_LIBRARIES})
else()
  target_link_libraries(AnglePicture vtkHybrid vtkWidgets)
endif()

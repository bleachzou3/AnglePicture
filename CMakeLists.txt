
cmake_minimum_required(VERSION 2.8)

PROJECT(AnglePicture)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})
set(Srcs AnglePicture.cpp ONB.cpp)
set(Hdrs ONB.h Vector3.h)
add_executable(AnglePicture MACOSX_BUNDLE ${Srcs} ${Hdrs})

if(VTK_LIBRARIES)
  target_link_libraries(AnglePicture ${VTK_LIBRARIES})
else()
  target_link_libraries(AnglePicture vtkHybrid vtkWidgets)
endif()

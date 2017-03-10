
cmake_minimum_required(VERSION 2.8)

PROJECT(AnglePicture)

find_package(VTK REQUIRED)
find_package(ITK REQUIRED)
include(${VTK_USE_FILE})
include(${ITK_USE_FILE})
include_directories(
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}  
)
set(CMAKE_AUTOMOC ON)
set(Srcs AnglePicture.cpp ONB.cpp InteractorStyleRollBall.cpp AnglePictureUtility.cpp)
set(Hdrs ONB.h Vector3.h InteractorStyleRollBall.h AnglePictureUtility.h)
set(LOGs log4cpp.conf )
set(LinkedLog log4cpp)

 source_group("Sources" FILES
   ${Srcs}
  )
   source_group("Headers" FILES
   ${Hdrs}
  )

MESSAGE(STATUS "{VTK_USE_FILE}" ${VTK_USE_FILE})
MESSAGE(STATUS "{ITK_USE_FILE}" ${VTK_USE_FILE})
MESSAGE(STATUS "{CMAKE_CURRENT_SOURCE_DIR}" ${CMAKE_CURRENT_SOURCE_DIR})
MESSAGE(STATUS "{CMAKE_CURRENT_BINARY_DIR}" ${CMAKE_CURRENT_BINARY_DIR})


include_directories(${CMAKE_CURRENT_SOURCE_DIR}/LOG_HOME/include ..)
link_directories(${CMAKE_CURRENT_SOURCE_DIR}/LOG_HOME/lib)  
add_executable(AnglePicture MACOSX_BUNDLE ${Srcs} ${Hdrs})

if(VTK_LIBRARIES)
  target_link_libraries(AnglePicture ${VTK_LIBRARIES} ${LinkedLog} )
  MESSAGE(STATUS "VTK_LIBRARIES" ${VTK_LIBRARIES})
else()
  target_link_libraries(AnglePicture vtkHybrid vtkWidgets)
endif()


if(ITK_LIBRARIES)
  target_link_libraries(AnglePicture ${ITK_LIBRARIES})
  MESSAGE(STATUS "ITK_LIBRARIES" ${ITK_LIBRARIES})
endif()


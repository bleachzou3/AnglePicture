
cmake_minimum_required(VERSION 2.8)

PROJECT(AnglePicture)

find_package(ITK REQUIRED NO_MODULE)



find_package(VTK REQUIRED NO_MODULE)
#set(VTK_DIR D:/VTK_H/lib/cmake/vtk-6.3)
#set (VTK_USE_FILE D:/VTK_H/lib/cmake/vtk-6.3/UseVTK.cmake)

#set(ITK_DIR E:/ITK_4.11/lib/cmake/ITK-4.11)


#set (ITK_USE_FILE E:/ITK_4.11/lib/cmake/ITK-4.11/UseITK.cmake)
include(${VTK_USE_FILE})
include(${ITK_USE_FILE})
include_directories(
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}  
)
set(CMAKE_AUTOMOC ON)
#set(Srcs AnglePicture.cpp ONB.cpp InteractorStyleRollBall.cpp AnglePictureUtility.cpp DisplayVoxelUtility.cpp AnglePictureSegmentation.cpp)
set(Srcs  ONB.cpp InteractorStyleRollBall.cpp AnglePictureUtility.cpp DisplayVoxelUtility.cpp AnglePictureSegmentation.cpp)
set(Hdrs ONB.h Vector3.h InteractorStyleRollBall.h AnglePictureUtility.h DisplayVoxelUtility.h MeanPoint.h MeanShiftAlgo.h Common.h CommonUtility.h MeanShiftAlgo.h)
set(LOGs log4cpp.conf )
set(LinkedLog log4cpp)

 source_group("Sources" FILES
   ${Srcs}
  )
   source_group("Headers" FILES
   ${Hdrs}
  )

MESSAGE(STATUS "{VTK_USE_FILE}" ${VTK_USE_FILE})
MESSAGE(STATUS "{ITK_USE_FILE}" ${ITK_USE_FILE})
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


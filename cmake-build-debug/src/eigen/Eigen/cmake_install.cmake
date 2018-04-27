# Install script for directory: /home/maytee/Documents/2.29/finiteVolumeSolver/src/eigen/Eigen

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Debug")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xDevelx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/eigen3/Eigen" TYPE FILE FILES
    "/home/maytee/Documents/2.29/finiteVolumeSolver/src/eigen/Eigen/Cholesky"
    "/home/maytee/Documents/2.29/finiteVolumeSolver/src/eigen/Eigen/CholmodSupport"
    "/home/maytee/Documents/2.29/finiteVolumeSolver/src/eigen/Eigen/Core"
    "/home/maytee/Documents/2.29/finiteVolumeSolver/src/eigen/Eigen/Dense"
    "/home/maytee/Documents/2.29/finiteVolumeSolver/src/eigen/Eigen/Eigen"
    "/home/maytee/Documents/2.29/finiteVolumeSolver/src/eigen/Eigen/Eigenvalues"
    "/home/maytee/Documents/2.29/finiteVolumeSolver/src/eigen/Eigen/Geometry"
    "/home/maytee/Documents/2.29/finiteVolumeSolver/src/eigen/Eigen/Householder"
    "/home/maytee/Documents/2.29/finiteVolumeSolver/src/eigen/Eigen/IterativeLinearSolvers"
    "/home/maytee/Documents/2.29/finiteVolumeSolver/src/eigen/Eigen/Jacobi"
    "/home/maytee/Documents/2.29/finiteVolumeSolver/src/eigen/Eigen/LU"
    "/home/maytee/Documents/2.29/finiteVolumeSolver/src/eigen/Eigen/MetisSupport"
    "/home/maytee/Documents/2.29/finiteVolumeSolver/src/eigen/Eigen/OrderingMethods"
    "/home/maytee/Documents/2.29/finiteVolumeSolver/src/eigen/Eigen/PaStiXSupport"
    "/home/maytee/Documents/2.29/finiteVolumeSolver/src/eigen/Eigen/PardisoSupport"
    "/home/maytee/Documents/2.29/finiteVolumeSolver/src/eigen/Eigen/QR"
    "/home/maytee/Documents/2.29/finiteVolumeSolver/src/eigen/Eigen/QtAlignedMalloc"
    "/home/maytee/Documents/2.29/finiteVolumeSolver/src/eigen/Eigen/SPQRSupport"
    "/home/maytee/Documents/2.29/finiteVolumeSolver/src/eigen/Eigen/SVD"
    "/home/maytee/Documents/2.29/finiteVolumeSolver/src/eigen/Eigen/Sparse"
    "/home/maytee/Documents/2.29/finiteVolumeSolver/src/eigen/Eigen/SparseCholesky"
    "/home/maytee/Documents/2.29/finiteVolumeSolver/src/eigen/Eigen/SparseCore"
    "/home/maytee/Documents/2.29/finiteVolumeSolver/src/eigen/Eigen/SparseLU"
    "/home/maytee/Documents/2.29/finiteVolumeSolver/src/eigen/Eigen/SparseQR"
    "/home/maytee/Documents/2.29/finiteVolumeSolver/src/eigen/Eigen/StdDeque"
    "/home/maytee/Documents/2.29/finiteVolumeSolver/src/eigen/Eigen/StdList"
    "/home/maytee/Documents/2.29/finiteVolumeSolver/src/eigen/Eigen/StdVector"
    "/home/maytee/Documents/2.29/finiteVolumeSolver/src/eigen/Eigen/SuperLUSupport"
    "/home/maytee/Documents/2.29/finiteVolumeSolver/src/eigen/Eigen/UmfPackSupport"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xDevelx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/eigen3/Eigen" TYPE DIRECTORY FILES "/home/maytee/Documents/2.29/finiteVolumeSolver/src/eigen/Eigen/src" FILES_MATCHING REGEX "/[^/]*\\.h$")
endif()


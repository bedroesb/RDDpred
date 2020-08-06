# Install script for directory: /data/project/breast_cancer_han/RDD_Project/RDDpred.Dir/Test.Dir/RDDpred_v1.0.Dir/ToolBox.Dir/bamtools-2.4.0/src

# Set the install prefix
IF(NOT DEFINED CMAKE_INSTALL_PREFIX)
  SET(CMAKE_INSTALL_PREFIX "/usr/local")
ENDIF(NOT DEFINED CMAKE_INSTALL_PREFIX)
STRING(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
IF(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  IF(BUILD_TYPE)
    STRING(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  ELSE(BUILD_TYPE)
    SET(CMAKE_INSTALL_CONFIG_NAME "Release")
  ENDIF(BUILD_TYPE)
  MESSAGE(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
ENDIF(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)

# Set the component getting installed.
IF(NOT CMAKE_INSTALL_COMPONENT)
  IF(COMPONENT)
    MESSAGE(STATUS "Install component: \"${COMPONENT}\"")
    SET(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  ELSE(COMPONENT)
    SET(CMAKE_INSTALL_COMPONENT)
  ENDIF(COMPONENT)
ENDIF(NOT CMAKE_INSTALL_COMPONENT)

# Install shared libraries without execute permission?
IF(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  SET(CMAKE_INSTALL_SO_NO_EXE "0")
ENDIF(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/bamtools/shared" TYPE FILE FILES "/data/project/breast_cancer_han/RDD_Project/RDDpred.Dir/Test.Dir/RDDpred_v1.0.Dir/ToolBox.Dir/bamtools-2.4.0/src/shared/bamtools_global.h")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  INCLUDE("/data/project/breast_cancer_han/RDD_Project/RDDpred.Dir/Test.Dir/RDDpred_v1.0.Dir/ToolBox.Dir/bamtools-2.4.0/build/src/api/cmake_install.cmake")
  INCLUDE("/data/project/breast_cancer_han/RDD_Project/RDDpred.Dir/Test.Dir/RDDpred_v1.0.Dir/ToolBox.Dir/bamtools-2.4.0/build/src/third_party/cmake_install.cmake")
  INCLUDE("/data/project/breast_cancer_han/RDD_Project/RDDpred.Dir/Test.Dir/RDDpred_v1.0.Dir/ToolBox.Dir/bamtools-2.4.0/build/src/toolkit/cmake_install.cmake")
  INCLUDE("/data/project/breast_cancer_han/RDD_Project/RDDpred.Dir/Test.Dir/RDDpred_v1.0.Dir/ToolBox.Dir/bamtools-2.4.0/build/src/utils/cmake_install.cmake")

ENDIF(NOT CMAKE_INSTALL_LOCAL_ONLY)


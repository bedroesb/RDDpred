# Install script for directory: /data/project/breast_cancer_han/RDD_Project/RDDpred.Dir/Test.Dir/RDDpred_v1.0.Dir/ToolBox.Dir/bamtools-2.4.0/src/api

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
  FOREACH(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/bamtools/libbamtools.so.2.4.0"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/bamtools/libbamtools.so"
      )
    IF(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      FILE(RPATH_CHECK
           FILE "${file}"
           RPATH "")
    ENDIF()
  ENDFOREACH()
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/bamtools" TYPE SHARED_LIBRARY FILES
    "/data/project/breast_cancer_han/RDD_Project/RDDpred.Dir/Test.Dir/RDDpred_v1.0.Dir/ToolBox.Dir/bamtools-2.4.0/lib/libbamtools.so.2.4.0"
    "/data/project/breast_cancer_han/RDD_Project/RDDpred.Dir/Test.Dir/RDDpred_v1.0.Dir/ToolBox.Dir/bamtools-2.4.0/lib/libbamtools.so"
    )
  FOREACH(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/bamtools/libbamtools.so.2.4.0"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/bamtools/libbamtools.so"
      )
    IF(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      IF(CMAKE_INSTALL_DO_STRIP)
        EXECUTE_PROCESS(COMMAND "/usr/bin/strip" "${file}")
      ENDIF(CMAKE_INSTALL_DO_STRIP)
    ENDIF()
  ENDFOREACH()
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/bamtools" TYPE STATIC_LIBRARY FILES "/data/project/breast_cancer_han/RDD_Project/RDDpred.Dir/Test.Dir/RDDpred_v1.0.Dir/ToolBox.Dir/bamtools-2.4.0/lib/libbamtools.a")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/bamtools/api" TYPE FILE FILES "/data/project/breast_cancer_han/RDD_Project/RDDpred.Dir/Test.Dir/RDDpred_v1.0.Dir/ToolBox.Dir/bamtools-2.4.0/src/api/api_global.h")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/bamtools/api" TYPE FILE FILES "/data/project/breast_cancer_han/RDD_Project/RDDpred.Dir/Test.Dir/RDDpred_v1.0.Dir/ToolBox.Dir/bamtools-2.4.0/src/api/BamAlgorithms.h")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/bamtools/api" TYPE FILE FILES "/data/project/breast_cancer_han/RDD_Project/RDDpred.Dir/Test.Dir/RDDpred_v1.0.Dir/ToolBox.Dir/bamtools-2.4.0/src/api/BamAlignment.h")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/bamtools/api" TYPE FILE FILES "/data/project/breast_cancer_han/RDD_Project/RDDpred.Dir/Test.Dir/RDDpred_v1.0.Dir/ToolBox.Dir/bamtools-2.4.0/src/api/BamAux.h")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/bamtools/api" TYPE FILE FILES "/data/project/breast_cancer_han/RDD_Project/RDDpred.Dir/Test.Dir/RDDpred_v1.0.Dir/ToolBox.Dir/bamtools-2.4.0/src/api/BamConstants.h")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/bamtools/api" TYPE FILE FILES "/data/project/breast_cancer_han/RDD_Project/RDDpred.Dir/Test.Dir/RDDpred_v1.0.Dir/ToolBox.Dir/bamtools-2.4.0/src/api/BamIndex.h")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/bamtools/api" TYPE FILE FILES "/data/project/breast_cancer_han/RDD_Project/RDDpred.Dir/Test.Dir/RDDpred_v1.0.Dir/ToolBox.Dir/bamtools-2.4.0/src/api/BamMultiReader.h")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/bamtools/api" TYPE FILE FILES "/data/project/breast_cancer_han/RDD_Project/RDDpred.Dir/Test.Dir/RDDpred_v1.0.Dir/ToolBox.Dir/bamtools-2.4.0/src/api/BamReader.h")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/bamtools/api" TYPE FILE FILES "/data/project/breast_cancer_han/RDD_Project/RDDpred.Dir/Test.Dir/RDDpred_v1.0.Dir/ToolBox.Dir/bamtools-2.4.0/src/api/BamWriter.h")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/bamtools/api" TYPE FILE FILES "/data/project/breast_cancer_han/RDD_Project/RDDpred.Dir/Test.Dir/RDDpred_v1.0.Dir/ToolBox.Dir/bamtools-2.4.0/src/api/IBamIODevice.h")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/bamtools/api" TYPE FILE FILES "/data/project/breast_cancer_han/RDD_Project/RDDpred.Dir/Test.Dir/RDDpred_v1.0.Dir/ToolBox.Dir/bamtools-2.4.0/src/api/SamConstants.h")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/bamtools/api" TYPE FILE FILES "/data/project/breast_cancer_han/RDD_Project/RDDpred.Dir/Test.Dir/RDDpred_v1.0.Dir/ToolBox.Dir/bamtools-2.4.0/src/api/SamHeader.h")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/bamtools/api" TYPE FILE FILES "/data/project/breast_cancer_han/RDD_Project/RDDpred.Dir/Test.Dir/RDDpred_v1.0.Dir/ToolBox.Dir/bamtools-2.4.0/src/api/SamProgram.h")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/bamtools/api" TYPE FILE FILES "/data/project/breast_cancer_han/RDD_Project/RDDpred.Dir/Test.Dir/RDDpred_v1.0.Dir/ToolBox.Dir/bamtools-2.4.0/src/api/SamProgramChain.h")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/bamtools/api" TYPE FILE FILES "/data/project/breast_cancer_han/RDD_Project/RDDpred.Dir/Test.Dir/RDDpred_v1.0.Dir/ToolBox.Dir/bamtools-2.4.0/src/api/SamReadGroup.h")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/bamtools/api" TYPE FILE FILES "/data/project/breast_cancer_han/RDD_Project/RDDpred.Dir/Test.Dir/RDDpred_v1.0.Dir/ToolBox.Dir/bamtools-2.4.0/src/api/SamReadGroupDictionary.h")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/bamtools/api" TYPE FILE FILES "/data/project/breast_cancer_han/RDD_Project/RDDpred.Dir/Test.Dir/RDDpred_v1.0.Dir/ToolBox.Dir/bamtools-2.4.0/src/api/SamSequence.h")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/bamtools/api" TYPE FILE FILES "/data/project/breast_cancer_han/RDD_Project/RDDpred.Dir/Test.Dir/RDDpred_v1.0.Dir/ToolBox.Dir/bamtools-2.4.0/src/api/SamSequenceDictionary.h")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/bamtools/api/algorithms" TYPE FILE FILES "/data/project/breast_cancer_han/RDD_Project/RDDpred.Dir/Test.Dir/RDDpred_v1.0.Dir/ToolBox.Dir/bamtools-2.4.0/src/api/algorithms/Sort.h")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  INCLUDE("/data/project/breast_cancer_han/RDD_Project/RDDpred.Dir/Test.Dir/RDDpred_v1.0.Dir/ToolBox.Dir/bamtools-2.4.0/build/src/api/internal/cmake_install.cmake")

ENDIF(NOT CMAKE_INSTALL_LOCAL_ONLY)


cmake_minimum_required(VERSION 3.16)
project(MetaGSCA_cpp VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ===========================================================================
# Dependencies:
#   - Eigen  (linear algebra, eigen decomposition)
#   - OpenMP (optional, for parallelism)
#   - nlopt  (nonlinear optimization for ZINB MLE)
# ===========================================================================
find_package(Eigen3 REQUIRED NO_MODULE)
find_package(OpenMP QUIET)

# On macOS with AppleClang, find_package(OpenMP) often fails its compile test
# even though libomp is installed.  Detect and force-enable manually.
if(NOT OpenMP_CXX_FOUND AND APPLE)
    find_path(OMP_INCLUDE_DIR omp.h
        HINTS /opt/homebrew/opt/libomp/include /usr/local/opt/libomp/include)
    find_library(OMP_LIBRARY omp
        HINTS /opt/homebrew/opt/libomp/lib /usr/local/opt/libomp/lib)
    if(OMP_INCLUDE_DIR AND OMP_LIBRARY)
        message(STATUS "Manually enabling OpenMP: ${OMP_LIBRARY}")
        set(OPENMP_MANUAL_FOUND TRUE)
    else()
        message(STATUS "OpenMP not found. Building single-threaded.")
    endif()
endif()

# nlopt: try multiple strategies to find it
find_package(NLopt QUIET CONFIG)
if(NOT NLopt_FOUND)
    # Homebrew on macOS often installs nlopt without CMake config files
    find_package(PkgConfig QUIET)
    if(PkgConfig_FOUND)
        pkg_check_modules(NLOPT IMPORTED_TARGET nlopt)
    endif()
endif()
if(NOT NLopt_FOUND AND NOT NLOPT_FOUND)
    # Last resort: search for the header and library manually
    find_path(NLOPT_INCLUDE_DIR nlopt.h
        HINTS /opt/homebrew/include /usr/local/include /usr/include)
    find_library(NLOPT_LIBRARY nlopt
        HINTS /opt/homebrew/lib /usr/local/lib /usr/lib)
    if(NLOPT_INCLUDE_DIR AND NLOPT_LIBRARY)
        set(NLOPT_MANUAL_FOUND TRUE)
        message(STATUS "Found NLopt manually: ${NLOPT_LIBRARY}")
    else()
        message(FATAL_ERROR "NLopt not found. Install with: brew install nlopt")
    endif()
endif()

# ===========================================================================
# Library: core MetaGSCA functions
# ===========================================================================
add_library(metagsca_core STATIC
    src/DataFrame.cpp
    src/gsnca_check.cpp
    src/gsnca_sd.cpp
    src/gsnca_filt.cpp
    src/zinb_model.cpp
    src/zinb_qr_residuals.cpp
    src/gsnca_stat_qr.cpp
    src/gsnca_perm_dc_qr.cpp
    src/gsnca_perm_ct_qr.cpp
    src/MetaGSCA.cpp
    src/tsv_io.cpp
    src/rng.cpp
)

target_include_directories(metagsca_core PUBLIC include)
target_link_libraries(metagsca_core PUBLIC Eigen3::Eigen)

if(NLopt_FOUND)
    target_link_libraries(metagsca_core PUBLIC NLopt::nlopt)
elseif(NLOPT_FOUND)
    target_link_libraries(metagsca_core PUBLIC PkgConfig::NLOPT)
elseif(NLOPT_MANUAL_FOUND)
    target_include_directories(metagsca_core PUBLIC ${NLOPT_INCLUDE_DIR})
    target_link_libraries(metagsca_core PUBLIC ${NLOPT_LIBRARY})
else()
    message(FATAL_ERROR "NLopt not found.")
endif()

if(OpenMP_CXX_FOUND)
    target_link_libraries(metagsca_core PUBLIC OpenMP::OpenMP_CXX)
    target_compile_definitions(metagsca_core PUBLIC HAS_OPENMP)
    message(STATUS "OpenMP enabled via find_package")
elseif(OPENMP_MANUAL_FOUND)
    target_include_directories(metagsca_core PUBLIC ${OMP_INCLUDE_DIR})
    target_compile_options(metagsca_core PUBLIC -Xpreprocessor -fopenmp)
    target_link_libraries(metagsca_core PUBLIC ${OMP_LIBRARY})
    target_compile_definitions(metagsca_core PUBLIC HAS_OPENMP)
    message(STATUS "OpenMP enabled manually for macOS")
else()
    message(STATUS "OpenMP NOT enabled — running single-threaded")
endif()

# ===========================================================================
# Executable: batch runner (replaces examples_all_kegg.R)
# ===========================================================================
add_executable(scgsca_run src/main.cpp)
target_link_libraries(scgsca_run PRIVATE metagsca_core)

include(CTest)
if(BUILD_TESTING)
    add_executable(test_ct_tail tests/test_ct_tail.cpp)
    target_include_directories(test_ct_tail PRIVATE include)
    add_test(NAME fixed_ct_tail COMMAND test_ct_tail)
endif()

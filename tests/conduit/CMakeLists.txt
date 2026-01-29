#=============================================================================
# Copyright (C) 2026 Commissariat a l'energie atomique et aux energies alternatives (CEA)
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# * Redistributions of source code must retain the above copyright
#   notice, this list of conditions and the following disclaimer.
# * Redistributions in binary form must reproduce the above copyright
#   notice, this list of conditions and the following disclaimer in the
#   documentation and/or other materials provided with the distribution.
# * Neither the name of CEA nor the names of its contributors may be used to
#   endorse or promote products derived from this software without specific
#   prior written permission.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLEq FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#=============================================================================

cmake_minimum_required(VERSION 3.22...4.2)
project(pdi_conduit_plugin LANGUAGES C CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Includes
include(CTest)
include(GNUInstallDirs)

# set(conduit_DIR "../../vendor/conduit/install-debug/lib/cmake/conduit/")
set(RUNTEST_DIR "${CMAKE_CURRENT_LIST_DIR}/../../cmake/runtest-dir")

# If Conduit is in a non-standard path, run cmake with:
# -Dconduit_DIR=/path/to/conduit/lib/cmake/conduit
find_package(conduit REQUIRED)
# find_package(PDI REQUIRED COMPONENTS plugins)

add_executable(conduit_demo conduit.cxx)
target_link_libraries(conduit_demo PRIVATE conduit::conduit)

add_executable(conduit_example example.cxx)
target_link_libraries(conduit_example conduit::conduit m)

file(COPY ${CMAKE_CURRENT_LIST_DIR}/config.yaml DESTINATION ${CMAKE_BINARY_DIR}/conduit/)
message(" ${CMAKE_BINARY_DIR}")

################################################################################
#
#
#    WashingtonDC Dreamcast Emulator
#    Copyright (C) 2016 snickerbockers
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#
################################################################################

cmake_minimum_required(VERSION 2.6)
project(WashingtonDC)
set(WashingtonDC_VERSION_MAJOR 0)
set(WashingtonDC_VERSION_MINOR 1)

option(ENABLE_SH4_ICACHE "emulate the sh4's instruction cache" ON)
option(ENABLE_SH4_OCACHE "emulate the sh4's operand cache" ON)
option(ENABLE_SH4_MMU "emulate the sh4's Memory Management Unit (MMU)" ON)

set(sh4_sources "${PROJECT_SOURCE_DIR}/src/hw/sh4/sh4.cpp"
                "${PROJECT_SOURCE_DIR}/src/hw/sh4/sh4_inst.cpp"
                "${PROJECT_SOURCE_DIR}/src/hw/sh4/sh4_mem.cpp"
                "${PROJECT_SOURCE_DIR}/src/hw/sh4/Memory.cpp")

if (ENABLE_SH4_ICACHE)
    add_definitions(-DENABLE_SH4_ICACHE)
    set(sh4_sources ${sh4_sources} "${PROJECT_SOURCE_DIR}/src/hw/sh4/Icache.cpp")
endif()

if (ENABLE_SH4_OCACHE)
    add_definitions(-DENABLE_SH4_OCACHE)
    set(sh4_sources ${sh4_sources} "${PROJECT_SOURCE_DIR}/src/hw/sh4/Ocache.cpp")
endif()

if (ENABLE_SH4_MMU)
    add_definitions(-DENABLE_SH4_MMU)
    set(sh4_sources ${sh4_sources} "${PROJECT_SOURCE_DIR}/src/hw/sh4/sh4_mmu.cpp")
endif()

set(sh4asm_sources "${PROJECT_SOURCE_DIR}/src/tool/sh4asm/Inst.cpp"
                   "${PROJECT_SOURCE_DIR}/src/tool/sh4asm/main.cpp"
                   "${PROJECT_SOURCE_DIR}/src/tool/sh4asm/sh4asm.cpp")

include_directories("${PROJECT_SOURCE_DIR}/src/"
                    "${PROJECT_SOURCE_DIR}/src/hw/sh4"
                    "${PROJECT_SOURCE_DIR}/src/common"
                    "${PROJECT_SOURCE_DIR}/external")

add_library(sh4asm ${sh4asm_sources})
add_library(sh4 ${sh4_sources})

add_executable(sh4asm_test "${PROJECT_SOURCE_DIR}/unit_tests/sh4asm_test.cpp")
target_link_libraries(sh4asm_test sh4asm)

add_executable(sh4_test "${PROJECT_SOURCE_DIR}/unit_tests/sh4_test.cpp")
target_link_libraries(sh4_test sh4)

add_executable(sh4inst_test "${PROJECT_SOURCE_DIR}/unit_tests/sh4inst_test.cpp")
target_link_libraries(sh4inst_test sh4 sh4asm)
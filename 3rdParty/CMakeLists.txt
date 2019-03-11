add_subdirectory("sleef" "sleef")



set(OBJ_CONV_SOURCES
        ${CMAKE_SOURCE_DIR}/3rdParty/objconv/cmdline.cpp
        ${CMAKE_SOURCE_DIR}/3rdParty/objconv/cmdline.h
        ${CMAKE_SOURCE_DIR}/3rdParty/objconv/cof2asm.cpp
        ${CMAKE_SOURCE_DIR}/3rdParty/objconv/cof2cof.cpp
        ${CMAKE_SOURCE_DIR}/3rdParty/objconv/cof2elf.cpp
        ${CMAKE_SOURCE_DIR}/3rdParty/objconv/cof2omf.cpp
        ${CMAKE_SOURCE_DIR}/3rdParty/objconv/coff.cpp
        ${CMAKE_SOURCE_DIR}/3rdParty/objconv/coff.h
        ${CMAKE_SOURCE_DIR}/3rdParty/objconv/containers.cpp
        ${CMAKE_SOURCE_DIR}/3rdParty/objconv/containers.h
        ${CMAKE_SOURCE_DIR}/3rdParty/objconv/converters.h
        ${CMAKE_SOURCE_DIR}/3rdParty/objconv/disasm.h
        ${CMAKE_SOURCE_DIR}/3rdParty/objconv/disasm1.cpp
        ${CMAKE_SOURCE_DIR}/3rdParty/objconv/disasm2.cpp
        ${CMAKE_SOURCE_DIR}/3rdParty/objconv/elf.cpp
        ${CMAKE_SOURCE_DIR}/3rdParty/objconv/elf.h
        ${CMAKE_SOURCE_DIR}/3rdParty/objconv/elf2asm.cpp
        ${CMAKE_SOURCE_DIR}/3rdParty/objconv/elf2cof.cpp
        ${CMAKE_SOURCE_DIR}/3rdParty/objconv/elf2elf.cpp
        ${CMAKE_SOURCE_DIR}/3rdParty/objconv/elf2mac.cpp
        ${CMAKE_SOURCE_DIR}/3rdParty/objconv/error.cpp
        ${CMAKE_SOURCE_DIR}/3rdParty/objconv/error.h
        ${CMAKE_SOURCE_DIR}/3rdParty/objconv/extras.zip
        ${CMAKE_SOURCE_DIR}/3rdParty/objconv/library.cpp
        ${CMAKE_SOURCE_DIR}/3rdParty/objconv/library.h
        ${CMAKE_SOURCE_DIR}/3rdParty/objconv/mac2asm.cpp
        ${CMAKE_SOURCE_DIR}/3rdParty/objconv/mac2elf.cpp
        ${CMAKE_SOURCE_DIR}/3rdParty/objconv/mac2mac.cpp
        ${CMAKE_SOURCE_DIR}/3rdParty/objconv/macho.cpp
        ${CMAKE_SOURCE_DIR}/3rdParty/objconv/macho.h
        ${CMAKE_SOURCE_DIR}/3rdParty/objconv/main.cpp
        ${CMAKE_SOURCE_DIR}/3rdParty/objconv/maindef.h
        ${CMAKE_SOURCE_DIR}/3rdParty/objconv/omf.cpp
        ${CMAKE_SOURCE_DIR}/3rdParty/objconv/omf.h
        ${CMAKE_SOURCE_DIR}/3rdParty/objconv/omf2asm.cpp
        ${CMAKE_SOURCE_DIR}/3rdParty/objconv/omf2cof.cpp
        ${CMAKE_SOURCE_DIR}/3rdParty/objconv/omfhash.cpp
        ${CMAKE_SOURCE_DIR}/3rdParty/objconv/opcodes.cpp
        ${CMAKE_SOURCE_DIR}/3rdParty/objconv/stdafx.cpp
        ${CMAKE_SOURCE_DIR}/3rdParty/objconv/stdafx.h)

add_executable( ${TARGET_OBJCONV} ${OBJ_CONV_SOURCES})

#
# enable nasm for asmlib
#

enable_language(ASM_NASM)

set(CMAKE_ASM_NASM_OBJECT_FORMAT elf64)
set(CMAKE_ASM_NASM_COMPILE_OBJECT "<CMAKE_ASM_NASM_COMPILER> <INCLUDES> \
    <FLAGS> -f ${CMAKE_ASM_NASM_OBJECT_FORMAT} -o <OBJECT> <SOURCE>")

set(ASM_NASM_FLAGS "${ASM_NASM_FLAGS} -DUNIX -Worphan-labels -I${CMAKE_SOURCE_DIR}/3rdParty/asmlib/")

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_ASM_NASM_FLAGS "${ASM_NASM_FLAGS} -g -Fdwarf")
else()
    set(CMAKE_ASM_NASM_FLAGS "${ASM_NASM_FLAGS}")
endif()

set(ASMLIB_SOURCES
        ${CMAKE_SOURCE_DIR}/3rdParty/asmlib/asm/instrset64.asm
        ${CMAKE_SOURCE_DIR}/3rdParty/asmlib/asm/procname64.asm
        ${CMAKE_SOURCE_DIR}/3rdParty/asmlib/asm/cpuid64.asm
        ${CMAKE_SOURCE_DIR}/3rdParty/asmlib/asm/rdtsc64.asm
        ${CMAKE_SOURCE_DIR}/3rdParty/asmlib/asm/round64.asm
        ${CMAKE_SOURCE_DIR}/3rdParty/asmlib/asm/memcpy64.asm
        ${CMAKE_SOURCE_DIR}/3rdParty/asmlib/asm/memmove64.asm
        ${CMAKE_SOURCE_DIR}/3rdParty/asmlib/asm/memset64.asm
        ${CMAKE_SOURCE_DIR}/3rdParty/asmlib/asm/memcmp64.asm
        ${CMAKE_SOURCE_DIR}/3rdParty/asmlib/asm/strlen64.asm
        ${CMAKE_SOURCE_DIR}/3rdParty/asmlib/asm/strcpy64.asm
        ${CMAKE_SOURCE_DIR}/3rdParty/asmlib/asm/strcat64.asm
        ${CMAKE_SOURCE_DIR}/3rdParty/asmlib/asm/strstr64.asm
        ${CMAKE_SOURCE_DIR}/3rdParty/asmlib/asm/strcmp64.asm
        ${CMAKE_SOURCE_DIR}/3rdParty/asmlib/asm/stricmp64.asm
        ${CMAKE_SOURCE_DIR}/3rdParty/asmlib/asm/strtouplow64.asm
        ${CMAKE_SOURCE_DIR}/3rdParty/asmlib/asm/substring64.asm
        ${CMAKE_SOURCE_DIR}/3rdParty/asmlib/asm/strspn64.asm
        ${CMAKE_SOURCE_DIR}/3rdParty/asmlib/asm/strcountutf864.asm
        ${CMAKE_SOURCE_DIR}/3rdParty/asmlib/asm/strcountset64.asm
        ${CMAKE_SOURCE_DIR}/3rdParty/asmlib/asm/divfixedi64.asm
        ${CMAKE_SOURCE_DIR}/3rdParty/asmlib/asm/divfixedv64.asm
        ${CMAKE_SOURCE_DIR}/3rdParty/asmlib/asm/popcount64.asm
        ${CMAKE_SOURCE_DIR}/3rdParty/asmlib/asm/physseed64.asm
        ${CMAKE_SOURCE_DIR}/3rdParty/asmlib/asm/mother64.asm
        ${CMAKE_SOURCE_DIR}/3rdParty/asmlib/asm/mersenne64.asm
        ${CMAKE_SOURCE_DIR}/3rdParty/asmlib/asm/sfmt64.asm
        ${CMAKE_SOURCE_DIR}/3rdParty/asmlib/asm/cputype64.asm
        ${CMAKE_SOURCE_DIR}/3rdParty/asmlib/asm/debugbreak64.asm
        ${CMAKE_SOURCE_DIR}/3rdParty/asmlib/asm/unalignedisfaster64.asm
        ${CMAKE_SOURCE_DIR}/3rdParty/asmlib/asm/cachesize64.asm)

add_library( ${TARGET_ASMLIB} STATIC ${ASMLIB_SOURCES})


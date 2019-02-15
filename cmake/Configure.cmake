
include(CheckCCompilerFlag)
include(CheckCSourceCompiles)
include(CheckTypeSize)


set(PRICER_SUPPORTED_EXTENSIONS AVX AVX2 AVX512F)

# Allow to define the Gcc/Clang here
# As we might compile the lib with MSVC, but generates bitcode with CLANG
# Intel vector extensions.
set(CLANG_FLAGS_ENABLE_AVX "-mavx")
set(CLANG_FLAGS_ENABLE_AVX2 "-mavx2;-mfma")
set(CLANG_FLAGS_ENABLE_AVX512F "-mavx512f")

# All variables storing compiler flags should be prefixed with FLAGS_
if (CMAKE_C_COMPILER_ID MATCHES "(GNU|Clang)")
    # Always compile sleef with -ffp-contract.
    set(FLAGS_STRICTMATH "-ffp-contract=off")
    set(FLAGS_FASTMATH "-ffast-math")

    # Without the options below, gcc generates calls to libm
    set(FLAGS_NO_ERRNO "-fno-math-errno -fno-trapping-math")

    # Intel vector extensions.
    foreach (SIMD ${PRICER_SUPPORTED_EXTENSIONS})
        set(FLAGS_ENABLE_${SIMD} ${CLANG_FLAGS_ENABLE_${SIMD}})
    endforeach ()

    # Warning flags.
    set(FLAGS_WALL "-Wall -Wno-unused -Wno-attributes -Wno-unused-result")
    if (CMAKE_C_COMPILER_ID MATCHES "GNU")
        # The following compiler option is needed to suppress the warning
        # "AVX vector return without AVX enabled changes the ABI" at
        # src/arch/helpervecext.h:88
        string(CONCAT FLAGS_WALL ${FLAGS_WALL} " -Wno-psabi")
    endif (CMAKE_C_COMPILER_ID MATCHES "GNU")
elseif (MSVC)
    # Intel vector extensions.
    set(FLAGS_ENABLE_AVX /D__SSE2__ /D__SSE3__ /D__SSE4_1__ /D__AVX__ /arch:AVX)
    set(FLAGS_ENABLE_AVX2 /D__SSE2__ /D__SSE3__ /D__SSE4_1__ /D__AVX__ /D__AVX2__ /arch:AVX2)
    set(FLAGS_ENABLE_AVX512F /D__SSE2__ /D__SSE3__ /D__SSE4_1__ /D__AVX__ /D__AVX2__ /D__AVX512F__ /arch:AVX2)
    set(FLAGS_NO_ERRNO "")
elseif (CMAKE_C_COMPILER_ID MATCHES "Intel")
    set(FLAGS_ENABLE_AVX "-mavx")
    set(FLAGS_ENABLE_AVX2 "-march=core-avx2")
    set(FLAGS_ENABLE_AVX512F "-xCOMMON-AVX512")
    set(FLAGS_STRICTMATH "-fp-model strict -Qoption,cpp,--extended_float_type")
    set(FLAGS_FASTMATH "-fp-model fast=2 -Qoption,cpp,--extended_float_type")
    set(FLAGS_WALL "-fmax-errors=3 -Wall -Wno-unused -Wno-attributes")
    set(FLAGS_NO_ERRNO "")
endif ()


option(DISABLE_AVX "Disable AVX" OFF)

if (NOT DISABLE_AVX)
    set(CMAKE_REQUIRED_FLAGS ${FLAGS_ENABLE_AVX})
    CHECK_C_SOURCE_COMPILES("
  #if defined(_MSC_VER)
  #include <intrin.h>
  #else
  #include <x86intrin.h>
  #endif
  int main() {
    __m256d r = _mm256_add_pd(_mm256_set1_pd(1), _mm256_set1_pd(2));
  }" COMPILER_SUPPORTS_AVX)
endif ()


option(DISABLE_AVX2 "Disable AVX2" OFF)

if (NOT DISABLE_AVX2)
    set(CMAKE_REQUIRED_FLAGS ${FLAGS_ENABLE_AVX2})
    CHECK_C_SOURCE_COMPILES("
  #if defined(_MSC_VER)
  #include <intrin.h>
  #else
  #include <x86intrin.h>
  #endif
  int main() {
    __m256i r = _mm256_abs_epi32(_mm256_set1_epi32(1)); }"
            COMPILER_SUPPORTS_AVX2)
endif ()


option(DISABLE_AVX512F "Disable AVX512F" OFF)

if (NOT DISABLE_AVX512F)
    set(CMAKE_REQUIRED_FLAGS ${FLAGS_ENABLE_AVX512F})
    CHECK_C_SOURCE_COMPILES("
  #if defined(_MSC_VER)
  #include <intrin.h>
  #else
  #include <x86intrin.h>
  #endif
  __m512 addConstant(__m512 arg) {
    return _mm512_add_ps(arg, _mm512_set1_ps(1.f));
  }
  int main() {
    __m512i a = _mm512_set1_epi32(1);
    __m256i ymm = _mm512_extracti64x4_epi64(a, 0);
    __mmask16 m = _mm512_cmp_epi32_mask(a, a, _MM_CMPINT_EQ);
    __m512i r = _mm512_andnot_si512(a, a); }"
            COMPILER_SUPPORTS_AVX512F)
endif ()


# Reset used flags
set(CMAKE_REQUIRED_FLAGS)
set(CMAKE_REQUIRED_LIBRARIES)


# keep the user informed

if (COMPILER_SUPPORTS_AVX)
    message(STATUS "AVX supported")
endif ()

if (COMPILER_SUPPORTS_AVX2)
    message(STATUS "AVX2 supported")
endif ()

if (COMPILER_SUPPORTS_AVX512F)
    message(STATUS "AVX512F supported")
endif ()

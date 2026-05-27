cmake_minimum_required(VERSION 3.10)
project(ZKP_Hardware_Bench)

# C++14 표준 사용 (libsnark 필수)
set(CMAKE_CXX_STANDARD 14)

# 🌟 [여기 추가됨!] 구형 C++ 문법 경고 무시 및 깐깐한 에러 검사 끄기 🌟
add_compile_options(-Wno-deprecated-declarations -Wno-error)

# 🚨 ARM 보드용 필수 최적화 및 에러 방지 플래그 🚨
set(CMAKE_BUILD_TYPE Release)
set(USE_ASM OFF CACHE BOOL "Disable x86 assembly for ARM" FORCE)
set(MULTICORE OFF CACHE BOOL "Disable multicore" FORCE)
set(WITH_PROCPS OFF CACHE BOOL "Disable memory profiling" FORCE)

# libsnark 하위 폴더를 라이브러리로 인식
add_subdirectory(libsnark)

# 내 코드를 'benchmark'라는 이름의 실행 파일로 만들기
add_executable(benchmark benchmark_sha256_groth16.cpp)

# curve
target_compile_definitions(benchmark PRIVATE CURVE_ALT_BN128)

# 컴파일러에게 헤더 파일(.hpp) 위치를 정확히 알려주기
target_include_directories(benchmark PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/libsnark
    ${CMAKE_CURRENT_SOURCE_DIR}/libsnark/depends/libff
)

# 실행 파일에 필요한 라이브러리들 연결
target_link_libraries(benchmark snark ff gmp)
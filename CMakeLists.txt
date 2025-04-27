# CMake 最低版本要求
cmake_minimum_required(VERSION 3.10)

# 设置项目名称和编程语言
project(xxn_hw1 C)  # 项目名称为 xxn_hw1，语言为 C

# 设置 C 标准（例如 C11）
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)  # 强制使用指定标准


# 包含头文件目录（如果头文件在项目根目录下的 inc 文件夹中）
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/inc)
# 添加可执行文件
add_executable(${PROJECT_NAME}
    ${CMAKE_CURRENT_SOURCE_DIR}/src/main.c
    ${CMAKE_CURRENT_SOURCE_DIR}/src/algebra.c
)
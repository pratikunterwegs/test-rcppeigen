find src -name "*.cpp" ! -name "RcppExports.cpp" -exec clang-format -style=google -i {} \;

clang-cl -fuse-ld=lld -flto=thin /D "NDEBUG" /D "_CONSOLE" /D "NOMINMAX" /D "_UNICODE" /D "UNICODE" -Xclang -fcxx-exceptions -Xclang -std=c++17 -Qunused-arguments -Xclang -ffast-math -Xclang -ffunction-sections -Xclang -Wno-deprecated-declarations -Xclang -Wno-unknown-pragmas -Xclang -Wno-ignored-pragmas -Xclang -Wno-unused-private-field -Xclang -Wno-inconsistent-dllimport -mmmx  -msse  -msse2 -msse3 -mssse3 -msse4.1 -msse4.2 -mavx -mavx2  -Xclang -Wno-unused-variable -Xclang -Wno-language-extension-token -I"z:\vc\x64\include" -Ox -MT integer_utils.cpp Source.cpp -link "z:\vc\x64\lib\libboost_random-clang60-mt-s-x64-1_66.lib" "z:\vc\x64\lib\libboost_system-clang60-mt-s-x64-1_66.lib" "kernel32.lib" "user32.lib" "winspool.lib" "comdlg32.lib" "advapi32.lib" "shell32.lib" "ole32.lib" "oleaut32.lib" "uuid.lib" "odbc32.lib" "odbccp32.lib"

:: -fuse-ld=lld -flto=thin  /LIBPATH:"z:\vc\x64\lib"

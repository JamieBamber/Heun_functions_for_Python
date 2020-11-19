#!/bin/bash
rm HeunC_lib.so HeunC.o 
g++ -shared -c -fPIC HeunC.cpp -o HeunC.o -std=c++14
g++ -shared -Wl,-soname,HeunC_lib.so -o HeunC_lib.so HeunC.o 
echo "made HeunC_lib.so"

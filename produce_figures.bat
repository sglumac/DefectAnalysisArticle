set REPO_DIR=%~dp0

del /q /s build
del /q /s out
mkdir build
mkdir out
cd build
cmake -D CMAKE_INSTALL_PREFIX="%REPO_DIR%/out" -D SUNDIALS="%REPO_DIR%/sundials" "%REPO_DIR%"
cmake --build . --target install
cd ..

python scripts\results_analysis.py

echo "Figures produced to the .\out folder..."
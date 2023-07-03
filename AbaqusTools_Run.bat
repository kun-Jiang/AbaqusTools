@echo off
setlocal

rem Get the path of current folder
set "SCRIPT_DIR=%~dp0"

rem Combine the folder path with python file's name
set "PYTHON_FILE=%SCRIPT_DIR%AbaqusTools.py"

rem Run the python program
python "%PYTHON_FILE%"

endlocal

@echo off
setlocal

rem 获取批处理脚本文件所在目录的路径
set "SCRIPT_DIR=%~dp0"

rem 在Python程序文件名前添加批处理脚本文件所在目录的路径
set "PYTHON_FILE=%SCRIPT_DIR%AbaqusTools.py"

rem 运行Python程序
python "%PYTHON_FILE%"

endlocal

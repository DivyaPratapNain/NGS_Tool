@echo off
REM Convert Windows path to WSL path
set WIN_PATH=%~dp0install_deps.py
for /f "usebackq tokens=*" %%i in (`wsl wslpath -a "%WIN_PATH%"`) do set WSL_PATH=%%i

REM Replace problematic characters (spaces, parentheses)
set ESCAPED_PATH=%WSL_PATH: =\ %
set ESCAPED_PATH=%ESCAPED_PATH:(=\(%
set ESCAPED_PATH=%ESCAPED_PATH:)=\)%

REM Run the Python script in WSL
wsl bash -c "python3 %ESCAPED_PATH%"
pause

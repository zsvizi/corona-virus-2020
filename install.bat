@echo off
title Install ImPro

echo ---------------------
echo ----   INSTALL   ----
echo ---------------------
echo.

echo Checking python version...
rem This line is to force everyone to use a similar version of Python
py --version | findstr /i "Python 3.7"
if %errorlevel% neq 0 (
    echo ERROR! Python 3.7.x is required
    exit /B %errorlevel%
)

if exist "venv" (
    echo Reinstalling environment...
    rmdir "venv" /S /Q
) else (
    echo Installing environment...
)
py -m pip install --user virtualenv
py -m venv venv

echo Activating environment...
call .\venv\Scripts\activate.bat

echo Installing packages...
py -m pip install --upgrade pip
pip install -r requirements.txt
if %errorlevel% neq 0 (
    echo ERROR! Something went wrong
    exit /B %errorlevel%
)

echo Successfully installed!
echo.

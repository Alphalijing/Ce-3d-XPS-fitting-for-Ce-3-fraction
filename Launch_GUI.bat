@echo off
chcp 65001 >nul
setlocal
pushd "%~dp0"
if exist "dist\Ce3d_XPS_Fitter_v2.0.1_EN\Ce3d_XPS_Fitter_v2.0.1_EN.exe" goto packaged
where py >nul 2>nul
if not errorlevel 1 goto python_launcher
where python >nul 2>nul
if not errorlevel 1 goto python_command
echo Python was not found. Install Python and the dependencies in requirements.txt, or use the packaged application folder.
set "ce3d_exit=1"
goto failed
:packaged
start "" "dist\Ce3d_XPS_Fitter_v2.0.1_EN\Ce3d_XPS_Fitter_v2.0.1_EN.exe"
set "ce3d_exit=%errorlevel%"
goto finish
:python_launcher
py -3 "ce3d_xps_fitter_gui.py"
set "ce3d_exit=%errorlevel%"
if not "%ce3d_exit%"=="0" goto failed
goto finish
:python_command
python "ce3d_xps_fitter_gui.py"
set "ce3d_exit=%errorlevel%"
if not "%ce3d_exit%"=="0" goto failed
goto finish
:failed
echo.
echo The application could not be started. Review the messages above. Install dependencies with python -m pip install -r requirements.txt.
pause
:finish
popd
exit /b %ce3d_exit%

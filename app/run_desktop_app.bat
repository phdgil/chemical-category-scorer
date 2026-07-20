@echo off
setlocal EnableExtensions
cd /d "%~dp0"

if exist "%~dp0launch_desktop_app.vbs" (
    start "" wscript.exe "%~dp0launch_desktop_app.vbs"
    exit /b 0
)

start "" pyw -3 "%~dp0launch_desktop_app.pyw"
exit /b 0

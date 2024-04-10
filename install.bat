@echo off
:: Check if WSL is available on the system
where wsl >nul 2>&1

if %ERRORLEVEL% neq 0 (
    echo WSL is not installed on this system. Please go to https://learn.microsoft.com/en-us/windows/wsl/install-manual#step-1---enable-the-windows-subsystem-for-linux to install 
    pause
    exit /b
)

:: Check if Ubuntu is already installed
wsl -l | findstr /C:"Ubuntu" >nul 2>&1
if %ERRORLEVEL% equ 0 (
    echo Ubuntu is already installed.
    goto skipUbuntuInstall
)

:: Set default WSL version to 2
wsl --set-default-version 2
if %ERRORLEVEL% neq 0 (
    echo Failed to set WSL default version to 2.
    pause
    exit /b
)

:: Install Ubuntu in WSL
wsl --install -d Ubuntu
if %ERRORLEVEL% neq 0 (
    echo Failed to install Ubuntu on WSL.
    pause
    exit /b
)

:skipUbuntuInstall
echo Continuing with the script...
pause
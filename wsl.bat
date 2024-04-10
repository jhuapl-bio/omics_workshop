@echo off
:: Check if WSL is installed
where wsl >nul 2>&1
if %ERRORLEVEL% neq 0 (
    echo WSL is not installed on this system. Please visit the link below to install:
    echo https://learn.microsoft.com/en-us/windows/wsl/install
    pause
    exit /b
)
:: Change directory to where the script is located if needed
:: cd /path/to/your/script

:: Execute the bash script in WSL
echo Running install.sh in WSL...
@REM wsl.exe -u root apt update
wsl.exe bash -c "echo yes"
wsl.exe bash -c "bash install.sh"



if %ERRORLEVEL% neq 0 (
    echo There was an error executing install.sh.
    pause
    exit /b
)

echo Bash script execution completed successfully.
pause
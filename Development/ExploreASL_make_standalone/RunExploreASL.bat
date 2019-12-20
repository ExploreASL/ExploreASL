REM # Some administration

@echo off
cls
set MCR_ROOT=cd;
dir /b ExploreASL*.exe > temp.temp 
REM # this puts the ExploreASL*.exe in temp.temp
REM # option /b avoids the header & footer of "dir"
set /p xASLpath=< temp.temp
REM # this puts the ExploreASL path in variable
del temp.temp
REM # this removes the temporary file

REM # Runs ExploreASL executable
%xASLpath%

REM # Remove temporary folder
dir /b ExploreASL*mcr > temp.temp 
set /p mcrPath=< temp.temp
del temp.temp
rmdir /s /q %mcrPath%
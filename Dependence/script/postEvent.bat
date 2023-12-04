set SolutionDir=%1
set SolutionName=%2
set ProjectDir=%3
set ProjectName=%4
set Platform=%5
set Configuration=%6
set TargetExt=%7


echo %ProjectName%|findstr "^Plugin">nul
if %errorlevel% equ 0 (
xcopy "%SolutionDir%\lib\%Platform%\%Configuration%\%ProjectName%%TargetExt%" "%SolutionDir%\bin\%Platform%\%Configuration%\plugins\" /Y
if "%Configuration%"=="Debug" (
	xcopy "%SolutionDir%\lib\%Platform%\%Configuration%\%ProjectName%.pdb" "%SolutionDir%\bin\%Platform%\%Configuration%\plugins\" /Y
)
) else (
	xcopy "%SolutionDir%\lib\%Platform%\%Configuration%\%ProjectName%%TargetExt%" "%SolutionDir%\bin\%Platform%\%Configuration%\" /Y
	if "%Configuration%"=="Debug" (
	xcopy "%SolutionDir%\lib\%Platform%\%Configuration%\%ProjectName%.pdb" "%SolutionDir%\bin\%Platform%\%Configuration%\" /Y
)
)

pause

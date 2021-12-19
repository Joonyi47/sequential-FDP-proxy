@REM+- OmniWorks Replacement History - gss_macros`pc:$convert.bat;4 
@REM       4*[100705] 20-MAR-2007 16:59:06 (GMT) thk 
@REM         "2006.3: Move convert to bin" 
@REM       3*[100704] 20-MAR-2007 16:33:57 (GMT) psandhu 
@REM         "2006.3 add lib path to convert.bat" 
@REM     2B1 [ 43465] 21-MAR-2003 14:32:54 (GMT) thk 
@REM         "2003a-2004a CODE SPLIT" 
@REM     2A1 [ 28765] 17-MAY-2002 15:33:17 (GMT) thk 
@REM         "2002a-2003a CODE SPLIT" 
@REM       2*[ 16252] 08-AUG-2001 08:40:28 (GMT) axd 
@REM         "2002a change PC macros so they can run execs in both the new and old directory structures" 
@REM     1A1 [ 10523] 06-MAR-2001 11:15:06 (GMT) jmb 
@REM         "2001a-2002a CODE SPLIT" 
@REM       1*[  2961] 27-JUN-2000 16:40:50 (GMT) SuperUser 
@REM         "2001a" 
@REM+- OmniWorks Replacement History - gss_macros`pc:$convert.bat;4 
@echo off
CALL $ECLRC.BAT
REM Handle any command line arguments
REM
CALL $ecltidy
set $$ecltar=
:NewArg
if "%1"=="" goto NoArgs
set $$ecltar=%$$ecltar% %1%
shift
goto NewArg
:NoArgs
CALL $ChkArgs %$$ecltar%
set $$ecltar=

REM Initailse key variables
REM
CALL $eclrc
if not exist %ECLHOME%\$eclrc.BAT SET ECLERR=TRUE
if NOT "%ECLERR%"=="TRUE" goto Continue
echo *
echo * ERROR - The macro $ECLRC.BAT is either set incorrectly, or can
echo *         not be found.  Please adjust $ECLRC.BAT.
echo *
goto End

:Continue
REM Check for valid versions
REM
:ChkVer
CALL $ChkVer CONVERT
if "%ECLERR%"=="TRUE" goto End

REM If (valid) data directory has been given, change to it
if "%ecl_data%"=="" goto NoData
if exist %ecl_data%\nul cd %ecl_data%
:NoData

REM Display banner
cls
echo *
echo * CONVERT VERSION %ecl_ver%
echo *

:Run
REM Check Config
CALL $chkcfg
set TEMP_PATH=%PATH%
set PATH=%ECLARCH%\%ecl_ver%\lib\%ecl_bindir%;%PATH%

%ecl_exe%
set PATH=%TEMP_PATH%
set TEMP_PATH=

if "%$$eclcfg%"=="TRUE" del ECL.CFG
set $$eclcfg=

:Tidy
if "%$$eclcfg%"=="TRUE" del ECL.CFG
set $$eclcfg=
if exist CONFIG del CONFIG
CALL $ecltidy
:End
CALL $ecltidy
:End

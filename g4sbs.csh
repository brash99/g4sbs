#!/bin/csh

#script to set up the environment for g4sbs:
setenv G4SBS /Users/brash/g4sbs

if( ! ${?PATH} ) then
    setenv PATH /Users/brash/g4sbs/bin
else
    setenv PATH /Users/brash/g4sbs/bin:${PATH}
endif

set OS=`uname -s`


if( "$OS" == "Darwin" ) then # Mac OS: set DYLD_LIBRARY_PATH to library directory:
    if(! ${?DYLD_LIBRARY_PATH}) then
	setenv DYLD_LIBRARY_PATH /Users/brash/g4sbs/lib
    else
	setenv DYLD_LIBRARY_PATH /Users/brash/g4sbs/lib:${DYLD_LIBRARY_PATH}
    endif
endif

# set LD_LIBRARY_PATH regardless of OS:
if( ! ${?LD_LIBRARY_PATH}) then
    setenv LD_LIBRARY_PATH /Users/brash/g4sbs/lib
else
    setenv LD_LIBRARY_PATH /Users/brash/g4sbs/lib:${LD_LIBRARY_PATH}
endif



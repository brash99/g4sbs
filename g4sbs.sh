#!/bin/sh

#script to set up the environment for g4sbs:
export G4SBS=/Users/brash/g4sbs

if test "x$PATH" = "x" ; then
    export PATH=/Users/brash/g4sbs/bin
else
    export PATH=/Users/brash/g4sbs/bin:$PATH
fi

OS=`uname -s`


if [ "$OS" = "Darwin" ]
then # Mac OS: set DYLD_LIBRARY_PATH to library directory:
    if test "x$DYLD_LIBRARY_PATH" = "x"; then
	export DYLD_LIBRARY_PATH=/Users/brash/g4sbs/lib
    else
	export DYLD_LIBRARY_PATH=/Users/brash/g4sbs/lib:$DYLD_LIBRARY_PATH
    fi
fi

# set LD_LIBRARY_PATH regardless of OS:
if test "x$LD_LIBRARY_PATH" = "x"; then
    export LD_LIBRARY_PATH=/Users/brash/g4sbs/lib
else
    export LD_LIBRARY_PATH=/Users/brash/g4sbs/lib:$LD_LIBRARY_PATH
fi



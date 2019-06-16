#!/bin/csh
#Isaac Mooney, WSU, June 2019

#this script is capable of taking in multiple root files (currently necessary, actually),
#applies a generic macro (specified by an argument passed to this file) to each one
#and writes each output to a similarly named file, to be hadded afterward

# command line arguments
if ( $# < "1") then
    echo 'Error: illegal number of parameters'
    exit #need a filename, so have to have at least one argument
endif
cd macros #necessary because of the path hard-coded into the "runroot" file
runroot $1 || exit #compiles the macro, where $1 is the first argument given to the script (the name of the file without the extension)
cd .. #back to the top-level directory

set ExecPath = `pwd`

set uscore = "_" #for easy concatenation of variables later
set execute = "./macros/bin/$1" #or 'root macros/bin/$1.cxx+' on the command line for root compiling via ACLiC

set base = out/file #this happens to be the shared basename of the two files in the test - change if your input files change!
set analysisTag = root_macro #for actual physics code, this could be e.g. "QA", or "unfolded"

if ( ! -d out/$1) then
mkdir -p out/$1 #making output directory
endif

if ( ! -d log/$1) then
mkdir -p log/$1 #making log directory
endif
 
foreach input ( ${base}* )

set OutBase = `basename $input | sed 's/.root//g'` #this removes the filetype so we can append to the filename
set OutBase = `basename $OutBase | sed 's:out/::g'` #this removes the out/ from the input since it is unnecessary
echo $OutBase #should just be "file1", "file2", etc.
set OutBase = "$OutBase$uscore$analysisTag" #appending what we did to produce this file, e.g. "unfolded"

set outLocation = out/$1/ #producing a separate folder for output of different steps of the analysis
set outName = ${OutBase}.root #same as input name but with the analysisTag appended

set inFiles = ${input}

set LogFile = log/${analysisTag}/${OutBase}.log
set ErrFile = log/${analysisTag}/${OutBase}.err

echo "Logging output to " $LogFile
echo "Logging errors to " $ErrFile

set arg = "$outLocation $outName $inFiles" 

#./macros/bin/root_macro out/root_macro/ file1_root_macro.root out/file1.root <- this is a (currently) working execution line on which to model the qsub

echo qsub -V -q erhiq -l mem=4GB -o $LogFile -e $ErrFile -N $1 -- ${ExecPath}/submit/qwrap.sh ${ExecPath} $execute $arg

qsub -V -q erhiq -l mem=4GB -o $LogFile -e $ErrFile -N $1 -- ${ExecPath}/submit/qwrap.sh ${ExecPath} $execute $arg

end #end of loop over input files

#!/bin/csh
#Isaac Mooney, WSU - June 2019

#this script is capable of taking in multiple root files (currently necessary, actually),
#applies a generic macro (specified by an argument passed to this file) to each one
#and writes each output to a similarly named file, to be hadded afterward

# command line arguments
    #1: the executable; this also doubles as the "analysisTag", i.e. the generic action taken on the files, e.g. "unfolded", or "closure"
    #2: the input type, e.g. 'data', or 'QA'
    #3: the collision species: pp, pA, or AA
    #4: an analysis-specific wildcard (not required). Currently being used in data to take either full or charged jets.

if ( $# < "1") then
    echo 'Error: illegal number of parameters'
    exit #need a filename, so have to have at least one argument
endif

#compile!
cd macros #necessary because of the path hard-coded into the "runroot" file
runroot $1 || exit #compiles the macro, where $1 is the first argument given to the script (the name of the file without the extension)
cd .. #back to the top-level directory

set base = ""

set analysisTag = $1 #this could be e.g. "unfolded" or "bin_drop"
set inputType = $2 #this could be e.g. "data" or "sim"
set species = $3 #i.e. pp, pA, AA
set tag = $4 #analysis-specific tag (currently for full v. charged jets)

#this block of conditionals will become more complicated when I have more macros than just hists.cxx. Will need to add in an "if $1 == x" capability 

echo 'Pulling files from out/$2!'
echo 'Species is $3!'
#setting input files based on command-line arguments
if (${species} == 'pp') then 
    set base = out/${inputType}/sum #Be careful! This also catches charged jets in its net, if the files exist. Need to think of how to distinguish...
else if (${species} == 'pA') then
    set base = out/${inputType}/pAu_2015_200_HT_ #Be careful! This also catches charged jets in its net, if the files exist. Need to think of how to distinguish...
else if (${species} == 'AA') then
    echo 'Error: AA is not ready yet! Be patient!'
    exit
else
    echo 'Error: unrecognized collision species!'
    exit
endif

set ExecPath = `pwd`

set uscore = "_" #for easy concatenation of variables later
set execute = "./macros/bin/$1" #or 'root macros/bin/$1.cxx+' on the command line for root compiling via ACLiC

if ( ! -d out/${inputType}/${analysisTag}) then #e.g. out/data/hists. Beautiful!
mkdir -p out/${inputType}/${analysisTag} #making output directory
endif

if ( ! -d log/${inputType}/${analysisTag}) then
mkdir -p log/${inputType}/${analysisTag} #making log directory
endif
 
foreach input ( ${base}*_${tag}* )

    set OutBase = `basename $input | sed 's/.root//g'` #this removes the filetype so we can append to the filename
    set OutBase = `basename $OutBase | sed 's:out/::g'` #this removes the out/ from the input since it is unnecessary
    echo $OutBase #should just be "sum1", "sum2", etc.
    set OutBase = "$OutBase$uscore$analysisTag" #appending what we did to produce this file, e.g. "unfolded"

    set outLocation = out/${inputType}/${analysisTag}/ #producing a separate folder for output of different steps of the analysis
    set outName = ${OutBase}.root

    set inFiles = ${input}

    set LogFile = log/${inputType}/${analysisTag}/${OutBase}.log
    set ErrFile = log/${inputType}/${analysisTag}/${OutBase}.err

    echo "Logging output to " $LogFile
    echo "Logging errors to " $ErrFile

    set arg = "$outLocation $outName $inFiles" 

    #./macros/bin/root_macro out/root_macro/ file1_root_macro.root out/file1.root <- this is a (currently) working execution line on which to model the qsub

    echo qsub -V -q erhiq -l mem=4GB -o $LogFile -e $ErrFile -N $1 -- ${ExecPath}/submit/qwrap.sh ${ExecPath} $execute $arg

    qsub -V -q erhiq -l mem=4GB -o $LogFile -e $ErrFile -N $1 -- ${ExecPath}/submit/qwrap.sh ${ExecPath} $execute $arg

end #end of loop over input files

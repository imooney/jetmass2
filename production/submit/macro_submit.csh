#!/bin/csh
#Isaac Mooney, WSU - June 2019

#this script is capable of taking in multiple root files,
#applies a generic macro (specified by an argument passed to this file) to each one
#and writes each output to a similarly named file, to be hadded afterward

# command line arguments
    #1: the executable; this also doubles as the "analysisTag", i.e. the generic action taken on the files, e.g. "hists"
    #2: the input type, e.g. 'pythia', or 'herwig'
    #3: the hadronic/partonic option given during production of the input files, i.e. undecayed, decayed, or partonicFS
    #4: the radius parameter for the jets, e.g. 04 for R = 0.4.

if ( $# < "4") then
    echo 'Error: illegal number of parameters'
    exit #need a filename, so have to have at least one argument
endif

#compile!
cd macros #necessary because of the path hard-coded into the "runroot" file
runroot $1 || exit #compiles the macro, where $1 is the first argument given to the script (the name of the file without the extension)
cd .. #back to the top-level directory

set base = ""

set analysisTag = $1 #this could be e.g. "hists"
set inputType = $2 #this could be e.g. "pythia" or "herwig"
set finalstate = $3 #i.e. undecayed, decayed, partonicFS
set radius = $4 #radius parameter

#this block of conditionals will become more complicated when I have more macros than just hists.cxx. Will need to add in an "if $1 == x" capability 

echo 'Pulling files from out/!'
echo 'MC is '$2
echo 'Final state is '$3
echo 'Jet radius is '$4

#setting input files based on command-line arguments
if (${inputType} == 'pythia') then 
    if (${analysisTag} == 'hists') then
	set base = out/star_mass_pythia8_IsaacsHepmcs_ #FSR_only_
    endif
    if (${analysisTag} == 'piKp_JEF') then
	set base = out/pythia8_
    endif
else if (${inputType} == 'herwig') then
    set base = out/star_mass_herwig7_IsaacsHepmcs_
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
 
foreach input ( ${base}*R${radius}_${finalstate}_HeavyIons.root )

    set OutBase = `basename $input | sed 's/.root//g'` #this removes the filetype so we can append to the filename
    set OutBase = `basename $OutBase | sed 's:out/::g'` #this removes the out/ from the input since it is unnecessary
    echo $OutBase #should just be "sum1", "sum2", etc.
#    set OutBase = "$OutBase$uscore$tag" #appending what we did to produce this file, e.g. "unfolded"

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

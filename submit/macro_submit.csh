#!/bin/csh
#Isaac Mooney, WSU - June 2019

#this script is capable of taking in multiple root files,
#applies a generic macro (specified by an argument passed to this file) to each one
#and writes each output to a similarly named file, to be hadded afterward

# command line arguments
    #1: the executable; this also doubles as the "analysisTag", i.e. the generic action taken on the files, e.g. "unfolded", or "closure"
    #2: the input type, e.g. 'data', 'QA', 'sim'
    #3: the collision species: pp, pA, or AA
    #4: the jet radius. Options: any number > 0 & <= 9.9; passed in without the decimal, e.g. "06" for 0.6.        
    #5: an analysis-specific wildcard (not required)

#something to think about: for e.g. bin_drop or stat_err_scaling, now that we have jet radius as a knob, it would be nice to pass this in to this code as well. Is this doable?

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
set radius = $4 #e.g. 04 for 0.4
set tag = $5 #analysis-specific tag (currently for full v. charged jets in data or matched v. unmatched in simulation)

echo 'Pulling files from out/'${inputType}'!'
echo 'Species is '${species}'!'
#setting input files based on command-line arguments
if (${species} == 'pp') then 
    if (${inputType} == 'data') then
	#set base = out/${inputType}/pp12Pico_pass #this is for HT! Comment this out and uncomment the next line for ppJP2 instead!
	set base = out/${inputType}/sum #Be careful! This also catches charged jets in its net, if the files exist. But later, using $tag makes it work fine.
    else if (${inputType} == 'sim') then
	set base = out/${inputType}/Cleanpp12Pico
    endif
else if (${species} == 'pA') then
    set base = out/${inputType}/pAu_2015_200_HT_ #Be careful! This also catches charged jets in its net, if the files exist. Need to think of how to distinguish...
else if (${species} == 'AA') then
    echo 'Error: AA is not ready yet! Be patient!'
    exit
else
    if (${analysisTag} != 'bin_drop' && ${analysisTag} != 'stat_err_scaling' && ${analysisTag} != 'unfold') then
	echo 'Error: unrecognized collision species!'
	exit
    endif
endif

set ExecPath = `pwd`

set uscore = "_" #for easy concatenation of variables later
set execute = "./macros/bin/${analysisTag}" #or 'root macros/bin/${analysisTag}.cxx+' on the command line for root compiling via ACLiC

if (${analysisTag} == 'bin_drop' || ${analysisTag} == 'stat_err_scaling' || ${analysisTag} == 'unfold') then

    if (! -d log/${analysisTag}) then
	mkdir -p log/${analysisTag}
    endif
    
    set LogFile = log/${analysisTag}/${analysisTag}.log
    set ErrFile = log/${analysisTag}/${analysisTag}.err	    
    
    echo qsub -V -q erhiq -l mem=4GB -o $LogFile -e $ErrFile -N ${analysisTag} -- ${ExecPath}/submit/qwrap.sh ${ExecPath} $execute
    qsub -V -q erhiq -l mem=4GB -o $LogFile -e $ErrFile -N ${analysisTag} -- ${ExecPath}/submit/qwrap.sh ${ExecPath} $execute

else #anything but procedures which don't need any args or input/output help
    
    if ( ! -d out/${inputType}/${analysisTag}) then #e.g. out/data/hists. Beautiful!
	mkdir -p out/${inputType}/${analysisTag} #making output directory
    endif

    if ( ! -d log/${inputType}/${analysisTag}) then
	mkdir -p log/${inputType}/${analysisTag} #making log directory
    endif

    foreach input ( ${base}*_${tag}_R${radius}* )
    
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
	
	set arg = "$outLocation $outName $inFiles $inputType" 
	
	#./macros/bin/root_macro out/root_macro/ file1_root_macro.root out/file1.root <- this is a (currently) working execution line on which to model the qsub
	
	echo qsub -V -q erhiq -l mem=4GB -o $LogFile -e $ErrFile -N ${analysisTag} -- ${ExecPath}/submit/qwrap.sh ${ExecPath} $execute $arg
	
	qsub -V -q erhiq -l mem=4GB -o $LogFile -e $ErrFile -N ${analysisTag} -- ${ExecPath}/submit/qwrap.sh ${ExecPath} $execute $arg
       
    end #end of loop over input files

endif

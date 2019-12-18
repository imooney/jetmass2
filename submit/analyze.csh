#!/bin/csh
#Isaac Mooney, WSU - June 2019
#This script will expand as I add in more of the analyses: pp, pA, AA
#It currently handles the QA and data (for all three systems). Next it will be expanded to handle the Pythia6 and Pythia6+Geant simulation on disk.

#args:
#1: the analysis type. Current options: QA, data, sim
#2: the trigger. Current options for pp: JP2, HT2, VPDMB(-nobsmd); for pA: JP2, HT2, BBCMB. For MC analysis, only JP2 is allowed at the moment.
#3: the species. Current options: pp, pA, AA
#4: the jet radius. Options: any number > 0 & <= 9.9; passed in without the decimal, e.g. "06" for 0.6.
#5: an analysis-specific wildcard (not required). Currently being used in data to take either full (full) or charged (ch) jets.
#6: second analysis-specific wildcard (not required). Currently being used in sim to either require matching (matched) - for responses - or not (unmatched) - for normal spectra

# Arguments       
if ( $# < "4") then
    echo 'Error: illegal number of parameters'
    exit
endif
    
set ExecPath = `pwd`
set arg = '' 

set analysisType = $1
set trigger = $2
set species = $3
set radius = $4
set wc1 = $5
set wc2 = $6

#~~~input checking and directory / file setup~~~#

if (${analysisType} == 'QA') then
    make bin/QA || exit
    set execute = './bin/QA'
    # Create the folder name for output
    set outFile = QA
    echo 'Running in QA mode!'
else if (${analysisType} == 'data') then
    make bin/data || exit
    set execute = './bin/data'
    # Create the folder name for output
    set outFile = data
    echo 'Running in data mode!'
else if (${analysisType} == 'sim') then
    make bin/sim || exit
    set execute = './bin/sim'
    # Create the folder name for output
    set outFile = sim
    echo 'Running in sim mode!'
else
    echo 'Error: unrecognized analysis type option'
    exit
endif

if (${analysisType} != 'sim') then
    if (${trigger} == 'JP2' && ${species} == 'pp') then
	set trigger = "ppJP2"
	set base = /tier2/home/groups/rhi/STAR/Data/ppJP2Run12/sum	
	echo "Running on the ppJP2-triggered data!"
    else if (${trigger} == 'HT2' && ${species} == 'pp') then
	set trigger = "ppHT2"
	set base = /tier2/home/groups/rhi/STAR/Data/ppHT2Run12/pp12Pico
	echo "Running on the ppHT2-triggered data!"
    else if (${trigger} == 'VPDMB' && ${species} == 'pp') then
	set trigger = "ppVPDMB"
	set base = /tier2/home/groups/rhi/STAR/Data/ppMBRun12/sum
	echo "Running on the ppMB-triggered data!"
    else if (${trigger} == 'JP2' && ${species} == 'pA') then
	set trigger = "pAuJP2"
	set base = /tier2/home/groups/rhi/STAR/Data/P16id/production_pAu200_2015/HT/pAu_2015_200_HT_1
	echo "Running on the pAuJP2-triggered data!"
    else if (${trigger} == 'HT2' && ${species} == 'pA') then
	set trigger = "pAuHT2"
	set base = /tier2/home/groups/rhi/STAR/Data/P16id/production_pAu200_2015/HT/pAu_2015_200_HT_1
	echo "Running on the pAuHT2-triggered data!"	
    else if (${trigger} == 'BBCMB' && ${species} == 'pA') then
	set trigger = "pAuBBCMB"
	set base = /tier2/home/groups/rhi/STAR/Data/P16id/production_pAu200_2015/MB/pAu_2015_200_MB_1
	echo "Running on the pAuBBCMB-triggered data!"
    else if (${species} == 'AA' || ${species} == 'AuAu') then
	echo "AA is not ready yet - be patient!"
    else
	echo "Error: unrecognized trigger and/or species option!"
	exit
    endif
else
    set base = /tier2/home/groups/rhi/STAR/Data/AddedEmbedPythiaRun12pp200/Cleanpp12Pico_
    echo "Running on the Pythia6 and Pythia6+Geant (JP2-triggered) simulation!"
endif

# Make output directories if they don't already exist            
if ( ! -d out/${outFile} ) then
mkdir -p out/${outFile}
endif

if ( ! -d log/${outFile} ) then
mkdir -p log/${outFile}
endif

echo ${base} #for debugging

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Now Submit jobs for each data file
foreach input ( ${base}* )

    #synthesize the output file base name from the input file base name and the options passed
    set OutBase = `basename $input | sed 's/.root//g'`
    set uscore = "_" #useful for chaining variables together
    if (${analysisType} != 'sim') then
	set OutBase = "$OutBase$uscore$trigger$uscore${wc1}${uscore}R${radius}"
    else
	set OutBase = "$OutBase$uscore${wc1}$uscore${wc2}${uscore}R${radius}"
    endif
    #make the output path and names
    set outLocation = out/${outFile}/
    set outName = ${OutBase}.root

    #input files
    set Files = ${input}
    
    #log files
    set LogFile = log/${outFile}/${OutBase}.log
    set ErrFile = log/${outFile}/${OutBase}.err

    echo "Logging output to " $LogFile
    echo "Logging errors to " $ErrFile

    set dummy = 'dummy' #this extra argument is for compatibility between QA and the other segments of the analysis (or for future expansion)

    #setting command-line args based on the analysis file being used:
    if (${analysisType} == 'QA') then
	set arg = "$outLocation $outName $trigger $dummy $Files"
    else if (${analysisType} == 'data') then
	set arg = "$outLocation $outName ${radius} $trigger ${wc1} $Files"
    else if (${analysisType} == 'sim') then
	set arg = "$outLocation $outName ${radius} ${wc1} ${wc2} $Files"
    endif

    echo "now submitting this script: "
    echo qsub -V -q erhiq -l mem=4GB -o $LogFile -e $ErrFile -N ${analysisType} -- ${ExecPath}/submit/qwrap.sh ${ExecPath} $execute $arg

    qsub -V -q erhiq -l mem=4GB -o $LogFile -e $ErrFile -N ${analysisType} -- ${ExecPath}/submit/qwrap.sh ${ExecPath} $execute $arg   
    #add back in a second: -q erhiq

end

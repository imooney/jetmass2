#!/bin/csh
#Isaac Mooney, WSU - June 2019
#This script will expand as I add in more of the analyses: pp, pA, AA
#It currently handles just the QA (for all three systems). Next it will be expanded to handle the data.

#args:
#1: the analysis type. Current options: QA, data
#2: the trigger. Current options for pp: JP2, VPDMB(-nobsmd); for pA: JP2, BBCMB
#3: the species. Current options: pp, pA, AA
#4: an analysis-specific wildcard (not required). Currently being used in data to take either full or charged jets.

set ExecPath = `pwd`
set arg = '' 

#will hopefully extend this later to be capable of running QA for any system: pp, pA, AA [shouldn't be too hard - "if $3 == 'pp' set base = ..." etc.]
if ($1 == 'QA') then
    make bin/QA || exit
    set execute = './bin/QA'
    # Create the folder name for output
    set outFile = QA
    echo 'Running in QA mode!'
else if ($1 == 'data') then
    make bin/data || exit
    set execute = './bin/data'
    # Create the folder name for output
    set outFile = data
    echo 'Running in data mode!'
else
    echo 'Error: unrecognized analysis type option'
    exit
endif

if ($2 == 'JP2' && $3 == 'pp') then
    set trigger = "ppJP2"
    set base = /nfs/rhi/STAR/Data/ppJP2Run12/sum	
    echo "Running on the ppJP2-triggered data!"
else if ($2 == 'VPDMB' && $3 == 'pp') then #NOTE: for data mode, a request for anything but JP2 for pp is ignored!!!
    set trigger = "ppVPDMB"
    set base = /nfs/rhi/STAR/Data/ppMBRun12/sum
    echo "Running on the ppMB-triggered data, if this is QA mode! Otherwise, reverting back to ppJP2!"
else if ($2 == 'JP2' && $3 == 'pA') then
    set trigger = "pAuJP2"
    set base = /nfs/rhi/STAR/Data/P16id/production_pAu200_2015/HT/pAu_2015_200_HT_1
    echo "Running on the pAuJP2-triggered data!"
else if ($2 == 'BBCMB' && $3 == 'pA') then #NOTE: for data mode, a request for anything but JP2 for pA is ignored!!!
    set trigger = "pAuBBCMB"
    set base = /nfs/rhi/STAR/Data/P16id/production_pAu200_2015/MB/pAu_2015_200_MB_1
    echo "Running on the pAuBBCMB-triggered data, if this is QA mode! Otherwise, reverting back to pAuJP2!"
else if ($3 == 'AA' || $3 == 'AuAu') then
    echo "AA is not ready yet - be patient!"
else
    echo "Error: unrecognized trigger and/or species option!"
    exit
endif

# Arguments       
if ( $# < "3") then
    echo 'Error: illegal number of parameters'
    exit
endif

# Make output directories if they don't already exist            
if ( ! -d out/${outFile} ) then
mkdir -p out/${outFile}
endif

if ( ! -d log/${outFile} ) then
mkdir -p log/${outFile}
endif

echo ${base} #for debugging

# Now Submit jobs for each data file
foreach input ( ${base}* )

    #synthesize the output file base name from the input file base name and the options passed
    set OutBase = `basename $input | sed 's/.root//g'`
    set uscore = "_" #useful for chaining variables together
    set OutBase = "$OutBase$uscore$trigger"
    
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
    if ($1 == 'QA') then
	set arg = "$outLocation $outName $trigger $dummy $Files"
    else if ($1 == 'data') then
	set arg = "$outLocation $outName $3 $4 $Files" #currently the only place that argument 4 shows up
    endif

    echo "now submitting this script: "
    echo qsub -V -q erhiq -l mem=4GB -o $LogFile -e $ErrFile -N $1 -- ${ExecPath}/submit/qwrap.sh ${ExecPath} $execute $arg

    qsub -V -q erhiq -l mem=4GB -o $LogFile -e $ErrFile -N $1 -- ${ExecPath}/submit/qwrap.sh ${ExecPath} $execute $arg   
    #add back in a second: -q erhiq

end

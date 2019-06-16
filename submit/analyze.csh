#!/bin/csh
#Isaac Mooney, WSU - June 2019
#This script will expand as I add in more of the analyses: pp, pA, AA
#It currently handles just the QA for pp. The next step is for it to handle pA QA:
    #a "system" string can be added as $3. This will delineate what "base" should be for file input and will be chained with $2 to set the trigger as well.
    #to test, though, will need also to work on the QA.cxx file to make it inclusive of the two.

set ExecPath = `pwd`
set arg = '' 

#will hopefully extend this later to be capable of running QA for any system: pp, pA, AA [shouldn't be too hard - "if $3 == 'pp' set base = ..." etc.]
if ($1 == 'QA') then
    make bin/QA || exit
    set execute = './bin/QA'
    # Create the folder name for output
    set outFile = QA
    echo 'Running in QA mode!'
else
    echo 'Error: unrecognized analysis type option'
    exit
endif

if ($2 == 'JP2' && $3 == 'pp') then
    set trigger = "ppJP2"
    set base = /nfs/rhi/STAR/Data/ppJP2Run12/sum	
    echo "Running on the ppJP2-triggered data!"
else if ($2 == 'VPDMB' && $3 == 'pp') then
    set trigger = "ppVPDMB"
    set base = /nfs/rhi/STAR/Data/ppMBRun12/sum
    echo "Running on the ppMB-triggered data!"
else if ($2 == 'JP2' && $3 == 'pA') then
    set trigger = "pAuJP2"
    set base = /nfs/rhi/STAR/Data/P16id/production_pAu200_2015/HT/pAu_2015_200_HT_1
    echo "Running on the pAuJP2-triggered data!"
else if ($2 == 'BBCMB' && $3 == 'pA') then
    set trigger = "pAuBBCMB"
    set base = /nfs/rhi/STAR/Data/P16id/production_pAu200_2015/MB/pAu_2015_200_MB_1
    echo "Running on the pAuBBCMB-triggered data!"
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

    set arg = "$outLocation $outName $trigger $dummy $Files"

    echo "now submitting this script: "
    echo qsub -V -q erhiq -l mem=4GB -o $LogFile -e $ErrFile -N $1 -- ${ExecPath}/submit/qwrap.sh ${ExecPath} $execute $arg

    qsub -V -q erhiq -l mem=4GB -o $LogFile -e $ErrFile -N $1 -- ${ExecPath}/submit/qwrap.sh ${ExecPath} $execute $arg   
    #add back in a second: -q erhiq

end

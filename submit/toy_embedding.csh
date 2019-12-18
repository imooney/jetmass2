#!/bin/csh
#Isaac Mooney, WSU - Nov. 2019

#args:
#1: the jet radius. Options: any number > 0 & <= 9.9; passed in without the decimal, e.g. "06" for 0.6.
#2: an analysis-specific wildcard (not required). Currently being used in data to take either full (full) or charged (ch) jets.

# Arguments       
if ( $# < "2") then
    echo 'Error: illegal number of parameters'
    exit
endif
    
set ExecPath = `pwd`
set arg = '' 

set radius = $1
set wc1 = $2

#~~~directory / file setup~~~#
make bin/toy_embedding || exit
set execute = './bin/toy_embedding'
# Create the folder name for output
set outFile = toy_embedding
echo 'Running toy embedding!'

set trigger_pp = "ppJP2"
#set base_pp = /tier2/home/groups/rhi/STAR/Data/ppJP2Run12/sum	
set base_pp = ~/jetmass2/lists/pp_JP2_filelist
echo "Running on the ppJP2-triggered data!"
	
set trigger_pAu = "pAuBBCMB"
#set base_pAu = /tier2/home/groups/rhi/STAR/Data/P16id/production_pAu200_2015/MB/pAu_2015_200_MB_1
set base_pAu = ~/jetmass2/lists/pAu_MB_filelist
echo "Running on the pAuBBCMB-triggered data!"

# Make output directories if they don't already exist            
if ( ! -d out/${outFile} ) then
mkdir -p out/${outFile}
endif

if ( ! -d log/${outFile} ) then
mkdir -p log/${outFile}
endif

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

set file_arr = ( 0 1 2 3 4 5 6 7 8 9 ) #for 10 lists of files ending in 0 through 9

# Now Submit jobs for each data file
foreach input ( ${file_arr} ) #${base_pp}* )

    #synthesize the output file base name from the input file base name and the options passed
    #set OutBase = `basename ${input} | sed 's/.txt//g'`
    set uscore = "_" #useful for chaining variables together
    set OutBase = "ppJP2embedpAuBBCMB$uscore${wc1}${uscore}R${radius}${uscore}${input}"
   #make the output path and names
    set outLocation = out/${outFile}/
    set outName = ${OutBase}.root

    #input files
    set Files_pp = ${base_pp}${input}'.txt'
    set Files_pAu = ${base_pAu}${input}'.txt'
    echo ${Files_pp}' '${Files_pAu}
    
    #log files
    set LogFile = log/${outFile}/${OutBase}.log
    set ErrFile = log/${outFile}/${OutBase}.err

    echo "Logging output to " $LogFile
    echo "Logging errors to " $ErrFile

    set dummy = 'dummy' #this extra argument is for compatibility between QA and the other segments of the analysis (or for future expansion)

    #setting command-line args based on the analysis file being used:
    set arg = "$outLocation $outName ${radius} ${wc1} ${Files_pp} ${Files_pAu}"

    echo "now submitting this script: "
    echo qsub -V -q erhiq -l mem=4GB -o $LogFile -e $ErrFile -N "toy_embed" -- ${ExecPath}/submit/qwrap.sh ${ExecPath} $execute $arg

    qsub -V -q erhiq -l mem=4GB -o $LogFile -e $ErrFile -N "toy_embed" -- ${ExecPath}/submit/qwrap.sh ${ExecPath} $execute $arg   
    #add back in a second: -q erhiq

end

#! /usr/bin/env csh

#  submit_rivet.sh
#  Isaac Mooney, WSU - July 2019

#arguments:
#   1: the MC which produced the input HepMCs. Options: P6, P8, H7 [ = Pythia6, Pythia8, or Herwig7]
#   2: whether the input HepMCs decayed the intermediate/long-lived hadrons, or whether, instead, hadronization was turned off. Options: decayed, undecayed, partonicFS
#   3: the desired jet radius. Format: e.g. 04 for R = 0.4.
#   4: whether this is a test or a real analysis run. Options: test, showtime [or anything else]


#checks that the proper number arguments was given
if ( $# < "4") then
    echo "Error: illegal number of parameters. Exiting."
    exit
endif

#build the library first
cd src
rivet-buildplugin ../submit/RivetSTAR_MASS.so STAR_MASS.cc -lfastjetcontribfragile `(root-config --cflags --evelibs)`
cd ..

set TYPE = $1
set DECAYS = $2
set RADIUS = $3
set TEST = $4
set nHardBins = 11
set nEvents = "_10000Events_"
set base = `pwd`
set ExecPath = ${base}/submit
set P6Path = /nfs/rhi/STAR/Data/RhicJEWEL
set P8Path = /nfs/rhi/STAR/Data/IsaacsPy8
set H7Path = /nfs/rhi/STAR/Data/RhicHerwig
#some default file path/name components
set whichdir = "woDecay" #default
set inpath = $P8Path #default
set Filename = "PYTHIA8" #default
set Fileend = "_1MEvents_200GeV.hepmc" #default

#Runs over the full Pythia files if not in debug mode
if ($TEST != "test") then
    set nEvents = "_1MEvents_"
endif

#Adjusts all necessary path/name components to match the desired input/output
if ($TYPE == "P6") then
    set inpath = $P6Path
    set whichdir = ''
    set Filename = "pp200GeV"
    set Fileend = "_1Mevt.hepmc"
else if ($TYPE == "P8") then
    #writing flags to the file
    echo "#Determines whether Herwig or Pythia HepMC files will be the input\ninputIsPythia = TRUE\n#Sets the jet radius parameter for clustering\njetRadius = "$RADIUS"\n#If hadronic final state, turns weak decays on or off; otherwise, partonic FS\ndecayFlag = "$DECAYS"\n" > submit/rivet_flags.txt
    set inpath = $P8Path
    set Filename = "PYTHIA8"
    if ($DECAYS == "decayed") then
	set whichdir = "wDecay"
        set Fileend = ${nEvents}"200GeV_decays_on.hepmc"
    else if ($DECAYS == "undecayed") then
        set whichdir = "woDecay"
        set Fileend = ${nEvents}"200GeV_decays_off.hepmc"
    else if ($DECAYS == "partonicFS") then
        set whichdir = "noHad"
        set Fileend = ${nEvents}"200GeV_decays_off.hepmc"
    else
        echo "Error: unrecognized option for hadronic decays / hadronization. Exiting."
        exit
    endif
else if ($TYPE == "H7") then
    #writing flags to the file
    echo "#Determines whether Herwig or Pythia HepMC files will be the input\ninputIsPythia = FALSE\n#Sets the jet radius parameter for clustering\njetRadius = "$RADIUS"\n#If hadronic final state, turns weak decays on or off; otherwise, partonic FS\ndecayFlag = "$DECAYS"\n" > submit/rivet_flags.txt
    set inpath = $H7Path
    if ($DECAYS == "decayed") then
        set whichdir = "wDecay"
        set Fileend = "-S76456.hepmc"
        set Filename = RHIC200_${whichdir}
    else if ($DECAYS == "undecayed") then
        set whichdir = "woDecay"
        set Fileend = "-S76456.hepmc"
        set Filename = RHIC200_${whichdir}
    else if ($DECAYS == "partonicFS") then
        set whichdir = "noHad"
        set Fileend = "-S876.hepmc"
        set Filename = "noHadro_rhic"
    else
        echo "Error: unrecognized option for hadronic decays / hadronization. Exiting."
	exit
    endif
else
    echo "Error: unable to interpret analysis type. Next time, try one of the following: P6, P8, H7. Exiting."
    exit
endif

@ i = 0
@ low = 0
@ high = 0

#WHILE LOOP TO SUBMIT ALL 11 JOBS
while ($i < $nHardBins)
    #the pT-hat range is part of the file names, so need to set the strings that appear there for each case
    if ($i == 0) then
        set low_alt = ''
        set high_alt = ''
        set analysis_string = ''
        set low = "_5" #this is use for the input name - different than the analysis class name
        set high = 10 #same as above line
    else if ($i > 0 && $i < 9) then
	@ low = ($i + 1)
	@ low = ($low * 5)
        @ high = ($low + 5)
	set low = "_"$low
	set analysis_string = "PTHAT"
	set low_alt = $low
	set high_alt = $high
    else if ($i == 9) then
        set low = "_50"
        set high = 60
	set low_alt = $low
	set high_alt = $high
    else if ($i == 10) then
        set low = "_60"
        set high = 80
	set low_alt = $low
	set high_alt = $high
    else
        echo "Error: unable to interpret pT hard bin edges. Exiting."
        exit
    endif
    
    #basic format:
    #qsub [basic args] -N ~nameOfJob~ -o ~nameOfLogFile~ -e ~nameOfErrorFile~ -- ~QwrapFileName~ ~pathToAnalysisSharedLibrary~ rivet -a ~analysisClassName~ --pwd ~inputFileName~ -H ~yodaOutputFileName~
    echo "DEBUG: "qsub -V -p 10 -q erhiq -l mem=4gb -W umask=0022 -N rivet_${TYPE}${DECAYS}_${i} -o log/rivet_${TYPE}${DECAYS}_${i}.log -e log/rivet_${TYPE}${DECAYS}_${i}.err -- ${ExecPath}/qwrap.sh ${ExecPath} rivet -a STAR_MASS${low_alt}${analysis_string}${high_alt} --pwd ${inpath}/${whichdir}/${Filename}${low}pthat${high}${Fileend} -H yoda/${TYPE}${low}pthat${high}${nEvents}${DECAYS}.yoda

    qsub -V -p 10 -q erhiq -l mem=4gb -W umask=0022 -N rivet_${TYPE}${DECAYS}_${i} -o log/rivet_${TYPE}${DECAYS}_${i}.log -e log/rivet_${TYPE}${DECAYS}_${i}.err -- ${ExecPath}/qwrap.sh ${ExecPath} rivet -a STAR_MASS${low_alt}${analysis_string}${high_alt} --pwd ${inpath}/${whichdir}/${Filename}${low}pthat${high}${Fileend} -H yoda/${TYPE}${low}pthat${high}${nEvents}${DECAYS}.yoda

    @ i ++ #increments the current job number
end
#END OF SUBMIT LOOP

#Writing a stock comment to the Rivet parameters file
echo "\n#Isaac Mooney, WSU - July 2019\n#Since Rivet does not have an obvious way to pass in command line arguments,\n#this file will serve the same purpose.\n#The submit/submit_rivet.sh submit script will change the flags depending on\n#the arguments passed to it just before runtime. The edited parameters will be read in\n#during runtime by the init() function in the submit/RivetSTAR_MASS.so shared library." >> submit/rivet_flags.txt

echo "Done."

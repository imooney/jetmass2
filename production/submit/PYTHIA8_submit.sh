#! /usr/bin/env csh

set ExecPath = `pwd`
set Exec = "./main42"
set Decay = "on"
# submit the jobs
# pthat 5-10
qsub -V -p 10 -q erhiq -l mem=4gb -W umask=0022 -N PYTHIA8_200GeV -o Log/PYTHIA8_5pthat10_200GeV.log -e Log/PYTHIA8_5pthat10_200GeV.err -- ${ExecPath}/qwrap.sh ${ExecPath} $Exec PYTHIA8_5pthat10_1MEvents_200GeV.cmnd PYTHIA8_5pthat10_MinBias_200GeV_decays_$Decay.hepmc $Decay

# pthat 10-15
qsub -V -p 10 -q erhiq -l mem=4gb -W umask=0022 -N PYTHIA8_200GeV -o Log/PYTHIA8_10pthat15_200GeV.log -e Log/PYTHIA8_10pthat15_200GeV.err -- ${ExecPath}/qwrap.sh ${ExecPath} $Exec PYTHIA8_10pthat15_1MEvents_200GeV.cmnd PYTHIA8_10pthat15_MinBias_200GeV_decays_$Decay.hepmc $Decay

# pthat 15-20
qsub -V -p 10 -q erhiq -l mem=4gb -W umask=0022 -N PYTHIA8_200GeV -o Log/PYTHIA8_15pthat20_200GeV.log -e Log/PYTHIA8_15pthat20_200GeV.err -- ${ExecPath}/qwrap.sh ${ExecPath} $Exec PYTHIA8_15pthat20_1MEvents_200GeV.cmnd PYTHIA8_15pthat20_MinBias_200GeV_decays_$Decay.hepmc $Decay

# pthat 20-25
qsub -V -p 10 -q erhiq -l mem=4gb -W umask=0022 -N PYTHIA8_200GeV -o Log/PYTHIA8_20pthat25_200GeV.log -e Log/PYTHIA8_20pthat25_200GeV.err -- ${ExecPath}/qwrap.sh ${ExecPath} $Exec PYTHIA8_20pthat25_1MEvents_200GeV.cmnd PYTHIA8_20pthat25_MinBias_200GeV_decays_$Decay.hepmc $Decay

# pthat 25-30
qsub -V -p 10 -q erhiq -l mem=4gb -W umask=0022 -N PYTHIA8_200GeV -o Log/PYTHIA8_25pthat30_200GeV.log -e Log/PYTHIA8_25pthat30_200GeV.err -- ${ExecPath}/qwrap.sh ${ExecPath} $Exec PYTHIA8_25pthat30_1MEvents_200GeV.cmnd PYTHIA8_25pthat30_MinBias_200GeV_decays_$Decay.hepmc $Decay

# pthat 30-35
qsub -V -p 10 -q erhiq -l mem=4gb -W umask=0022 -N PYTHIA8_200GeV -o Log/PYTHIA8_30pthat35_200GeV.log -e Log/PYTHIA8_30pthat35_200GeV.err -- ${ExecPath}/qwrap.sh ${ExecPath} $Exec PYTHIA8_30pthat35_1MEvents_200GeV.cmnd PYTHIA8_30pthat35_MinBias_200GeV_decays_$Decay.hepmc $Decay

# pthat 35-40
qsub -V -p 10 -q erhiq -l mem=4gb -W umask=0022 -N PYTHIA8_200GeV -o Log/PYTHIA8_35pthat40_200GeV.log -e Log/PYTHIA8_35pthat40_200GeV.err -- ${ExecPath}/qwrap.sh ${ExecPath} $Exec PYTHIA8_35pthat40_1MEvents_200GeV.cmnd PYTHIA8_35pthat40_MinBias_200GeV_decays_$Decay.hepmc $Decay

# pthat 40-45
qsub -V -p 10 -q erhiq -l mem=4gb -W umask=0022 -N PYTHIA8_200GeV -o Log/PYTHIA8_40pthat45_200GeV.log -e Log/PYTHIA8_40pthat45_200GeV.err -- ${ExecPath}/qwrap.sh ${ExecPath} $Exec PYTHIA8_40pthat45_1MEvents_200GeV.cmnd PYTHIA8_40pthat45_MinBias_200GeV_decays_$Decay.hepmc $Decay

# pthat 45-50
qsub -V -p 10 -q erhiq -l mem=4gb -W umask=0022 -N PYTHIA8_200GeV -o Log/PYTHIA8_45pthat50_200GeV.log -e Log/PYTHIA8_45pthat50_200GeV.err -- ${ExecPath}/qwrap.sh ${ExecPath} $Exec PYTHIA8_45pthat50_1MEvents_200GeV.cmnd PYTHIA8_45pthat50_MinBias_200GeV_decays_$Decay.hepmc $Decay

# pthat 50-60
qsub -V -p 10 -q erhiq -l mem=4gb -W umask=0022 -N PYTHIA8_200GeV -o Log/PYTHIA8_50pthat60_200GeV.log -e Log/PYTHIA8_50pthat60_200GeV.err -- ${ExecPath}/qwrap.sh ${ExecPath} $Exec PYTHIA8_50pthat60_1MEvents_200GeV.cmnd PYTHIA8_50pthat60_MinBias_200GeV_decays_$Decay.hepmc $Decay

# pthat 60-80
qsub -V -p 10 -q erhiq -l mem=4gb -W umask=0022 -N PYTHIA8_200GeV -o Log/PYTHIA8_60pthat80_200GeV.log -e Log/PYTHIA8_60pthat80_200GeV.err -- ${ExecPath}/qwrap.sh ${ExecPath} $Exec PYTHIA8_60pthat80_1MEvents_200GeV.cmnd PYTHIA8_60pthat80_MinBias_200GeV_decays_$Decay.hepmc $Decay

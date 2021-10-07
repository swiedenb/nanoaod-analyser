import os

with open('rpv.txt', 'r') as fp:
    for line in fp:
        mass = line.split("M-")[-1].split("_")[0]
        print (line)
        with open('run_RPV_' + mass + '.json', 'w') as op:
            op.write('{\n')
            op.write('"dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/user/chschule/NanoAOD/NanoAODv7_Post/2018/RPVresonantToMuTau_M-' + mass + '_LLE_LQD-001_TuneCP5_13TeV-calchep-pythia8/ : cfg/mc_2018.cfg"\n')
            op.write('}')

import fileinput
import glob

files=glob.glob("*.json")
for filename in files:
    with fileinput.FileInput(filename, inplace=True, backup='.bak') as file:
        for line in file:
            print(line.replace("dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/", "davs://grid-webdav.physik.rwth-aachen.de:2889/"), end='')

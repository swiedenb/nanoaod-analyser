import fileinput
import glob

files=glob.glob("*.json")
for filename in files:
    with fileinput.FileInput(filename, inplace=True, backup='.bak') as file:
        for line in file:
            print(line.replace("signal_2017.cfg", "signal_2017_zprime.cfg"), end='')

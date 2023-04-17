import fileinput
import glob

files=glob.glob("*.json")
for filename in files:
    with fileinput.FileInput(filename, inplace=True, backup='.bak') as file:
        for line in file:
            print(line.replace("mc_2018.cfg", "signal_2017.cfg").replace("chschule", "swiedenb").replace("2018", "2017"),  end='')

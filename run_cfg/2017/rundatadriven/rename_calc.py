import fileinput
import glob

files=glob.glob("*.json")
for filename in files:
       # Read in the file
    with open(filename, 'r') as file :
      filedata = file.read()

    # Replace the target string
    filedata = filedata.replace('calc', 'run')

    # Write the file out again
    with open(filename, 'w') as file:
      file.write(filedata) 
#    with fileinput.FileInput(filename, inplace=True, backup='.bak') as file:
#        for line in file:
#            line.replace(".json", "_calc_datadriven.json")

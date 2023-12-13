import json
import glob
import subprocess

files=glob.glob("bak/*.json")

for filename in files:
    with open(filename,'r') as f:
        injson = json.load(f)

    folder = injson.keys()[0] 
    cfg = injson[injson.keys()[0]]
    a = subprocess.Popen('gfal-ls ' + folder, shell=True, stdout=subprocess.PIPE)
    filenames = a.stdout.read()
    filenames = filenames.split("\n")
    if filenames[-1] == '':del filenames[-1]
    sample = folder.split('/')[-2]
    jsonout = {
                "sample":sample,
                "config": cfg,
                "folder": folder,
                "files" : filenames
              }

    print filename.split('/')[-1]
    with open(filename.split('/')[-1],'w') as json_file:
        json.dump(jsonout, json_file,indent = 4, sort_keys=True)

import glob
import subprocess

listdir=glob.glob("*bak")
for oldname in listdir:
    try:
        newname=oldname.replace(".bak","")
        subprocess.call(["mv",oldname,newname])
    except:
        pass


#import glob
#import subprocess
#
#listdir=glob.glob("*")
#for oldname in listdir:
    #try:
        #newname=oldname.replace("Num","")
        #subprocess.call(["mv",oldname,newname])
    #except:
        #pass

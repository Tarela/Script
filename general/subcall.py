import subprocess
def sp(cmd):
    a=subprocess.Popen(cmd,stdout=subprocess.PIPE,shell='TRUE')
    ac = a.communicate()
    return ac

def sperr(cmd):
    a=subprocess.Popen(cmd,stderr=subprocess.PIPE,shell='TRUE')
    ac = a.communicate()
    return ac 

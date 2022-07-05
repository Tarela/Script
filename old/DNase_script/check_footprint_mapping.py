'''
Created on 2011-5-24

@author: undergraduate
'''
#!/usr/bin/env python
#Time-stamp:<2011-05-24 Sheng'en>
"""
Description:

"""

# ------------------------------------
# Python Modual
# ------------------------------------

import os,sys,re
from optparse import OptionParser
import logging
import string

# ------------------------------------
# error and warning
# ------------------------------------


# ------------------------------------
# Misc functions
# ------------------------------------

def








# ------------------------------------
# Main function
# ------------------------------------

def main():
    usage = "python %prog <-p deletepath> "
    description = """delete all unuse pipeline result in given path"""

    optparser = OptionParser(version="%prog 0.72",description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")

#========core options=============
    optparser.add_option("-i","--inputbed",dest="inputbed",type="str",
                         help="")
    optparser.add_option("-o","--outputbed",dest="outputbed",type="str",
                         help="")
    optparser.add_option("-l","--length",dest="length",type="",
                         help="")


    (options,args) = optparser.parse_args()

    inputbed = options.inputbed
    outputbed = options.outputbed
    length = options.length

    if not inputbed or not outputbed:
        optparser.print_help()
        sys.exit(1)
    



if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)


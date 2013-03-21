# ===========================================================================

import string

# ======================================================================

class Configuration(object):

    def __init__(self,configfile):
        self.file = configfile
        self.parameters = {}
        self.read()
        self.convert()
        return

    # ------------------------------------------------------------------
    # Read in values by keyword and populate the parameters dictionary.

    def read(self):

        thisfile = open(self.file)
        # Go through file, reading it line by line:
        for line in thisfile:
            if line[0] == '#': continue
            if line[0:2] == '\n': continue
            line.strip()
            line = line.split('#')[0]
            line = ''.join(line.split())
            line = line.split(':')
            name = line[0]
            value = line[1]
            self.parameters[name] = value

        thisfile.close()

        return

    # ------------------------------------------------------------------
    # Convert string values into floats where necessary.

    def convert(self):

        for key in self.parameters.keys():
            try: 
                self.parameters[key] = float(self.parameters[key])
                # print "Converted "+key+" to floating point:",self.parameters[key]
            except ValueError:
                pass
        
        return

# ======================================================================


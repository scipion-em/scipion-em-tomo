
from os.path import join

import pyworkflow as pw
import tests.DataSet

class DataSet(tests.DataSet):

    _datasetDict = {}  # store all created datasets

    def __init__(self, name, folder, files, url=None):
        print "init?"
        """ 
        Params:
            
        #filesDict is dict with key, value pairs for each file
        """
        self._datasetDict[name] = self
        self.folder = folder
        self.path = join(pw.Config.SCIPION_TESTS, folder)
        self.filesDict = files
        self.url = url
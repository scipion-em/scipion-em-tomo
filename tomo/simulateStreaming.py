# **************************************************************************
# *
# * Authors:    Alberto Garcia Mena (alberto.garcia@cnb.csic.es)
# *
# * [1] SciLifeLab, Stockholm University
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import time
import sys
import os
import shutil
'''
Generates the mrc files in a specific path simulating the microscope 
adquisition rate. Also generates a mdoc with the information of each mrc file 
generated at the same rate that the files.
INPUTS
    input pathOrg: path with all the files for the simulation. Needs a set of mrc 
    files endwith '_X.X.mrc' (X.X is the angle of the tilt), a general mdoc file with
    the information of all the tilts
    fileNameMdoc: Name of the general mdoc file 
    timeSleep: Time waiting until the next file is generated
    
The delay (timeTrigger) until starts to generate the files is 20 segs
'''


timeTrigger = 2

class simulateStreaming():
    def __int__(self):
        pass

    def simulateStreaming(self, pathOrg, fileNameMdoc, pathStreaming, timeSleep):
        self.numtilts = 200
        self.pathStreaming = pathStreaming
        self.fileNameMdoc = fileNameMdoc
        self.pathMdoc = os.path.join(pathOrg, fileNameMdoc)
        pathFolderData, mdocName = os.path.split(self.pathMdoc)
        # copied all files in streamingFolder. The first the mdoc, after that
        # added files with a sleeptime
        if os.path.isdir(pathStreaming):
            shutil.rmtree(pathStreaming)
        os.makedirs(pathStreaming, exist_ok=True)
        self.parseMdoc(self.pathMdoc)
        self.copyMdoc()#just the header
        print('sleeping to launch protocol in Scipion {} s...'.format(timeSleep))
        time.sleep(timeTrigger)

        for item in self.angleList:
            if item != None:
                for nameFile in os.listdir(pathFolderData):
                    index = nameFile.rfind('_') + 1
                    if nameFile[-3:] == 'mrc' and nameFile[index:-4] == str(item):
                            print('path:', nameFile)
                            shutil.copyfile(os.path.join(pathFolderData, nameFile),
                                            os.path.join(pathStreaming, nameFile))
                            #print('self.angleList.index(item): {}, item: {}'.format(
                            #    self.angleList.index(item), item))
                            self.copyMdoc(self.angleList.index(item))
                            print('sleeping  {} s...'.format(timeSleep))

                            time.sleep(timeSleep)

        mrcBool = False
        mrcMainFile = self.fileNameMdoc.split('mdoc')[0]
        mrcMainFilePath = os.path.join(pathFolderData, mrcMainFile)
        if os.path.isfile(mrcMainFilePath):
            mrcBool = True
        else:
            mrcMainFile = mrcMainFile + '.mrcs'
            mrcMainFilePath = os.path.join(pathFolderData, mrcMainFile)
            if os.path.isfile(mrcMainFilePath):
                mrcBool = True
                mrcMainFilePath = os.path.join(pathFolderData, mrcMainFile)
            else:
                print('mrc main file not found')

        if mrcBool == True:
            mrcMainFileStreamingPath = os.path.join(self.pathStreaming,
                                                    mrcMainFile)
            shutil.copyfile(mrcMainFilePath, mrcMainFileStreamingPath)
        else:
            print('mrc main file not found')

    def copyMdoc(self, item=None):
        if item == None: #header
            self.pathmdocStreaming = os.path.join(self.pathStreaming, self.fileNameMdoc )
            with open(self.pathmdocStreaming, 'w') as f:
                for line in self.mainText:
                    if line != None:
                        f.write(line)
        else:
            with open(self.pathmdocStreaming, 'a') as f:
                for line in self.tiltList[item]:
                    if line != None:
                        f.write(line)

    def parseMdoc(self, pathMdoc):
        header = True
        self.mainText = []
        self.tiltList = [None] * self.numtilts
        self.angleList = [None] * self.numtilts
        index = 0
        print('pathMdoc:', pathMdoc, '\n')
        with open(pathMdoc, 'r') as f:
            #myline = f.readline()
            # while myline:
            for line in f:
                if line.startswith('[ZValue'):
                    index += 1
                    header = False
                    self.tiltList[index] = line
                elif header == True:
                    if line != None:
                        self.mainText.append(line)
                else:
                    if line[:9] == 'TiltAngle':
                        # print('index: {}, float(line[12:]): {}'.
                        #       format(index, float(line[12:])))
                        self.angleList[index] = round(float(line[12:]), 1)
                    self.tiltList[index] = self.tiltList[index] + line

            self.tiltList.remove(None)
            self.angleList.remove(None)

if __name__ == '__main__':
    args = sys.argv[1:]
    if len(args) >= 2:
        pathOrg = args[0]
        fileNameMdoc = args[1]
        pathStreaming = args[2]
        timeSleep = int(args[3])
        print('\n', pathOrg, '\n', fileNameMdoc, '\n', pathStreaming,  '\n', timeSleep, '\n')
    else:
        print('Args mandatories: pathOrg, fileNameMdoc, timeSleep')

    simulateStreamingC = simulateStreaming()
    simulateStreamingC.simulateStreaming(pathOrg, fileNameMdoc, pathStreaming, timeSleep)
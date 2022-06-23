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

class simulateStreaming():
    def __int__(self):
        pass

    def simulateStreaming(self, pathOrg, fileNameMdoc, pathStreaming, timeSleep):
        self.numtilts = 200
        self.pathStreaming = pathStreaming
        pathMdoc = os.path.join(pathOrg, fileNameMdoc)
        pathFolderData, mdocName = os.path.split(pathMdoc)
        #print('pathMdoc: ', pathMdoc)
        # copied all files in streamingFolder. The first the mdoc, after that
        # added files with a sleeptime
        os.makedirs(pathStreaming, exist_ok=True)
        self.parseMdoc(pathMdoc)

        #shutil.copyfile(pathMdoc, os.path.join(pathStreaming, mdocName))
        self.copyMdoc()#just the header
        #set streaming path for protocolfilesPath.set(pathStreaming)
        #print('pathStreaming: ', pathStreaming)
        #print('pathStreaming updated: ', glob(join(pathOrg, self.filesPattern.get()))[0])

        print('sleeping to launch protocol in Scipion {} s...'.format(timeSleep))
        time.sleep(30)

        for item in self.angleList:
            for nameFile in os.listdir(pathFolderData):
                #headPath, nameFile = os.path.split(path)
                if nameFile[-3:] == 'mrc' and float(nameFile[-7:-4]) == item:
                        print('path:', nameFile)
                        shutil.copyfile(os.path.join(pathFolderData, nameFile),
                                        os.path.join(pathStreaming, nameFile))
                        self.copyMdoc(self.angleList.index(item))
                        print('sleeping  {} s...'.format(timeSleep))

                        time.sleep(timeSleep)

    def copyMdoc(self, item=None):
        if item == None: #header
            self.pathmdocStreaming = os.path.join(self.pathStreaming, "stack10.mdoc")
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
# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
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

import os
import pwem
from .constants import (NAPARI_ENV_ACTIVATION, NAPARI_ACTIVATION_CMD,
                        getNaparyEnvName, V0_3_11, NAPARI_DEF_VER)

__version__ = '3.4.0'
_logo = "icon.png"
_references = []


class Plugin(pwem.Plugin):
    @classmethod
    def _defineVariables(cls):
        cls._defineVar(NAPARI_ENV_ACTIVATION, NAPARI_ACTIVATION_CMD)

    @classmethod
    def getEnviron(cls):
        return None

    @classmethod
    def getDependencies(cls):
        """ Return a list of dependencies. Include conda if
            activation command was not found. """
        condaActivationCmd = cls.getCondaActivationCmd()
        neededProgs = ['wget']
        if not condaActivationCmd:
            neededProgs.append('conda')

        return neededProgs

    @classmethod
    def addNapariPackage(cls, env, version, default=False):
        ENV_NAME = getNaparyEnvName(version)
        NAPARI_INSTALLED = f"napari_{version}_installed"
        installCmd = [cls.getCondaActivationCmd(),
                      f'conda create -y -n {ENV_NAME} -c conda-forge',
                      'python=3.10 napari=0.4.17 pyqt pip &&',
                      f'conda activate {ENV_NAME} &&',
                      f'pip install napari-tomotwin napari-boxmanager=={version}']

        # Flag installation finished
        installCmd.append(f'&& touch {NAPARI_INSTALLED}')

        napari_commands = [(" ".join(installCmd), NAPARI_INSTALLED)]

        envPath = os.environ.get('PATH', "")
        # keep path since conda likely in there
        installEnvVars = {'PATH': envPath} if envPath else None
        env.addPackage(f'napari', version=version,
                       tar='void.tgz',
                       commands=napari_commands,
                       neededProgs=cls.getDependencies(),
                       default=default,
                       vars=installEnvVars)

    @classmethod
    def defineBinaries(cls, env):
        cls.addNapariPackage(env, V0_3_11)
        cls.addNapariPackage(env, NAPARI_DEF_VER, default=True)

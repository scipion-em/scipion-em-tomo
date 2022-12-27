# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Estrella Fernandez Gimenez (me.fernandez@cnb.csic.es)
# *
# *  BCU, Centro Nacional de Biotecnologia, CSIC
# *
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

from pyworkflow.viewer import DESKTOP_TKINTER
from pyworkflow.protocol.params import EnumParam
from tomo.protocols import ProtSplitEvenOddTomoSet
from pwem.viewers.viewer_base import EmProtocolViewer


class SplitEvenOddViewer(EmProtocolViewer):
    """ Viewer for protocol_split_evenodd_subtomos (or tomos)."""
    _label = 'viewer Split Even/Odd'
    _targets = [ProtSplitEvenOddTomoSet]
    _environments = [DESKTOP_TKINTER]

    def _defineParams(self, form):
        form.addSection(label='Viewer')
        form.addParam('set', EnumParam, default=0, choices=['odd set', 'even set'], display=EnumParam.DISPLAY_HLIST,
                      label='Show')

    def _getVisualizeDict(self):
        return {'set': self._showVolumeSlices}

    def _showVolumeSlices(self, param=None):
        if self.getEnumText('set') == 'odd set':
            return [self.objectView(self.protocol.outputset_odd)]
        else:
            return [self.objectView(self.protocol.outputset_even)]

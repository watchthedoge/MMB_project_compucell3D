# -*- coding: utf-8 -*-

from __future__ import division
from PlayerPython import *
import CompuCellSetup
from PySteppables import *
from datetime import datetime
import CompuCell
import sys
import numpy as np
import random
from tempfile import TemporaryFile
from PySteppablesExamples import MitosisSteppableBase

# ---------------- GLOBALS ----------------
growthratewt = 6.3
growthratescrb = 3.4
stiffness_kd = 0.9
stiffness_wt = 2.0
CI = 0.1

minvol = 1100
maxvol = 2200
relaxtime = 200


APOP_SCALE = 0.1

adderlist = []

# ---------------- INITIALIZATION ----------------
class ConstraintInitializerSteppableAdder(SteppableBasePy):
    def __init__(self, _simulator, _frequency=1):
        SteppableBasePy.__init__(self, _simulator, _frequency)

    def start(self):
        for cell in self.cellList:
            self.moveCell(cell, [random.randint(-5,5), random.randint(-5,5), 0])

        cellarr = list(range(110))
        random.shuffle(cellarr)
        kdcells = cellarr[:20]

        for cell in self.cellList:
            if cell.id in kdcells:
                cell.type = 2

        for cell in self.cellList:
            cell.targetVolume = random.normalvariate(1800, 500)
            if cell.type == 1:
                cell.lambdaVolume = stiffness_wt
            else:
                cell.lambdaVolume = stiffness_kd
            adderlist.append([cell.id, 1800, cell.type])

# ---------------- GROWTH ----------------
class GrowthSteppableLinear(SteppableBasePy):
    def __init__(self,_simulator,_frequency=10):
        SteppableBasePy.__init__(self,_simulator,_frequency)

    def step(self,mcs):
        if mcs <= relaxtime:
            return

        for cell in self.cellList:
            if cell.type in (1,3,5,7):
                gr = growthratewt
            elif cell.type in (2,4,6,8):
                gr = growthratescrb
            else:
                continue

            delta = max(0, round(random.normalvariate(gr,2.5)))
            cell.targetVolume += delta * np.exp(-(CI/mcs)*((cell.volume-cell.targetVolume)**2))

# ---------------- MITOSIS ----------------
class MitosisSteppableAdder(MitosisSteppableBase):
    def __init__(self,_simulator,_frequency=10):
        MitosisSteppableBase.__init__(self,_simulator,_frequency)

    def step(self,mcs):
        if mcs <= relaxtime:
            return

        for cell in list(self.cellList):
            for x in adderlist:
                if x[0] == cell.id and cell.volume - x[1] > (1800 - mcs/60):
                    self.divideCellAlongMinorAxis(cell)
                    adderlist.append([self.childCell.id, self.childCell.targetVolume])

    def updateAttributes(self):
        self.parentCell.targetVolume /= 2
        self.cloneParent2Child()

# ---------------- DENSITY APOPTOSIS ----------------
class DeathSteppable(SteppableBasePy):
    def __init__(self,_simulator,_frequency=10):
        SteppableBasePy.__init__(self,_simulator,_frequency)

    def start(self):
        global LocalDensity, P_apo, cumulitiveapop
        LocalDensity = []
        P_apo = []
        cumulitiveapop = []

    def step(self,mcs):
        if mcs <= relaxtime:
            return

        apopwt = apopscrb = 0

        for cell in self.cellList:
            if cell.type in (9,10) or cell.volume < 50:
                continue

            dens = 1.0 / cell.volume
            for neighbor,_ in self.getCellNeighborDataList(cell):
                if neighbor:
                    dens += 1.0 / neighbor.volume

            if cell.type in (2,4,6,8):
                p = 0.00072194 / (1 + np.exp(-509.4*(dens*3 - 0.0067)))
                if APOP_SCALE * p > random.random():
                    cell.type = 10
                    cell.targetVolume = 0
                    cell.lambdaVolume = 2
                    cell.dict["apt"] = mcs
                    apopscrb += 1

            elif cell.type in (1,3,5,7):
                p = 0.00033074 / (1 + np.exp(-235.8*(dens*3 - 0.0152)))
                if APOP_SCALE * p > random.random():
                    cell.type = 9
                    cell.targetVolume = 0
                    cell.lambdaVolume = 2
                    cell.dict["apt"] = mcs
                    apopwt += 1

            LocalDensity.append([cell.id, cell.type, (mcs-relaxtime)/10.0, dens])
            P_apo.append([cell.id, cell.type, (mcs-relaxtime)/10.0, dens])

        cumulitiveapop.append([(mcs-relaxtime)/10.0, apopscrb, apopwt])


class CellMotilitySteppable(SteppableBasePy):
    def __init__(self,_simulator,_frequency=10):
        SteppableBasePy.__init__(self,_simulator,_frequency)

    def start(self):
        # iterating over all cells in simulation
        for cell in self.cellList:
            break
            cell.lambdaVecX = 10.1 * random.uniform(-1.0,1.0)
            cell.lambdaVecY = 10.1 * random.uniform(-1.0,1.0)

    def step(self, mcs):

        # warm-up: no motility at all
        if mcs < relaxtime:
            for cell in self.cellList:
                cell.lambdaVecX = 0.0
                cell.lambdaVecY = 0.0
            return

        # ramp from 0 -> 1 over RAMP_MCS
        RAMP_MCS = 500.0
        ramp = min(1.0, (mcs - relaxtime) / RAMP_MCS)

        for cell in self.cellList:
            amp = ramp * 10.1
            cell.lambdaVecX = amp * random.uniform(-1.0, 1.0)
            cell.lambdaVecY = amp * random.uniform(-1.0, 1.0)


class neighbourdata(SteppableBasePy):
    def __init__(self,_simulator,_frequency=10):
        SteppableBasePy.__init__(self,_simulator,_frequency)

    def start(self):
        self.pW = self.addNewPlotWindow(
            _title='CellCount v Time',
            _xAxisTitle='Time(Frames)',
            _yAxisTitle='Cell Count',
            _xScaleType='linear',
            _yScaleType='linear')

        self.pW.addPlot('cellcount WT', _style='Dots', _color='green', _size=5)
        self.pW.addPlot('cellcount scrb', _style='Dots', _color='red', _size=5)

    def step(self,mcs):
        cellcount_wt = 0
        cellcount_scrb = 0

        for cell in self.cellList:
            if cell.type in (1,3,5,7):
                cellcount_wt += 1
            elif cell.type in (2,4,6,8):
                cellcount_scrb += 1

        time = (mcs - relaxtime) / float(10)
        self.pW.addDataPoint('cellcount WT', time, cellcount_wt)
        self.pW.addDataPoint('cellcount scrb', time, cellcount_scrb)

class tracking(SteppableBasePy):
    def __init__(self,_simulator,_frequency=10):
        SteppableBasePy.__init__(self,_simulator,_frequency)

    def start(self):
        global trackingfile
        trackingfile = []

    def step(self,mcs):
        for cell in self.cellList:
            if cell.type == 9 or cell.type == 10:
                state = 1
            else:
                state = 0

            time = (mcs - relaxtime) / float(10)
            trackingfile.append([
                cell.xCOM,
                cell.yCOM,
                int(time),
                int(cell.id),
                int(cell.type),
                state
            ])

    def finish(self):
        np.savetxt(
            'C:\\Users\\Dgradeci\\Dropbox\\CC3D Workspace\\data\\tracking-%s'
            % datetime.now().strftime('%H-%M-%m-%d'),
            trackingfile
        )


# ---------------- CLEANUP ----------------
class cleanup(SteppableBasePy):
    def __init__(self,_simulator,_frequency=1000):
        SteppableBasePy.__init__(self,_simulator,_frequency)
        self.extrusions = []

    def step(self,mcs):
        if mcs <= relaxtime:
            return

        to_delete = []
        for cell in self.cellList:
            if cell.type in (9,10):
                apt = cell.dict.get("apt", None)
                if apt and mcs - apt >= 5:
                    self.extrusions.append([cell.id, cell.type, (mcs-relaxtime)/10.0])
                    to_delete.append(cell)

        for cell in to_delete:
            self.deleteCell(cell)

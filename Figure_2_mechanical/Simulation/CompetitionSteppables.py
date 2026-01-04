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

# SUPER KILLER (type 3)
growthrate_super = 9.0
stiffness_super = 3.5

stiffness_kd = 0.9
stiffness_wt = 2.0
CI = 0.1

minvol = 1100
maxvol = 2200
relaxtime = 200

# apoptosis enabled scaling
APOP_SCALE = 0.1

# --- type IDs ---
WT_TYPE = 1
SCRB_TYPE = 2
SUPER_TYPE = 3

APOP_WT = 9
APOP_SCRB = 10
APOP_SUPER = 11  

USE_PIF_TYPES = True

adderlist = []

# ---------------- INITIALIZATION ----------------
class ConstraintInitializerSteppableAdder(SteppableBasePy):
    def __init__(self, _simulator, _frequency=1):
        SteppableBasePy.__init__(self, _simulator, _frequency)

    def start(self):
        # small random displacement
        for cell in self.cellList:
            self.moveCell(cell, [random.randint(-5,5), random.randint(-5,5), 0])


        if not USE_PIF_TYPES:
            cellarr = list(range(110))
            random.shuffle(cellarr)
            kdcells = cellarr[:20]
            supercells = cellarr[20:30]

            for cell in self.cellList:
                if cell.id in kdcells:
                    cell.type = SCRB_TYPE
                elif cell.id in supercells:
                    cell.type = SUPER_TYPE
                else:
                    cell.type = WT_TYPE

        # initialize mechanical parameters based on type
        for cell in self.cellList:
            cell.targetVolume = random.normalvariate(1800, 500)

            if cell.type == WT_TYPE:
                cell.lambdaVolume = stiffness_wt
            elif cell.type == SCRB_TYPE:
                cell.lambdaVolume = stiffness_kd
            elif cell.type == SUPER_TYPE:
                cell.lambdaVolume = stiffness_super
            else:
                # fallback
                cell.lambdaVolume = stiffness_wt

            adderlist.append([cell.id, 1800, cell.type])

# ---------------- GROWTH ----------------
class GrowthSteppableLinear(SteppableBasePy):
    def __init__(self,_simulator,_frequency=10):
        SteppableBasePy.__init__(self,_simulator,_frequency)

    def step(self,mcs):
        if mcs <= relaxtime:
            return

        for cell in self.cellList:
            if cell.type in (APOP_WT, APOP_SCRB, APOP_SUPER):
                continue

            if cell.type == WT_TYPE:
                gr = growthratewt
            elif cell.type == SCRB_TYPE:
                gr = growthratescrb
            elif cell.type == SUPER_TYPE:
                gr = growthrate_super
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
            if cell.type in (APOP_WT, APOP_SCRB, APOP_SUPER):
                continue

            for x in adderlist:
                if x[0] == cell.id and cell.volume - x[1] > (1800 - mcs/60.0):
                    self.divideCellAlongMinorAxis(cell)
                   
                    adderlist.append([self.childCell.id, self.childCell.targetVolume, self.childCell.type])
                    break

    def updateAttributes(self):
        self.parentCell.targetVolume /= 2
        self.cloneParent2Child()

        # preserve lineage/type on division
        self.childCell.type = self.parentCell.type

        # preserve stiffness by type
        if self.childCell.type == WT_TYPE:
            self.childCell.lambdaVolume = stiffness_wt
            self.parentCell.lambdaVolume = stiffness_wt
        elif self.childCell.type == SCRB_TYPE:
            self.childCell.lambdaVolume = stiffness_kd
            self.parentCell.lambdaVolume = stiffness_kd
        elif self.childCell.type == SUPER_TYPE:
            self.childCell.lambdaVolume = stiffness_super
            self.parentCell.lambdaVolume = stiffness_super

# ---------------- DENSITY APOPTOSIS (MECHANICAL) ----------------
class DeathSteppable(SteppableBasePy):
    def __init__(self,_simulator,_frequency=10):
        SteppableBasePy.__init__(self,_simulator,_frequency)

        # super is most resistant to density-driven death 
        self.super_pmax = 0.00010
        self.super_k = 180.0
        self.super_thr = 0.030  # higher threshold > harder to kill

    def start(self):
        global LocalDensity, P_apo, cumulitiveapop
        LocalDensity = []
        P_apo = []
        cumulitiveapop = []

    def step(self,mcs):
        if mcs <= relaxtime:
            return

        apopwt = 0
        apopscrb = 0
        apopsuper = 0

        for cell in self.cellList:
            if cell.type in (APOP_WT, APOP_SCRB, APOP_SUPER) or cell.volume < 50:
                continue

            dens = 1.0 / cell.volume
            for neighbor,_ in self.getCellNeighborDataList(cell):
                if neighbor and neighbor.type not in (APOP_WT, APOP_SCRB, APOP_SUPER):
                    dens += 1.0 / neighbor.volume

            if cell.type == SCRB_TYPE:
                p = 0.00072194 / (1.0 + np.exp(-509.4*(dens*3.0 - 0.0067)))
                if APOP_SCALE * p > random.random():
                    cell.type = APOP_SCRB
                    cell.targetVolume = 0
                    cell.lambdaVolume = 2
                    cell.dict["apt"] = mcs
                    apopscrb += 1

            elif cell.type == WT_TYPE:
                p = 0.00033074 / (1.0 + np.exp(-235.8*(dens*3.0 - 0.0152)))
                if APOP_SCALE * p > random.random():
                    cell.type = APOP_WT
                    cell.targetVolume = 0
                    cell.lambdaVolume = 2
                    cell.dict["apt"] = mcs
                    apopwt += 1

            elif cell.type == SUPER_TYPE:
                # density-driven death for super 
                p = self.super_pmax / (1.0 + np.exp(-self.super_k*(dens*3.0 - self.super_thr)))
                if APOP_SCALE * p > random.random():
                    cell.type = APOP_SUPER
                    cell.targetVolume = 0
                    cell.lambdaVolume = 2
                    cell.dict["apt"] = mcs
                    apopsuper += 1

            t = (mcs-relaxtime)/10.0
            LocalDensity.append([cell.id, cell.type, t, dens])
            P_apo.append([cell.id, cell.type, t, dens])

        cumulitiveapop.append([(mcs-relaxtime)/10.0, apopscrb, apopwt, apopsuper])


# ---------------- MOTILITY ----------------
class CellMotilitySteppable(SteppableBasePy):
    def __init__(self,_simulator,_frequency=10):
        SteppableBasePy.__init__(self,_simulator,_frequency)

    def start(self):
        for cell in self.cellList:
            break
            cell.lambdaVecX = 10.1 * random.uniform(-1.0,1.0)
            cell.lambdaVecY = 10.1 * random.uniform(-1.0,1.0)

    def step(self, mcs):
        if mcs < relaxtime:
            for cell in self.cellList:
                cell.lambdaVecX = 0.0
                cell.lambdaVecY = 0.0
            return

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
        self.pW.addPlot('cellcount super', _style='Dots', _color='blue', _size=5)

    def step(self,mcs):
        cellcount_wt = 0
        cellcount_scrb = 0
        cellcount_super = 0

        for cell in self.cellList:
            if cell.type == WT_TYPE:
                cellcount_wt += 1
            elif cell.type == SCRB_TYPE:
                cellcount_scrb += 1
            elif cell.type == SUPER_TYPE:
                cellcount_super += 1

        time = (mcs - relaxtime) / float(10)
        self.pW.addDataPoint('cellcount WT', time, cellcount_wt)
        self.pW.addDataPoint('cellcount scrb', time, cellcount_scrb)
        self.pW.addDataPoint('cellcount super', time, cellcount_super)


class tracking(SteppableBasePy):
    def __init__(self,_simulator,_frequency=10):
        SteppableBasePy.__init__(self,_simulator,_frequency)

    def start(self):
        global trackingfile
        trackingfile = []

    def step(self,mcs):
        for cell in self.cellList:
            state = 1 if cell.type in (APOP_WT, APOP_SCRB, APOP_SUPER) else 0
            time = (mcs - relaxtime) / float(10)
            trackingfile.append([cell.xCOM, cell.yCOM, int(time), int(cell.id), int(cell.type), state])

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
            if cell.type in (APOP_WT, APOP_SCRB, APOP_SUPER):
                apt = cell.dict.get("apt", None)
                if apt is not None and (mcs - apt) >= 5:
                    self.extrusions.append([cell.id, cell.type, (mcs-relaxtime)/10.0])
                    to_delete.append(cell)

        for cell in to_delete:
            self.deleteCell(cell)

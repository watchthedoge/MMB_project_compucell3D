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

# SUPER-WINNER (type 3)
growthrate_super = 9.0        # higher than WT
stiffness_super = 3.5         # stiffer than WT

# Contact-competition weights
SUPER_KILLS_WT = 1.5
SUPER_KILLS_SCRB = 3.0
WT_KILLS_SUPER = 0.2
SCRB_KILLS_SUPER = 0.1

minvol = 1100
maxvol = 2200
relaxtime = 200

# scale death probabilities
APOP_SCALE = 0.1

adderlist = []

# living types
WT_TYPE = 1
SCRB_TYPE = 2
SUPER_TYPE = 3

# apoptotic types
APOP_WT = 9
APOP_SCRB = 10
APOP_SUPER = 11


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


class ConstraintInitializerSteppableAdder(SteppableBasePy):
    def __init__(self, _simulator, _frequency=1):
        SteppableBasePy.__init__(self, _simulator, _frequency)

    def start(self):

        # small random displacement
        for cell in self.cellList:
            self.moveCell(cell, [random.randint(-5,5), random.randint(-5,5), 0])

        for cell in self.cellList:
            # target volume initialization
            cell.targetVolume = random.normalvariate(1800, 500)

            # stiffness based ONLY on PIF-defined type
            if cell.type == WT_TYPE:
                cell.lambdaVolume = stiffness_wt

            elif cell.type == SCRB_TYPE:
                cell.lambdaVolume = stiffness_kd

            elif cell.type == SUPER_TYPE:
                cell.lambdaVolume = stiffness_super

            else:
                # fallback safety
                cell.lambdaVolume = stiffness_wt

            # preserve adder logic
            adderlist.append([cell.id, cell.targetVolume, cell.type])


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
                if x[0] == cell.id and cell.volume - x[1] > (1800 - mcs/60):
                    self.divideCellAlongMinorAxis(cell)
                    adderlist.append([self.childCell.id, self.childCell.targetVolume, self.childCell.type])
                    break

    def updateAttributes(self):
        self.parentCell.targetVolume /= 2
        self.cloneParent2Child()
        # preserve lineage/type (SUPER stays SUPER)
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


# ---------------- CONTACT-BASED COMPETITION (SAFE) ----------------
# Uses getCellNeighborDataList + commonSurfaceArea
class DeathSteppable(SteppableBasePy):
    def __init__(self,_simulator,_frequency=10):
        SteppableBasePy.__init__(self,_simulator,_frequency)

        # base death strengths
        self.base_scrb = 3.0 * 0.000416
        self.base_wt = 10.0 * 0.0000416
        # super is hard to kill
        self.base_super = 0.00002

    def start(self):
        global perimeterarray
        global P_apo
        global cumulitiveapop
        cumulitiveapop = []
        perimeterarray = []
        P_apo = []

    def step(self,mcs):
        if mcs <= relaxtime:
            return

        apopwt = 0
        apopscrb = 0
        apopsuper = 0
        deathcount = 0

        for cell in self.cellList:
            if cell.type in (APOP_WT, APOP_SCRB, APOP_SUPER):
                continue
            if cell.volume < 50:
                continue
            if cell.type not in (WT_TYPE, SCRB_TYPE, SUPER_TYPE):
                continue

            total_contact = 0.0
            weighted_opp_contact = 0.0

            for neighbor, commonSurfaceArea in self.getCellNeighborDataList(cell):
                if not neighbor:
                    continue
                if neighbor.type in (APOP_WT, APOP_SCRB, APOP_SUPER):
                    continue
                if neighbor.type not in (WT_TYPE, SCRB_TYPE, SUPER_TYPE):
                    continue

                total_contact += commonSurfaceArea

                # weighted opponent contact
                if cell.type == WT_TYPE:
                    # WT dies from SCRB and SUPER (SUPER stronger)
                    if neighbor.type == SCRB_TYPE:
                        weighted_opp_contact += commonSurfaceArea
                    elif neighbor.type == SUPER_TYPE:
                        weighted_opp_contact += SUPER_KILLS_WT * commonSurfaceArea

                elif cell.type == SCRB_TYPE:
                    # SCRB dies from WT and SUPER (SUPER strongest)
                    if neighbor.type == WT_TYPE:
                        weighted_opp_contact += commonSurfaceArea
                    elif neighbor.type == SUPER_TYPE:
                        weighted_opp_contact += SUPER_KILLS_SCRB * commonSurfaceArea

                elif cell.type == SUPER_TYPE:
                    # SUPER dies weakly from WT and SCRB
                    if neighbor.type == WT_TYPE:
                        weighted_opp_contact += WT_KILLS_SUPER * commonSurfaceArea
                    elif neighbor.type == SCRB_TYPE:
                        weighted_opp_contact += SCRB_KILLS_SUPER * commonSurfaceArea

            if total_contact > 0.0:
                contact_frac = weighted_opp_contact / total_contact
            else:
                contact_frac = 0.0

            # death decision (APOP_SCALE applied)
            if cell.type == SCRB_TYPE:
                p = APOP_SCALE * (self.base_scrb * contact_frac)
                if p > random.random():
                    cell.type = APOP_SCRB
                    cell.dict["apt"] = mcs
                    cell.targetVolume = 0
                    cell.lambdaVolume = 2
                    deathcount += 1
                    apopscrb += 1

            elif cell.type == WT_TYPE:
                p = APOP_SCALE * (self.base_wt * contact_frac)
                if p > random.random():
                    cell.type = APOP_WT
                    cell.dict["apt"] = mcs
                    cell.targetVolume = 0
                    cell.lambdaVolume = 2
                    deathcount += 1
                    apopwt += 1

            elif cell.type == SUPER_TYPE:
                p = APOP_SCALE * (self.base_super * contact_frac)
                if p > random.random():
                    cell.type = APOP_SUPER
                    cell.dict["apt"] = mcs
                    cell.targetVolume = 0
                    cell.lambdaVolume = 2
                    deathcount += 1
                    apopsuper += 1

            time = (mcs - relaxtime) / float(10)
            P_apo.append([cell.id, cell.type, time, contact_frac, deathcount])
            perimeterarray.append([cell.id, cell.type, time, contact_frac])

        time = (mcs - relaxtime) / float(10)
        cumulitiveapop.append([time, apopscrb, apopwt, apopsuper])

    def finish(self):
        np.savetxt('C:\\Users\\Dgradeci\\Dropbox\\CC3D Workspace\\data\\perimeterarray-%s'
                   % datetime.now().strftime('%H-%M-%m-%d'), perimeterarray)
        np.savetxt('C:\\Users\\Dgradeci\\Dropbox\\CC3D Workspace\\data\\P_apo-%s'
                   % datetime.now().strftime('%H-%M-%m-%d'), P_apo)
        np.savetxt('C:\\Users\\Dgradeci\\Dropbox\\CC3D Workspace\\data\\cumulitiveapop-%s'
                   % datetime.now().strftime('%H-%M-%m-%d'), cumulitiveapop)
        return



class DeathSteppablePerimiter(DeathSteppable):
    def __init__(self,_simulator,_frequency=10):
        DeathSteppable.__init__(self,_simulator,_frequency)


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
        np.savetxt('C:\\Users\\Dgradeci\\Dropbox\\CC3D Workspace\\data\\tracking-%s'
                   % datetime.now().strftime('%H-%M-%m-%d'), trackingfile)


class neighbourdata(SteppableBasePy):
    def __init__(self,_simulator,_frequency=10):
        SteppableBasePy.__init__(self,_simulator,_frequency)

    def start(self):
        self.pW = self.addNewPlotWindow(
            _title='CellCount v Time',
            _xAxisTitle='Time(Frames)',
            _yAxisTitle='Cell Count',
            _xScaleType='linear',
            _yScaleType='linear'
        )
        self.pW.addPlot('cellcount WT', _style='Dots', _color='green', _size=5)
        self.pW.addPlot('cellcount scrb', _style='Dots', _color='red', _size=5)
        self.pW.addPlot('cellcount super', _style='Dots', _color='blue', _size=5)

    def step(self,mcs):
        wt = 0
        scrb = 0
        sup = 0

        for cell in self.cellList:
            if cell.type == WT_TYPE:
                wt += 1
            elif cell.type == SCRB_TYPE:
                scrb += 1
            elif cell.type == SUPER_TYPE:
                sup += 1

        time = (mcs - relaxtime) / float(10)
        self.pW.addDataPoint('cellcount WT', time, wt)
        self.pW.addDataPoint('cellcount scrb', time, scrb)
        self.pW.addDataPoint('cellcount super', time, sup)


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
                # more stable delay
                if apt is not None and (mcs - apt) >= 10:
                    self.extrusions.append([cell.id, cell.type, (mcs-relaxtime)/10.0])
                    to_delete.append(cell)

        for cell in to_delete:
            self.deleteCell(cell)

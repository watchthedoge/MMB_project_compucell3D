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

# 🔧 FIX: apoptosis enabled
APOP_SCALE = 0.1

adderlist = []


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

# ---------------- CONTACT-BASED COMPETITION (SAFE) ----------------
class DeathSteppable(SteppableBasePy):
    #  Apoptosis:
    def __init__(self,_simulator,_frequency=10):
        SteppableBasePy.__init__(self,_simulator,_frequency)

    def start(self):
       global perimeterarray
       global P_apo
       global cumulitiveapop
       cumulitiveapop=[]
       perimeterarray=[]
       P_apo=[]

    def step(self,mcs):
        apopwt=0
        apopscrb=0
        deathcount=0
        cellcount=0 
        typecell=0

        if mcs>relaxtime:
            for cell in self.cellList: 

                if cell.type==2 or cell.type==4 or cell.type==6 or cell.type==8:
                    typecell=2
                else:
                    typecell=1

                neighbourpixel=[]
                perimeterpercentage=0

                if cell.type!=9 and cell.type!=10:

                    for pixel in self.getCopyOfCellBoundaryPixels(cell):
                        a=self.point3DToNumpy(pixel)
                        centrecell = self.cellField[a[0],a[1],a[2]]

                        leftpixel = self.cellField[a[0]-1,a[1],a[2]]
                        if hasattr(leftpixel, "id"):
                            if centrecell.id!=leftpixel.id:
                                neighbourpixel.append(leftpixel.type)
                        else:
                            neighbourpixel.append(0)

                        rightpixel = self.cellField[a[0]+1,a[1],a[2]]
                        if hasattr(rightpixel, "id"):
                            if centrecell.id!=rightpixel.id:
                                neighbourpixel.append(rightpixel.type)
                        else:
                            neighbourpixel.append(0)

                        uppixel = self.cellField[a[0],a[1]+1,a[2]]
                        if hasattr(uppixel, "id"):
                            if centrecell.id!=uppixel.id:
                                neighbourpixel.append(uppixel.type)
                        else:
                            neighbourpixel.append(0)

                        downpixel = self.cellField[a[0],a[1]-1,a[2]]
                        if hasattr(downpixel, "id"):
                            if centrecell.id!=downpixel.id:
                                neighbourpixel.append(downpixel.type)
                        else:
                            neighbourpixel.append(0)

                if typecell==2 and len(neighbourpixel)>0:
                    perimeterpercentage = (
                        neighbourpixel.count(1)
                        + neighbourpixel.count(3)
                        + neighbourpixel.count(5)
                        + neighbourpixel.count(7)
                    ) / float(len(neighbourpixel))

                    ### APOP_SCALE applied
                    if APOP_SCALE * (3*(0.000416*perimeterpercentage)) > random.random():
                        cell.type=10
                        cell.targetVolume=0
                        cell.lambdaVolume=2
                        deathcount+=1
                        apopscrb+=1

                elif typecell==1 and len(neighbourpixel)>0:
                    perimeterpercentage = (
                        neighbourpixel.count(2)
                        + neighbourpixel.count(4)
                        + neighbourpixel.count(6)
                        + neighbourpixel.count(8)
                    ) / float(len(neighbourpixel))

                    ### APOP_SCALE applied
                    if APOP_SCALE * (10*(0.0000416*perimeterpercentage)) > random.random():
                        cell.type=9
                        cell.targetVolume=0
                        cell.lambdaVolume=2
                        deathcount+=1
                        apopwt+=1

                time = (mcs-relaxtime)/float(10)
                P_apo.append([cell.id,cell.type,time,perimeterpercentage,deathcount])

            time = (mcs-relaxtime)/float(10)
            perimeterarray.append([cell.id,cell.type,time,perimeterpercentage])
            cumulitiveapop.append([time,apopscrb,apopwt])
    def finish(self):
        np.savetxt('C:\Users\Dgradeci\Dropbox\CC3D Workspace\data\perimeterarray-%s'%datetime.now().strftime('%H-%M-%m-%d'), perimeterarray)   
        np.savetxt('C:\Users\Dgradeci\Dropbox\CC3D Workspace\data\P_apo-%s'%datetime.now().strftime('%H-%M-%m-%d'),P_apo)
        np.savetxt('C:\Users\Dgradeci\Dropbox\CC3D Workspace\data\cumulitiveapop-%s'%datetime.now().strftime('%H-%M-%m-%d'),cumulitiveapop)
        return


# You can keep this class name so Competition.py imports don't break,
# but internally it uses the SAME safe contact-fraction approach.
class DeathSteppablePerimiter(SteppableBasePy):
    def __init__(self,_simulator,_frequency=10):
        SteppableBasePy.__init__(self,_simulator,_frequency)

    def start(self):
        global perimeterarray
        global P_apo
        global cumulitiveapop
        cumulitiveapop=[]
        perimeterarray=[]
        P_apo=[]

    def step(self,mcs):
        if mcs <= relaxtime:
            return

        apopwt = 0
        apopscrb = 0
        deathcount = 0

        for cell in self.cellList:
            if cell.type == 9 or cell.type == 10:
                continue
            if cell.volume < 50:
                continue

            is_scrb = (cell.type in (2,4,6,8))
            is_wt   = (cell.type in (1,3,5,7))
            if not (is_scrb or is_wt):
                continue

            total_contact = 0.0
            opp_contact = 0.0

            for neighbor, commonSurfaceArea in self.getCellNeighborDataList(cell):
                if not neighbor:
                    continue
                if neighbor.type == 9 or neighbor.type == 10:
                    continue

                total_contact += commonSurfaceArea

                if is_scrb and (neighbor.type in (1,3,5,7)):
                    opp_contact += commonSurfaceArea
                elif is_wt and (neighbor.type in (2,4,6,8)):
                    opp_contact += commonSurfaceArea

            if total_contact > 0.0:
                perimeterpercentage = opp_contact / total_contact
            else:
                perimeterpercentage = 0.0

            # same competition rules as your original perimeter steppable
            if is_scrb:
                if 3.0 * (0.000416 * perimeterpercentage) > random.random():
                    cell.type = 10
                    cell.dict["apt"] = mcs
                    cell.targetVolume = 0
                    cell.lambdaVolume = 2
                    deathcount += 1
                    apopscrb += 1

            elif is_wt:
                if 10.0 * (0.0000416 * perimeterpercentage) > random.random():
                    cell.type = 9
                    cell.dict["apt"] = mcs
                    cell.targetVolume = 0
                    cell.lambdaVolume = 2
                    deathcount += 1
                    apopwt += 1

            time = (mcs - relaxtime) / float(10)
            P_apo.append([cell.id, cell.type, time, perimeterpercentage, deathcount])
            perimeterarray.append([cell.id, cell.type, time, perimeterpercentage])

        time = (mcs - relaxtime) / float(10)
        cumulitiveapop.append([time, apopscrb, apopwt])

    def finish(self):
        np.savetxt('C:\\Users\\Dgradeci\\Dropbox\\CC3D Workspace\\data\\perimeterarray-%s'
                   % datetime.now().strftime('%H-%M-%m-%d'), perimeterarray)
        np.savetxt('C:\\Users\\Dgradeci\\Dropbox\\CC3D Workspace\\data\\P_apo-%s'
                   % datetime.now().strftime('%H-%M-%m-%d'), P_apo)
        np.savetxt('C:\\Users\\Dgradeci\\Dropbox\\CC3D Workspace\\data\\cumulitiveapop-%s'
                   % datetime.now().strftime('%H-%M-%m-%d'), cumulitiveapop)
        return




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

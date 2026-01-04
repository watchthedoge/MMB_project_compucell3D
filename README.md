# MMB_project_compucell3D
This repository contains a **CompuCell3D (CC3D)** implementation of cell competition models used to study **biochemical and mechanical competition**, including extensions to **three competing cell types** and **apoptosis-dependent dynamics**.

The project reproduces and extends published cell-competition models and was developed for a modeling methods course project.


Based on the work of  Gradeci D, Bove A, Vallardi G, Lowe AR, Banerjee S, Charras G. Cell-
scale biophysical determinants of cell competition in epithelia. eLife.
2021;10:e61011. doi:10.7554/eLife.61011
---

## 1. Software Requirements

This project was developed and tested with:

- **CompuCell3D 3.7.x – 3.9.x (64-bit)**
- **Python 2.7** (bundled with CC3D)
- **Operating System**: Windows 11  
  (Linux should work)


⚠️ **Important**  
Do **not** use Python 3. This code relies on CC3D’s Python 2.7 API.

---


## 2. File Roles 
### **Competition.xml**
Defines the **physical model**:
- Lattice size and dimensionality
- Cell types and IDs
- Plugins (Volume, Surface, NeighborTracker, ExternalPotential, etc.)
- Contact energies (mechanical interactions)


---

### **Competition.piff**
Defines the **initial condition**:
- Cell IDs, types, positions, and sizes
- Used to seed winner/loser/killer configurations
- Controls spatial setup (clusters, isolated cells, invasion fronts)


---

### **Competition.py**
Acts as the **simulation entry point**:
- Initializes the CC3D simulation
- Registers steppables
- Selects which biological mechanisms are active (growth, apoptosis, motility, etc.)

---

### **CompetitionSteppables.py**
Contains **all biological rules**, implemented as steppables:
- Growth and target volume dynamics
- Mitosis and inheritance
- Apoptosis (density-based, perimeter-based, scaled by `APOP_SCALE`)
- Mechanical motility
- Cell tracking and cleanup of apoptotic cells


---

### **ParameterScanSpecs.xml**
(Optional) Defines **parameter sweeps**:
- Used for batch simulations
- Can vary growth rates, stiffness, apoptosis scale, etc.

Not required for single simulations.

---

## 3. Running the Simulation

1. Open **CompuCell3D Player**
2. Load `Competition.cc3d` (put in the same directory as the "Simulation" folder from one of the experiments)
3. Click **Run**
4. View the lattice in real time

---

## 4. Key Model Parameters

These are typically defined in `CompetitionSteppables.py`:

- `growthratewt`, `growthratescrb`, `growthrate_killer`
- `stiffness_wt`, `stiffness_kd`, `stiffness_killer`
- `APOP_SCALE`
  - `0` → no apoptosis
  - `0.1` → weak apoptosis
  - `1` → strong apoptosis
- `relaxtime` → mechanical equilibration time before biology activates

---

## 5. Output

During simulation, the model can generate:
- Cell counts over time
- Apoptosis events
- Cell tracking data

Outputs are written to the `data/` directory (paths may need adjustment on non-Windows systems).

---

## 6. Known Limitations

- Python 2.7 only
- Long simulations can crash if apoptosis cleanup is disabled
- Very dense tissues may jam and suppress division
- Results are qualitative and model-dependent

---

## 7. Reproducibility Notes

- Results depend on random seeds (not fixed by default)
- Exact lattice size and plugin settings matter
- Apoptotic cells must be removed to avoid lattice saturation

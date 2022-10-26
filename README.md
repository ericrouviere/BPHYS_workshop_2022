# Minimal Modeling and Simulation
BPHYS workshop, Fall 2022, Eric Rouviere

---
## Contents
1. `notebooks/` contains the jupyter notebooks we will be working through.
2. `src/` contains the source code of much of the workshop. 
---
## Installation of software

### **1. Install Julia**
1. Check type of CPU.
2. Download from https://julialang.org/downloads/
3. Install by clicking on the downloaded file and dragging to applications.

### **4. Installing the nessicary Julia packages**
1. Start the julia app. This should open up the terminal with julia running. 
3. Type `]`. This gets you to the package manager
4. Enter `add PyPlot, LaTeXStrings, BenchmarkTools`
5. press the delete key to go back to julia
6. Enter `using PyPlot, LaTeXStrings, BenchmarkTools` 

### **2. Install text and notebook editor**
**Option 1: VS Code**
1. Download from https://code.visualstudio.com/download
2. Install by clicking on the downloaded file and dragging to applications.
3. Start VS Code.
4. Inside VS Code, go to the Extensions view by clicking View on the top menu bar and then selecting Extensions.
5. In the Extensions view, search for the term "julia" in the Marketplace search box, then select the Julia extension (julialang.language-julia) and select the Install button.
6. Install the "Jupyter" extension in the same way.
7. Quit VS Code

**Option 2: Any editor and IJulia**
1. Install your favorite editor (I leave this up to you. VS code works.)
2. In a julia terminal enter `] add IJulia`
3. Press the delete key to get back to julia
4. Enter `using IJulia`
5. When prompted to install conda, agree. 
6. Enter `notebook()` to start Jupyter. 


### **5. Installing git**
1. Got to https://git-scm.com/
2. click "Download for XXX" and follow installer instructions.
3. In the terminal, move to a directory where you want to put the repo (e.g. `cd ~/Desktop`).
4. Enter `git clone git@github.com:ericrouviere/BPHYS_workshop_2022.git`

To keep your up repo up to date with changes I make, cd into `BPHYS_workshop_2022/` and enter `git pull`.


### **6. Restart computer**
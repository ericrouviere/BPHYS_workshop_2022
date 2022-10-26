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
4. Add julia to your PATH by adding the following line to the `~/.bash_profile` file.

    ```export PATH="$PATH:/Applications/Julia-1.8.app/Contents/Resources/julia/bin/"```
5. In the command line enter `source ~/.bash_profile`

### **2. Installing the nessicary Julia packages**
1. Start the julia app. This should open up the terminal with julia running. 
3. Type `]`. This gets you to the package manager
4. Enter `add PyPlot, IJulia, LaTeXStrings, BenchmarkTools`
5. Press the delete key to go back to julia.
6. Enter `using PyPlot, IJulia, LaTeXStrings, BenchmarkTools`.
7. Enter `jupyterlab()` to start Jupyter IDE. When prompted to install conda, agree. 

### **4. Installing git**
1. Got to https://git-scm.com/
2. click "Download for XXX" and follow installer instructions.
3. In the terminal, move to a directory where you want to put the repo (e.g. `cd ~/Desktop`).
4. Enter `git clone https://github.com/ericrouviere/BPHYS_workshop_2022.git`

To keep your up repo up to date with changes I make, cd into `BPHYS_workshop_2022/` and enter `git pull`
# Deployment guide for GUI verification

## Step 1: Create a new conda environment

We’ll use Python 3 (recommended):

```bash
conda create -n smg2s-gui python=3.10
conda activate smg2s-gui
```

## Step 2: Install required Python packages

Install the additional modules and dependencies:

```bash
conda install numpy scipy matplotlib pandas tabulate
```

## Step 3: GUI dependencies

Since the program uses matplotlib with GUI backends, it needs Tkinter.

Install Tk support in conda:
```bash
conda install tk
```

On Linux (Ubuntu/Debian), if Tk is still missing, also run:
```bash
sudo apt-get update
sudo apt-get install python3-tk
```
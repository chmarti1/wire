#!/usr/bin/python3
import wire
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import os

# This is an example script that uses the wire module to create more
# sophisticated plots.  Since plotting publication quality plots is 
# such an iterative process, the intent is for users to be able to edit
# this script as needed to achieve the results they want.  Parameters
# that will need to be edited are commented with a ===> above them to
# draw attention.

# ===> List the source coefficient files and their meta data
# ===> The STASH directory is a kind of working directory where the results
# ===> of prior WSOLVE runs can be "stashed".  
stashdir = 'stash/run2'
files = ['000.wcf', '001.wcf', '002.wcf', '003.wcf', '004.wcf']
titles = ['z = 1.0mm', 'z = 1.5mm', 'z = 2.0mm', 'z = 2.5mm', 'z = 3.0mm']
unit_length = 'mm'

# ===>  Create the plot with subplots
# ===>  If users want to change the number of plots, they will also want to 
# ===>  change the subplot arrangement and figure size
nrows = 3           # Number of rows in the subplot array
ncols = 2           # Number of columns in the subplot array
figsize = (6.5,9)   # Width, Height (in)
fig,axs = plt.subplots(nrows,ncols,sharex=True,sharey=True,figsize=figsize)

# ===>  These options deal with how the plot appears
# ===>  The current measurement is negative, but many colormaps produce
# ===>  A more meaningful result if the maximum is positive.  These options
# ===>  allow you to invert the signal or choose a different color map
# ===>  See https://matplotlib.org/stable/tutorials/colors/colormaps.html
invert = True       # Should the negative current be plotted as positive signal?
cmap = 'inferno'


# Next, we'll load all of the image data from the files listed
# Users probably won't need to edit this part at all
row = 0
col = 0
vmin = float('inf')
vmax = float('-inf')
pcms = []
index = 0
for fn,title in zip(files, titles):
    filename = os.path.join(stashdir, fn)
    
    # Get the row and column
    row = index // ncols
    col = index % ncols
    
    # Load the file and generate x,y,I data
    x,y,I = wire.WireCoefficients(filename).grid()
    # If we're inverting the colors axis, do it now
    if invert:
        I = -I
    
    # Detect the data range - we'll use these later
    vmin = min(vmin, I.min())
    vmax = max(vmax, I.max())
    
    ax = axs[row,col]
    
    # Generate the plot
    # Stash the pcolormesh for later
    pcms.append(ax.pcolormesh(x,y,I,cmap=cmap))
    # Apply the label
    if title:
        ax.text(x.min()+1,y.max()-1,title,color='w',ha='left',va='bottom',zorder=1000)
        
        # ax.set_title(title)
    # Clear the inner ticks and apply a unit
    
    index += 1

# Remove the remaining axes
for index in range(index,ncols*nrows):
    row = index // ncols
    col = index % ncols
    axs[row,col].remove()

# Clip the max/min based on the inversion setting    
if invert:
    vmin = 0
    print('scale:', vmax)
else:
    vmax = 0
    print('scale:', vmin)
# Re-normalize all of the plots to share the same color scale
norm = colors.Normalize(vmin,vmax, clip=True)
for pcm in pcms:
    pcm.set(norm=norm)
    
# Set unit length labels
for row in range(nrows):
    axs[row,0].set_ylabel(unit_length)
for col in range(ncols):
    axs[nrows-1,col].set_xlabel(unit_length)
    
# Make everything snug
fig.tight_layout()

# Finish
filename = os.path.join(stashdir,'myplot.png')
fig.savefig(filename)

#!/usr/bin/env python

'''
   Created: 03//2022: First attempt.
   Author : Matthew Mumpower {matthew@mumpower.net}
   Description: Common analysis routines for PRISM output.
   Intended: Python 2.7.9
   Modified: 12/21/2018: Nicole Vassh {nvassh@gmail.com}
   Modified: 06//2023: Lauren Harewood {harewol@gmail.com}
'''

# ========================

# Matplotlib
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import MaxNLocator
from matplotlib.mlab import griddata
import matplotlib.gridspec as gridspec
from matplotlib.widgets import Slider
import os
import os.path


# Local Namespace
from output_lauren import AbTime #to read PRISM time output
from nzplane import drawBox,drawBoxes

#from SnContour_simp_4Lauren import plot
# ========================
abmin = 10**(-6)

blue = (float(64)/float(255),float(107)/float(255),float(255)/float(255)) # blue
lightblue = (float(64)/float(255),float(202)/float(255),float(255)/float(255)) #light blue
teal = (float(100)/float(255),float(255)/float(255),float(213)/float(255)) # teal
lightgreen = (float(135)/float(255),float(255)/float(255),float(128)/float(255)) # light green
greenyellow = (float(196)/float(255),float(255)/float(255),float(54)/float(255)) # green yellow
yellow = (float(255)/float(255),float(219)/float(255),float(64)/float(255)) # yellow
lightorange= (float(255)/float(255),float(166)/float(255),float(64)/float(255)) # light orange
orange = (float(255)/float(255),float(108)/float(255),float(64)/float(255)) # orange
red = (float(204)/float(255),float(21)/float(255),float(10)/float(255)) # red
darkred = (float(128)/float(255),0,0) # dark red

# Generate colormapping
colorarr = [blue,lightblue,teal,lightgreen,greenyellow,yellow,lightorange,orange,red,darkred]
cmap = colors.ListedColormap(colorarr)
bounds = [abmin,abmin*10,abmin*100,abmin*1000,abmin*(10**4),abmin*(10**5),abmin*(10**6),abmin*(10**7),abmin*(10**8),abmin*(10**9)]
boundsl = ['$10^{%s}$'%int(np.log10(abmin)),'$10^{%s}$'%int(np.log10(abmin*10.)),'$10^{%s}$'%int(np.log10(abmin*100.)),'$10^{%s}$'%int(np.log10(abmin*(10**3))),'$10^{%s}$'%int(np.log10(abmin*(10**4))),'$10^{%s}$'%int(np.log10(abmin*(10**5))),'$10^{%s}$'%int(np.log10(abmin*(10**6))),'$10^{%s}$'%int(np.log10(abmin*(10**7))),'$10^{%s}$'%int(np.log10(abmin*(10**8))),'$10^{%s}$'%int(np.log10(abmin*(10**9)))]
norm = colors.BoundaryNorm(bounds,cmap.N)

cmapBlack = mpl.colors.ListedColormap([(0,0,0)],name='M-Black')
cmapGreen = mpl.colors.ListedColormap([(0,1,0)],name='N-Green')
cmapRed = mpl.colors.ListedColormap([(1,0,0)],name='N-Red')

def cfnc(z):
  if z > abmin and z <= abmin*10.:
    c1 = colorarr[0]
  elif z > abmin*10. and z <= abmin*100.:
    c1 = colorarr[1]
  elif z > abmin*100. and z <= abmin*1000.:
    c1 = colorarr[2]
  elif z > abmin*1000.: 
#  and z <= abmin*(10**4):
    c1 = colorarr[3]
    '''
  elif z > abmin*(10**4) and z <= abmin*(10**5):
   c1 = colorarr[4]
  elif z > abmin*(10**5) and z <= abmin*(10**6):
    c1 = colorarr[5]
  elif z > abmin*(10**6) and z <= abmin*(10**7):
    c1 = colorarr[6]
  elif z > abmin*(10**7) and z <= abmin*(10**8):
    c1 = colorarr[7]
  elif z > abmin*(10**8) and z <= abmin*(10**9):
    c1 = colorarr[8]
  elif z > abmin*(10**9):
    c1 = colorarr[9]
    '''
  else:
    c1 = "none"
  return c1


class NZPlaneSliderFigure(object):
  '''
     Draw abundance population the NZ-plane and Y(A) controlled by a slider.
  '''

  def __init__(self, t_file=None, path_file=None, scale=1., **kw):
    '''
       Initialize the file and read the data
    '''
    self.time = AbTime(t_file)

    # Load solar data
    sol = open("arnould2007_ssdata.dat",'r')
    #sol = open("/Users/nicolevassh/research/reverseengineering_old/masses/dzmp/data/solar_sne08_EB0075.dat",'r')
    
    self.solarA, self.solarY, self.solarYl, self.solarYh = [],[],[],[]
    # Skip header line
    sol.readline()
    # Read the rest of file
    for line in sol:
      line = line.split()
      Z = int(line[0])
      A = int(line[1])
      Y  = float(line[2])
      Yl = Y - float(line[3])
      Yh = float(line[4]) - Y

      self.solarA.append(A)
      self.solarY.append(Y/scale)
      self.solarYl.append(Yl/scale)
      self.solarYh.append(Yh/scale)

    sol.close()


    # Construct figure
    self.fig1 = plt.figure(num=None,figsize=(14,14),facecolor='w',edgecolor='k')

    # Set Axes & Resize Subplots
    self.gs = gridspec.GridSpec(1,1)#,height_ratios=[0.05,1])
    #self.ax = plt.subplot(self.gs[0,0])
    self.bx = plt.subplot(self.gs[0,0])
   # self.cx = plt.subplot(self.gs[1,0]) # Slider

    # Set NZ-plane bounds (fixed for now)
    xmin = 94
    xmax = 114
    ymin = 46
    ymax = 68

    # Set Axis Labels
    #self.ax.set_ylabel('Y(A)')
    #self.ax.set_xlabel('A')
    #self.ax.set_yscale('log')
    self.bx.set_ylabel('Z',fontsize=20)
    self.bx.set_xlabel('N',fontsize=20)

    # Set Axis Bounds
    #self.ax.axis([120,205,0.0000001,0.001])
    self.bx.axis([xmin+0.5,xmax+0.5,ymin-.5,ymax+0.5])
    self.bx.tick_params(axis='both',which='major', labelsize=14)
    self.bx.xaxis.set_major_locator(MaxNLocator(integer=True))
  def gen_movie_images(self):
    '''
       Generate images for the movie.
    '''
    #for m in sorted(self.time.data.keys(), reverse=True):
    print 'it read this'
    for m in self.time.data.keys():
      print("Frame: %s of %s" % (m,len(self.time.data)))
      fn = "new_xft_images/withContour/nlod/03/ee/betadCont/highres/pic%04d.png" % m
      if not os.path.isfile(fn):
        self.update_movie(m)
        plt.tight_layout()
        plt.savefig(fn,dpi=400)

  def th_info(self, m):
    '''
      Outputs time, temp, density
    '''
    self.t.set_text("t = %.6e" % (self.time.data[m].t))
    self.t9.set_text("t9 = %.6e" % (self.time.data[m].T9))
    self.rho.set_text("rho = %.6e" % (self.time.data[m].rho))


  def update_movie(self, m):
    '''
       Dynamically update the figure based off current timestep (m)
    '''
    #self.plotYa(m)
    self.reDrawY(m)
    #self.th_info(m)
    
  def plotYa(self, m):
    '''
       Plot Y(A) in the top panel.
    '''
    try:
      A = self.time.data[m].Ya.keys()
      Y = self.time.data[m].Ya.values()
      #self.Yline.set_xdata(A)
      #self.Yline.set_ydata(Y)
    except:
      print 'step not present'
  

  def reDrawY(self, m):
    '''
       Redraw each box in the bottom panel.
    '''
    for l in self.NZboxes:
      l.remove()
    keys = self.time.data[m].Y.keys()
    
    N = [n[0] for n in keys]
    Z = [n[1] for n in keys]
    
    Y = [self.time.data[m].Y[n,z] for (n,z) in keys] #if self.time.data[m].Y[n,z] >= 1e-6]
    ''' 
    N = [n for (n,z) in keys if self.time.data[m].Y[n,z] >= 1e-6]
    Z = [n for (n,z) in keys if self.time.data[m].Y[n,z] >= 1e-6]
    '''
    self.NZboxes = []
    for z,n,y in zip(Z,N,Y):
      self.NZboxes.append(drawBox(ax=self.bx, y=z, x=n, fc=cfnc(y), ec='none', zorder=5))

  def plot(self):
    '''
       Generate the plot on the screen.
    '''


    # Top Panel -----
    # Plot solar data
   # self.ax.errorbar(self.solarA, self.solarY, yerr=[self.solarYl, self.solarYh],color='k',fmt='.', ecolor='k', capsize=2.2, capthick=1.2, elinewidth=1.2)

    # Initialize the Y-line
    #self.Yline = self.ax.plot(0, 0.000001, ls='-',color='red')[0]
    #self.plotYa(3)

    # Bottom Panel -----
    # Plot initial abundances
    self.NZboxes = []
    self.reDrawY(3)
    
    self.t = self.bx.text(168,0.002,"t = %.6e" % (self.time.data[3].t))
    self.t9 = self.bx.text(168,0.001,"t9 = %.6e" % (self.time.data[3].T9))
    self.rho = self.bx.text(168,0.0005,"rho = %.6e" % (self.time.data[3].rho))
    
    ax2 = make_axes_locatable(self.bx)
    cax = ax2.append_axes("right",size="3%",pad=0.05)
    # Set colorbar
    cb1 = mpl.colorbar.ColorbarBase(cax,cmap=cmap,norm=norm,ticks=bounds,boundaries=bounds,spacing='uniform',orientation='vertical')
    cb1.set_label('Abundance',**{'size':20,'labelpad':6})
    cb1.set_ticklabels(boundsl)
    #cb1.cax.tick_params(labelsize=14)

    # Plot stable isotopes
    stable = open("/afs/crc.nd.edu/group/tna-shared/NDI/data/nubase2016/stable_nubase2016.dat",'r')
    for line in stable:
      line = line.split()
      z = int(line[0])
      a = int(line[1])
      drawBox(self.bx, y=z, x=a-z, fc='k', ec='black')
    stable.close()

    self.gen_movie_images()

    #plt.tight_layout()
   #plt.show()


#fig = NZPlaneSliderFigure("/Users/nicolevassh/Dowloads/YVT-10GK5GK-000032.dat", scale=2000.)
#fig = NZPlaneSliderFigure("output/newtraj/yesReacLib_ytime_04.txt", scale=2000.)
fig = NZPlaneSliderFigure("output/newtraj/xft_new_ytime_nlod_03_gef.txt", scale=2000.)
#fig = NZPlaneSliderFigure("output/newtraj/newNubaseFiles/output_gef_ytime_07.txt", scale=2000.)
#fig = NZPlaneSliderFigure("output/newtraj/ReacLibONLY_ytime_01.txt", scale=2000.)

#below will find the region within the FRDM2012 one neutron dripline
#ext = open("/afs/crc.nd.edu/group/tna-shared/PRISM/prism/input/extent/full_reaclib.dat",'r')
ext = open("ame2020_ext.dat",'r')
#frdm2012 = open("input/nuclear/NDI/s1n_frdm2012.dat",'r')
frdm2012 = open("input/nuclear/NDI/s1n_xft_new_nlod_me.dat",'r')
frdm2012_bmd = open("input/nuclear/NDI/bmd_beoh350-gtff-sdn_frdm2012.bin.dat")
#frdm2012_bmd = open("input/nuclear/NDI/bmd_beoh350-mkt_frdm2012.bin.dat")
#frdm2012_bmd = open("input/nuclear/NDI/bmd_beoh350-gtff-nesext_frdm2012.bin.dat")

innet = {}
for line in ext:
  line = line.split()
  Z = int(line[0])
  Nmin = int(line[1])
  Nmax = int(line[2])
  for N in range(Nmin,Nmax+1):
    innet[Z,N] = 1
ext.close()

ratez=[]
bmdr, bmd_nuc ={},{}
numnuc = int(frdm2012_bmd.readline().split()[0])
for i in range(0,numnuc):
    Z = int(frdm2012_bmd.readline().split()[0])
    N = int(frdm2012_bmd.readline().split()[0])
    bmd_nuc[Z,N] = 1
    frdm2012_bmd.readline()
    frdm2012_bmd.readline()
    rate = float(frdm2012_bmd.readline().split()[0])
    ratez.append(rate)
    frdm2012_bmd.readline()
    if str(rate)=="NaN":
        rate=0
    else:
        bmdr[Z,N] = rate

#weighted_rate = [r *p for r,p in zip(rate,


nuctype = 'ee'
print 'nuctype', nuctype
snlevels = [0.0,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2,2.25,2.5,2.75,3,3.25,3.5,3.75,4,4.25,4.5,4.75,5,5.25,5.5,5.75,6,6.25,6.5]
#bmd_levels = np.linspace(np.nanmin(min(zi)),np.log10(max(ratez)))

bzee,bnee,bee = [],[],[]
bzeo,bneo,beo = [],[],[]
bzoe,bnoe,boe = [],[],[]
bzoo,bnoo,boo = [],[],[]
z,n,bmd = [],[],[]

for (Z,N) in bmdr.keys():
    if N % 2 == 0 and Z % 2 == 0:
        bee.append(bmdr[Z,N])
        bnee.append(N)
        bzee.append(Z)
    elif N % 2 == 0 and Z % 2 != 0:
        beo.append(bmdr[Z,N])
        bneo.append(N)
        bzeo.append(Z)
    elif N % 2 != 0 and Z % 2 == 0:
        boe.append(bmdr[Z,N])
        bnoe.append(N)
        bzoe.append(Z)
    elif N % 2 != 0 and Z % 2 != 0:
        boo.append(bmdr[Z,N])
        bnoo.append(N)
        bzoo.append(Z)
print 'bee', bee, 'bzee', bzee, 'bnee'
print 'beo',beo,'bzeo',bzeo,'bneo'
print 'boo', boo,'bzoo', 'bnoo',bnoo
print 'boe',boe, 'bzoe', bzoe,'bnoe', bnoe


if nuctype == 'ee':
    bmd = bee
    z = bzee
    n = bnee
elif nuctype == 'eo': #evenN-oddZ
    bmd = beo
    z = bzeo
    n = bneo
elif nuctype == 'oe': #oddN- evenZ
    bmd = boe
    z = bzoe
    n = bnoe
elif nuctype == 'oo':
    bmd = boo
    z = bzoo
    n = bnoo
else:
    print 'nuctype not defined'



k = np.linspace(min(n), max(n),num=1+int((max(n)-min(n))/2))
l = np.linspace(min(z), max(z),num=1+int((max(z)-min(z))/2))
xi = k
yi = l
zi = griddata(n,z,bmd,xi,yi,interp='linear')

bmd_valid = np.array([v for v in bmd if v>0.0])
bmd_levels = np.logspace(np.log10(np.min(bmd_valid)),np.log10(np.max(bmd_valid)))
print 'zi shape:', np.shape(zi)
print 'zi min/max:', np.nanmin(zi), np.nanmax(zi)

'''
points = np.array(zip(n,z))
values = np.array(bmd)
grid_z = griddata(points, values, (k[None,:], l[:,None]), method='linear')
'''
fig.bx.contour(xi, yi, zi, levels=bmd_levels, cmap=cmapRed, zorder=200,lw=10, alpha=.2)

x,y,z=[],[],[]
sepn,sepnpos,frdmnuc,frdmdrip = {}, {},{},{}
nnuc = int(frdm2012.readline().split()[0])
for i in range(0,nnuc):
  Z = int(frdm2012.readline().split()[0])
  N = int(frdm2012.readline().split()[0])
  frdmnuc[Z,N] = 1
  Zd = frdm2012.readline()
  Nd = frdm2012.readline()
  sn = float(frdm2012.readline().split()[0])


  if sn == "NaN":
      sn=0

  elif sn>0:
      sepnpos[Z,N]=sn
      frdmdrip[Z,N] = 1
      sepn[Z,N]=sn
  else:
      sepn[Z,N]=sn



zee,nee,snee = [],[],[]
zeo,neo,sneo = [],[],[]
zoe,noe,snoe = [],[],[]
zoo,noo,snoo = [],[],[]
z,n,sn = [],[],[]




for (Z,N) in sepn.keys():
    print 'Z',Z ,'N',N
    print 'N%2', N%2, 'Z%2', Z%2
    if N % 2 == 0 and Z % 2 == 0:
        snee.append(sepn[Z,N])
        nee.append(N)
        zee.append(Z)
    elif N % 2 == 0 and Z % 2 != 0:
        sneo.append(sepn[Z,N])
        neo.append(N)
        zeo.append(Z)
    elif N % 2 != 0 and Z % 2 == 0:
        snoe.append(sepn[Z,N])
        noe.append(N)
        zoe.append(Z)
    elif N % 2 != 0 and Z % 2 != 0:
        snoo.append(sepn[Z,N])
        noo.append(N)
        zoo.append(Z)
print 'snee', snee, 'zee', zee, 'nee'
print 'sneo',sneo,'zeo',zeo,'neo'
print 'snoo', snoo,'zoo', 'noo',noo
print 'snoe',snoe, 'zoe', zoe,'noe', noe
if nuctype == 'ee':
    sn = snee
    z = zee
    n = nee
elif nuctype == 'eo': #evenN-oddZ
    sn = sneo
    z = zeo
    n = neo
elif nuctype == 'oe': #oddN- evenZ
    sn = snoe
    z = zoe
    n = noe
elif nuctype == 'oo':
    sn = snoo
    z = zoo
    n = noo
else:
    print 'nuctype not defined'
x = np.linspace(min(n), max(n),num=1+int((max(n)-min(n))/2))
y = np.linspace(min(z), max(z),num=1+int((max(z)-min(z))/2))

z = griddata(n, z, sn, x, y, interp='linear')
fig.bx.contour(x, y, z, levels=snlevels, cmap=cmapBlack, zorder=100,lw=10, alpha=.2)

'''
for (Z,N) in sepnpos.keys():
    drawboxes(ax=fig.bx,x=N,y=Z,fc='light grey', ec='none',union=True,zorder=-5)
    '''
frdm2012.close()

frdm2012_bmd.close()



m1 = {}
for (Z,N) in innet:
  if frdmnuc.get((Z,N),-1000)==-1000 and N<20:
    m1[Z,N] = 1
  else:
    if frdmdrip.get((Z,N),-1000)!=-1000:
      m1[Z,N] = 1

Na,Za=[],[]
for (Z,N) in sorted(m1):
  Na.append(N)
  Za.append(Z)



# Plot mass data
drawBoxes(ax=fig.bx, N=Na, Z=Za, fc='#dedee2', ec='#bbbbc0', union=True, zorder=-5)
# Plot mass data
#drawBoxes(ax=fig_1.bx, N_1=Na_1, Z_1=Za_1, fc='#dedee2', ec='#bbbbc0', union=True, zorder=-5)




#fig_1.plot()
fig.plot()

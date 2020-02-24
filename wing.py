import numpy as np
import numpy.linalg as lg
import scipy.linalg as slg
import scipy.linalg.lapack as slapack
import scipy.sparse.linalg as splg
from math import *
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import scipy.sparse as sps
import time as tm

import toolkit
from utils import *

class wing:
    def __init__(self, sld): #wing initializer
        self.sld=sld
    def patchcompose(self, afls=['n4412'], chords=[1.0], ypos=[0.0, 1.0], xpos=[0.0], incidences=[0.0], \
        gammas=[0.0], xdisc=20, ydisc=20, yspacing=np.array([]), xstrat=lambda x: (np.sin(pi*x-pi/2)+1)/2, zbase=0.0, \
            prevlines={}):
        #patch composer based on wing's sectional parameters. The wing is meant to be a guideline for other abuted
        #surfaces, therefore this function does not accept prevlines dictionary as input.
        nsects=len(ypos)
        afls=trimlist(nsects, afls)
        chords=trimlist(nsects, chords)
        xpos=trimlist(nsects, xpos)
        incidences=trimlist(nsects, incidences)
        gammas=trimlist(nsects, gammas)
        if len(yspacing)==0:
            yspacing=np.linspace(ypos[0], ypos[-1], 20)
        
        zpos=[zbase]
        for i in range(1, nsects):
            zpos+=[zpos[-1]+tan(gammas[i])*(ypos[i]-ypos[i-1])]
        zspacing=np.interp(yspacing, np.array(ypos), np.array(zpos))
        xspacing=np.interp(yspacing, np.array(ypos), np.array(xpos))
        gspacing=np.interp(yspacing, np.array(ypos), np.array(gammas))
        cspacing=np.interp(yspacing, np.array(ypos), np.array(chords))

        afl=np.zeros((nsects, 2*xdisc-1, 2))
        discafl=np.zeros((len(yspacing), 2*xdisc-1, 2))
        for i in range(nsects):
            afl[i, :, :]=read_afl(afls[i], ext_append=True, header_lines=1, disc=xdisc, \
                strategy=xstrat, remove_TE_gap=True, extra_intra=False, \
                    incidence=incidences[i], inverse=False)
        for i in range(np.size(afl, 1)):
            discafl[:, i, 0]=np.interp(yspacing, np.array(ypos), afl[:, i, 0])
            discafl[:, i, 1]=np.interp(yspacing, np.array(ypos), afl[:, i, 1])
        
        pointpos=np.zeros((len(yspacing), np.size(afl, 1), 3))
        for i in range(len(yspacing)):
            pointpos[i, :, :]=wing_afl_positprocess(discafl[i, :, :], gamma=gspacing[i], xpos=xspacing[i], \
                ypos=yspacing[i], zpos=zspacing[i], c=cspacing[i])
        pts=[]
        for i in range(np.size(afl, 1)):
            pts+=[[]]
            for j in range(len(yspacing)-1, -1, -1):
                pts[-1]+=[pointpos[j, i, :]]
        
        horzlines, vertlines, paninds, sldpts=self.sld.addpatch(pts, [1], prevlines=prevlines)
        self.paninds=paninds
        self.horzlines=horzlines
        self.vertlines=vertlines
        wakecombs=[[paninds[0][i], paninds[-1][i]] for i in range(len(paninds[0]))]
        return horzlines, vertlines, paninds, sldpts, wakecombs
    def plotpressure(self):
        for i in range(len(self.paninds[0])):
            Cps=[]
            xs=[]
            for j in range(len(self.paninds)):
                Cps+=[self.sld.Cps[self.paninds[j][i]]]
                xs+=[self.sld.panels[self.paninds[j][i]].colpoint[0]]
            cloc=max(xs)
            xs=[x/cloc for x in xs]
            plt.plot(xs, Cps)
        plt.xlabel('x/c')
        plt.ylabel('Cp')
        plt.show()
    def plotpressure3D(self, xlim=[], ylim=[], zlim=[]):
        fig=plt.figure()
        ax=plt.axes(projection='3d')

        xsl=[]
        ysl=[]
        Cpsl=[]
        xsu=[]
        ysu=[]
        Cpsu=[]

        panindsl=[]
        panindsu=[]

        if self.sld.solavailable:
            for panlist in self.paninds:
                for p in panlist:
                    if self.sld.panels[p].nvector[2]<0:
                        xsl+=[self.sld.panels[p].colpoint[0]]
                        ysl+=[self.sld.panels[p].colpoint[1]]
                        Cpsl+=[self.sld.Cps[p]]
                        panindsl+=[p]
                    else:
                        xsu+=[self.sld.panels[p].colpoint[0]]
                        ysu+=[self.sld.panels[p].colpoint[1]]
                        Cpsu+=[self.sld.Cps[p]]
                        panindsu+=[p]
            ax.scatter3D(xsl, ysl, Cpsl, 'blue')
            ax.scatter3D(xsu, ysu, Cpsu, 'red')
        if len(xlim)!=0:
            ax.set_xlim3d(xlim[0], xlim[1])
        if len(ylim)!=0:
            ax.set_ylim3d(ylim[0], ylim[1])
        if len(zlim)!=0:
            ax.set_zlim3d(zlim[0], zlim[1])
        plt.xlabel('x')
        plt.ylabel('y')
        plt.show()
    def plotpatch(self, xlim=[], ylim=[], zlim=[]):
        fig=plt.figure()
        ax=plt.axes(projection='3d')

        for hllist in self.horzlines:
            for l in hllist:
                ax.plot3D(self.sld.lines[l, 0, :], self.sld.lines[l, 1, :], self.sld.lines[l, 2, :], 'gray')
        for vllist in self.vertlines:
            for l in vllist:
                ax.plot3D(self.sld.lines[l, 0, :], self.sld.lines[l, 1, :], self.sld.lines[l, 2, :], 'gray')
        if len(xlim)!=0:
            ax.set_xlim3d(xlim[0], xlim[1])
        if len(ylim)!=0:
            ax.set_ylim3d(ylim[0], ylim[1])
        if len(zlim)!=0:
            ax.set_zlim3d(zlim[0], zlim[1])
        plt.xlabel('x')
        plt.ylabel('y')
        plt.show()
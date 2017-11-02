# -*- coding: utf-8 -*-
"""
Created on Fri Mar 24 13:38:58 2017

@author: Jan Maka
"""
from matplotlib.dates import epoch2num
from matplotlib.patches import Rectangle
import numpy as np
import matplotlib.pyplot as plt
import os
import scipy.stats as st

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def plotfsec(self,fols, ops, to_file = True):
    x = np.arange(1,ts+1)
    plt.stem(x,self.fols/np.sum(self.fols)*100)
    plt.axis([0,max_ts+1,0,15])
    plt.xlabel('1 s time bin')
    plt.ylabel('Percetage of positive f')
    plt.savefig('../Results/rdn_fols%s'%max_ts+'.png')
    plt.show()
    plt.stem(x,self.ops/np.sum(self.ops)*100)
    plt.axis([0,max_ts+1,0,15])
    plt.xlabel('1 s time bin')
    plt.ylabel('Percetage of negative f')
    plt.savefig('../Results/rdn_ops%s'%max_ts+'.png')
    plt.show()

def autolabel(rects,ax):
        """
        Attach a text label above each bar displaying its height
        """
        for rect in rects:
            height = rect.get_height()
            ax.text(rect.get_x() + rect.get_width()/2., 1.05*height,
                    '%d' % int(height),
                    ha='center', va='bottom')

def barplot(stats, names, groups,colors, directory = "Barplots",name = "",ylab = ""):
    if not os.path.exists('../Results/'+directory):
        os.makedirs('../Results/'+directory)
    N = len(groups)
    ind = np.arange(N)
    width = 0.35
    fig, ax = plt.subplots()
    rects = []
    for i,key in enumerate(names.keys()):
        mean = [np.mean(stats[key][group]) for group in groups]
        err = [st.sem(stats[key][group]) for group in groups]
        rects.append(ax.bar(ind+i*width, mean, width, color=colors[key], yerr=err))
        #autolabel(rects[-1],ax)
    # add some text for labels, title and axes ticks
    ax.set_ylabel(ylab)
    ax.set_title('Group comparison')
    ax.set_xticks(ind + width / 2)
    ax.set_xticklabels(groups)
    ax.set_xlabel('Parameter')
    
    ax.legend([rects[i][0] for i in range(len(rects))], [key for key in names.keys()])
    
    #autolabel(rects1)
    #autolabel(rects2)
    plt.savefig('../Results/%s/%s.png'%(directory, name))
    plt.show()

def scaling(p):
    if p<0.05 and p>0:
        return 2.5-p*50
    elif p>-0.05 and p<0:
        return -2.5-p*50
    else:
        return 0

def forceAspect(ax,aspect=1):
    im = ax.get_images()
    extent =  im[0].get_extent()
    ax.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/aspect)

def createRasterPlots(FAM,IPP,names,scalefactor,to_file = True,directory = 'RasterPlots'):
    if not os.path.exists('../Results/'+directory):
        os.makedirs('../Results/'+directory)
    for key in names.keys():
            for exp in range(len(names[key])):
                fig = plt.figure(figsize=(12,12))
                ax = fig.add_subplot(111, aspect='equal')
                plt.suptitle('%s'%(names[key][exp]), fontsize=14, fontweight='bold')
                n_s,n_l,n_f = FAM[key][exp].shape
                for s in range(n_s):
                    plt.text(0.06+s*0.125, 1.025,ExperimentConfigFile('../RawData/'+names[key][exp]).sections()[s],
                     horizontalalignment='center',
                     verticalalignment='center',
                     fontsize=10,
                     transform = ax.transAxes)
                    #MakeRelationGraph(FAM[key][exp][s,:,:],IPP[key][exp][s,:,:],exp,s,key,directory,scalefactor)
                    _FAM = FAM[key][exp][s,:,:]
                    _IPP = IPP[key][exp][s,:,:]
                    n_l,n_f = _FAM.shape
                    print n_l, n_f
                    pair_labels = []
                    pos = 0
                    for i in range(n_l):
                        for j in range(i,n_f):
                            if i!=j and abs(_FAM[i,j])<0.05 and _FAM[i,j]>0:
                                ax.add_patch(patches.Rectangle((
                                        s, -1*pos),1 , 1,facecolor=(1,0,0,_IPP[i,j]*0.5/scalefactor)))  
                            elif i!=j and abs(_FAM[i,j])<0.05 and _FAM[i,j]<0:
                                ax.add_patch(patches.Rectangle((
                                        s, -1*pos),1 , 1,facecolor=(0,0,1,_IPP[i,j]*0.5/scalefactor)))
                            if i!=j and abs(_FAM[j,i])<0.05 and _FAM[j,i]>0:
                                ax.add_patch(patches.Rectangle((
                                        s, -1*pos),1 , 1,facecolor=(1,0,0,_IPP[j,i]*0.5/scalefactor))) 
                            elif i!=j and abs(_FAM[j,i])<0.05 and _FAM[j,i]<0:
                                ax.add_patch(patches.Rectangle((
                                        s, -1*pos),1 , 1,facecolor=(0,0,1,_IPP[j,i]*0.5/scalefactor)))
                            if i!=j:
                                #ax.add_patch(patches.Rectangle((0, -1*pos+1),8,1,facecolor="black",fill=False))
                                pair_labels.append(str(i+1)+'|'+str(j+1))
                                #pair_labels.append(str(j+1)+'|'+str(i+1))
                                pos+=1
                for i in range(8-n_s):
                    ax.add_patch(patches.Rectangle((
                                        n_s+i, -pos+1),1, pos,facecolor="lightgrey"))
                plt.axis([0,8,-pos+1,1])
                ax.set_aspect('auto')
                ax.xaxis.grid()
                ax.xaxis.set_ticklabels([])
                ax.get_yaxis().set_ticks([-1*i+0.5 for i in range(pos)])
                ax.set_yticklabels(pair_labels)
                plt.xlabel("session")
                plt.ylabel("following strength in pair")
                plt.savefig('../Results/%s/%s.png'%(directory,names[key][exp]))
                #plt.show()
                plt.close(fig)   

def createRasterPlotsSUM(FAM,IPP,names,scalefactor,to_file = True,directory = 'RasterPlotsSUM'):
    if not os.path.exists('../Results/'+directory):
        os.makedirs('../Results/'+directory)
    for key in names.keys():
            fig = plt.figure(figsize=(12,12))
            ax = fig.add_subplot(111, aspect='equal')
            plt.suptitle(key, fontsize=14, fontweight='bold')
            exp_pos = 0
            for exp in range(len(names[key])):
                n_s,n_l,n_f = FAM[key][exp].shape
                for s in range(n_s):
                    plt.text(0.06+s*0.125, 1.025,ExperimentConfigFile('../RawData/'+names[key][exp]).sections()[s],
                     horizontalalignment='center',
                     verticalalignment='center',
                     fontsize=10,
                     transform = ax.transAxes)
                    #MakeRelationGraph(FAM[key][exp][s,:,:],IPP[key][exp][s,:,:],exp,s,key,directory,scalefactor)
                    _FAM = FAM[key][exp][s,:,:]
                    _IPP = IPP[key][exp][s,:,:]
                    n_l,n_f = _FAM.shape
                    print n_l, n_f
                    pair_labels = []
                    pos = exp_pos+0
                    for i in range(n_l):
                        for j in range(i,n_f):
                            if i!=j and abs(_FAM[i,j])<0.05 and _FAM[i,j]>0:
                                ax.add_patch(patches.Rectangle((
                                        s, -1*pos),1 , 1,facecolor=(1,0,0,_IPP[i,j]*0.5/scalefactor)))  
                            elif i!=j and abs(_FAM[i,j])<0.05 and _FAM[i,j]<0:
                                ax.add_patch(patches.Rectangle((
                                        s, -1*pos),1 , 1,facecolor=(0,0,1,_IPP[i,j]*0.5/scalefactor)))
                            if i!=j and abs(_FAM[j,i])<0.05 and _FAM[j,i]>0:
                                ax.add_patch(patches.Rectangle((
                                        s, -1*pos),1 , 1,facecolor=(1,0,0,_IPP[j,i]*0.5/scalefactor))) 
                            elif i!=j and abs(_FAM[j,i])<0.05 and _FAM[j,i]<0:
                                ax.add_patch(patches.Rectangle((
                                        s, -1*pos),1 , 1,facecolor=(0,0,1,_IPP[j,i]*0.5/scalefactor)))
                            if i!=j:
                                #ax.add_patch(patches.Rectangle((0, -1*pos+1),8,1,facecolor="black",fill=False))
                                pair_labels.append(str(i+1)+'|'+str(j+1))
                                #pair_labels.append(str(j+1)+'|'+str(i+1))
                                pos+=1
                for i in range(8-n_s):
                    ax.add_patch(patches.Rectangle((
                                        n_s+i, -pos+1),1, pos,facecolor="lightgrey"))
                pos -=exp_pos
                exp_pos += pos
            plt.axis([0,8,-exp_pos+1,1])
            ax.set_aspect('auto')
            ax.xaxis.grid()
            ax.xaxis.set_ticklabels([])
            #ax.get_yaxis().set_ticks([-1*i+0.5 for i in range(pos)])
            ax.set_yticklabels([])
            plt.xlabel("session")
            #plt.ylabel("following strength in pair")
            plt.savefig('../Results/%s/%s.png'%(directory,key))
            #plt.show()
            plt.close(fig) 


def CreateRelationGraphs(FAM,IPP,names,scalefactor,to_file = True,directory = 'InteractionsGraphs'):
    if not os.path.exists('../Results/'+directory):
        os.makedirs('../Results/'+directory)
    for key in names.keys():
        for exp in range(len(names[key])):
            n_s,n_l,n_f = FAM[key][exp].shape
            for s in range(n_s):
                MakeRelationGraph(FAM[key][exp][s,:,:],IPP[key][exp][s,:,:],exp,s,key,directory,scalefactor,names)
    
def MakeRelationGraph(FAM,IPP,exp,s,key,directory,scalefactor,names, fig = None, xy0 = 0):
    print exp
    n_l,n_f = FAM.shape
    power = []
    conn = []
    for i in range(n_l):
        for j in range(n_f):
            if i!=j and abs(FAM[i,j])<0.05:
                power.append(IPP[i,j])
                conn.append([IPP[i,j]*FAM[i,j]/abs(FAM[i,j]),i+1,j+1])
    #scalefactor = np.max(power)/50
    G = nx.MultiDiGraph(multiedges=True, sparse=True)
    for i in range(len(conn)):
        G.add_edges_from([(conn[i][1],conn[i][2])], weight=conn[i][0])
    edge_colors = [conn[i][0] for i in range(len(conn))]
    size = 10
    pos=nx.circular_layout(G)
    for key in pos.keys():
        pos[key]+=np.array([2,0])
    node_labels = {node:node for node in G.nodes()}  
    if not fig:
        fig = plt.figure(figsize=(10*size,10*size))
    plt.suptitle('%s'%(names[key][exp]), fontsize=14, fontweight='bold')
    ax = fig.add_subplot(111, aspect='equal')
    #fig.suptitle(self.path, fontsize=14*size, fontweight='bold')
    nx.draw_networkx_labels(G, pos, labels=node_labels,font_size=120)
    cmap = mcol.LinearSegmentedColormap.from_list(name='red_white_blue', 
                                             colors =[(0, 0, 1), 
                                                      (1, 1., 1), 
                                                      (1, 0, 0)],
                                             N=20-1,
                                             )
    nx.draw(G,pos,ax, node_size=500*size**2,node_color= 'grey',edge_color=edge_colors,edge_cmap=cmap,width = 0,arrows=False, edge_style='dashed',zorder=10)
    for c in conn:
        p = patches.FancyArrowPatch(pos[c[2]],pos[c[1]],connectionstyle='arc3, rad=-0.3',arrowstyle="simple",shrinkA=12*size, shrinkB=12*size,mutation_scale=size*c[0]/scalefactor, color = cmap(c[0]+0.5),zorder=1,alpha=0.5)
        ax.add_patch(p)
    plt.savefig('../Results/%s/%s_exp%s_s%s_graph.png'%(directory,key,exp,s))
    #plt.show()
    plt.close(fig)    


def create_group_graph(FAM,IPP,names,scalefactor,to_file = True,directory = 'InteractionsGraphs'):
    if not os.path.exists('../Results/'+directory):
        os.makedirs('../Results/'+directory)
    max_s = 0
    for key in names.keys():
        for exp in range(len(names[key])):
            n_s,n_l,n_f = FAM[key][exp].shape
            if n_s>max_s:
                max_s = n_s
    figs = []
    size = 10
    for s in range(max_s):
        fig = plt.figure(figsize=(10*size,10*size)) 
        figs.append(fig)       
    for key in names.keys():
        for exp in range(len(names[key])):
            n_s,n_l,n_f = FAM[key][exp].shape
            for s in range(n_s):
                MakeRelationGraph(FAM[key][exp][s,:,:],IPP[key][exp][s,:,:],exp,s,key,directory,scalefactor)

def plotphist(data,names,colors,to_file = True,directory = 'Interactions',vrange = [], prange = [],nbins = 24):
    if not os.path.exists('../Results/'+directory):
        os.makedirs('../Results/'+directory)
    stats = {}
    for key in data.keys():
        stats[key] = {}
        stats[key][directory] = []
        stats[key]["mean"] = []
        stats[key]["median"] = []
        for i in range(len(data[key])):
            x = data[key][i].flatten()
            stats[key][directory]+=list(x[x!=0])
            stats[key]["mean"].append(np.mean(x[x!=0]))
            stats[key]["median"].append(np.median(x[x!=0]))
            plt.suptitle(names[key][i], fontsize=14, fontweight='bold')
            plt.hist(x, bins = np.arange(vrange[0],vrange[1],(vrange[1]-vrange[0])*1./nbins),normed=True)
            plt.axvline(np.mean(x), color=colors[key], linestyle='solid', linewidth=2)
            plt.axvline(np.median(x), color=colors[key], linestyle='dashed', linewidth=2)
            plt.ylim(prange)
            plt.savefig('../Results/'+directory+'/'+names[key][i]+'.png')
            plt.show()
    for key in data.keys():
        for i in range(len(data[key])):
            x = data[key][i].flatten()
            plt.subplot(2,1,1)
            plt.axvline(np.mean(x[x!=0]), color=colors[key], linestyle='solid', linewidth=2)
            plt.subplot(2,1,2)
            plt.axvline(np.median(x[x!=0]), color=colors[key], linestyle='dashed', linewidth=2)
    plt.savefig('../Results/'+directory+'/'+"WT_KO_RD"+'.png')
    plt.show()
    plt.show()
    return stats    
#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 26 15:51:36 2017

@author: jmaka
"""

    
def plot_tf(self):
    t1 = self.cf.gettime('SNIFF 1 dark')[0]-self.tstart
    t2 = self.cf.gettime('SNIFF 1 light')[1]-self.tstart
    t2find = []
    for i in range(self.lm):
        sig = self.sd[t1*self.fs:t2*self.fs,i]
        if len(np.where(sig==smells[path]['soc'])[0])>0:
            t2find.append(np.where(sig==smells[path]['soc'])[0][0]/self.fs)
        else:
            t2find.append(-1)
    print np.sum(self.f_sum[:2,:],axis = 0),
    plt.plot(np.sum(self.f_sum[:-1,:],axis = 0),t2find,'o')
    plt.xlabel('F parameter')
    plt.ylabel('time to find social [s]')
    plt.savefig(self.directory+'/'+'_tf_'+'ts%s'%self.treshold+'.png')
    plt.show()
    
    plt.hist(x,bins=np.arange(0,3.5,0.1))
    plt.ylim(0, 1400)
    plt.savefig('../Results/global_data_hist_%s.png'%(ts))
    plt.show()
    print np.mean(x[x>0]), np.std(x[x>0])
    i_arr = []
    for el in iexps:
        i_arr+=el
    x = np.array(i_arr)
    plt.hist(x,bins=np.arange(0,20,1))
    plt.savefig('../Results/global_exp_interhist_%s.png'%(ts))
    plt.show()
    
    plt.hist(x[x>0],bins=np.arange(0,20,1))
    plt.savefig('../Results/global_exp_interhist_no0_%s.png'%(ts))
    plt.show()




def plot_graphs(self, nr=''):
    for s in range(len(self.phases)):
        self.plot_graph(s, nr=nr)
        
        
def plot_fsums(self,base_name = 'following_sum',nr=''):
    plt.figure(figsize=(15,15))
    plt.suptitle(self.path, fontsize=14, fontweight='bold')
    for s in range(len(self.phases)):
        plt.subplot(len(self.phases),1,s+1)
        plt.title('Session %s'%(s+1))
        self.plot_fsum(s)
    plt.savefig(self.directory+'/'+'following_sum_%s_%s%s.png'%(s,self.treshold,nr))
    plt.show()
    
    plt.figure(figsize=(15,15))
    plt.suptitle(self.path, fontsize=14, fontweight='bold')
    for s in range(len(self.phases)):
        plt.subplot(len(self.phases),1,s+1)
        plt.title('Session %s'%(s+1))
        self.plot_fsum2(s)
    plt.savefig(self.directory+'/'+'all_sum_%s_%s%s.png'%(s,self.treshold,nr))
    plt.show()
    
def plot_fsum(self,s):
    for i in range(self.lm):
        data = []
        mouse_nr = []
        for j in range(self.lm):
            if i!=j:
                if self.f[i,j,s] > 1:
                    data.append(self.f[i,j,s]-1)
                    mouse_nr.append(j)
        self.f_sum[s,i] = np.sum(data)
    nr = np.arange(1,self.lm+1)
    #print nr, f_sum
    plt.stem(nr, self.f_sum[s,:])
    plt.axis([0,self.lm+1.5,0,5])
    plt.xlabel('Mouse number')
    plt.ylabel('Sum of f parameter')
    plt.xticks(nr)
    
def plot_fsum2(self,s):
    for i in range(self.lm):
        data = []
        mouse_nr = []
        for j in range(self.lm):
            if i!=j:
                if self.f[i,j,s] > 0:
                    data.append(self.f[i,j,s]-1)
                    mouse_nr.append(j)
        self.f_sum[s,i] = np.sum(data)
    nr = np.arange(1,self.lm+1)
    #print nr, f_sum
    plt.stem(nr, self.f_sum[s,:])
    plt.axis([0,self.lm+1.5,-5,5])
    plt.xlabel('Mouse number')
    plt.ylabel('Sum of f parameter')
    plt.xticks(nr)
    

def plot_graph(self,s, base_name = 'following_graph', nr=''):
    pairs = []       
    for i in range(self.lm):
        for j in range(self.lm):
            if i!=j:
                pairs.append([self.f[i][j][s],i+1,j+1])
    pairs.sort(key=lambda x: -x[0])
    conn = [(np.round((x[0]),2),x[1],x[2]) for x in pairs]
    G = nx.MultiDiGraph(multiedges=True, sparse=True)
    for i in range(len(conn)):
        G.add_edges_from([(conn[i][1],conn[i][2])], weight=conn[i][0])
    edge_colors = [conn[i][0] for i in range(len(conn))]
    size = 10
    pos=nx.circular_layout(G)
    node_labels = {node:node for node in G.nodes()}     
    fig = plt.figure(figsize=(10*size,10*size)) 
    ax = fig.add_subplot(111, aspect='equal')
    fig.suptitle(self.path, fontsize=14*size, fontweight='bold')
    nx.draw_networkx_labels(G, pos, labels=node_labels,font_size=120)
    cmap = mcol.LinearSegmentedColormap.from_list(name='red_white_blue', 
                                             colors =[(0, 0, 1), 
                                                      (1, 1., 1), 
                                                      (1, 0, 0)],
                                             N=20-1,
                                             )
    nx.draw(G,pos,ax, node_size=500*size**2,node_color= 'grey',edge_color=edge_colors,edge_cmap=cmap,width = 0,arrows=False, edge_style='dashed',zorder=10)
    for c in conn:
        if abs(c[0]-1)>0.01:
            #print c[0]
            p = patches.FancyArrowPatch(pos[c[2]],pos[c[1]],connectionstyle='arc3, rad=-0.3',arrowstyle="simple",shrinkA=12*size, shrinkB=12*size,mutation_scale=20*size*abs(c[0]), color = cmap(c[0]+0.5),zorder=1,alpha=0.5)
            ax.add_patch(p)
    plt.savefig(self.directory+'/'+base_name+'_%s_ts%s%s.png'%(s,self.treshold,nr))
    plt.close(fig)
    
def plotfhist(self,ts,to_file = True,nr = ''):
    plt.suptitle(self.path, fontsize=14, fontweight='bold')
    x = self.f[:,:,:].flatten()
    print 'Neutral%s'%len(x[x==0])
    plt.hist(x[x!=0])#,bins=np.arange(-5,5,0.1))
    plt.savefig(self.directory+'/fhist_%s%s.png'%(ts,nr))
    plt.show()



def fvalue(self, (m1,m2),t1, t2,s ):
        follow_stat, detected_idx = self.fpatterns((m1,m2),t1, t2)
        #print follow_stat
        #print follow_stat
        '''t = np.arange(-3,3,1.0/self.fs)
        for i in range(len(detected_idx[0])):
            s = detected_idx[1][i]
            plt.plot(t,self.sd[s-3*self.fs:s+3*self.fs,m1]-0.05,'ro')
            plt.plot(t,self.sd[s-3*self.fs:s+3*self.fs,m2]+0.05,'bo')
            plt.axis([-3.1,3.1,-0.5,9.5])
            plt.show()
        #print follow_stat'''
        print '#################%s%s##################'%(m1,m2)
        gpos = 0
        gneg = 0
        gp = []
        pp = [0]
        for key in follow_stat.keys():
            pos,neg,p = follow_stat[key]
            gpos+=pos
            gneg+=neg
            gp.append(p)
            pf,pa =round(prob(int(neg), 1-p, int(pos+neg)),3), round(prob(int(pos), p, int(pos+neg)),3)
#            if pf <0.05:
#                pp.append(2.5-pf*50)
#                print '\x1b[6;30;41m', key, follow_stat[key], pf,pa,'\x1b[0m'
#            elif pa <0.05:
#                pp.append(-2.5+pa*50)
#                print '\x1b[6;30;44m', key, follow_stat[key], pf,pa,'\x1b[0m'
#            else:
#                print key, follow_stat[key], pf,pa
        gpmean,gpstd = round(np.mean(gp),3), round(np.std(gp),3)
        pf,pa =round(prob(int(gneg), 1-gpmean, int(gpos+gneg)),3), round(prob(int(gpos), gpmean, int(gpos+gneg)),3)
        if pf <0.05:
            pp = 2.5-pf*50
            print '\x1b[6;30;41m',"global", [gpos,gneg,gpmean], pf,pa,'\x1b[0m'
        elif pa <0.05:
            pp = -2.5+pa*50
            print '\x1b[6;30;44m', "global", [gpos,gneg,gpmean], pf,pa,'\x1b[0m'
        else:
            pp = 0
            print "global", [gpos,gneg,gpmean], pf,pa
        print gpstd
        return np.sum(pp), gpos+gneg


def calculatestats():
    f_arr = []
    for i in range(len(fexps)):
        f_arr+=fexps[i]
        el = np.array(fexps[i])
        x = el[el!=0]
        fstats.append([len(el[el==0])]+[round(f(x),5) for f in [np.mean, np.std, st.skew, st.kurtosis]])
        istats.append([round(f(iexps[i]),5) for f in [np.mean, np.std, st.skew, st.kurtosis]])
        print fstats[-1]    
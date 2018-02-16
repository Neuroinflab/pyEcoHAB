import numpy as np
class MaximumEntropy:

    def __init__(self,data):
        """ 
        data in the form of Experiment.sd (dictionary with mice positions zeros removed).
        """
        self.data = data
        self.mice = data.keys()
        self.mice.pop("time")
        self.lm = len(self.mice)
        self.calculate_all_states()
        
    def maximum_entropy_distribution(self,positions,alpha,beta,Z):
        mice = [x for x in range(self.mice)]
        linear = alpha[mice,positions-1].sum()
        fij = np.zeros((self.lm,self.lm,8,8))
        for m1 in mice:
            new_mice = mice.copy()
            new_mice.remove(m1)
            fij[m1,new_mice,positions[m1],positions[new_mice]] = 1
            
        correlations += 0.5*beta*fij.sum()
        return np.exp(linear+correlations)/Z

    def maximum_entropy_distribution_array(self,positions,alpha,beta,Z):
        mice = [x for x in range(self.mice)]
        new_alpha = np.ones((self.lm,8,positions.shape[1]))*alpha[:,:,None]
        linear = alpha[mice,positions-1].sum(axis=(0,1))
        correlations = np.zeros((positions.shape[1],))
        for i, positions in enumerate(self.all_states):
            fij = np.zeros((self.lm,self.lm,8,8))
            for m1 in mice:
                new_mice = mice.copy()
                new_mice.remove(m1)
                fij[m1,new_mice,positions[m1],positions[new_mice]] = 1
            
            correlations[i] += 0.5*beta*fij.sum()
        return np.exp(linear+correlations)/Z
    
    def calculate_all_states(self):
        state_no = np.linspace(0,8**self.lm-1,8**self.lm,dtype=int)
        basis = [8**i for i in range(self.lm)][::-1]
        self.all_states = np.zeros((8**self.lm,self.lm),dtype=int)
        for j,n in enumerate(basis):
            self.all_states[:,j] = state_no//n
            state_no = state_no - self.all_states[:,j]*n
        self.all_states +=1

    def calculate_Z(self,alpha,beta):
        
        Z = self.maximum_entropy_distribution_array(self.all_states,alpha,beta,1)
        return Z.sum()
        
class InformationTheoryMethods:
    def remove_zeros(self):
        indices = set() #set of indices, where position is zero
        for mouse in self.mice:
            zeros = np.where(self.sd[mouse] == 0)[0]
            if len(zeros):
                indices.update(set(zeros))
        indices = sorted(list(indices))
        mask = np.ones(len(self.sd['time']), dtype=bool)
        mask[indices] = False
        for key in self.sd.keys():
            self.sd[key] = self.sd[key][mask]

    def get_mask(self,signal, starttime,endtime):
        arr = signal['time']
        idcs = np.where((arr >= starttime) & (arr < endtime))[0]
        if len(idcs) >= 2:
            return (idcs[0], idcs[-1] + 1)
        if len(idcs) == 1:
            return (idcs[0], idcs[0] + 1)
        return (0,0)
         
    def single_state_probability(self,mouse,starttime,endtime):
        result = np.zeros((8,))
        _mask_slice = self.get_mask(self.sd,starttime,endtime)
        for i in range(1,9):
            result[i-1] += len(np.where(self.sd[mouse][_mask_slice[0]:_mask_slice[1]] == i)[0])
        if len(self.sd[mouse][_mask_slice[0]:_mask_slice[1]]):
            return result/len(self.sd[mouse][_mask_slice[0]:_mask_slice[1]])
        return []

    def vector_state_probability(self,signal,starttime,endtime):
        shape = int(8**self.lm)
        result = np.empty((shape,))
        _mask_slice = self.get_mask(signal,starttime,endtime)
        length = _mask_slice[1]-_mask_slice[0]
        new_signal = np.empty((length,self.lm),dtype=np.int)
        for i,mouse in enumerate(self.mice):
            new_signal[:,i] = np.array(signal[mouse][_mask_slice[0]:_mask_slice[1]],dtype=int)
            
        basis = np.array([8**i for i, mouse in enumerate(self.mice)],dtype=np.int)
        print(length)
        for row in new_signal:
            indx = ((row-1)*basis).sum()
            result[indx] += 1
        
        return result/length
    
    def generate_independent_data_sequences(self,starttime,endtime):
        probabilities = {}
        for mouse in self.mice:
            probabilities[mouse] = self.single_state_probability(mouse,starttime,endtime)
            if not len(probabilities[mouse]):
                return []
        print(probabilities)
        mask = self.get_mask(self.sd,starttime,endtime)
        length = mask[1]-mask[0]
        surrogate_data = {}
        positions = [int(x) for x in range(1,9)]
        for mouse in self.mice:
            surrogate_data[mouse] = numpy.random.choice(positions,length,p=probabilities[mouse])

        print('generated surrogate data')
        surrogate_data['time'] = np.linspace(starttime,endtime,length)
        return surrogate_data


    def independent_data_comparison(self,starttime,endtime,title=''):
        sur_data = self.generate_independent_data_sequences(starttime,endtime)
        if not len(sur_data):
            return
        independent_states = self.vector_state_probability(sur_data,starttime,endtime)
        states = self.vector_state_probability(self.sd,starttime,endtime)
        non_zero = np.where(np.array(independent_states)!=0)[0]
        independent_states = independent_states[non_zero]
        
        independent_states = np.sort(independent_states)[::-1]
        non_zero_states =  np.where(np.array(states)!=0)[0]
        states = states[non_zero_states]
        states = np.sort(states)[::-1]
        
        return states, independent_states
    
    def plot_independent_data_comparison(self,states, independent_states):
        
        rank = np.linspace(len(non_zero),0,len(non_zero))
        rank_states = np.linspace(len(non_zero_states),0,len(non_zero_states))
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        ax.semilogy(rank,independent_states[::-1],'k',label='Independent model')
        ax.semilogy(rank_states,states[::-1],'r', label='Empirical')
        ax.set_xlabel('Location pattern rank')
        ax.set_ylabel('Location pattern probability')
        if title:
            ax.set_title(title)
        ax.legend(loc=1)
        plt.show()
        
    def plot_heat_maps(self,result,name,xlabels=None,ylabels=None,subdirectory=None,vmax=None,vmin=None,xticks=None,yticks=None):
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        cax = ax.imshow(result,interpolation='none',aspect='auto',cmap="viridis",origin="lower")#,extent=[1,8,1,8])
        cbar = fig.colorbar(cax)
        if not xlabels:
            xlabels = self.mice
        if not ylabels:
            ylabels = self.mice
        ax.get_yaxis().set_ticks([i for i,x in enumerate(ylabels)])
        ax.get_xaxis().set_ticks([i for i,x in enumerate(xlabels)])
        ax.set_xticklabels(xlabels)
        ax.set_yticklabels(ylabels)
        
        if subdirectory:
            dir_name = utils.check_directory(self.directory,subdirectory)
            new_name = os.path.join(dir_name,name)
        else:
            new_name = os.path.join(self.directory,name)
        fig.savefig(new_name+'.png',transparent=False, bbox_inches=None, pad_inches=2,frameon=None)
    def mutual_probability(self,mouse1,mouse2,starttime,endtime):
        
        result = np.zeros((8,8))
        _mask_slice = self.get_mask(self.sd,starttime,endtime)
        m1_sig = self.sd[mouse1][_mask_slice[0]:_mask_slice[1]]
        m2_sig = self.sd[mouse2][_mask_slice[0]:_mask_slice[1]]
        
        for  i,x in enumerate(m1_sig):
            result[x-1,m2_sig[i]-1] += 1
            
        normalization = result.sum()
       
        if normalization:
            return result/normalization
        return []

    def mutual_information_mouse_pair(self,mouse1,mouse2,starttime,endtime):
        "Iij = sumxy(pij(x,y)*log2(pij(x,y)/(pi(x)*pj(y))))"
        pm1 = self.single_state_probability(mouse1,starttime,endtime)
        pm2 = self.single_state_probability(mouse2,starttime,endtime)
        pm1pm2 = self.mutual_probability(mouse1,mouse2,starttime,endtime)
        
        Iij = 0
        if len(pm1) and len(pm2):
            for i in range(8):
                for j in range(8):
                    #print(i+1,j+1,pm1pm2[i,j])
                    if pm1pm2[i,j]:
                        Iij += pm1pm2[i,j]*np.log(pm1pm2[i,j]/(pm1[i]*pm2[j]))
        return Iij

    def mutual_information(self,starttime,endtime):
        
        result = np.zeros((self.lm,self.lm))
        
        for i,mouse1 in enumerate(self.mice):
            for j,mouse2 in enumerate(self.mice):
                if i > j:
                    result[i,j] = self.mutual_information_mouse_pair(mouse1,mouse2,starttime,endtime)
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)

        cax = ax.imshow(result,interpolation='none',aspect='auto',cmap="viridis",origin="lower")#,extent=[1,8,1,8])
        cbar = fig.colorbar(cax)
        
        ax.get_yaxis().set_ticks([i for i,x in enumerate(self.mice)])
        ax.get_xaxis().set_ticks([i for i,x in enumerate(self.mice)])
        ax.set_xticklabels(self.mice)
        ax.set_yticklabels(self.mice)
   
        if subdirectory:
            dir_name = utils.check_directory(self.directory,subdirectory)
            new_name = os.path.join(dir_name,name)
        else:
            new_name = os.path.join(self.directory,name)
        fig.savefig(new_name+'.png',transparent=False, bbox_inches=None, pad_inches=2,frameon=None)
    def mutual_probability(self,mouse1,mouse2,starttime,endtime):
        
        result = np.zeros((8,8))
        _mask_slice = self.get_mask(self.sd,starttime,endtime)
        m1_sig = self.sd[mouse1][_mask_slice[0]:_mask_slice[1]]
        m2_sig = self.sd[mouse2][_mask_slice[0]:_mask_slice[1]]
        
        for  i,x in enumerate(m1_sig):
            result[x-1,m2_sig[i]-1] += 1
            
        normalization = result.sum()
       
        if normalization:
            return result/normalization
        return []

    def mutual_information_mouse_pair(self,mouse1,mouse2,starttime,endtime):
        "Iij = sumxy(pij(x,y)*log2(pij(x,y)/(pi(x)*pj(y))))"
        pm1 = self.single_state_probability(mouse1,starttime,endtime)
        pm2 = self.single_state_probability(mouse2,starttime,endtime)
        pm1pm2 = self.mutual_probability(mouse1,mouse2,starttime,endtime)
        
        Iij = 0
        if len(pm1) and len(pm2):
            for i in range(8):
                for j in range(8):
                    #print(i+1,j+1,pm1pm2[i,j])
                    if pm1pm2[i,j]:
                        Iij += pm1pm2[i,j]*np.log(pm1pm2[i,j]/(pm1[i]*pm2[j]))
        return Iij

    def mutual_information(self,starttime,endtime):
        
        result = np.zeros((self.lm,self.lm))
        
        for i,mouse1 in enumerate(self.mice):
            for j,mouse2 in enumerate(self.mice):
                if i > j:
                    result[i,j] = self.mutual_information_mouse_pair(mouse1,mouse2,starttime,endtime)
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)

        cax = ax.imshow(result,interpolation='none',aspect='auto',cmap="viridis",origin="lower")#,extent=[1,8,1,8])
        cbar = fig.colorbar(cax)
        ax.get_yaxis().set_ticks([i for i,x in enumerate(self.mice)])
        ax.get_xaxis().set_ticks([i for i,x in enumerate(self.mice)])
        ax.set_xticklabels(self.mice)
        ax.set_yticklabels(self.mice)
  
        for label in ax.xaxis.get_ticklabels():
            label.set_rotation(90)
        fig.subplots_adjust(left=0.2)
        fig.subplots_adjust(bottom=0.3)
        plt.show()
        return result

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
        

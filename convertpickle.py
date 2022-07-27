import pickle

with open('Pd.pkl', 'rb') as fin : 
    Pd = pickle.load(fin)
pickle.dump(Pd, open("Pd.pkl","wb"), protocol=2) 
 
with open('Pgmax.pkl', 'rb') as fin : 
    Pgmax = pickle.load(fin)
pickle.dump(Pgmax, open("Pgmax.pkl","wb"), protocol=2) 
    
with open('Pgmin.pkl', 'rb') as fin : 
    Pgmin = pickle.load(fin)   
pickle.dump(Pgmin, open("Pgmin.pkl","wb"), protocol=2)     
    
with open('Gen2Bus.pkl', 'rb') as fin : 
    Gen2Bus = pickle.load(fin)
pickle.dump(Gen2Bus, open("Gen2Bus.pkl","wb"), protocol=2)    
    
with open('GenVarCost0.pkl', 'rb') as fin : 
    GenVarCost0 = pickle.load(fin)
pickle.dump(GenVarCost0, open("GenVarCost0.pkl","wb"), protocol=2)    
    
with open('GenVarCost.pkl', 'rb') as fin : 
    GenVarCost = pickle.load(fin)
pickle.dump(GenVarCost, open("GenVarCost.pkl","wb"), protocol=2)     
    
with open('Bij.pkl', 'rb') as fin : 
    Bij = pickle.load(fin)
pickle.dump(Bij, open("Bij.pkl","wb"), protocol=2)     
    
with open('MVAlim.pkl', 'rb') as fin : 
    MVAlim = pickle.load(fin)
pickle.dump(MVAlim, open("MVAlim.pkl","wb"), protocol=2)     
    
with open('A.pkl', 'rb') as fin : 
    A = pickle.load(fin)
pickle.dump(A, open("A.pkl","wb"), protocol=2)    

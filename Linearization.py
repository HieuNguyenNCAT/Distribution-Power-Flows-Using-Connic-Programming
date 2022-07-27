from math import sin, cos, tan, pi, sqrt, log, acos, ceil 
import cvxpy as cvx
import scipy.io as spio #Reading Matlab mat data as Python dictionary
# to run the following commands, please follow  https://www.instructables.com/id/Call-MATLAB-Script-and-Function-From-Python/
import matlab.engine
eng = matlab.engine.start_matlab()

def data_import(case):
    ret = eng.Read_RadialNet_Matpower(case)  #case33bw, case123  # these command is equivalent to run Read_Matpower_Radial(case33bw) in Matlab
    data = spio.loadmat('Radial_Netdata', squeeze_me=True) #taken from Matlab command Read_Matpower_Radial(Matpower_case_format);
    #data contains the following parameters
    #'Sbase', 'Pgmax', 'Pgmin', 'Qgmax', 'Qgmin', 'GenVarCost0', 'GenVarCost', 'G', 'B', 'MVAlim', 'connex', 'Pd', 'Qd', 'Vmax', 'Vmin','Gen2Bus'
    slack_bus =0 #case18ieee, case141, case33bw
    if case == "case18ieee":
        scale_Sbase(10,data)
    if case == "case123":
        slack_bus =113
    r = data['r']
    x = data['x']
    #slack_bus =55 #case 56 sackbus is 56 (Python index count from 0)
    #slack_bus =113
    I = len(data['Pd'])
    G = data['G']
    B = -data['B']
    MVAlim =data['MVAlim']
    status = data['status']
    print('the following parameters are set up: data, slack_bus, r, x, I, G, B, MVAlim, status')
    return [data, slack_bus, r, x, I, G, B, MVAlim, status]
  
def scale_Sbase(scale,data):  
    #change Sbase of testcase
    data['Sbase']= scale*data['Sbase']
    data['r'] = scale*data['r']
    data['x'] = scale*data['x']
    data['gs'] = data['gs']/scale
    data['bs'] = data['bs']/scale
    
def connex(i,j,data): #return 1 of (i,j) or (j,i) is a "connected" line
    temp = False
    for line in range(data['L']):
        z1= (i== data['from_bus'][line]-1) and (j== data['to_bus'][line]-1)
        z2= (j== data['from_bus'][line]-1) and (i== data['to_bus'][line]-1)
        z3= (data['status'][line]==1) #some lines are open switches
        z= (z1 or z2) and z3
        temp = temp or z
    return temp 
    
def rotatedCONE(x1,x2,x3,x4):
    #x1^2 + x2^2 <= x3*x4
    # this is equivalent to x1^2 +x2^2 + 1/4*(x3-x4)^2 <= 1/4*(x3+x4)^2
    return  [cvx.norm2(cvx.bmat( [[x1, x2, 0.5*(x3-x4)]])) <= 0.5*(x3+x4)] 
    
def rotatedCONE2(x1,x2,x3):
    #x1^2  <= x3*x4
    # this is equivalent to x1^2 +1/4*(x2-x3)^2 <= 1/4*(x2+x3)^2
    return  [cvx.norm2(cvx.bmat( [[x1, 0.5*(x2-x3)]])) <= 0.5*(x2+x3)]     

def linearSOCP(x1,x2,x3,approximation_para):
    #x1^2+x2^2 <= x3^2
    #upgrate to have auxiliary variable installation inside 
    v= approximation_para
    xi = cvx.Variable( approximation_para+1,nonneg=True)
    eta =cvx.Variable( approximation_para+1 ,nonneg=True) 
    const=[]
    const += [xi[v] <= x3 ] #(a[n]+b[n]) ]
    const += [eta[v] <= xi[v]*tan(pi/(2**(v+1)))]
    for k in range(1,v+1):
        const += [xi[k]== xi[k-1]*cos(pi/(2**(k+1))) + eta[k-1]*sin(pi/(2**(k+1)))]
        const += [eta[k]>= -xi[k-1]*sin(pi/(2**(k+1))) + eta[k-1]*cos(pi/(2**(k+1)))]
        const += [eta[k]>= xi[k-1]*sin(pi/(2**(k+1))) - eta[k-1]*cos(pi/(2**(k+1)))]              
    const += [xi[0] >=  x1 ]
    const += [eta[0] >= x2 ] 
    const += [eta[0] >= -x2 ]
    return const

def linearized_rotated_Cone(x1,x2,x3,x4,approximation_para):
    #x1^2 + x2^2 <= x3*x4
    x1x2 = cvx.Variable( 1, nonneg=True)
    const=[]
    const+=linearSOCP(x1,x2,x1x2,approximation_para)
    const+=linearSOCP(x1x2,0.5*(x3-x4), 0.5*(x3+x4), approximation_para)
    return const

# ------------------- Old version No longer used ------------------------------------------------
def linearSOCP_old(x1,x2,x3,xi,eta,approximation_para):
    #x1^2+x2^2 <= x3^2
    # this function requires declaration of auxiliary variables beforehand (xi,eta)
    v= approximation_para
    const=[]
    const += [xi[v] <= x3 ] #(a[n]+b[n]) ]
    const += [eta[v] <= xi[v]*tan(pi/(2**(v+1)))]
    for k in range(1,v+1):
        const += [xi[k]== xi[k-1]*cos(pi/(2**(k+1))) + eta[k-1]*sin(pi/(2**(k+1)))]
        const += [eta[k]>= -xi[k-1]*sin(pi/(2**(k+1))) + eta[k-1]*cos(pi/(2**(k+1)))]
        const += [eta[k]>= xi[k-1]*sin(pi/(2**(k+1))) - eta[k-1]*cos(pi/(2**(k+1)))]              
    const += [xi[0] >=  x1 ]
    const += [eta[0] >= x2 ] 
    const += [eta[0] >= -x2 ]
    return const
    
def arcsec(x):
    return acos(1/x)
def calculate_polyhedral_size(delta):
    #print('calculate the number of constraints v given the approximation error delta')
    return  ceil(log( pi/(2*arcsec(1+delta)))/log(2))
    
def find_Pythegorean(lower_theta, upper_theta):
    m_s = cvx.Variable(1, integer=True)
    n_s = cvx.Variable(1, integer=True)
    const =[]
    const += [n_s>=1]
    const += [sqrt(1+ sin(lower_theta))*n_s <= m_s*sqrt(1- sin(lower_theta))]
    const += [sqrt(1+ sin(upper_theta))*n_s >= m_s*sqrt(1- sin(upper_theta))]
    problem1 = cvx.Problem(cvx.Minimize(cvx.square(m_s)+cvx.square(n_s)), const)
    problem1.solve(solver='CPLEX',verbose=False, cplex_params={"mip.tolerances.absmipgap": 1e-12})
    
    m_c = cvx.Variable(1, integer=True)
    n_c = cvx.Variable(1, integer=True)
    const =[]
    const += [n_c >= 1]
    const += [sqrt(1+ cos(lower_theta))*n_c >= m_c*sqrt(1- cos(lower_theta))]
    const += [sqrt(1+ cos(upper_theta))*n_c <= m_c*sqrt(1- cos(upper_theta))]
    problem2 = cvx.Problem(cvx.Minimize(cvx.square(m_c)+cvx.square(n_c)), const)
    problem2.solve(solver='CPLEX',verbose=False, cplex_params={"mip.tolerances.absmipgap": 1e-12})
    if (lower_theta == pi/2):
        return [ int(2*m_c.value[0]*n_c.value[0]),  int(m_c.value[0]**2-n_c.value[0]**2),  int(m_c.value[0]**2+n_c.value[0]**2)]
    elif (m_s.value[0]**2+n_s.value[0]**2)  < m_c.value[0]**2+n_c.value[0]**2:
        return [ int(m_s.value[0]**2-n_s.value[0]**2),  int(2*m_s.value[0]*n_s.value[0]),  int(m_s.value[0]**2+n_s.value[0]**2)]
    else:
        return [ int(2*m_c.value[0]*n_c.value[0]),  int(m_c.value[0]**2-n_c.value[0]**2),  int(m_c.value[0]**2+n_c.value[0]**2)]
#print(find_Pythegorean(pi/4, 1.1*pi/4))

def generate_rational_coeff(delta):
    v = calculate_polyhedral_size(delta)
    k_v = 2*((2/pi*2*arcsec(1+delta) )**(1/v))
    k_v =(1-1e-6)*k_v
    a={}
    b={}
    c={}
    theta = pi/4
    for i in range(1,v+1):
        [a[i], b[i], c[i]] = find_Pythegorean(theta, k_v*theta)
        if b[i] ==0:
            theta =acos(0)
        else:
            theta =arcsec(c[i]/b[i])
        theta = pi/(2**(i+1))
    return [a, b, c, v]
           
def Rational_linearSOCP(x1,x2,x3,coeff):
    #x1^2+x2^2 <= x3^2
    [a_coff, b_coff, c_coff, v] = coeff 
    xi = cvx.Variable( v+1,nonneg=True)
    eta =cvx.Variable( v+1,nonneg=True) 
    const=[]
    const += [xi[v] <= x3 ] 
    const += [b_coff[v]*eta[v] <= xi[v]*a_coff[v]]
    for k in range(1,v+1):
        const += [c_coff[k]*xi[k]== xi[k-1]*b_coff[k] + eta[k-1]*a_coff[k]]
        const += [c_coff[k]*eta[k]>= -xi[k-1]*a_coff[k] + eta[k-1]*b_coff[k]]
        const += [c_coff[k]*eta[k]>= xi[k-1]*a_coff[k] - eta[k-1]*b_coff[k]]              
    const += [xi[0] >=  x1 ]
    const += [eta[0] >= x2 ] 
    const += [eta[0] >= -x2 ]
    return const
    
def Rational_linearized_rotated_Cone(x1,x2,x3,x4,delta):
    #x1^2 + x2^2 <= x3*x4
    x1x2 = cvx.Variable( 1, nonneg=True)
    const=[]
    const+=Rational_linearSOCP(x1,x2,x1x2,delta)
    const+=Rational_linearSOCP(x1x2, 0.5*(x3-x4), 0.5*(x3+x4), delta)
    return const
    
def approximation_error_of_BenTal(v):
    return (1/cos(pi/2**(v+1))-1)
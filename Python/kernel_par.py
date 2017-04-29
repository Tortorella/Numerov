import potentials as pot
import numpy as np
class params:

# Physical params
  m        = 1   # Particle mass  
  Lambda   = 0   # Cut-off
  C        = 0   # Interaction parameter
  hbaron2m = 1   # Kinetic constant

# System parm
  E           = -20   # energy
  Estep       = 1     # initial energy step  
  Max_dist    = 10    # Max distance
  Precision   = 0.1   # Precision in distance
  E_precision = 0.001 # Precision in energy
  l           = 0     # Angular momentum
  nodi_voluti = 0     # Energy level needed

# Will be inizialied lately
  Step        = 0     
  Rstart      = 0
  
# discretization parm
  N_of_bins   = 10000     # Discretization array
  Match       = 5000      # Where to match the two wavefunctions

# to be saved parms
  done       = 0    # Starts not done
  Last_Estep = 1    # Starts increasing the energy

# Potentials
  potenziale = pot.pionless

# From bin to distance
  def Distance(self,i):
      return (i)*self.Step

  #def __init__(self):


  def setup(self):
      self.Step       = self.Max_dist/self.N_of_bins  # Simulation step
      self.Rstart     = self.Distance(0)              # Initial step init





    
# kernel
def numerov(P, WF):

  # Local variables
  max_WF     = WF[0]
  nodi       = 0
  Range      = P.Max_dist
  Step       = P.Step
  Binni      = P.N_of_bins

  # From 0 to Match
  #  First two points
  R     = P.Distance(0)
  f0    = 0.
  f1    = F(P, P.Rstart-Step)
  WF[0] = 0.
  WF[1] = Step
  
  for i in range(2,P.Match+2):
    # new WF point: 
    R   = P.Distance(i)
    f2  = f1
    f1  = f0
    f0  = F(P,R)
    WF[i]= (  (12.-10.*f1)*WF[i-1]-f2*WF[i-2]  )/f0
    if WF[i]*WF[i-1] < 0:
      nodi=nodi+1
    if np.abs(WF[i])>max_WF:
      max_WF=np.abs(WF[i])
      
  # Preparing the matching point for the revers run    
  WF_Match   = WF[P.Match]
  WF_Match_1 = WF[P.Match+1]

  # From inf to Match
  #  Last two points
  R     = P.Distance(Binni)
  f0    = F(P,R)
  f1    = F(P,R-Step)
  WF[Binni-1]   = 0.001
  WF[Binni-2] = 0.5 * (12. - f0 * 10.) * WF[Binni-1] / f1;
  
  for i in range(Binni-3,P.Match-1,-1):
    # new WF point:
    R   = P.Distance(i)
    f2  = f1
    f1  = f0
    f0  = F(P,R)
    WF[i]= (  (12.-10.*f1)*WF[i+1]-f2*WF[i+2]  )/f0
    if WF[i]*WF[i+1] < 0:
      nodi=nodi+1
    if np.abs(WF[i])>max_WF:
      max_WF=np.abs(WF[i])
  
  # Rescaling vector in order to get the match
  WF[P.Match-1:] = WF[P.Match-1:] * (WF_Match/WF[P.Match])



  # Changing energy according the energy (Bisection algorithm)
  if ((WF[P.Match+1]-WF_Match_1) < 0):
  #if (nodi <= nodi_voluti): #aumento l'energia
    if P.Last_Estep < 0:
      P.Estep      = P.Estep/2.
      P.Last_Estep = 1
    P.E = P.E+P.Estep
  else:                    # Diminuisco l'energia
    if (P.Last_Estep > 0):
      P.Estep      = P.Estep/2.
      P.Last_Estep = - 1
    P.E = P.E-P.Estep
    #print 'Step= ', Step, 'N= ',nodi,'dE= ', '%.1E' %P.Estep,'E= ', '%2.8E' %P.E  , '       Match: ', '%.1E' %(WF[P.Match+1]-WF_Match_1)

  # Is the job done?
  if (np.abs(P.Estep) < P.E_precision):
    P.done=1
    return



def G(P,x):
    return (P.E-P.potenziale(x))/P.hbaron2m  #-l*(l+1)/x**2
def F(P,x):
    return 1. + P.Step**2*2*G(P,x)/12.


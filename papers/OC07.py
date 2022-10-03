from numpy import *

def y_a(St):
  """
  solves exactly for y_a (Mathematica)
  """
  a = -1./(1+St**-1)
  #the minimum value of a is -0.5
  #corresponding to St=1 particles
  if type(St)==ndarray:
    a[a<-0.5] = -0.5
  else:
    if a<-0.5: a = -0.5

  out = (3 - sqrt(3)*sqrt(11 - (4*(28 + 45*a))/
             (-530 + 27*a*(-32 + 9*a) + 
                9*sqrt(5636 + a*(21760 + 3*a*(7612 + 3*a*(424 + 81*a)))))**
              0.3333333333333333 + 
            2*(-530 + 27*a*(-32 + 9*a) + 
                9*sqrt(5636 + a*(21760 + 3*a*(7612 + 3*a*(424 + 81*a)))))**
              0.3333333333333333) + 
         sqrt(6)*sqrt(11 + (56 + 90*a)/
             (-530 + 27*a*(-32 + 9*a) + 
                9*sqrt(5636 + a*(21760 + 3*a*(7612 + 3*a*(424 + 81*a)))))**
              0.3333333333333333 - 
            (-530 + 27*a*(-32 + 9*a) + 
               9*sqrt(5636 + a*(21760 + 3*a*(7612 + 3*a*(424 + 81*a)))))**
             0.3333333333333333 + (9*sqrt(3)*(1 - 4*a))/
             sqrt(11 - (4*(28 + 45*a))/
                (-530 + 27*a*(-32 + 9*a) + 
                   9*sqrt(5636 + a*(21760 + 3*a*(7612 + 3*a*(424 + 81*a)))))**
                 0.3333333333333333 + 
               2*(-530 + 27*a*(-32 + 9*a) + 
                   9*sqrt(5636 + a*(21760 + 3*a*(7612 + 3*a*(424 + 81*a)))))**
                 0.3333333333333333)))/12.

  return out


def Dv_turb(St_a, St_b, Re=inf, doExact=True):
  """
    The Ormel and Cuzzi expressions for relative velocities between particles. We
    follow Eq. (3.15).  Results are expressed in terms of the gas velocity v_g

    y*=1.6 is assumed (for St<1); for St>1 y*=1

    Input:
      St_a,St_b :: Stokes number particles 1,2
      Re        :: Reynolds number of the flow (default: infinite)
  """
  if type(St_a)==ndarray:
    St1 = where(St_a>=St_b,St_a,St_b)
    St2 = where(St_a<St_b, St_a,St_b)
  else:
    St1 = max(St_a,St_b)
    St2 = min(St_a,St_b)

  #try
  if doExact:  St12a = y_a(St1) *St1
  else: St12a = 1.6*St1

  #but..
  if type(St_a)==ndarray:
    ia = St12a>=1
    ib = St12a<=Re**-0.5
    St12a[ia] = 1
    St12a[ib] = Re**-0.5
  else:
    if St12a>=1.0: St12a=1
    if St12a<=Re**-0.5: St12a=Re**-0.5

  #this is just equation (17) of OC07
  D2v1 = (St1-St2)/(St1+St2) *( St1**2/(St12a+St1) -St1**2/(1+St1) 
                               -St2**2/(St12a+St2) +St2**2/(1+St2))

  #below formula is somewhat different... WHY???
  #[13.07.09] and wrong I think
  #D2v2 = (St12a -Re**-0.5) *(St1*(St12a+Re**-0.5) +St12a*Re**-0.5)\
  #      /( (St2+St12a)*(St2+Re**-0.5) )\
  #      +(St12a -Re**-0.5) *(St2*(St12a+Re**-0.5) +St12a*Re**-0.5)\
  #      /( (St2+St12a)*(St2+Re**-0.5) )


  #this is just equation (18) of OC07
  D2v2 = 2*(St12a -Re**-0.5)\
        +St1**2/(St1+St12a) -St1**2/(St1+Re**-0.5)\
        +St2**2/(St2+St12a) -St2**2/(St2+Re**-0.5)
       #+St1**2 *(Re**-0.5-St12a) /( (St1+St12a)*(St1+Re**-0.5) )\
       #+St2**2 *(Re**-0.5-St12a) /( (St2+St12a)*(St2+Re**-0.5) )


  if type(St_a)==ndarray:
    D2v2[D2v2<0] = 0
  else:
    if D2v2<0: D2v2=0

  return sqrt(D2v1+D2v2)

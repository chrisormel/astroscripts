import numpy as np

def epsilon_plateau (St, xp): 
    """
    returns \epsilon in the plateau limit
    See HO23, eqs. 23 & Table 3.

    xp -- planet radius normalized by Hill radius
    """
    ps = 3.8/St
    pb = np.minimum(1, 2.36 *xp**0.5)
    pi = 0.610


    nom = (1-ps)*(1-pb)*pi
    den = 1 -(1-ps)*(1-pb)*(1-pi)

    return 1 -nom/den


def epsilon_balset (St, qpl, eta, rpl):
    """
    return the ballistic epsilon, following Liu & Ormel (2018)

    input:
        qpl :planet mass normalized to stellar mass
        rpl :planet radius normalized to stellar radius
        eta :pressure gradient parameter
    """
    dv = eta/(1+5.7*qpl*St/eta**3) +0.52*(qpl*St)**(1/3)
    epsset =  0.32 *np.sqrt(qpl*dv /(St*eta**2))

    epsbal = rpl/(2*np.pi*eta*St) *np.sqrt(2*qpl/rpl +dv**2)

    return epsbal, epsset


def e_eq (qp, eta, b=1.9):
    """
    resonance equilibrium eccentricity, eq.28 of HO23
    """
    eeq = np.sqrt(1.5*eta*b*(qp/3)**(1/3))
    return eeq


def delp_HN90 (p1, p2, b):
    """
    change of eccentricity vector due to encounter (Hasegawa & Nakazwa 1990)
    All in Hill units
    """
    R1 = 0.747
    R2 = 2.38
    R4 = 3.46
    R5 = 0.147
    R6 = -1.86
    R7 = 1.73

    delp1 = -9/b**3 *(( 2/9 -R2/2)*p2 -R4*p1*p2/b )
    delp2 = -9/b**2 *(-R1 -(2/9 +R2/2)*p1/b +R6*p1**2/b**2 +R7*p2**2/b**2)

    return delp1, delp2


def b_kick (qp, eta, eeq, b=1.9):
    """
    Hasegawa expressions with p1=0

    Hasegawa uses Hill units (e.g, b=1.9)
    but we convert back/forth to orbital units
    """
    hM = (qp/3)**(1/3)

    #convert eccentricity vector to Hill units
    #p1 = -2*eeq**2 /(qp/3)**(1/3)
    #p2 = eeq /(qp/3)**(1/3)
    p1 = 0.
    p2 = np.sqrt(1.5*eta*b /hM)

    delp1, delp2 = delp_HN90 (p1, p2, b)

    #convert back to orbital units
    e2 = ((p1+delp1)**2 +(p2+delp2)**2) *hM**2
    dele2 = e2 -eeq**2

    #and the corresponding change in \Delta b
    #(Eq.10 of HO23)
    bnew2 = 4/3 *dele2 + (b*hM)**2
    delb = np.sqrt(bnew2) -b*hM

    return delb


def St_res_weak (qp, eta, b=1.9):
    Awd = 2.7

    #equilibrium eccentricity for j+1:j resonance
    #(independent of j)
    eeq = e_eq(qp, eta)

    #orbital units (not Hill units as in paper)
    bkick = b_kick (qp, eta, eeq)   
    borb = b*(qp/3)**(1/3)

    #this is eq.27 of HO23
    Stres_wk = Awd *8*np.pi*eta /(3*borb*bkick)

    return Stres_wk


def epsilon_HO (Starr, qp, eta, rcap, reduceRes=False):
    """
    rpl     :planet radius in terms of orbital radius (semi-major axis)
    rcap    :capture radius in terms of orbital radius
    """

    #capture radius in terms of Hill radius
    Rcap_to_Rh = rcap /(qp/3)**(1/3)

    epsplat = epsilon_plateau (Starr, Rcap_to_Rh)
    epsbal1, epsset1 = epsilon_balset (1.0, qp, eta, rcap)
    epsbal, epsset = epsilon_balset (Starr, qp, eta, rcap)

    #expressions to obtain critical Stokes numbers, Eq. (24,25,27,31) of HO23
    Stplat = 5.56 *eta/qp**(2/3)    
    Stcrit = 168 *(qp/3e-6)**(-1/3)
    Stres_str = 36.0 *(eta/1e-3) *(qp/3e-6)**(-2/3) #the strong version
    Stres_wk = St_res_weak (qp, eta)

    #decide to use weak or strong limit
    #Stres = np.where(Stcrit>Stres_str, Stres_wk, Stres_str)
    Stres = np.where(Stcrit>Stres_wk, Stres_str, Stres_wk)

    #print('crit strong weak', Stcrit, Stres_str, Stres_wk)

    epsarr = epsplat.copy()
    ires = Starr>Stres

    if reduceRes:
        #here we simply assing the plateau value
        #rather than interatively
        #decreasing the mass/size of particles in some way
        epsarr[ires] = epsilon_plateau (Stres[ires], Rcap_to_Rh[ires])
    else:
        epsarr[ires] = 0.0 #put zero for simplicity

    ii = Starr<Stplat

    #epsarr[ii] = np.sqrt(epsset1[ii]**2 +epsbal1[ii]**2 *Starr[ii]**2)
    epsarr[ii] = np.maximum(epsset1, epsbal1 *Starr[ii])


    ii = Starr<1.0
    epsarr[ii] = epsset[ii]

    return epsarr, Stplat, Stcrit, Stres



import numpy as np
#job finding probability
def jf(theta, nu_L=None, eta_L=None, A=None):
    if nu_L is None:
        f = A*theta**(1-eta_L)
        f = np.minimum(f, 1)
        return f
    else:
        return theta/(1+theta**( nu_L))**(1/ nu_L)
    
    
#vacancy filling probability
def vf(theta,  nu_L=None, eta_L=None, A=None):
    if nu_L is None:
        q = A*theta**(-eta_L)
        q = np.minimum(q, 1)
        return q
    else:
        return 1/(1+theta**nu_L)**(1/ nu_L)


def jf_invert(f,  nu_L=None, eta_L=None, A=None):
    # theta as a function of the job finding rate
    f = np.maximum(f, 0)
    f = np.minimum(f, 0.999)
    if nu_L is None:
        return (f/A)**(1/(1-eta_L))
    else:
        x = f**nu_L/(1-f**nu_L)
        return x**(1/ nu_L)


def theta_invert(q,  nu_L=None, eta_L=None, A=None):
    #q = np.asarray(q)
    q = np.maximum(q, 1e-10)
    q = np.minimum(q, 1)
    if nu_L is None:
        return (A/q)**(1/eta_L)
    else:
        return (1/q**( nu_L) - 1)**(1/ nu_L)

    
    
def eta_L_fun(theta,  nu_L):
    # Elasticity of matching with respect to job seekers, which equals
    # negative elasticity of vacancy filling rate with respect to tightnss
    return theta**nu_L/(1+theta**nu_L)


def epsi_f_theta_fun(theta,  nu_L):
    # Elasticity of job finding rate with respect to tightness
    return 1/(1+theta**nu_L)
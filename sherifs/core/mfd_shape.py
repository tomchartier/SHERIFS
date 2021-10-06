 # -*- coding: utf-8 -*-

"""SHERIFS
Seismic Hazard and Earthquake Rates In Fault Systems

Populates the magnitude bins with the faults and scenarios that can generate these magnitudes.

@author: Thomas Chartier
"""
import numpy as np
import scipy
from sherifs.core.seismic_moment import mag_to_M0

def GR(mfd_param, bin_mag):
    b_value = mfd_param['b_value']
    Beta =  b_value * np.log(10)
    p_MFD = [] #probability distribution
    for mag in bin_mag :
        p_i = (((np.exp(-Beta * (mag - 0.05 - bin_mag[0]))) / ( 1. - np.exp(-Beta * (bin_mag[-1] - bin_mag[0]))))
        - ((np.exp(-Beta * (mag + 0.05 - bin_mag[0]))) / ( 1. - np.exp(-Beta * (bin_mag[-1] - bin_mag[0])))))
        p_MFD.append(p_i)

    return p_MFD

def tapered_GR(mfd_param,bin_mag):
    b_value = mfd_param['b_value']
    Beta =  b_value * 2./3.
    M_corner = mfd_param['M_corner']
    M0_c = mag_to_M0(M_corner)
    M_t = bin_mag[0]-0.05
    M0_t = mag_to_M0(M_t)

    p_MFD = [] #probability distribution
    for mag in bin_mag :
        M0_min = mag_to_M0(mag-0.05)
        M0_max = mag_to_M0(mag+0.05)
        p_i = (((M0_t/M0_min)**Beta*np.exp((M0_t-M0_min)/M0_c))
        - ((M0_t/M0_max)**Beta*np.exp((M0_t-M0_max)/M0_c)))
        p_MFD.append(p_i)

    return p_MFD

def user_defined(mfd_param,bin_mag):

    return p_MFD

def double_GR(mfd_param,bin_mag):
    b_value = mfd_param['b_value']
    #Mrupt = mfd_param['Mrupt']
    Mrupt = mfd_param["Mrupt"]
    Beta =  b_value * np.log(10)
    print('use of the MFD "double GR" with Mrupt = '+str(Mrupt))
    p_MFD = [] #probability distribution
    for mag in [float(round(i,1)) for i in bin_mag] :
        if mag < Mrupt-0.2:
            p_i = (((np.exp(-Beta * (mag - 0.05 - bin_mag[0]))) / ( 1. - np.exp(-Beta * (Mrupt - bin_mag[0]))))
            - ((np.exp(-Beta * (mag + 0.05 - bin_mag[0]))) / ( 1. - np.exp(-Beta * (Mrupt - bin_mag[0])))))
        elif mag == Mrupt :
            p_i = p_MFD[-5]
        elif mag > Mrupt:
            p_i = p_MFD[-1] * (10**(-b_value*0.1))
        else :
            p_i = min(p_MFD)/3.
        p_MFD.append(p_i)

    return p_MFD


def YC(mfd_param,bin_mag,Mmax):
    b_value = mfd_param['b_value']
    beta = b_value * np.log(10)

    p_MFD = [] #probability distribution
    for mag in bin_mag :
        if mag < (Mmax - 0.5) :
            pi = 0.1 * ((beta*np.exp(-beta*(mag -bin_mag[0])))/(1.-np.exp(-beta*(Mmax-bin_mag[0]))))
        elif mag < (Mmax + 0.1) :
            pi = 0.1*((beta*(np.exp(-beta*(Mmax-1.5-bin_mag[0]))))/(1-np.exp(-beta*(Mmax-bin_mag[0]))))
        p_MFD.append(pi)
    return p_MFD

def YC_modified(mfd_param,bin_mag,Mmax):
    #Youngs and coppersmith 1985, see Makrup et al 2015 for simplified formulation
    #this equation has been modified in order to have a stable moment for the last three bins of magnitude and keep the
    # target stable for the GR part of the distribution
    # these value are good for fitting the not declestered catalog for the greater marmara region
    # with a b value of 1.15
    b_value = mfd_param['b_value']
    Mf = mfd_param['Mf'] #Mmax with the best a priori fit, value for the Marmara sea
    size_of_bump = mfd_param['size_of_bump'] #in order to modify the respective size of the two parts of the
    beta = b_value * np.log(10)
    '''
    if Mmax <= 7.8 :
        Mf= 7.5 #Mmax with the best a priori fit, value for the Marmara sea
        size_of_bump = 0.25 #in order to modify the respective size of the two parts of the
    else :
        Mf= 7.8 #Mmax with the best a priori fit, value for the Marmara sea
        size_of_bump = 0.5 #in order to modify the respective size of the two parts of the
    '''
    #Mf= 7.5
    #size_of_bump = 0.25

    p_MFD = [] #probability distribution
    for mag in bin_mag :
        if mag < (Mmax - 0.8) :
            pi = 0.1 * (((beta*np.exp(-beta*(mag -bin_mag[0])))/(1.-np.exp(-beta*(Mmax-bin_mag[0]))))
            *np.exp(1./b_value * (Mmax-Mf)))
        elif mag < (Mmax + 0.09) :
            pi = 0.1*(((beta*(np.exp(-beta*(Mmax-1.5-bin_mag[0]))))/(1-np.exp(-beta*(Mmax-bin_mag[0]))))
            *scipy.stats.norm(0, 2).pdf((mag-(Mmax-0.4))*5.)/scipy.stats.norm(0, 2).pdf(0.))*size_of_bump
        elif mag > (Mmax + 0.09) :
            pi = 0.
        p_MFD.append(pi)
    return p_MFD

def YC_marmara(mfd_param,bin_mag,Mmax):
    #Youngs and coppersmith 1985, see Makrup et al 2015 for simplified formulation
    #this equation has been modified in order to have a stable moment for the last three bins of magnitude and keep the
    # target stable for the GR part of the distribution
    # these value are good for fitting the not declestered catalog for the greater marmara region
    # with a b value of 1.15
    b_value = mfd_param['b_value']
    Mf = mfd_param['Mf'] #Mmax with the best a priori fit, value for the Marmara sea
    size_of_bump = mfd_param['size_of_bump'] #in order to modify the respective size of the two parts of the
    beta = b_value * np.log(10)
    '''
    if Mmax <= 7.8 :
        Mf= 7.5 #Mmax with the best a priori fit, value for the Marmara sea
        size_of_bump = 0.25 #in order to modify the respective size of the two parts of the
    else :
        Mf= 7.8 #Mmax with the best a priori fit, value for the Marmara sea
        size_of_bump = 0.5 #in order to modify the respective size of the two parts of the
    '''
    #Mf= 7.5
    #size_of_bump = 0.25

    p_MFD = [] #probability distribution
    for mag in bin_mag :
        if mag < (Mmax - 0.8) :
            pi = 0.1 * (((beta*np.exp(-beta*(mag -bin_mag[0])))/(1.-np.exp(-beta*(Mmax-bin_mag[0]))))
            *np.exp(1./b_value * (Mmax-Mf)))
        elif mag < (Mmax + 0.09) :
            pi = 0.1*(((beta*(np.exp(-beta*(Mmax-1.5-bin_mag[0]))))/(1-np.exp(-beta*(Mmax-bin_mag[0]))))
            *scipy.stats.norm(0, 2).pdf((mag-(Mmax-0.4))*5.)/scipy.stats.norm(0, 2).pdf(0.))*size_of_bump
        elif mag > (Mmax + 0.09) :
            pi = 0.
        p_MFD.append(pi)
    return p_MFD

def user_defined(mfd_param,bin_mag):#Youngs and coppersmith 1985, see Makrup et al 2015 for simplified formulation
    #this  equation has been modified in order to have a stable moment for the last three bins of magnitude and keep the
    # target stable for the GR part of the distribution

    b_value = mfd_param['b_value']
    Mf = mfd_param['Mf'] #Mmax with the best a priori fit
    size_of_bump = mfd_param['size_of_bump'] #in order to modify the respective size of the two parts of the distribution

    beta = b_value * np.log(10)
    #Mf= 6.4 #Mmax with the best a priori fit
    #size_of_bump = 1.

    self.calculation_log_file.write('\nUsing YC_modified with Mf='+ str(Mf)+' and size_of_bump='+str(size_of_bump)+'\n')

    p_MFD = [] #probability distribution
    for mag in bin_mag :
        if mag < (Mmax - 0.5) :
            pi = 0.1 * (((beta*np.exp(-beta*(mag -bin_mag[0])))/(1.-np.exp(-beta*(Mmax-bin_mag[0]))))
            *np.exp(1./b_value * (Mmax-Mf)))
        elif mag < (Mmax + 0.09) :
            pi = 0.1*(((beta*(np.exp(-beta*(Mmax-1.5-bin_mag[0]))))/(1-np.exp(-beta*(Mmax-bin_mag[0]))))
            *scipy.stats.norm(0, 1).pdf((mag-(Mmax-0.25))*5.)/scipy.stats.norm(0, 1).pdf(0.))*size_of_bump
        elif mag > (Mmax + 0.09) :
            pi = 0.
        p_MFD.append(pi)

    return p_MFD

def UCERF_DV(bin_mag):
    mfd_ucerf_old = [0.0802928535,0.006255187,0.0811119499,0.0158967773,0.0852733969,
                 0.1133026657,0.0851402625,0.0811713962,0.017219035,0.011174207,
                 0.0092182702,0.0123154789,0.0066646595,0.0160229518,0.0328559211,
                 0.0630301509,0.2736356088,0.008430495,0.0009887329]

    mfd_ucerf_guessed = [0.1138501463,0.0088694563,0.1150115728,0.0225406165,0.1209122392,0.1606559551,
                 0.1207234632,0.1150958638,0.0244154936,0.015844313,0.0130709194,0.0032216881,
                 0.0041427486,0.0104758837,0.0226150741,0.0223162261,0.1019376482,0.0040203001,0.0002803921]

    mfd_ucerf = [0.,0.1225479995,0.0095470667,0.1237981981,0.0181124989,
                 0.1286539225,0.1709131163,0.1226901936,0.072599789,0.0215865667,
                 0.0120345436,0.0086226427,0.0115716064,0.0049446472,0.0100869533,
                 0.0209929462,0.0310850941,0.107063831,0.0028465709,0.0003018133]
    #without the two last bin (very low because derived from a mean mfd)
    mfd_ucerf = [0.,0.1225479995,0.0095470667,0.1237981981,0.0181124989,
                 0.1286539225,0.1709131163,0.1226901936,0.072599789,0.0215865667,
                 0.0120345436,0.0086226427,0.0115716064,0.0049446472,0.0100869533,
                 0.0209929462,0.0310850941,0.107063831]
    #mfd for the mean branch of UCERF
    #taken from the binnin in 0.05 but mean between two bins in order to smooth the curve a bit
    mfd_ucerf = [0.,0.1170913873,0.1180649338,0.0434295744,0.0431719421,
                 0.0279491815,0.1776244068,0.1772698507,0.0686700644,
                 0.0460331714,0.0073358841,0.008953384,0.0109225241,
                 0.0065276199,0.0022473075,0.0059243291,0.0092550551,
                 0.0662428769,0.0627078821,0.0005786247]



    bin_mfd_ucerf = np.linspace(6.0,7.9,len(mfd_ucerf))

    interp_ucerf_mfd = interp1d(bin_mfd_ucerf,mfd_ucerf)

    p_MFD = [] #probability distribution
    for mag in bin_mag :
        p_MFD.append(interp_ucerf_mfd(mag))

    return p_MFD

def incrementally_defined(bin_mag,inc_rates):
    mfd_ucerf_old = [0.0802928535,0.006255187,0.0811119499,0.0158967773,0.0852733969,
                 0.1133026657,0.0851402625,0.0811713962,0.017219035,0.011174207,
                 0.0092182702,0.0123154789,0.0066646595,0.0160229518,0.0328559211,
                 0.0630301509,0.2736356088,0.008430495,0.0009887329]

    mfd_ucerf_guessed = [0.1138501463,0.0088694563,0.1150115728,0.0225406165,0.1209122392,0.1606559551,
                 0.1207234632,0.1150958638,0.0244154936,0.015844313,0.0130709194,0.0032216881,
                 0.0041427486,0.0104758837,0.0226150741,0.0223162261,0.1019376482,0.0040203001,0.0002803921]

    mfd_ucerf = [0.,0.1225479995,0.0095470667,0.1237981981,0.0181124989,
                 0.1286539225,0.1709131163,0.1226901936,0.072599789,0.0215865667,
                 0.0120345436,0.0086226427,0.0115716064,0.0049446472,0.0100869533,
                 0.0209929462,0.0310850941,0.107063831,0.0028465709,0.0003018133]
    #without the two last bin (very low because derived from a mean mfd)
    mfd_ucerf = [0.,0.1225479995,0.0095470667,0.1237981981,0.0181124989,
                 0.1286539225,0.1709131163,0.1226901936,0.072599789,0.0215865667,
                 0.0120345436,0.0086226427,0.0115716064,0.0049446472,0.0100869533,
                 0.0209929462,0.0310850941,0.107063831]
    #mfd for the mean branch of UCERF
    #taken from the binnin in 0.05 but mean between two bins in order to smooth the curve a bit
    mfd_ucerf = [0.,0.1170913873,0.1180649338,0.0434295744,0.0431719421,
                 0.0279491815,0.1776244068,0.1772698507,0.0686700644,
                 0.0460331714,0.0073358841,0.008953384,0.0109225241,
                 0.0065276199,0.0022473075,0.0059243291,0.0092550551,
                 0.0662428769,0.0627078821,0.0005786247]



    bin_mfd_ucerf = np.linspace(6.0,7.9,len(mfd_ucerf))

    interp_ucerf_mfd = interp1d(bin_mfd_ucerf,mfd_ucerf)

    p_MFD = [] #probability distribution
    for mag in bin_mag :
        p_MFD.append(interp_ucerf_mfd(mag))

    return p_MFD

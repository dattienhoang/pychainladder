# -*- coding: utf-8 -*-
"""
Created on Thu Feb 23 13:52:12 2017

@author: basaleh, dahoang, yocho
"""

import numpy as np
import pandas as pd

# something to convert an incremental triangle to a cumulative triangle
def incre2cumul(incre_tri):
    if type(incre_tri) != np.ndarray:
        incre_tri = np.asarray(incre_tri)
    cumul_tri = np.copy(incre_tri)
    for i in range(len(incre_tri)):
        cumul_tri[i] = np.cumsum(incre_tri[i])
    return cumul_tri

# something to convert a cumulative triangle to an incremental triangle
def cumul2incre(cumul_tri):
    if type(cumul_tri) != np.ndarray:
        cumul_tri = np.asarray(cumul_tri)
    incre_tri = np.copy(cumul_tri)
    for i in range(len(cumul_tri)):
        for j in range(len(cumul_tri)):
            if j == 0:
                incre_tri[i][j] = cumul_tri[i][j]
            if j != 0:
                incre_tri[i][j] = cumul_tri[i][j] - cumul_tri[i][j-1]
    return incre_tri

# basic chain ladder technique
def chain_ladder(tri):
    if type(tri) != np.ndarray:
        tri = np.asarray(tri)
    dy_k = np.zeros(len(tri[0]))
    # take a sum of each column
    dy_k = np.nansum(tri, axis=0)
    # calculate the development pattern factors
    for i in reversed(range(len(dy_k))):
        if i != 0:
            dy_k[i] = dy_k[i] / dy_k[i-1]
    dy_k[0] = np.nan
    # now apply the development pattern factors to forecast the triangle
    for i in range(len(tri)):
        w = np.where(np.isnan(tri[i]) == True) #check the format of this... output is indices right?
        if len(w) > 0:
            for j in w[0]:
                tri[i][j] = tri[i][j-1] * dy_k[j]
    # now return the squared triangle
    return tri

# Cape Cod Method...specifically, the Clark Method
def capecod(tri, vol, cumul=True, devincre=12):
    from scipy.optimize import curve_fit
    # define ways to compute G(x), two options:
    # using a loglogistic
    def loglogistic(x, omega, theta):
        return (x**(omega)) / ((x**(omega)) + (theta**(omega)))
    # using a Weibull function
    def weibull(x, omega, theta):
        return 1. - np.exp(-1 * ((x / theta)**omega))
    # NOTE: for now, Weibull. If using Loglogistic, data needs to be truncated
    
    # now some checks...do later!
    # check to make sure inputs are numpy arrays, if not, then convert them!
    if type(tri) != np.ndarray or type(vol) != np.ndarray:
        tri = np.asarray(tri)
        vol = np.asarray(vol)
    # make sure that volume and rows match up!
    if len(tri) != len(vol):
        print 'your exposure base and your triangle dimensions do not match!'
        return
    # if no development increment is provided, assume 12
    
    # satisfy data assumptions for cape cod
    # perform math as if triangle is cumulative, conversion if needed
    # also make a copy of the triangle for ELR
    ELR = np.copy(tri)
    mu  = np.copy(tri)
    if cumul == False:
        tri = incre2cumul(tri)
    #    ELR = incre2cumul(ELR)
    #else:
    #    ELR = cumul2incre(ELR)
    
    # normalize each row by its max observed payment
    Ntri = np.copy(tri)
    for i in range(len(Ntri)):
        Ntri[i] = Ntri[i] / np.nanmax(Ntri[i])
    
    # prepare independant variable
    dy = []
    dev = devincre
    for i in range(len(Ntri[0,:])):
        dy.append([dev] * np.nansum(np.isfinite(Ntri[:,i])))
        dev += devincre
    # now prepare format the dependant variable
    y = []
    for i in range(len(Ntri[0,:])):
        y.append(Ntri[np.isfinite(Ntri[:,i]),i].tolist())
    # do curve fitting, estimate the omega and theta parameters respectively
    dy = np.asarray(sum(dy,[]))
    y = np.asarray(sum(y,[]))
    popt, pcov = curve_fit(weibull, dy, y, method='lm')#method='lm', 
    #popt, pcov = curve_fit(loglogistic, dy, y, method='lm')#method='lm', 
    #print loglogistic(120, popt[0], popt[1])
    #popt = [1.274717, 1.829030]
    print '[omega, theta]: ', popt
    # popt = [omega, theta]
    
    # quick check to observe the results of Weilbull curve fitting
    plot = 1
    if plot == 1:
        import matplotlib.pyplot as plt
        trialX = np.linspace(0, np.max(dy), 1000)
        y1 = loglogistic(trialX, 1.477251, 21.4625) 
        y2 = weibull(trialX, 1.441024, 22.3671)
        y3 = weibull(trialX, 1.274717, 1.829030 *12)# R ybc
        y4 = weibull(trialX, popt[0], popt[1])      # dth
        plt.figure()
        plt.title('Growth Functions of Test Triangle')
        plt.xlabel('Loss Development (mo)')
        plt.ylabel('Cumulative Payments (%)')
        plt.plot(trialX, y1, label = 'LDF textbook loglogistic', linestyle='--', color='r')
        plt.plot(trialX, y2, label = 'clark textbook loglogistic')
        plt.plot(trialX, y3, label = 'R weibull', linestyle='--', color='g')
        plt.plot(trialX, y4, label = 'python weibull', linestyle='--', color='m')
        plt.legend(loc='lower right')
        plt.show()
    
    # now figure out the ELR using the MLE approximation
    ELR *= 0.
    for i in range(len(ELR)):
        ELR[i] += vol[i]
    dev = devincre
    for i in range(len(ELR[0,:])):
        #ELR[:,i] *= weibull(dev, popt[0], popt[1]) - weibull(dev - devincre, popt[0], popt[1])
        ELR[:,i] *= weibull((dev + dev - devincre)/2, popt[0], popt[1])
        dev += devincre
    ELR = np.nansum(cumul2incre(tri)) / np.nansum(np.diagonal(np.fliplr(ELR)))
    #ELR = 0.553
    print 'ELR: ', ELR
    
    # finally, do the cape cod calculation
    mu *= 0.
    mu[np.where(np.isnan(mu)==True)] = 0.
    for i in range(len(mu)):
        mu[i] += vol[i] * ELR
    dev = devincre
    for i in range(len(mu[0,:])):
        #mu[i] *=  weibull(dev, popt[0], popt[1]) - weibull(dev - devincre, popt[0], popt[1])
        mu[:,i] *=  weibull((dev + dev - devincre)/2, popt[0], popt[1])
        print 'i: ', i,  weibull((dev + dev - devincre)/2, popt[0], popt[1])
        dev += devincre
    print weibull(60, popt[0], popt[1])
    print weibull(120, popt[0], popt[1])
    # if the input triangle was incremental, return it as such
    if cumul == False:
        mu = cumul2incre(mu)
    return mu

# read triangle from csv as numpy matrix


# tri = pd.read_csv('M:\PC Analytics\ARMS\ARMS MCMC\YB_R\R_triangles_1.csv').as_matrix()
tri = pd.read_csv('M:\PC Analytics\ARMS\ARMS MCMC\YB_R\R_triangles_1.csv')
tri = tri.drop(tri.columns[[0]],1)
tri = tri.as_matrix()

# 1. Testing incre2cumul and cumul2incre
test = cumul2incre(tri)
tri - incre2cumul(test)
df_yb = pd.DataFrame(test)
df_yb.to_csv('M:\PC Analytics\ARMS\ARMS MCMC\YB_R\R_triangles_1_python.csv')

tri = pd.read_csv('C:\Users\dhoang\Documents\GitHub\ARMS\Data\clark.csv')
tri = tri.as_matrix()
test1 = capecod(tri, [5000,5200,5400,5600,5800], cumul=True)
for i in range(len(test1)):
    print np.sum(test1[i])
print test1
print 'ultimtes at 60 mo'
for i in test1[:,4]: print i
print 'total:', np.sum(test1[:,4])






import numpy as np
import matplotlib.pyplot as plt


#b) Cross section
#hc2 = 0.389e9     #hc2 = 38809e4  (hbar^2c^2, Gev^-2 = pb = 2.568e9)
hc2 = 2.56810e9

m = 4.18      # mass b-quark [GeV]
M = 0.10565838 #mass mu [GeV]
mZ = 91.1876    # mass Z-boson [GeV]

e = 0.31345
ee = e**2
alpha = 1/137.
gZ = 0.718082

# coupling constants for Z (see tab. 15.1 on p. 422 in Thomson for weak mixing angle)
cV_mu = -0.04
cA_mu = -0.5
cV_b = -0.35
cA_b = -0.5

cMu_plus = (cV_mu**2 + cA_mu**2)
cMu_minus = (cV_mu**2 - cA_mu**2)
cb_plus = (cV_b**2 + cA_b**2)
cb_minus = (cV_b**2 - cA_b**2)


def QED(E_cm, c):

    if(isinstance(c, (list, tuple, np.ndarray))==1 and isinstance(E_cm, (list, tuple, np.ndarray))==0):
        print "Calculating QED differential cross-section as function of cos(theta)...\n"

        s = E_cm**2         #E_cm = 2E
        E_square = (0.5*E_cm)**2

        p = np.sqrt(E_square - M**2)
        pp = np.sqrt(E_square - m**2)
        p1p2 = E_square + p**2
        p3p4 = E_square + pp**2

        ksi = 3 * (2*np.pi) * 1/float(64*np.pi**2*s) * (pp/float(p)) * hc2

        ds_QED = np.zeros(len(c))

        for i in xrange(len(c)):
            p1p3 = p2p4 = E_square - p*pp*c[i]
            p1p4 = p2p3 = E_square + p*pp*c[i]

#            a[i] = (24*np.pi*alpha**2)/float(9*s)*np.sqrt((s - m**2)/float(s - M**2)) * ((1 + c[i])**2 + (1 - c[i]**2)*(-m**2/float(s) - M**2/float(s)) + M**2*m**2/float(s**2)*c[i]**2)
#            a[i] = np.pi*alpha**2/float(6*s)*np.sqrt(s-m**2)/float(np.sqrt(s-M**2)) * (1 + c[i]**2 + (1 - c[i]**2)*(m**2/float(s) + M**2/float(s)) + m**2*M**2/float(s)*c[i]**2)
            ds_QED[i] = ksi*ee**2/float(9*s**2) * 8*(p1p4*p2p3 + p1p3*p2p4 + m**2*p1p2 + M**2*p3p4 + 2*m**2*M**2)

        return ds_QED #(??)

    elif(isinstance(E_cm, (list, tuple, np.ndarray))==1 and isinstance(c, (list, tuple, np.ndarray))==0):
        print "Calculating QED cross-section as function of s..."

        b = np.zeros(len(E_cm))
        for i in xrange(len(E_cm)):
            b[i] = np.pi*alpha**2/3.*np.sqrt((s[i]-m**2)/float(s[i]-M**2)) * (4/3.*(1+m**2/float(s[i]) + M**2/float(s[i])) + M**2*m**2/float(s[i]))
        return b*hc2 # *2(??)

def Z(s, theta):
    ksi = gZ**4/float(64*np.pi)
    if(isinstance(theta, (list, tuple, np.ndarray))==1 and isinstance(s, (list, tuple, np.ndarray))==0):
        print "Calculating EW differential cross-section as function of cos(theta)...\n"

        a = np.zeros(len(theta))
        for i in xrange(len(theta)):
            b[i] = ksi/float((s-mZ**2)**2) * np.sqrt((s-m**2)/float(s-M**2)) * ( 2*(cMu_plus)*(cb_plus)*(s*(1+np.cos(theta[i])**2) - (m**2 + M**2)*np.cos(theta[i])**2 + 1/float(s)*m**2*M**2*np.cos(theta[i])**2)\
            + m**2*cMu_plus*cb_minus * (2-M**2/float(s))\
            + M**2*cMu_minus*cb_plus * (2-m**2/float(s))\
            + 2/float(s)*M**2*m**2*cMu_minus*cb_minus\
            + 16*cV_mu*cA_mu*cV_b*cA_b * np.sqrt(1-1/float(s)*(m**2 + M**2) + 1/float(s**2)*m**2*M**2)*np.cos(theta[i]) )
        return a*hc2 # *2(??)

    elif(isinstance(s, (list, tuple, np.ndarray))==1 and isinstance(theta, (list, tuple, np.ndarray))==0):
        print "Calculating EW cross-section as function of s..."


E_cm = 10.
s1 = E_cm**2  #GeV
s2 = np.linspace(s1, s1*10, 10000)
c1 = np.linspace(-1.0, 1.0, 100)
c2 = 0.

#QED
plt.plot(c1, QED(s1, c1), 'r-')
plt.title(r'QED, $\frac{d\sigma}{d(cos\theta)}$')
plt.grid('on')
plt.show()

plt.plot(np.sqrt(s2), QED(s2, c2), 'b-')
plt.title('QED')
plt.show()
"""
#Electroweak
#print Z(s, theta)
plt.plot(np.cos(theta1), Z(s1, theta1), 'b-')
plt.title('EW')
plt.show()
"""

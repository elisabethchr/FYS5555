import numpy as np
import matplotlib.pyplot as plt

hc2 = 2.56810e-9    # conversion factor GeV^-2 -> pb

# masses
m = 4.18      # mass b-quark [GeV]
M = 0.10566   # mass mu [GeV]
mZ = 91.1876  # mass Z-boson [GeV]
mH = 125.18   # mass Higgs [GeV]
mW = 80.379   # mass W-bosons [GeV]

# fine-structure and coupling constants
sinW_square = 0.23146; sinW = np.sqrt(sinW_square)
cosW_square = 1 - sinW_square; cosW = np.sqrt(cosW_square)

alpha = 1/137.
e = np.sqrt(4*np.pi*alpha)
ee = e**2
gZ = e/(sinW*cosW)
gW = e/sinW

# decay widths
widthH = 0.013      #0.0061744
widthZ = 2.4952

# coupling constants for Z (see tab. 15.1 on p. 423 in Thomson for weak mixing angle)
cV_mu = -0.04
cA_mu = -0.5
cV_b = -0.35
cA_b = -0.5

cMu_minus = (cV_mu**2 - cA_mu**2)
cMu_plus = (cV_mu**2 + cA_mu**2)
cb_plus = (cV_b**2 + cA_b**2)
cb_minus = (cV_b**2 - cA_b**2)


def CS(E_cm, c):
    if(isinstance(c, (list, tuple, np.ndarray))==1 and isinstance(E_cm, (list, tuple, np.ndarray))==0):
        print "Calculating differential cross-section as function of cos(theta)...\n"

        s = E_cm**2         #E_cm = 2E
        E_square = (0.5*E_cm)**2

        p = np.sqrt(E_square - M**2)*np.ones(np.size(c))
        pp = np.sqrt(E_square - m**2)*np.ones(np.size(c))
        p1p2 = E_square + p**2
        p3p4 = E_square + pp**2
        p1p3 = p2p4 = E_square + p*pp*c
        p1p4 = p2p3 = E_square - p*pp*c

        M_A = np.zeros(len(c))
        M_Z = np.zeros(len(c))
        M_AZ = np.zeros(len(c))
        M_H = np.zeros(len(c))

        ksi = 3 * (2*np.pi)/(64.*np.pi**2*s) * (pp[0]/(p[0])) / hc2


        for i in xrange(len(c)):
            M_H[i] = (m**2*M**2*gW**4)/(4*(s-mH**2-1j*mH*widthH)**2*mW**4) * (p1p2[i]*p3p4[i] - m**2*p1p2[i] - M**2*p3p4[i] + m**2*M**2)

            M_A[i] = (4*np.pi*alpha)**2/float(9*s**2) * 8*(p1p4[i]*p2p3[i] + p1p3[i]*p2p4[i] + m**2*p1p2[i] + M**2*p3p4[i] + 2*m**2*M**2)

            M_Z[i] = gZ**4/(2*((s-mZ**2)**2 + mZ**2*widthZ**2)) * (cMu_plus*cb_plus*(p1p4[i]*p2p3[i] + p1p3[i]*p2p4[i])\
                        + m**2*cMu_plus*cb_minus*p1p2[i]\
                        + M**2*cMu_minus*cb_plus*p3p4[i]\
                        - 4*cV_mu*cA_mu*cV_b*cA_b*(p1p4[i]*p2p3[i] - p1p3[i]*p2p4[i])\
                        + 2*M**2*m**2*cb_plus*cb_minus)
            # M_AZ[i] =

        M_H *= ksi; M_A *= ksi; M_Z *= ksi
        return M_H, M_A, M_Z


E_cm = 150.
s1 = E_cm**2  #GeV
s2 = np.linspace(s1, s1*10, 1e4)
c1 = np.linspace(-1.0, 1.0, 1e3)
c2 = 0.

#Higgs
plt.plot(c1, CS(E_cm, c1)[0], 'r-')
plt.legend([r'$M_H^2$'])
plt.xlabel(r'$\cos\theta$')
plt.ylabel(r'$\frac{d\sigma}{d(\cos\theta})$', size=14)
plt.title(r'$\sqrt{s} = %i$ GeV' %E_cm)
plt.tight_layout()
plt.grid('on')
plt.show()

#QED
plt.plot(c1, CS(E_cm, c1)[1], 'r-')
plt.legend([r'$M_\gamma^2$'])
plt.xlabel(r'$\cos\theta$')
plt.ylabel(r'$\frac{d\sigma}{d(\cos\theta})$', size=14)
plt.title(r'$\sqrt{s} = %i$ GeV' %E_cm)
plt.tight_layout()
plt.grid('on')
plt.show()

#Electroweak
plt.plot(c1, CS(E_cm, c1)[2], 'r-')
plt.legend([r'$M_Z^2$'])
plt.xlabel(r'$\cos\theta$')
plt.ylabel(r'$\frac{d\sigma}{d(\cos\theta})$', size=14)
plt.title(r'$\sqrt{s} = %i$ GeV' %E_cm)
plt.tight_layout()
plt.grid('on')
plt.show()

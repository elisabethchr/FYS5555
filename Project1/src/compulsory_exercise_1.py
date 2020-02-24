import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolor

#-----------------------------
#    import comphep data
#-----------------------------
def importFile(filename, rows):
    data = np.loadtxt(filename, skiprows = rows)
#    a, b = np.shape(data)
    var = data[:, 0]
    dat = data[:, 1]
    return var, dat

A10theta, A10dsigma = importFile("../data/comphep_dsigma_A_10.txt", 3)
A150theta, A150dsigma = importFile("../data/comphep_dsigma_A_150.txt", 3)
H150theta, H150dsigma = importFile("../data/comphep_dsigma_H_150.txt", 3)
Z150theta, Z150dsigma = importFile("../data/comphep_dsigma_Z_150.txt", 3)
tot150theta, tot150dsigma = importFile("../data/comphep_dsigma_total_150.txt", 3)


#-----------------------------
#     define variables
#-----------------------------
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


#--------------------------------------------------------
#   calculate differential cross-sections for A, Z, H
#--------------------------------------------------------
def CS(E_cm, c):
    if(isinstance(c, (list, tuple, np.ndarray))==1 and isinstance(E_cm, (list, tuple, np.ndarray))==0):
        print "\nCalculating differential cross-section as function of cos(theta)...\n"
        a = c

    elif(isinstance(c, (list, tuple, np.ndarray))==0 and isinstance(E_cm, (list, tuple, np.ndarray))==1):
        print "\nCalculating total cross-section as function of center of mass energy...\n"
        a = E_cm

    s = E_cm**2         #E_cm = 2E
    E_square = (0.5*E_cm)**2

    p = np.sqrt(E_square - M**2)*np.ones(np.size(a))
    pp = np.sqrt(E_square - m**2)*np.ones(np.size(a))
    p1p2 = E_square + p**2
    p3p4 = E_square + pp**2
    p1p3 = p2p4 = E_square + p*pp*c
    p1p4 = p2p3 = E_square - p*pp*c

    M_A = np.zeros(len(a))
    M_Z = np.zeros(len(a))
    M_AZ = np.zeros(len(a))
    M_H = np.zeros(len(a))

    ksi = 3 * (2*np.pi)/(64.*np.pi**2*s) * (pp/p) / hc2 * np.ones(len(a))

    # amplitudes
    M_H = (m**2*M**2*gW**4)/(4*(s-mH**2-1j*mH*widthH)**2*mW**4) * (p1p2*p3p4 - m**2*p1p2 - M**2*p3p4 + m**2*M**2)

    M_A = (4*np.pi*alpha)**2/float(9*s**2) * 8*(p1p4*p2p3 + p1p3*p2p4 + m**2*p1p2 + M**2*p3p4 + 2*m**2*M**2)

    M_Z = gZ**4/(2*((s-mZ**2)**2 + mZ**2*widthZ**2)) * (cMu_plus*cb_plus*(p1p4*p2p3 + p1p3*p2p4)\
                + m**2*cMu_plus*cb_minus*p1p2\
                + M**2*cMu_minus*cb_plus*p3p4\
                - 4*cV_mu*cA_mu*cV_b*cA_b*(p1p4*p2p3 - p1p3*p2p4)\
                + 2*M**2*m**2*cMu_minus*cb_minus)   #-4cVcA or +4cVcA??

    M_AZ = (4*np.pi*alpha)*gZ**2/float(3*s*(s-mZ**2)) * 2*(cV_mu*cV_b*(p1p4*p2p3 + p1p3*p2p4)\
                + cV_mu*cV_b*m**2*p1p2\
                + cV_mu*cV_b*M**2*p3p4\
                - cA_mu*cA_b*(p1p4*p2p3 - p1p3*p2p4)\
                +2*m**2*M**2)

    M_H *= ksi; M_A *= ksi; M_Z *= ksi; M_AZ*=ksi
    return M_H, M_A, M_Z, M_AZ


E_cm = 150.
s1 = E_cm**2  #GeV
s2 = np.linspace(s1, s1*10, 1e4)
c1 = np.linspace(-1.0, 1.0, 1e3)
c2 = 0.

#-----------------------------
#         plotting
#-----------------------------
M_H = CS(E_cm, c1)[0]
M_A = CS(E_cm, c1)[1]
M_Z = CS(E_cm, c1)[2]
M_AZ = CS(E_cm, c1)[3]
M_tot = M_A + M_Z + 2*M_AZ

# Higgs
plt.plot(c1, M_H, 'r--')
plt.plot(H150theta, H150dsigma, 'b-')
plt.legend([r'$|M_H|^2$', r'Comphep $|M_H|^2$'])
plt.xlabel(r'$\cos\theta$')
plt.ylabel(r'$\frac{d\sigma}{d(\cos\theta})$', size=14)
plt.title(r'$\sqrt{s} = %i$ GeV' %E_cm)
plt.tight_layout()
plt.grid('on')
plt.show()

# QED
plt.plot(c1, M_A, 'r--')
plt.plot(A150theta, A150dsigma, 'b-')
plt.legend([r'$|M_\gamma|^2$', r'Comphep $|M_\gamma|^2$'])
plt.xlabel(r'$\cos\theta$')
plt.ylabel(r'$\frac{d\sigma}{d(\cos\theta})$', size=14)
plt.title(r'$\sqrt{s} = %i$ GeV' %E_cm)
plt.tight_layout()
plt.grid('on')
plt.show()

# Electroweak
plt.plot(c1, M_Z, 'r--')
plt.plot(Z150theta, Z150dsigma, 'b-')
plt.legend([r'$|M_Z|^2$', r'Comphep $|M_Z|^2$'])
plt.xlabel(r'$\cos\theta$')
plt.ylabel(r'$\frac{d\sigma}{d(\cos\theta})$', size=14)
plt.title(r'$\sqrt{s} = %i$ GeV' %E_cm)
plt.tight_layout()
plt.grid('on')
plt.show()

# Interference term QED, EW
plt.plot(c1, M_AZ, 'r--')
#plt.plot(Z150theta, Z150dsigma, 'b-')
plt.legend([r'$|M_{AZ}|^2$', r'Comphep $|M_Z|^2$'])
plt.xlabel(r'$\cos\theta$')
plt.ylabel(r'$\frac{d\sigma}{d(\cos\theta})$', size=14)
plt.title(r'$\sqrt{s} = %i$ GeV' %E_cm)
plt.tight_layout()
plt.grid('on')
plt.show()


# Total differential cross-section
plt.plot(tot150theta, tot150dsigma, 'r-')
plt.plot(c1, M_tot, 'b--')
plt.plot(A150theta, A150dsigma, 'y--')
plt.plot(Z150theta, Z150dsigma, 'g--')
plt.plot(c1, M_AZ, 'm--')
plt.plot(H150theta, H150dsigma, 'c--')
plt.legend([r'Comphep $|M_{fi}|^2$', r'$|M_{fi}|^2$', r'$|M_A|^2$', r'$|M_Z|^2$', r'$|M_{AZ}|^2$', r'$|M_H|^2$'])
plt.xlabel(r'$\cos\theta$')
plt.ylabel(r'$\frac{d\sigma}{d(\cos\theta})$', size=14)
plt.title(r'$\sqrt{s} = %i$ GeV' %E_cm)
plt.tight_layout()
plt.grid('on')
plt.show()

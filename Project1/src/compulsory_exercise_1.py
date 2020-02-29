import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolor
from scipy.integrate import simps
from scipy import integrate

#-----------------------------
#    import comphep data
#-----------------------------
def importFile(filename, rows):
    data = np.loadtxt(filename, skiprows = rows)
#    a, b = np.shape(data)
    var = data[:, 0]
    dat = data[:, 1]
    return var, dat
# differential cross-sections
A10theta, A10dsigma = importFile("../data/comphep_dsigma_A_10.txt", 3)
A150theta, A150dsigma = importFile("../data/comphep_dsigma_A_150.txt", 3)
H150theta, H150dsigma = importFile("../data/comphep_dsigma_H_150.txt", 3)
Z150theta, Z150dsigma = importFile("../data/comphep_dsigma_Z_150.txt", 3)
tot150theta, tot150dsigma = importFile("../data/comphep_dsigma_total_150.txt", 3)
# total cross-sections
totE_lin, totSigma_lin = importFile("../data/comphep_sigma_10_150_lin.txt", 3)
totE_log, totSigma_log = importFile("../data/comphep_sigma_10_150_log.txt", 3)
# forward-backward asymmetry
asymE_bb, asymSigma_bb = importFile("../data/comphep_Asym_tot_bB_10_150_lin.txt", 3)
asymE_cc, asymSigma_cc = importFile("../data/comphep_Asym_tot_cC_10_150_lin.txt", 3)
asymE_ee, asymSigma_ee = importFile("../data/comphep_Asym_tot_eE_10_150_lin.txt", 3)
#Zprime
ZpE_log, Zsigma_log = importFile("../data/comphep_sigma_Zprime_log.txt", 3)
dSigma_ZpE, ZpdSigma = importFile("../data/comphep_dsigma_Zprime_3000.txt", 3)
asymE_Zp, asymSigma_Zp = importFile("../data/comphep_Asym_Zprime.txt", 3)

#-----------------------------
#     define variables
#-----------------------------
hc2 = 2.56810e-9    # conversion factor GeV^-2 -> pb

# masses
m = [4.18, 1.275, 0.5109989e-3]      # mass b-quark, c, e- [GeV]
M = 0.10566   # mass mu [GeV]
mZ = 91.1876  # mass Z-boson [GeV]
mH = 125.18   # mass Higgs [GeV]
mW = 80.379   # mass W-bosons [GeV]

# fine-structure and coupling constants
sinW_square = 0.23146; sinW = np.sqrt(sinW_square)
#print sinW
cosW_square = 1 - sinW_square; cosW = np.sqrt(cosW_square)

alpha = 1/137.
e = np.sqrt(4*np.pi*alpha)
ee = e**2
gZ = e/(sinW*cosW)
gW = e/sinW

# decay widths
#widthH = 0.013      #0.0061744
widthH = 0.00407
widthZ = 2.4952

# coupling constants for Z (see tab. 15.1 on p. 423 in Thomson for weak mixing angle)
cV_mu = -0.04
cA_mu = -0.5
cV_b = [-0.35, 0.19, -0.04] # b-quark, c-quark, e
cA_b = [-0.5, 0.5, -0.5]    # b-quark, c-quark, e

cMu_minus = (np.square(cV_mu) - np.square(cA_mu))
cMu_plus = (np.square(cV_mu) + np.square(cA_mu))
cb_plus = (np.square(cV_b) + np.square(cA_b))
cb_minus = (np.square(cV_b) - np.square(cA_b))


# Forwards-backwards asymmetry
def AFB(dCS, Ecm):
    Bw = np.trapz(dCS[:int(len(Ecm)/2.)]).real
    Fw = np.trapz(dCS[(int(len(Ecm)/2.)):]).real
    if Fw==0 and Bw==0:
        print Fw, Bw
        Fw = 1e-9
        Bw = 1e-9
    AFB = (Fw - Bw)/float(Fw + Bw)
    return AFB

# Differential cross-section
def dCS(E_cm, c, m, cV_b, cA_b, cb_plus, cb_minus):
    if(isinstance(c, (list, tuple, np.ndarray))==1 and isinstance(E_cm, (list, tuple, np.ndarray))==0):
#        print "\nCalculating differential cross-section as function of cos(theta)...\n"
        a = c

    elif(isinstance(c, (list, tuple, np.ndarray))==0 and isinstance(E_cm, (list, tuple, np.ndarray))==1):
#        print "\nCalculating total cross-section as function of center of mass energy...\n"
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

    ksi = (2*np.pi)/(64.*np.pi**2*s) * (pp/p) / hc2 * np.ones(len(a))

    # amplitudes
    M_H = (m**2*M**2*gW**4)/(4*(s-mH**2-1j*mH*widthH)**2*mW**4) * (p1p2*p3p4 - m**2*p1p2 - M**2*p3p4 + m**2*M**2)

    M_A = (4*np.pi*alpha)**2/float(s**2) * 8*(p1p4*p2p3 + p1p3*p2p4 + m**2*p1p2 + M**2*p3p4 + 2*m**2*M**2)
    M_Z = gZ**4/(2*((s-mZ**2)**2 + mZ**2*widthZ**2)) * (cMu_plus*cb_plus*(p1p4*p2p3 + p1p3*p2p4)\
                + m**2*cMu_plus*cb_minus*p1p2\
                + M**2*cMu_minus*cb_plus*p3p4\
                - 4*cV_mu*cA_mu*cV_b*cA_b*(p1p4*p2p3 - p1p3*p2p4)\
                + 2*M**2*m**2*cMu_minus*cb_minus)   #-4cVcA or +4cVcA??

    M_AZ = (4*np.pi*alpha)*gZ**2/float(s*(s-mZ**2)) * 2*(cV_mu*cV_b*(p1p4*p2p3 + p1p3*p2p4)\
                + cV_mu*cV_b*m**2*p1p2\
                + cV_mu*cV_b*M**2*p3p4\
                - cA_mu*cA_b*(p1p4*p2p3 - p1p3*p2p4)\
                + 2*m**2*M**2*cV_mu*cV_b)

    if m == 4.18:
        M_A *=1/9.      #(charge of b)^2
        ksi *= 3        #colour charge rr, gg, bb
        M_AZ *= 1/3.    #charge of b
    elif m == 1.275:
        M_A *= 4/9.     #(charge of c)^2
        ksi *= 3        #colout charges rr, gg, bb
        M_AZ *= -2/3.   #charge of c

    M_H *= ksi; M_A *= ksi; M_Z *= ksi; M_AZ*=ksi
    return M_H, M_A, M_Z, M_AZ

#------------------------------------------------------
#    Total cross-section calculation and plotting
#------------------------------------------------------
n = 1001
E_cm = np.linspace(10, 150, n)
c = np.linspace(-1.0, 1.0, n)

CS_H = np.zeros((len(m), len(E_cm)))
CS_A = np.zeros((len(m), len(E_cm)))
CS_Z = np.zeros((len(m), len(E_cm)))
CS_AZ = np.zeros((len(m), len(E_cm)))
CS_tot = np.zeros((len(m), len(E_cm)))
#def dCS(E_cm, c, m, cV_b, cA_b, cb_plus, cb_minus):

for j in range(len(m)):
    for i in range(len(E_cm)):
        CS_H[j, i] = simps(dCS(E_cm[i], c, m[j], cV_b[j], cA_b[j], cb_plus[j], cb_minus[j])[0], c)
        CS_A[j, i] = simps(dCS(E_cm[i], c, m[j], cV_b[j], cA_b[j], cb_plus[j], cb_minus[j])[1], c)
        CS_Z[j, i] = simps(dCS(E_cm[i], c, m[j], cV_b[j], cA_b[j], cb_plus[j], cb_minus[j])[2], c)
        CS_AZ[j, i] = simps(dCS(E_cm[i], c, m[j], cV_b[j], cA_b[j], cb_plus[j], cb_minus[j])[3], c)
        CS_tot[j, i] = CS_A[j, i] + CS_Z[j, i] + 2*CS_AZ[j, i]
    print m[j]
print CS_tot[0, -1]

plt.plot(E_cm, CS_tot[0, :], 'b.')
plt.plot(E_cm, CS_tot[1, :], 'r.')
plt.plot(E_cm, CS_tot[2, :], 'g.')
plt.legend([r'$\sigma_{\mu^+\mu^-\rightarrow\overline{b}b}$', r'$\sigma_{\mu^+\mu^-\rightarrow\overline{c}c}$', r'$\sigma_{\mu^+\mu^-\rightarrow e^+e^-}$'])
plt.grid('on')
plt.show()

#linear
plt.plot(totE_lin, totSigma_lin, 'r.')
plt.plot(E_cm, CS_tot[0, :], 'b--')
plt.plot(E_cm, CS_H[0, :], 'y--')
plt.plot(E_cm, CS_A[0, :], 'g--')
plt.plot(E_cm, CS_Z[0, :], 'm--')
plt.xlabel(r'$\sqrt{s}$', size=12)
plt.ylabel(r'$\sigma$ [pb]')
#plt.xlim(10, 150)
plt.legend([r'Numerical $\sigma_{\mu^+\mu^-\rightarrow\overline{b}b}$', r'Analytical $\sigma_{\mu^+\mu^-\rightarrow\overline{b}b}$', r'Higgs', 'QED', r'Z'])
plt.grid('on')
plt.tight_layout()
#plt.savefig('../data/sigma_tot_10_150_lin.pdf')
plt.show()

#logarithmic
plt.plot(totE_log, totSigma_log, 'r.')
plt.loglog(E_cm, CS_tot[0, :], 'b--')
plt.loglog(E_cm, CS_A[0, :], 'g--')
plt.loglog(E_cm, CS_Z[0, :], 'y--')
plt.xlabel(r'$\sqrt{s}$', size=12)
plt.ylabel(r'$\sigma$ [pb]')
#plt.xlim(10, 150)
plt.legend([r'Numerical $\sigma_{\mu^+\mu^-\rightarrow\overline{b}b}$', r'Analytical $\sigma_{\mu^+\mu^-\rightarrow\overline{b}b}$', r'QED', r'Z'])
plt.grid('on')
plt.tight_layout()
#plt.savefig('../data/sigma_tot_10_150_log.pdf')
plt.show()


#-----------------------------------------
#     Forwards-Backwards asymmetry
#------------------------------------------
n = 100
E_cm = np.linspace(10, 150, n)
c = np.linspace(-1.0, 1.0, n)

CS_H = np.zeros(len(E_cm))
CS_A = np.zeros(len(E_cm))
CS_Z = np.zeros(len(E_cm))
CS_AZ = np.zeros(len(E_cm))
CS_tot = np.zeros(len(E_cm))
#def dCS(E_cm, c, m, cV_b, cA_b, cb_plus, cb_minus):

AFB_H = np.zeros((len(m), n))
AFB_A = np.zeros((len(m), n))
AFB_Z = np.zeros((len(m), n))
AFB_AZ = np.zeros((len(m), n))
AFB_tot = np.zeros((len(m), n))
AFB_test = np.zeros((len(m), n))

for j in range(len(m)):
    for i in range(len(E_cm)):
        dCS_H = dCS(E_cm[i], c, m[j], cV_b[j], cA_b[j], cb_plus[j], cb_minus[j])[0]
        dCS_A = dCS(E_cm[i], c, m[j], cV_b[j], cA_b[j], cb_plus[j], cb_minus[j])[1]
        dCS_Z = dCS(E_cm[i], c, m[j], cV_b[j], cA_b[j], cb_plus[j], cb_minus[j])[2]
        dCS_AZ = dCS(E_cm[i], c, m[j], cV_b[j], cA_b[j], cb_plus[j], cb_minus[j])[3]
        dCS_test = dCS_Z + dCS_A
        dCS_tot = dCS_A + dCS_Z + 2*dCS_AZ

        AFB_H[j, i] = AFB(dCS_H, E_cm)
        AFB_A[j, i] = AFB(dCS_A, E_cm)
        AFB_Z[j, i] = AFB(dCS_Z, E_cm)
        AFB_AZ[j, i] = AFB(dCS_AZ, E_cm)
        AFB_test[j, i] = AFB(dCS_test, E_cm)
        AFB_tot[j, i] = AFB(dCS_tot, E_cm)

plt.plot(E_cm, AFB_H[0, :], 'b--')
plt.plot(E_cm, AFB_A[0, :], 'r--')
plt.plot(E_cm, AFB_Z[0, :], 'g--')
#plt.plot(E_cm, AFB_AZ[0, :], 'y--')
plt.plot(E_cm, AFB_test[0, :], 'k--')
plt.xlabel(r'$\sqrt{s}[GeV]$')
plt.ylabel(r'$\sigma[pb]$')
plt.title('Foward-Backward Asymmetry')
plt.legend(['Higgs', 'QED', 'Z', 'QED + Z'])
plt.tight_layout()
plt.grid('on')
plt.savefig('../data/Asym_AZH.pdf')
plt.show()

plt.plot(E_cm, AFB_tot[0, :], 'b.')
plt.plot(E_cm, AFB_tot[1, :], 'r.')
plt.plot(E_cm, AFB_tot[2, :], 'g.')
plt.plot(asymE_bb, asymSigma_bb, 'b--')
plt.plot(asymE_cc, asymSigma_cc, 'r--')
plt.plot(asymE_ee, asymSigma_ee, 'g--')
plt.xlabel(r'$\sqrt{s}[GeV]$')
plt.ylabel(r'$\sigma [pb]$')
plt.title('Forward-Backward Asymmetry')
plt.legend([r'Analytic $\sigma_{\mu^+\mu^-\rightarrow\overline{b}b}$', r'Analytic $\sigma_{\mu^+\mu^-\rightarrow\overline{c}c}$', r'Analytic $\sigma_{\mu^+\mu^-\rightarrow e^+e^-}$',\
            r'Numerical $\sigma_{\mu^+\mu^-\rightarrow\overline{b}b}$', r'Numerical $\sigma_{\mu^+\mu^-\rightarrow\overline{c}c}$', r'Numerical $\sigma_{\mu^+\mu^-\rightarrow e^+e^-}$'])
plt.grid('on')
plt.savefig('../data/Asym_tot_150.pdf')
plt.show()

"""
plt.plot(E_cm, AFB_tot, 'r.')
plt.xlabel(r'$\sqrt{s}$', size=14)
plt.ylabel(r'$\sigma_{FB}[pb]$', size=14)
plt.title('Foward-Backward Asymmetry')
plt.legend(['\sigma_{tot}', 'A', 'Z', 'AZ'])
plt.grid('on')
plt.show()

"""
"""
#---------------------
#   Zprime plotting
#---------------------
plt.plot(dSigma_ZpE, ZpdSigma, 'r--')
plt.xlabel(r'$\cos\theta$')
plt.ylabel(r'$\frac{d\sigma}{d(\cos\theta)}[pb]$')
plt.legend([r'$\mu^+\mu^-\longrightarrow Z^\prime \longrightarrow b\overline{b}$'])
plt.title(r'$\sqrt{s}=3$TeV')
plt.grid('on')
plt.tight_layout()
plt.savefig('../data/diffCS_Zprime_3000.pdf')
plt.show()

plt.plot(ZpE_log, Zsigma_log, 'r--')
plt.xlabel(r'$\sqrt{s}[GeV]$')
plt.ylabel(r'$\sigma[pb]$')
plt.legend([r'$\mu^+\mu^-\longrightarrow Z^\prime \longrightarrow b\overline{b}$'])
plt.title(r'Total cross-section')
plt.yscale('log')
plt.grid('on')
plt.tight_layout()
plt.savefig('../data/CS_Zprime_3000.pdf')
plt.show()

plt.plot(asymE_Zp, asymSigma_Zp, 'r--')
plt.xlabel(r'$\sqrt{s}[GeV]$')
plt.ylabel(r'$\sigma[pb]$')
plt.legend([r'$\mu^+\mu^-\longrightarrow Z^\prime \longrightarrow b\overline{b}$'])
plt.title(r'Forward-Backwards Asymmetry')
plt.grid('on')
plt.tight_layout()
plt.savefig('../data/Asym_Zprime_3000.pdf')
plt.show()
"""
"""
#---------------------------------------------
#    Differential cross-section plotting
#---------------------------------------------
E_cm = 150.
dCS_H = dCS(E_cm, c, m[0], cV_b[0], cA_b[0], cb_plus[0], cb_minus[0])[0]
dCS_A = dCS(E_cm, c, m[0], cV_b[0], cA_b[0], cb_plus[0], cb_minus[0])[1]
dCS_Z = dCS(E_cm, c, m[0], cV_b[0], cA_b[0], cb_plus[0], cb_minus[0])[2]
dCS_AZ = dCS(E_cm, c, m[0], cV_b[0], cA_b[0], cb_plus[0], cb_minus[0])[3]
dCS_tot = dCS_A + dCS_Z + 2*dCS_AZ



# Higgs
plt.plot(c, dCS_H, 'r--')
plt.plot(H150theta, H150dsigma, 'b-')
plt.legend([r'$|M_H|^2$', r'Comphep $|M_H|^2$'])
plt.xlabel(r'$\cos\theta$')
plt.ylabel(r'$\frac{d\sigma}{d(\cos\theta)} [pb]$', size=14)
plt.title(r'$\sqrt{s} = %i$ GeV' %E_cm)
plt.tight_layout()
plt.grid('on')
#plt.savefig('../data/diffCS_H_150.pdf')
plt.show()

# QED
plt.plot(c, dCS_A, 'r--')
plt.plot(A150theta, A150dsigma, 'b-')
plt.legend([r'$|M_\gamma|^2$', r'Comphep $|M_\gamma|^2$'])
plt.xlabel(r'$\cos\theta$')
plt.ylabel(r'$\frac{d\sigma}{d(\cos\theta)} [pb]$', size=14)
plt.title(r'$\sqrt{s} = %i$ GeV' %E_cm)
plt.tight_layout()
plt.grid('on')
#plt.savefig('../data/diffCS_A_150.pdf')
plt.show()

# Electroweak
plt.plot(c, dCS_Z, 'r--')
plt.plot(Z150theta, Z150dsigma, 'b-')
plt.legend([r'$|M_Z|^2$', r'Comphep $|M_Z|^2$'])
plt.xlabel(r'$\cos\theta$')
plt.ylabel(r'$\frac{d\sigma}{d(\cos\theta)} [pb]$', size=14)
plt.title(r'$\sqrt{s} = %i$ GeV' %E_cm)
plt.tight_layout()
plt.grid('on')
#plt.savefig('../data/diffCS_Z_150.pdf')
plt.show()

# Interference term QED, EW
plt.plot(c, dCS_AZ, 'r--')
#plt.plot(Z150theta, Z150dsigma, 'b-')
plt.legend([r'$|M_{\gamma Z}|^2$', r'Comphep $|M_Z|^2$'])
plt.xlabel(r'$\cos\theta$')
plt.ylabel(r'$\frac{d\sigma}{d(\cos\theta)} [pb]$', size=14)
plt.title(r'$\sqrt{s} = %i$ GeV' %E_cm)
plt.tight_layout()
plt.grid('on')
#plt.savefig('../data/diffCS_AZ_150.pdf')
plt.show()

# Differential cross-sections, excluding interference terms
plt.plot(c, dCS_A, 'g--')
plt.plot(c, dCS_Z,'m--')
plt.plot(c, dCS_H, 'y--')
plt.xlabel(r'$\cos\theta$')
plt.ylabel(r'\frac{d\sigma}{d\cos\theta}[pb]')
plt.legend([r'$|M_\gamma|$', r'$|M_H|$', r'$|M_Z|$'])
plt.title(r'$\sqrt{s} = %i$ GeV' %E_cm)
plt.tight_layout()
plt.grid('on')
plt.savefig('../data/diffCS_H_A_Z.pdf')
plt.show()

# Total differential cross-section
plt.plot(tot150theta, tot150dsigma, 'r-')
plt.plot(c, dCS_tot, 'b--')
plt.plot(c, dCS_A, 'g--')
plt.plot(c, dCS_Z,'m--')
#plt.plot(A150theta, A150dsigma, 'y--')
#plt.plot(Z150theta, Z150dsigma, 'g--')
plt.plot(c, dCS_AZ, 'c--')
plt.plot(c, dCS_H, 'y--')
#plt.plot(H150theta, H150dsigma, 'c--')
plt.legend([r'Comphep $|M_{tot}|^2$', r'$|M_{tot}|^2$', r'$|M_{\gamma}|^2$', r'$|M_Z|^2$', r'$|M_{\gamma Z}|^2$', r'$|M_H|^2$'])
plt.xlabel(r'$\cos\theta$')
plt.ylabel(r'$\frac{d\sigma}{d(\cos\theta)} [pb]$', size=14)
plt.title(r'$\sqrt{s} = %i$ GeV' %E_cm)
plt.tight_layout()
plt.grid('on')
#plt.savefig('../data/diffCS_tot_150.pdf')
plt.show()
"""

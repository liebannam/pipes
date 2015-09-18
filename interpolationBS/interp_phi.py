import math
from scipy import integrate, optimize
import numpy as np
import matplotlib.pyplot as plt
import numpy. polynomial.chebyshev as cheb
from matplotlib import rc


iphi = lambda x: np.sqrt(9.8 / 8.) * (1. - np.cos(x)) / \
    np.sqrt((x - np.sin(x)) * np.sin(x / 2.))


def getChebNodes(N):
    return -np.cos(np.array(range(N + 1)) * np.pi / N)


def phi(t):
    tol = 1e-5  # ~ (machine epsilon)^(1/3)
    c = np.sqrt(9.8 * 3. / 8.)
    s = len(t)
    y = np.zeros(len(t))
    for i in range(s):
        if t[i] <= tol:
            y[i] = c * t[i]
        else:
            yt, err = integrate.quad(iphi, tol, t[i], epsabs=1e-10)
            y[i] = c * tol + yt
    return y


def getTheta(x):
    return [optimize.ridder(lambda x: 1. / 8. * (x - np.sin(x)) - xi, 0, 2 * np.pi) for xi in x]


def fX(a, p, a0, af):
    return 2 * pow((a - a0) / (af - a0), p) - 1


def fA(x, p, a0, af):
    return (af - a0) * pow((x + 1.) / 2., 1. / p) + a0


def main():
    M = 3 * pow(2, 5)
    N = 15
    # good values
    a1 = np.pi / 8.
    a2 = np.pi / 4.

    alphas = np.array([i / float(M) for i in range(1, M + 3)])
    r = np.zeros((len(alphas), 2))
    err = np.zeros((len(alphas), 2))
    A1 = np.linspace(0, a1, 50)
    T1 = getTheta(A1)
    Phitrue1 = phi(T1)
    A2 = np.linspace(a1, a2, 50)
    T2 = getTheta(A2)
    Phitrue2 = phi(T2)
    x = getChebNodes(N)

    for k in range(len(alphas)):
        p = alphas[k]
        print x
        #ax1 = [np.pi/8.*pow((xi+1.)/2,1./p) for xi in x]
        ax1 = fA(x, p, 0, a1)
        ax2 = fA(x, p, a2, a1)
        #ax2 = [np.pi/8.+np.pi/8.*pow((xi+1.)/2,1./p) for xi in x]
        theta1 = getTheta(ax1)
        theta2 = getTheta(ax2)
        # print "ax2"
        #print (ax2-1./8.*(theta2-np.sin(theta2)))
        # print "Theta"
        # print theta
        #print (1./8*(theta-np.sin(theta))-ax)
        phi1 = phi(theta1)
        phi2 = phi(theta2)
        fp1 = cheb.chebfit(x, phi1, N)
        fp2 = cheb.chebfit(x, phi2, N)
        print "Chebyshev Coeffs"
        print fp2
        r[k, 0] = abs(fp1[-1])
        r[k, 1] = abs(fp2[-1])
        err[k, 0] = np.linalg.norm(
            Phitrue1 - cheb.chebval(fX(A1, p, 0, a1), fp1))
        err[k, 1] = np.linalg.norm(
            Phitrue2 - cheb.chebval(fX(A2, p, a2, a1), fp2))
    print "Nth coefficient"
    print r
    print "Error"
    print err
    # print alphas
    p1 = 1. / 3.
    p2 = 5. / 12.
    ax1 = fA(x, p1, 0., a1)
    ax2 = fA(x, p2, a2, a1)
    theta1 = getTheta(ax1)
    theta2 = getTheta(ax2)
    phi1 = phi(theta1)
    phi2 = phi(theta2)
    fp1 = cheb.chebfit(x, phi1, N)
    fp2 = cheb.chebfit(x, phi2, N)
    print fp1

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')

    plt.subplots_adjust(hspace=0.4)
    plt.subplot(211)
    plt.semilogy(alphas, r)
    plt.legend([r"$\phi_1(A)$", r"$\phi_2(A)$"])
    xlab = ['1/6', '1/4', '1/3', '5/12', '1/2', '2/3', '3/4', '5/6', '1']
    ticks = [1. / 6., 1. / 4., 1. / 3., 5. / 12.,
             1. / 2., 2. / 3., 3. / 4., 5. / 6., 1.]
    plt.xticks(ticks, xlab)
    plt.grid(True)
    plt.title(
        r'$\phi(x) \approx \tilde{\phi}(x) =\sum_k a_k T_k(cx^{\alpha}-1)$')
    plt.ylabel('coefficient decay')
    # plt.ylabel(r'a_N')
    plt.xlabel(r'$\alpha$')
    plt.subplot(212)
    plt.xlabel(r'$\alpha$')
    # plt.ylabel(r'$\phi_i(x)-\phi_i^N(x)$')
    plt.semilogy(alphas, err)
    plt.xticks(ticks, xlab)
    plt.grid(True)
    # plt.title('Error')
    plt.ylabel('Error')
    plt.legend([r"$\phi_1(A)$", r"$\phi_2(A)$"])
    f = open("blah.txt", "w")
    MM = 200
    t = np.linspace(0, 2 * np.pi, MM)
    aa = 1. / 8. * (t - np.sin(t))
    phit = phi(t)
    phit1 = cheb.chebval(fX(aa[0:MM / 2], p1, 0., a1), fp1)
    phit2 = cheb.chebval(fX(aa[MM / 2:], p2, a2, a1), fp2)
    for k in range(len(t)):
        if k < MM / 2:
            ha = phit1[k]
        else:
            ha = phit2[k - MM / 2]
        f.write("%s   %s   %s\n" % (aa[k], phit[k] - ha, ha))
    f.close()
    plt.show()
    print phi([2 * np.pi])
if __name__ == "__main__":
    main()

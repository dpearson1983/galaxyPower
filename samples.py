import numpy as np
import matplotlib.pyplot as plt

N_s = 20000

def sinc(x):
    if (x != 0):
        return np.sin(x)/x
    else:
        return 1.0

A = 2.0
c = 3.0
w = 1.0
phi = np.pi/7.0
theta = np.linspace(0, 2*np.pi, 30)
x = A*np.sin(w*theta + phi) + c
np.savetxt("truth.dat", np.column_stack((theta,x)))

C = np.zeros((30,30))

for i in range(0,30):
    for j in range(0,30):
        if (i == j):
            C[i][j] = 0.10*x[i]
            print(C[i][j])
        else:
            C[i][j] = 0.10*np.sqrt(x[i]*x[j])*sinc(10.0*abs(theta[i]-theta[j]))

np.savetxt("covariance.dat", C)
L = np.linalg.cholesky(C)

np.savetxt("CholeskyDecomp.dat",L)

x_avg = np.zeros(30)
S = np.zeros((30,30))
S_t = np.zeros((30,30))
for i in range(0,N_s):
    R = np.random.normal(0.0, 1.0, 30)
    y = np.dot(L,R)
    x_i = x + y
    x_avg += x_i/float(N_s)
    filename = "sample" + str(i+1) + ".dat"
    np.savetxt(filename, np.column_stack((theta, x_i)))

errorbars = np.zeros(30)
for k in range(0,N_s):
    filename = "sample" + str(k+1) + ".dat"
    x_i = np.loadtxt(filename)
    for i in range(0,30):
        for j in range(0,30):
            S[i][j] += (x_i[i][1] - x_avg[i])*(x_i[j][1] - x_avg[j])
            S_t[i][j] += (x_i[i][1] - x_avg[i])*(x_i[j][1] - x_avg[j])/float(N_s-1)
        errorbars[i] = np.sqrt(S[i][i]/float(N_s-1))

S /= float(N_s - 1)
np.savetxt("sampleCovariance.dat", S)
np.savetxt("sampleCovariance2.dat", S_t)
np.savetxt("sampleAverage.dat", np.column_stack((theta,x_avg,errorbars)))

Q = (C-S)

delta = np.zeros(30*30);
dummy = np.zeros(30*30);
for i in range(0,30):
    for j in range(0,30):
        index = j + 30*i
        delta[index] = Q[i][j]
        dummy[index] = index+1

np.savetxt("covarDiffLin.dat", np.column_stack((dummy,delta)))
chisq = 0
count = 0
for i in range(0,30):
    for j in range(i,30):
        chisq += (S[i][j] - C[i][j])**2
        count += 1

print(chisq)

np.savetxt("covarianceDiff.dat", Q)

#plt.plot(theta, x, "o")
#plt.plot(theta, x_avg, "x")
#plt.show()

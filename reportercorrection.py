
import numpy as np
import scipy.optimize
import sys

labelnames = """
126C
127N
127C
128N
128C
129N
129C
130N
130C
131N
131C
""".split()
corrections = np.zeros(shape=(len(labelnames),9))
for i,l in enumerate(open(sys.argv[1])):
    if l.startswith('>'):
        continue
    sl = l.split()
    label = sl[0]
    correction = list(map(float,sl[1:]))
    # print(label,correction)
    correction = [ correction[0], 0.0, correction[1], 0.0, 100.0, 0.0, correction[2], 0.0, correction[3] ]
    corrections[labelnames.index(label),:] = correction

# print(corrections)

corrections /= corrections.sum(axis=1,keepdims=True)

# print(corrections)

correctionmat = np.zeros(shape=(len(labelnames),len(labelnames)+2))
for i in range(len(labelnames)):
    for k in range(0,9):
        j = i+(-4+k)
        if j >= 0 and j < correctionmat.shape[1]:
            correctionmat[i,j] = corrections[i,k]
correctionmat = correctionmat.transpose()

samples = 5000+10000*np.random.random((11,))
samples[2] = 0
samples[7] = 0

print(" ".join(map(lambda f: "%.2e"%f,samples)))

obs = np.dot(correctionmat,samples)

print(" ".join(map(lambda f: "%.2e"%f,obs)))
obs[obs==min(obs[[2,7]])] = 0
print(" ".join(map(lambda f: "%.2e"%f,obs)))

samples1,rnorm = scipy.optimize.nnls(correctionmat,obs)

print(" ".join(map(lambda f: "%.2e"%f,samples1)))



#! /usr/bin/env python3

import json
import numpy as np
import matplotlib.pyplot as plt

with open('_korali_result/latest') as f:
    d = json.load(f)

var ='InfectedDetailed'
percentages = [0.8, 0.9, 0.99]
ns = 100

N  = len(d['Samples'])
Ny = len( d['Samples'][0][var] )


x = range(Ny)

y = np.zeros((N,Ny))
s = np.zeros((N,Ny))
for k in range(N):
    y[k,:] = d['Samples'][k][var]
    s[k,:] = d['Samples'][k]['Sigma']

samples = np.zeros((N*ns,Ny))

for k in range(Ny):
    m = y[:,k]
    r = s[:,k]
    yy = [ np.random.normal(m,r) for _ in range(ns) ]
    samples[:,k] = np.asarray(yy).flatten()

mean   = np.zeros((Ny,1))
median = np.zeros((Ny,1))
for k in range(Ny):
  median[k] = np.quantile( samples[:,k],0.5)
  mean[k] = np.mean( samples[:,k] )

fig, ax = plt.subplots(1, 1)

for p in np.sort(percentages)[::-1]:
  q1 = np.zeros((Ny,));
  q2 = np.zeros((Ny,));
  for k in range(Ny):
    q1[k] = np.quantile( samples[:,k],0.5-p/2)
    q2[k] = np.quantile( samples[:,k],0.5+p/2)
  ax.fill_between( x, q1 , q2,  alpha=0.5, label=f' {100*p:.1f}% credible interval' )

ax.plot( x, mean, '-', lw=2, label='Mean', color='black')
ax.plot( x, median, '--', lw=2, label='Median', color='black')

data = [1, 2, 12, 22, 32, 45, 56, 87, 120, 181]
days = [1, 10, 20, 30, 40, 50, 60, 70, 80, 90]
ax.plot(days, data, 'or', label='Reference Data')

ax.legend(loc='upper left')
ax.grid()
ax.set_xlim(left=x[1])

plt.show()
plt.savefig('propagation.jpg')
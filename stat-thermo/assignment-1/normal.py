import numpy as np

#imput parameters
x_mean=1.
standard_deviation=2.
sample_size=10000

# numpy function to sample the normal distribution
x_norm=np.random.normal(x_mean,standard_deviation,sample_size)
#
f=open('x_normal.dat','w')
for x in x_norm:
    f.write(str(x)+'\n')
f.close()
# histogram of x
nbins=100
(h,x)=np.histogram(x_norm,nbins)
fhisto=open('x_histo.dat','w')
for i in range(len(h)):
    fhisto.write(str(x[i])+' '+str(h[i])+' '+'\n')
fhisto.close()

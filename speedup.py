
# Example output
# Creates a python script that will plot the scaling data

#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt

d=np.array([
[ 8 ,   0.0010 ],
[ 6 ,   0.0010 ],
[ 4 ,   0.0010 ],
[ 2 ,   0.0030 ],
[ 1,    0.0050 ]
])

plt.plot(d[:,0],d[:,0]/d[-1,0],"o",color="k", label="Ideal")
plt.plot(d[:,0],d[-1,1]/d[:,1],"x",color="r", label="CAF")

plt.xlabel("Number of Processors")
plt.ylabel("Speedup")
plt.tight_layout( )
plt.savefig('speedup.pdf')

plt.clf( )
plt.plot(d[:,0], d[:,1], "x",color="k", label="CAF")
plt.xlabel("Number of Processors")
plt.ylabel("Run time (s)")
plt.legend( )
plt.tight_layout( )
plt.savefig('runtime.pdf')

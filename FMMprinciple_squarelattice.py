import numpy as np
from matplotlib import pyplot as plt
# np.seterr('raise') # used this to track down a problem when dealing with np.inf

# parameters
M = 20 # number rows
N = 20 # number columns
h = 0.1 # spacing between lattice sites

# function to test whether point (tuple t) is in domain
def indomain(t):
    if (t[0] >= 0) and (t[0] < M) and (t[1] >= 0) and (t[1] < N):
        return True
    else:
        return False

# compute arrival time for point (tuple t)
def computeT(t):
    i = t[0]
    j = t[1]
    # set algorithm, comment out the appr
    # algo = 'KimmelSethian1996' # mostly based on https://escholarship.org/uc/item/7kx079v5
    # algo = 'KimmelSethian1998' # mostly based on https://doi.org/10.1073/pnas.95.15.8431
    algo = 'KimmelSethian1998'
    if algo == 'KimmelSethian1996':
        # find a and b, use fact that we can compare to np.inf
        a = np.inf
        for di in [-1,1]:
            if indomain((i+di,j)):
                a=np.min([a,T[i+di,j]])
        b = np.inf
        for dj in [-1,1]:
            if indomain((i,j+dj)):
                b=np.min([b,T[i,j+dj]])
        if (h/F[i,j] > np.abs(a-b)):
            return (a+b+np.sqrt(2*(h/F[i,j])**2-(a-b)**2))/2.0
        else:
            return h/F[i,j]+np.min([a,b])
    elif algo == 'KimmelSethian1998':
        minT = np.inf
        # loop over all pairs of neighbours
        for di in [-1,1]:
            for dj in [-1,1]:
                if indomain((i+di,j)) and indomain((i,j+dj)):
                    # slightly different variable names to avoid confusion with reference
                    Talpha = T[i+di,j]
                    Tgamma = T[i,j+dj]
                    # use short-ciruit evaluation in the following line,
                    # https://www.pythoninformer.com/python-language/intermediate-python/short-circuit-evaluation/
                    if (np.isfinite(Talpha) and np.isfinite(Tgamma)) and (h/F[i,j] > np.abs(Talpha-Tgamma)):
                        auxT = (Talpha+Tgamma+np.sqrt(2*(h/F[i,j])**2-(Talpha-Tgamma)**2))/2.0
                    else:
                        auxT = h/F[i,j]+np.min([Talpha,Tgamma])
                    minT = np.min([minT,auxT])
        return minT
    else:
        print('Invalid algorithm.')
        return np.nan

# set up arrays
F = np.full((M, N), 1.0) # array / matrix of speeds
T = np.full((M, N), np.inf) # array / matrix of arrival times
# here is the place to set the landscape
# F...

# initial conditions
# alive / accepted
alivedict = {}
# narrow
narrowdict = {}
# set the points with T=0, they go into narrowdict and T and will be set to accepted later
# (see https://escholarship.org/uc/item/7kx079v5)
for i in [0]:
#for i in range(M):
    T[i,0] = 0.0
    narrowdict[(i,0)] = T[i,0]

# test computeT function
# computeT((0,1))
print(narrowdict)

# marching
while bool(narrowdict): # while there are points in narrow band
    # find point from narrow dictionary that has lowest arrival time and act on it
    # https://stackoverflow.com/questions/3282823/get-the-key-corresponding-to-the-minimum-value-within-a-dictionary
    nminidx = min(narrowdict, key=narrowdict.get)
    i = nminidx[0]
    j = nminidx[1]
    alivedict[nminidx] = narrowdict[nminidx]
    T[i,j] = alivedict[nminidx]
    narrowdict.pop(nminidx)
    print('New point added (i, j, T[i,j]): ('+str(i)+', '+str(j)+', '+str(T[i,j])+')')
    # identify neighbors and update times
    neighbors = []
    for d in [(-1,0), (+1,0), (0,-1), (0,+1)]:
        di = d[0]
        dj = d[1]
        # only update if in domain and not alive / accepted yet
        if indomain((i+di,j+dj)) and (not ((i+di,j+dj) in alivedict)):
            neighbors.append((i+di,j+dj))
    for neighbor in neighbors:
        ni = neighbor[0]
        nj = neighbor[1]
        narrowdict[(ni,nj)]=np.min([T[ni,nj],computeT((ni,nj))])
    # now also update T, couldn't do this before
    for neighbor in neighbors:
        T[neighbor[0],neighbor[1]]=narrowdict[(neighbor[0],neighbor[1])]

# looking at the output
plt.imshow(T)
print(T[-1,0])
print(T[:,-1])
print(T[-1,-1])
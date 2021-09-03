import numpy as np
A=input()
[n , inp]=A.split('\n')
n=int(n)
inp=inp.split()
inp=list(map(int , inp))
counter=0;
out=np.ones([n])*-1
h=np.count_nonzero(out==-1)
while(h>2):
    counter+=1
    hap=0
    j=1
    for i in range(n):
        if out[i]!=-1:
            continue
        if i==n-1:
            if j>2:
                cc=0
                ii=i
                while(cc<j):
                    if out[ii-cc]==-1:
                        out[ii-cc]=counter
                        cc+=1
                    else:
                        ii=ii-1
                h=np.count_nonzero(out==-1)
                hap=1
            else:
                continue
        for k in range(i+1 , n):
            if out[k]==-1:
                b=inp[k]
                break
            else:
                b=-1
        if b==-1: 
            if j>2:
                cc=0
                ii=i
                while(cc<j):
                    if out[ii-cc]==-1:
                        out[ii-cc]=counter
                        cc+=1
                    else:
                        ii=ii-1
                    h=np.count_nonzero(out==-1)
                    hap=1
            continue
        if inp[i]==b:
            j+=1
        elif j>2:
            cc=0
            ii=i
            while(cc<j):
                if out[ii-cc]==-1:
                    out[ii-cc]=counter
                    cc+=1
                else:
                    ii=ii-1
            h=np.count_nonzero(out==-1)
            hap=1
            j=1
        else:
            j=1
    if hap==0:
        break
print(out)
            
        
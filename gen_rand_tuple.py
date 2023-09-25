import numpy as np
import random
import sys

if __name__ == "__main__":
    if (len(sys.argv)==3):
        n=int(sys.argv[1])
        n_max=int(sys.argv[2])
        X=random.sample(range(1,n_max+1),k=n)
        for i in X :
            print(i)

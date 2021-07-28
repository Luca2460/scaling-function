from types import MethodType
import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt
from math import sqrt
import json
import pandas as pd



# with open('cu2oseo3_scaling_H2000_Oe.txt') as f:
#     data = f.read()
#     f.close()

# print(data[0])
# print('***********\n')
# print(data[1], data[2], data[3])

data = pd.read_csv('cu2oseo3_scaling_H2000_Oe.txt', sep="\t", header=None)

x = list(data[1])
y = list(data[0])

plt.figure()
for xi, yi in zip(x, y):
    plt.scatter(xi, yi)
plt.xscale('log')

test = input("\n ************")
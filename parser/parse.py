import re

def rreplace(s, old, new):
    return (s[::-1].replace(old[::-1],new[::-1], 1))[::-1]

# Open a file: file
file = open('gustavo',mode='r')

# read all lines at once
all_of_it = file.read().replace("[*", "", 1)
tmp = all_of_it[:-2] 
print(tmp)


import re
#A, p, l, k, EC, j-inv, ratio, point, image
def rreplace(s, old, new):
    return (s[::-1].replace(old[::-1],new[::-1], 1))[::-1]

def get_parameters(st):
    res = re.sub(r'EllipticCurve(.*?)', '', st)
    broken_lines = res.split(",")
    A = int(broken_lines[0])
    p = int(broken_lines[1])
    l = int(broken_lines[2])
    k = int(broken_lines[3])
    j_inv = int(broken_lines[9])
    ratio = int(broken_lines[10])
    param = (A, p, l, k, j_inv, ratio)
    return param

# Open a file: file
file = open('gustavo',mode='r')

# read all lines at once
all_of_it = file.read().replace("[*", "", 1)
all_of_it = all_of_it.replace(" ", "")

tmp = all_of_it[:-2]
pat = r'\[\*|\*\]'
match = re.split(pat,tmp)
parameters = []
for e in match:
    new_str = e.replace("\n", "").replace("\\", "")
    if len(new_str) > 2:
        param = get_parameters(new_str)
        parameters.append(param)

to_save = open("params.txt", "w")
for param in parameters:
    str_f = str(param)+"\n"
    to_save.write(str_f)
    print(param)
to_save.close()

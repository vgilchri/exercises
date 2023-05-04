import re

parameters = []
with open("database_2", "r") as file:
    line = file.readline()
    while line:
        if "k is" in line:
            t = []
            k = int(re.findall(r'\b\d+\b', line)[0])
            t.append(k)
            for i in range(3):
                line = file.readline().strip()
                v = int(re.findall(r'\b\d+\b', line)[0])
                t.append(v)
            parameters.append(t)
        line = file.readline()

to_save = open("params.txt", "w")
for param in parameters:
    str_f = str(param)+"\n"
    to_save.write(str_f)
    print(param)
to_save.close()

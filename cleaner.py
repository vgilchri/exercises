import re
file = open('database',mode='r')
all_of_it = file.read().replace("problem", "")
all_of_it = all_of_it.replace("arisedsome", "")
all_of_it = all_of_it.replace("error", "")
all_of_it = all_of_it.replace("some", "")
all_of_it = all_of_it.replace("arised", "")
filtered = "\n".join([ll.rstrip() for ll in all_of_it.splitlines() if ll.strip()])


to_save = open("database_2", "w")

to_save.write(filtered)
to_save.close()

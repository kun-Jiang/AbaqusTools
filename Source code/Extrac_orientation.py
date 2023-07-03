import os

orientation_file = open('orientation.txt', 'r')
orientation_file_lines = orientation_file.readlines()[1:]
# print(orientation_file_lines)
Phase = []
phi1  = []
phi2  = []
phi3  = []
for line in orientation_file_lines:
    line = line.strip('\n')
    line = line.split('\t')

    Phase.append(line[1])
    phi1.append(float(line[2]))
    phi2.append(float(line[3]))
    phi3.append(float(line[4]))
print(phi1)
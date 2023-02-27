import os, sys, fnmatch

# for parametro in sys.argv:
#     print(parametro)
path = os.listdir(sys.argv[1])
fastas = fnmatch.filter(path, '*.fasta')
print(fastas)

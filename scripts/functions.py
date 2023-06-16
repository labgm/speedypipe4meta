from file import File
import csv

def read_file(path):

    # Open the file report.txt
    tsv_file = open(path)
    read_tsv = csv.reader(tsv_file, delimiter="\t")

    cont = 0 # Controle de linhas
    montadores = [] # Lista de Objetos (Montadores)
    for row in read_tsv:
        if cont == 0:
            objects = len(row)-1
            for i in range(objects):
                montadores.append(File(row[i+1]))

        elif cont == 13:
            for i in range(objects):
                montadores[i].contigs = float(row[i+1].replace("\n",""))

        elif cont == 15:
            for i in range(objects):
                montadores[i].total_length = float(row[i+1].replace("\n",""))

        elif cont == 17:
            for i in range(objects):
                montadores[i].n50 = float(row[i+1].replace("\n",""))

        elif cont == 18:
            for i in range(objects):
                montadores[i].n75 = float(row[i+1].replace("\n",""))

        elif cont == 19:
            for i in range(objects):
                montadores[i].l50 = float(row[i+1].replace("\n",""))

        elif cont == 20:
            for i in range(objects):
                montadores[i].l75 = float(row[i+1].replace("\n",""))
        cont += 1
    return montadores

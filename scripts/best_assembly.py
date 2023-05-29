import argparse
from functions import read_file
from itertools import groupby
import subprocess
import os

parser = argparse.ArgumentParser()
parser.add_argument('--input', type=str, default='', help='Input File')
parser.add_argument('--a1', type=str, default='N50', help='Input Argument 01')
parser.add_argument('--a2', type=str, default='', help='Input Argument 02')
parser.add_argument('--a3', type=str, default='', help='Input Argument 03')
parser.add_argument('--a4', type=str, default='', help='Input Argument 04')
parser.add_argument('--a5', type=str, default='', help='Input Argument 05')
parser.add_argument('--a6', type=str, default='', help='Input Argument 06')
parser.add_argument('--a7', type=str, default='', help='Input Argument 07')
args = parser.parse_args()


alist = read_file(args.input)

def N50():
    alist.sort(key=lambda x: x.n50, reverse=True)
    colocacao = 1
    score = len(alist)
    for file, elementos in groupby(alist, key=lambda x: x.n50):
        for elemento in elementos:
            elemento.score = elemento.score + score
        colocacao += 1
        score -= 1

def L50():
    alist.sort(key=lambda x: x.l50, reverse=False)
    colocacao = 1
    score = len(alist)
    for file, elementos in groupby(alist, key=lambda x: x.l50):
        for elemento in elementos:
            elemento.score = elemento.score + score
        colocacao += 1
        score -= 1

def contigs():
    alist.sort(key=lambda x: x.contigs, reverse=False)
    colocacao = 1
    score = len(alist)
    for file, elementos in groupby(alist, key=lambda x: x.contigs):
        for elemento in elementos:
            elemento.score = elemento.score + score
        colocacao += 1
        score -= 1

def largest():
    alist.sort(key=lambda x: x.largest_contigs, reverse=False)
    colocacao = 1
    score = len(alist)
    for file, elementos in groupby(alist, key=lambda x: x.largest_contigs):
        for elemento in elementos:
            elemento.score = elemento.score + score
        colocacao += 1
        score -= 1

def total():
    alist.sort(key=lambda x: x.total_length, reverse=False)
    colocacao = 1
    score = len(alist)
    for file, elementos in groupby(alist, key=lambda x: x.total_length):
        for elemento in elementos:
            elemento.score = elemento.score + score
        colocacao += 1
        score -= 1

def N75():
    alist.sort(key=lambda x: x.n75, reverse=True)
    colocacao = 1
    score = len(alist)
    for file, elementos in groupby(alist, key=lambda x: x.n75):
        for elemento in elementos:
            elemento.score = elemento.score + score
        colocacao += 1
        score -= 1

def L75():
    alist.sort(key=lambda x: x.l75, reverse=False)
    colocacao = 1
    score = len(alist)
    for file, elementos in groupby(alist, key=lambda x: x.l75):
        for elemento in elementos:
            elemento.score = elemento.score + score
        colocacao += 1
        score -= 1

# Pegar Valor dos Argumentos e executa como função
locals()[args.a1]()
if args.a2:
    locals()[args.a2]()
if args.a3:
    locals()[args.a3]()
if args.a4:
    locals()[args.a4]()
if args.a5:
    locals()[args.a5]()
if args.a6:
    locals()[args.a6]()
if args.a7:
    locals()[args.a7]()

# Ordenação dos Montadores(alist) do maior para o menor score
alist.sort(key=lambda x: x.score, reverse=True)

# if(os.path.isdir(alist[0].nome)):
print('Melhor montagem: ', alist[0].nome.upper())
best_assembly = alist[0].nome.upper()

path = args.input.replace('/metaquast/combined_reference/report.tsv','')
os.mkdir(f'{path}/assembly/best_assembly')
if best_assembly == 'IDBA':
    os.system(f'cp {path}/assembly/idba/scaffold.fa {path}/assembly/best_assembly/final.contigs.fa')
    print('Melhor montagem é o IDBA')
if best_assembly == 'SPADES':
    os.system(f'cp {path}/assembly/metaspades/scaffolds.fasta {path}/assembly/best_assembly/final.contigs.fa')
    print('Melhor montagem é o SPADES')     
if best_assembly == 'MEGAHIT':
    print('Melhor montagem é o MEGAHIT')
    os.system(f'cp {path}/assembly/megahit/final.contigs.fa {path}/assembly/best_assembly/final.contigs.fa')
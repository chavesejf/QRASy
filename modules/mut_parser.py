#!/bin/python3
# -----------------------------
# Written by Elton J. F. Chaves
# E-mail: chavesejf@gmail.com
# -----------------------------
class FileReader:
    def __init__(self, filename):
        self.filename = filename
        self.res = {}

    def read_file(self):
        with open(self.filename, 'r') as file:
            for line in file:
                chain, resnum = line.strip().split(':')

                if chain not in self.res:
                    self.res[chain] = []
                self.res[chain].append(int(resnum))

    def get_data(self):
        return self.res
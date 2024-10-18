#!/bin/python3
# -----------------------------
# Written by Elton J. F. Chaves
# E-mail: chavesejf@gmail.com
# -----------------------------
class PDBParser:
    def __init__(self, pdbfile):
        self.pdbfile = pdbfile
        self.atoms = []
        
    def parse(self):
        try:
            with open(self.pdbfile, 'r') as pdb:
                for line in pdb:
                    if line.startswith('ATOM') or line.startswith('HETATM'):
                        atom_data = self._parse_atom_line(line)
                        self.atoms.append(atom_data)
        except FileNotFoundError:
            print(f"Error: File '{self.file_path}' not found.")

    def _parse_atom_line(self, line):   
        atom_data = {
            'record_name': line[0:6].strip(),
            'atom_index': int(line[6:11].strip()),
            'atom_name': line[12:16].strip(),
            'residue_model': line[16].strip(),
            'residue_name': line[17:20].strip(),
            'chain_id': line[21].strip(),
            'residue_number': int(line[22:26].strip()),
            'x': float(line[30:38].strip()),
            'y': float(line[38:46].strip()),
            'z': float(line[46:54].strip()),
            'occupancy': 0.00,
            'bfactor': 0.00,
            'atom_type': str(line[66:].strip())
        }
        if atom_data['atom_type'] == '':
            if atom_data['atom_name'][0] ==  'N':  atom_data['atom_type'] =  'N'
            if atom_data['atom_name'][0] ==  'C':  atom_data['atom_type'] =  'C'
            if atom_data['atom_name'][0] ==  'O':  atom_data['atom_type'] =  'O'
            if atom_data['atom_name'][0] ==  'S':  atom_data['atom_type'] =  'S'
            if atom_data['atom_name'][0] ==  'H':  atom_data['atom_type'] =  'H'
        return atom_data

    def get_atoms(self):
        return self.atoms

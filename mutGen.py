#!/bin/python3    

# Written by Elton J. F. Chaves
# -----------------------------
import numpy as np
import re, argparse, os, copy, glob, subprocess

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
            if atom_data['residue_name'] in ligands: atom_data['record_name'] = 'HETATM'
            if atom_data['residue_name'] in ions: atom_data['record_name'] = 'HETATM'
            if atom_data['atom_name'][0] == 'N': atom_data['atom_type'] = 'N'
            if atom_data['atom_name'][0] == 'C': atom_data['atom_type'] = 'C'
            if atom_data['atom_name'][0] == 'O': atom_data['atom_type'] = 'O'
            if atom_data['atom_name'][0] == 'S': atom_data['atom_type'] = 'S'
            if atom_data['atom_name'][0] == 'H': atom_data['atom_type'] = 'H'
        return atom_data

    def get_atoms(self):
        return self.atoms

class FileReader:
    def __init__(self, filename):
        self.filename = filename
        self.res = {}

    def read_file(self):
        with open(self.filename, 'r') as file:
            for line in file:
                key, value = line.strip().split(':')
                
                if key not in self.res:
                    self.res[key] = []
                self.res[key].append([int(value[:-1]),value[-1]])

    def get_data(self):
        return self.res

def pre_processing(pdb, partner1, partner2, output_dir):
    """
    :parameters:
    ------------
    pdb_file - 
    partner1 -
    partner2 -

    :returns:
    ---------
    """

    # ---
    print_infos(message='pre-processing protocol', type='info')    
    
    # ---
    pdb_parser = PDBParser(pdb)
    pdb_parser.parse()
    _atoms = pdb_parser.get_atoms()

    # informa o índice e nome da estrutura
    # ------------------------------------
    print_infos(message=f'{os.path.basename(pdb)}', type='structure')

    # cria diretório para outputs
    # ---------------------------
    output_name = os.path.basename(pdb[:-4]).lower()
    output_dir = f'{output_dir}/outputs/{output_name}'

    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    # verifica se as cadeias de partner1 e partner2 estão no pdb
    # ----------------------------------------------------------
    chains = find_chains(pdb)
    count = 0
    for chain in chains:
        if chain in partner1:
            count += 1
        if chain in partner2:
            count += 1
        if partner2 in ligands:
            count += 1
    
    if count <= 1:
        print_infos(message=f'provided chain_id(s) not found', type='warning')
        print_end()
    else:
        str_partner1 = ''.join(str(x) for x in partner1)
        str_partner2 = ''.join(str(x) for x in partner2)
        print_infos(message=f'partners = {str_partner1}_{str_partner2}', type='info')

    # extrai as cadeias do complexo de acordo com o conteúdo das flags 'partner1' e 'partner2'
    # ----------------------------------------------------------------------------------------
    _partner1 = []
    _partner2 = []
    _ions1 = []
    _ions2 = []
    for atom in _atoms:
        if atom['chain_id'] in partner1 and atom['record_name'] == 'ATOM':
            if atom['residue_model'] in ['A', '']:
                _partner1.append(atom)
    
        if partner2 in ligands:
            add_hydrogens_to_lig = True
            if atom['residue_name'] in partner2 and atom['record_name'] == 'HETATM':
                if atom['residue_model'] in ['A', '']:
                    _partner2.append(atom)
        else:
            if atom['chain_id'] in partner2 and atom['record_name'] == 'ATOM':
                if atom['residue_model'] in ['A', '']:
                    _partner2.append(atom)

    # busca por ions na estrutura
    # ---------------------------
    for atom in _atoms:
        if atom['chain_id'] in partner1 and atom['record_name'] == 'HETATM':
            if atom['residue_name'] in ions:
                _ions1.append(atom)
        if atom['chain_id'] in partner2 and atom['record_name'] == 'HETATM':
            if atom['residue_name'] in ions:
                _ions2.append(atom)
    print_infos(message=f'{len(_ions1) + len(_ions2)} ions found.', type='info')
    
    # adiciona hidrogênios no ligante
    # -------------------------------
    if add_hydrogens_to_lig is True:
        ligand  = f"{output_dir}/{output_name}_lig.pdb"
        ligandH = f"{output_dir}/{output_name}_lig_withH.pdb"
        write_pdb(_partner2, outfile=f'{ligand}')
        subprocess.run(f'obabel -ipdb {ligand} -opdb -O {ligandH} -p > /dev/null 2>&1', stdout=subprocess.PIPE, shell=True)
        pdb_parser = PDBParser(ligandH)
        pdb_parser.parse()
        _partner2 = pdb_parser.get_atoms()

    # ---
    partners = _partner1 + _ions1 + _partner2 + _ions2

    # verifica se existem gaps na estrutura de 'partner1' e 'partner2'
    # ----------------------------------------------------------------
    ngaps, igaps = find_gaps(partners)
    if ngaps > 0:
        print_infos(message=f'{ngaps} gap(s) found', type='warning')
        for gap_info in igaps:
            print_infos(message=f'{gap_info}', type='warning')
            continue
        print_end()
    else:
        print_infos(message=f'{ngaps} gap(s)', type='info')

    # se não existir parâmetros de mutação, ativa o reconhecimento automático de resíduos na 
    # interface e cria arquivo c/ os parâmetros de mutação
    # --------------------------------------------------------------------------------------
    if mutant_list is None:
        print_infos(message=f'enabling automatic recognition of interface residues', type='info')
        interface1 = find_atoms_closest_to_protein(_partner2, _partner1, dist_cutoff=int_dist_cutoff)
        
        # ---
        if partner2 in ligands:
            interface2 = find_atoms_closest_to_protein(_partner1, _partner2, dist_cutoff=int_dist_cutoff)

        # escreve arquivo com parâmetros de mutação
        mutant_file = f'{output_dir}/{output_name}.mutants.list'
        with open(mutant_file, 'w') as f:
            for chain in interface1.keys():
                for resnum in interface1[chain]:
                    f.write(f'{chain}:{resnum}A\n')
            if partner2 not in ligands:
                for chain in interface2.keys():
                    for resnum in interface2[chain]:
                        f.write(f'{chain}:{resnum}A\n')
    else:
        mutant_file = mutant_list
    
    # carrega arquivo com parâmetros de mutação
    # -----------------------------------------
    mutants_reader = FileReader(mutant_file)
    mutants_reader.read_file()
    mutants = mutants_reader.get_data()
    
    # escreve arquivo .pdb do complexo selvagem
    # -----------------------------------------
    output_dir_wt = f'{output_dir}/{output_name}_WT'
    if not os.path.isdir(output_dir_wt):
        os.makedirs(output_dir_wt)
    write_pdb(partners, outfile=f'{output_dir_wt}/{output_name}_wt.pdb')
    wt_resids = write_resids(partners, outfile=f'{output_dir}/{output_name}.resids.txt')
    wild_type = [f'{output_dir_wt}/{output_name}_wt.pdb']
    
    # insere mutações (edita arquivo .pdb)
    # ------------------------------------
    pdbs = []
    for chain in mutants.keys():
        for resnum in mutants[chain]:
            mutant = []
            _output_name = f'{output_name}_{chain}_{resnum[0]}{resnum[1]}'
            _output_dir  = f'{output_dir}/{_output_name}'
            wild_type_structure = copy.deepcopy(partners)
            if not os.path.isdir(_output_dir):
                os.makedirs(_output_dir)
            
            count = 0
            for atom in wild_type_structure:
                if atom['chain_id'] == chain and atom['residue_number'] == resnum[0]:
                    resname = aminoacids1(atom['residue_name'])
                    resmut  = aminoacids2(resnum[1])
                    outfile = f'{_output_dir}/{output_name}_{chain}_{resname}{resnum[0]}{resnum[1]}.pdb'
                    if atom['atom_name'] in backbone:
                        atom['residue_name'] = resmut
                        count += 1
                    else:
                        continue
                mutant.append(atom)
            if count == 0:
                print_infos(message='the mutant does not correspond to any residue', type='info')
                print_end()
            else:
                write_pdb(mutant, outfile)
                pdbs.append(outfile)

    return wild_type, pdbs, wt_resids

def post_processing(pdb_files, partner1, partner2, wt_resids):
    """
    :parameters:
    ------------
    pdb_files
    partner1
    partner2
    
    :returns: 
    ---------
    """
    str_partner1 = ''.join(str(x) for x in partner1)
    str_partner2 = ''.join(str(x) for x in partner2)

    # ---
    print_infos(message='post-processing protocol', type='info')
    
    pdbs = []
    for mol, pdb in enumerate(pdb_files):
        os.chdir(submit_dir)

        # informa nome da estrutura
        # -------------------------
        print_infos(message=f'{os.path.basename(pdb)}', type='structure')

        # define caminho e prefixo para escrever outputs
        # ----------------------------------------------
        output_name = os.path.basename(pdb[:-4])
        output_dir  = os.path.dirname(pdb)

        # insere mutação com o tLeap
        # --------------------------                
        prefix = f'{output_dir}/{output_name}'
        
        # cria diretório de outputs
        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)
        
        # renumera os resíduos da estrutura (primeiro resíduo = 1)
        # obs.: este passo é necessário para modelagem das ligações dissulfeto no tleap
        _pdb = renumber_residues(pdb, outfile=f'{prefix}_00.pdb')
        
        # atualiza a variável '_atoms' com a estrutura que contém os resíduos renumerados
        pdb_parser = PDBParser(_pdb)
        pdb_parser.parse()
        _atoms = pdb_parser.get_atoms()

        # ---
        coords = []
        for atom in _atoms:
            coords.append(atom)
        _pdb = write_pdb(coords, outfile=f'{prefix}_01.pdb')
                
        # realiza busca por ligações dissulfeto e cria um novo arquivo .pdb
        _pdb, disulfides = find_disulfides(_pdb, outfile=f'{prefix}_02.pdb')
        conf = f'{prefix}_leap.conf'
        conf = write_tleap_conf(_pdb, conf, disulfides)
        subprocess.run(f"tleap -f {conf} > {prefix}_leap.log", stdout=subprocess.PIPE, shell=True)
        
        # regenera a numeração dos resíduos usando a estrutura wild-type como referência
        # ------------------------------------------------------------------------------
        # cria dicionário com a numeração dos resíduos de cada cadeia
        _pdb = glob.glob(f'{output_dir}/*_02.pdb')[0]
        pdb_parser = PDBParser(_pdb)
        pdb_parser.parse()
        _atoms = pdb_parser.get_atoms()
        
        ref = {}
        with open(wt_resids, 'r') as f:
            atoms = f.readlines()
            chain_id = atoms[0].split()[1]
            first_res = int(atoms[0].split()[2])
            ref[chain_id] = []
            for atom in atoms:
                atom = atom.split()
                cur_chain = atom[1]
                resnum = int(atom[2])
                if cur_chain != chain_id:
                    ref[cur_chain] = []
                ref[cur_chain].append(resnum)
                chain_id = cur_chain

        # reescreve a numeração dos resíduos de acordo com a estrutura wild_type
        seen_indices = set()
        for cur_chain in ref.keys():
            endat = len(ref[cur_chain])
            count = -1
            for atom in _atoms:
                if atom['atom_index'] in seen_indices:
                    continue
                else:
                    if atom['atom_name'] == 'N' \
                    or atom['residue_name'] in ions \
                    or atom['residue_name'] in ligands:
                        count += 1
                    if count <= endat - 1:
                        resnum = ref[cur_chain][count]
                        atom['chain_id'] = cur_chain
                        atom['residue_number'] = resnum
                        seen_indices.add(atom['atom_index'])
                    else:
                        break
        
        _pdb = write_pdb(_atoms, outfile=f'{output_dir}/{output_name}_03.pdb')
        pdbs.append(_pdb)

        # mostra total de ligações dissulfeto na tela
        if _pdb.__contains__('wt'):
            print_infos(message=f'{len(disulfides)} disulfide bond(s) found', type='info')

    return pdbs

def aminoacids1(resname):
    aa = {
        'ALA': 'A',
        'CYS': 'C',
        'ASP': 'D',
        'GLU': 'E', 'GLH': 'Q',
        'PHE': 'F',
        'GLY': 'G',
        'HIS': 'H', 'HID': 'H', 'HIE': 'H', 'HIP': 'H',
        'ILE': 'I',
        'LYS': 'K',
        'LEU': 'L',
        'MET': 'M',
        'ASN': 'N', 'ASH': 'N',
        'PRO': 'P',
        'GLN': 'Q',
        'ARG': 'R',
        'SER': 'S',
        'THR': 'T',
        'VAL': 'V',
        'TYR': 'Y',
        'TRP': 'W'
    }
    onelettercode = aa[resname]
    return onelettercode

def aminoacids2(resname):
    aa = {
        'A': 'ALA',
        'C': 'CYS',
        'D': 'ASP',
        'E': 'GLU',
        'F': 'PHE',
        'G': 'GLY',
        'H': 'HIS',
        'I': 'ILE',
        'K': 'LYS',
        'L': 'LEU',
        'M': 'MET',
        'N': 'ASN',
        'P': 'PRO',
        'Q': 'GLN',
        'R': 'ARG',
        'S': 'SER',
        'T': 'THR',
        'V': 'VAL',
        'Y': 'TYR',
        'W': 'TRP'
    }
    threelettercode = aa[resname]
    return threelettercode

def write_tleap_conf(pdb, conf, disulfides):
    with open(conf, 'w') as f:
        f.write("source leaprc.protein.ff14SBonlysc\n")
        f.write("source leaprc.water.tip3p\n")
        f.write(f"protein = loadPdb {pdb}\n")
        if disulfides:
            for bond in disulfides:
                f.write(f'{bond}\n')
        f.write(f"savePdb protein {pdb[:-4]}.pdb\n")
        f.write(f"quit")
    return conf

def write_resids(partners, outfile):
    with open(outfile, 'w') as f:
        for atom in partners:
            if atom['record_name'] == 'ATOM'   and atom['atom_name'] == 'CA':
                resname = atom['residue_name']
                resnum  = atom['residue_number']
                chain   = atom['chain_id']
                f.write(f'{resname} {chain} {resnum}\n')
    return outfile

def write_pdb(partners, outfile):
    with open(outfile, 'w') as f:
        chain_id = partners[0]['chain_id']
        for atom in partners:
            cur_chain = atom['chain_id']
            if atom['residue_name'] in ions:
                f.write('TER\n')
            elif chain_id != cur_chain:
                f.write('TER\n')
                chain_id = cur_chain
            f.write(
                f"{atom['record_name']:<6s}"
                f"{atom['atom_index']:5d} "
                f"{atom['atom_name']:^4s} "
                f"{atom['residue_name']:<3s} "
                f"{atom['chain_id']:1s}"
                f"{atom['residue_number']:4d}    "
                f"{atom['x']:8.3f}{atom['y']:8.3f}{atom['z']:8.3f}"
                f"{atom['occupancy']:6.2f}{atom['bfactor']:6.2f}           "
                f"{atom['atom_type']}\n"
                )
    return outfile

def renumber_residues(pdb, outfile):
    pdb_parser = PDBParser(pdb)
    pdb_parser.parse()
    _atoms = pdb_parser.get_atoms()
    count = 0
    for atom in _atoms:
        if atom['atom_name'] == 'N':
            count += 1
        if atom['record_name'] == 'HETATM' and atom['residue_name'] in ions:
            count += 1
        atom['residue_number'] = count
    write_pdb(_atoms, outfile)
    return outfile

def find_disulfides(pdb, outfile):
    cys_sg_atoms = []
    cys_residues = []
    j_coord_list = []
    ssbonds      = []
    pdb_parser = PDBParser(pdb)
    pdb_parser.parse()
    _atoms = pdb_parser.get_atoms()
    for atom in _atoms:
        if atom['residue_name'] == 'CYS' and atom['atom_name'] == 'SG':
            cys_sg_atoms.append(atom)
    for i in cys_sg_atoms:
        for j in cys_sg_atoms:
            i_resnum = i['residue_number']
            j_resnum = j['residue_number']
            if i_resnum not in j_coord_list:
                d = calc_distance(i, j)
            else:
                continue
            if d > 0 and d < 3.5 and i_resnum:
                j_coord_list.append(j_resnum)
                cys_residues.append(i_resnum)
                cys_residues.append(j_resnum)
                ssbonds.append(f'bond protein.{i_resnum}.SG protein.{j_resnum}.SG')
    for atom in _atoms:
        if atom['residue_number'] in cys_residues:
            atom['residue_name'] = 'CYX'
    write_pdb(_atoms, outfile)
    return outfile, ssbonds

def find_chains(pdbfile):
    chains = set()
    with open(pdbfile, 'r') as file:
        for line in file:
            if line.startswith('ATOM'):
                chain_id = line[21]
                chains.add(chain_id)
    return chains

def find_ligands(pdbfile):
    ligands = set()
    with open(pdbfile, 'r') as file:
        for line in file:
            if line.startswith('HETATM'):
                hetatm = line[17:20]
                ligands.add(hetatm)
    return ligands

def find_atoms_closest_to_protein(atoms1, atoms2, dist_cutoff):
    near_residues = {}
    seen_indices = set()
    for i in atoms1:
        for j in atoms2:
            if i['record_name'] in ['ATOM', 'HETATM']:
                d = calc_distance(i, j)
                if d <= dist_cutoff and j['atom_index'] not in seen_indices:
                    chain  = j['chain_id']
                    resnum = j['residue_number']
                    if chain not in near_residues:
                        near_residues[chain] = [resnum]
                    elif resnum not in near_residues[chain]:
                        near_residues[chain].append(resnum)
                    seen_indices.add(j['atom_index'])
    return near_residues

def find_gaps(partners):
    ca_atoms = [atom for atom in partners if atom['atom_name'] == 'CA' and atom['residue_name'] not in ions]
    ca_total = len(ca_atoms)
    gap_info = []
    ngaps = 0
    for n in range(ca_total - 1):
        i = ca_atoms[n]
        j = ca_atoms[n + 1]
        ichain = i['chain_id']
        jchain = j['chain_id']
        d = calc_distance(i, j)
        if d >= gap_dist_cutoff:
            if ichain == jchain:
                gap = \
                'gap = ' + \
                i['chain_id'] + ':' + i['residue_name'] + str(i['residue_number']) + \
                ' <--> ' + \
                j['chain_id'] + ':' + j['residue_name'] + str(j['residue_number'])
                gap_info.append(gap)
                ngaps += 1
    return ngaps, gap_info

def calc_distance(atom1, atom2):
    i = np.array([atom1['x'], atom1['y'], atom1['z']])
    j = np.array([atom2['x'], atom2['y'], atom2['z']])
    return np.linalg.norm(i - j)

def isdir(path):
    if os.path.isdir(path):
        return os.path.abspath(path)
    else:
        print_infos(message=f'path not found -> {path}', type='info')
        print_end()

def isfile(file):
    if os.path.isfile(file):
        return os.path.abspath(file)
    else:
        print_infos()

def print_end():
    exit(f'\n --- End process --- \n')

def print_infos(message, type):
    if type == 'info':
        print(f'      info: {message}')
    if type == 'warning':
        print(f'   warning: {message}')
    if type == 'structure':
        print(f' structure: {message}')
    if type == 'protocol':
        print(f'  protocol: {message}')

def header():
    print('')
    print(' =============================================')
    print('           MutGen - Mutant Generator          ')
    print('')   
    print('  Author: Chaves, EJF')
    print(' Contact: chavesejf@gmail.com')
    print(' =============================================')
    print('')

if (__name__ == "__main__"):   
    # define variável que armazena o diretório de submissão
    submit_dir = os.getcwd()

    # imprime o cabeçalho na tela
    header()

    # configura os argumentos do script
    # ---------------------------------
    parser = argparse.ArgumentParser()
    mandatory = parser.add_argument_group('mandatory arguments')
 
    # obrigatórios
    mandatory.add_argument('--ipdb', nargs=1, type=isfile, required=True, metavar='',
    help='str | input file in the PDB format')
    mandatory.add_argument('--partner1', nargs=1, required=True, metavar='',
    help='str | chain ID (e.g.: receptor)')
    mandatory.add_argument('--partner2', nargs=1, required=True, metavar='',
    help='str | chain ID (e.g.: ligand)')

    # opcionais
    parser.add_argument('--odir', nargs=1, type=isdir, default=[f'{submit_dir}'], metavar='',
    help=f'str | folder path to save the output files (default={submit_dir})')
    parser.add_argument('--mut_list', nargs=1, type=isfile, default=[None], metavar='',
    help='str | ...')
    parser.add_argument('--int_dist_cutoff', nargs=1, type=float, default=[4.5], metavar='',
    help=f'float | ')
    
    # finaliza configuração dos argumentos e variáveis
    # ------------------------------------------------
    args             = parser.parse_args()
    pdb_file         = args.ipdb[0]
    ligands          = find_ligands(pdb_file)
    partner1         = list(args.partner1)[0]
    partner2         = list(args.partner2)[0] if args.partner2[0] not in ligands else args.partner2[0]
    mutant_list      = args.mut_list[0]
    output_dir       = args.odir[0]
    int_dist_cutoff  = args.int_dist_cutoff[0]
    gap_dist_cutoff  = 16
    backbone         = ['N', 'CA', 'C', 'O']
    ions             = ['MG','CA','NA','CL','FE','K','ZN','MN']
    separate_chains  = True

    # pré-processamento
    # =================
    # a. 
    # b. 
    # -----------------
    wild_type, pdbs, wt_resids = pre_processing(pdb_file, partner1, partner2, output_dir)
    
    # pós-processamento
    # =================
    if len(wild_type) == 1 and len(pdbs) > 0:
        wild_type = post_processing(wild_type, partner1, partner2, wt_resids)
        pdbs = post_processing(pdbs, partner1, partner2, wt_resids)
    else:
        print_infos(message='nothing to do', type='info')
        print_end()
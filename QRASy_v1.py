#!/bin/python3

# -----------------------------
# Written by Elton J. F. Chaves
# E-mail: chavesejf@gmail.com
# -----------------------------
from modules.pdb_parser import *
from modules.mut_parser import *
from shutil import which
import numpy as np
import re
import os
import copy
import time
import glob
import argparse
import subprocess

def pre_processing(pdb, partner1, partner2, output_name, output_dir):
    """
    """
    # ---
    print_infos(message='pre-processing protocol', type='info')    
    
    # ---
    pdb_parser = PDBParser(pdb)
    pdb_parser.parse()
    atoms = pdb_parser.get_atoms()

    # informa o índice e nome da estrutura na tela
    # --------------------------------------------
    print_infos(message=f'{os.path.basename(pdb)}', type='structure')

    # cria diretório para outputs
    # ---------------------------
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    # verifica se as cadeias de partner1 e partner2 existem no arquivo .pdb
    # ---------------------------------------------------------------------
    chains = find_chains(pdb)
    count  = 0
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
    add_hydrogens_to_lig = False
    _partner1 = []
    _partner2 = []
    _ions1 = []
    _ions2 = []
    for atom in atoms:
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
    for atom in atoms:
        if atom['chain_id'] in partner1 and atom['record_name'] == 'HETATM':
            if atom['residue_name'] in ions:
                _ions1.append(atom)
        if partner2 in ligands:
            pass
        else:
            if atom['chain_id'] in partner2 and atom['record_name'] == 'HETATM':
                if atom['residue_name'] in ions:
                    _ions2.append(atom)

    print_infos(message=f'{len(_ions1) + len(_ions2)} ion(s) found', type='info')

    # escreve arquivo .pdb do complexo selvagem
    # -----------------------------------------
    partners  = _partner1 +  _partner2 + _ions1 + _ions2
    wt_resids = write_resids(partners, outfile=f'{output_dir}/{output_name}.resids.txt')
    
    if not os.path.isdir(f'{output_dir}/{output_name}_WT'):
        os.makedirs(f'{output_dir}/{output_name}_WT')
    
    if partner2 in ligands:
        partners = _partner1 + _ions1 + _ions2
    
    wt_struct = write_pdb(partners, outfile=f'{output_dir}/{output_name}_WT/{output_name}_wt.pdb')
    
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
        print_infos(message=f'{ngaps} gap(s) found', type='info')
    
    # se não existir parâmetros de mutação, ativa o reconhecimento automático de resíduos na 
    # interface e cria arquivo c/ os parâmetros de mutação
    # --------------------------------------------------------------------------------------
    if mutant_list is None:
        print_infos(message=f'mutant_list is False', type='info')
        print_infos(message=f'enabling automatic recognition of interface residues', type='info')
        interface1 = find_atoms_closest_to_protein(_partner2, _partner1, dist_cutoff=int_dist_cutoff)
        interface2 = find_atoms_closest_to_protein(_partner1, _partner2, dist_cutoff=int_dist_cutoff)

        # escreve arquivo com parâmetros de mutação
        mutant_file = f'{output_dir}/{output_name}.mutants.list'
        with open(mutant_file, 'w') as f:
            for chain in interface1.keys():
                for resnum in interface1[chain]:
                    f.write(f'{chain}:{resnum}\n')
            if partner2 not in ligands:
                for chain in interface2.keys():
                    for resnum in interface2[chain]:
                        f.write(f'{chain}:{resnum}\n')
    else:
        mutant_file = mutant_list
    
    # carrega arquivo com parâmetros de mutação
    # -----------------------------------------
    mutants_reader = FileReader(mutant_file)
    mutants_reader.read_file()
    mutants = mutants_reader.get_data()
    
    # insere mutações (edita arquivo .pdb)
    # ------------------------------------
    pdbs = []
    for chain in mutants.keys():
        for resnum in mutants[chain]:
            mutant = []
            _output_name = f'{output_name}_{chain}_{resnum}A'
            _output_dir  = f'{output_dir}/{_output_name}'
            wild_type_structure = copy.deepcopy(partners)
            
            if not os.path.isdir(_output_dir):
                os.makedirs(_output_dir)
            
            count = 0
            for atom in wild_type_structure:
                if atom['chain_id'] == chain and atom['residue_number'] == resnum:
                    resname = aminoacids_1lettercode(atom['residue_name'])
                    resmut  = aminoacids_3lettercode(resname)
                    outfile = f'{_output_dir}/{output_name}_{chain}_{resname}{resnum}A.pdb'
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
    
    # adiciona hidrogênios no ligante
    # -------------------------------
    if add_hydrogens_to_lig:
        ligH = openbabel(_partner2, lig_name=f"{output_dir}/{output_name}_lig")
        ligH_bcc, ligH_frcmod = antechamber(ligH, lig_net_charge)
    else:
        ligH_bcc    = None
        ligH_frcmod = None
    return wt_struct, wt_resids, pdbs, ligH_bcc, ligH_frcmod

def post_processing(pdb_files, partner1, partner2, wt_resids, ligH_bcc, ligH_frcmod):
    """
    """
    str_partner1 = ''.join(str(x) for x in partner1)
    str_partner2 = ''.join(str(x) for x in partner2)
    
    # ---
    print_infos(message='post-processing protocol', type='info')
    
    pdbs = []
    for mol, pdb in enumerate(pdb_files):
        os.chdir(submit_dir)

        # informa nome da estrutura
        print_infos(message=f'{os.path.basename(pdb)}', type='structure')

        # define caminho e prefixo para escrever outputs
        output_name = os.path.basename(pdb[:-4])
        output_dir  = os.path.dirname(pdb)
        prefix      = f'{output_dir}/{output_name}'
        
        # cria diretório de outputs
        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)
        
        # renumera resíduos da estrutura a partir do 1
        _pdb = renumber_residues(pdb, outfile=f'{prefix}_00.pdb')
        
        # atualiza a variável '_atoms' com a estrutura que contém os resíduos renumerados e escreve um novo arquivo .pdb
        pdb_parser = PDBParser(_pdb)
        pdb_parser.parse()
        _atoms = pdb_parser.get_atoms()
        coords = []
        for atom in _atoms:
            coords.append(atom)
        _pdb = write_pdb(coords, outfile=f'{prefix}_01.pdb')
            
        # modelagem de ligações dissulfeto
        _pdb, disulfides = find_disulfides(_pdb, outfile=f'{prefix}_02.pdb')

        # mostra total de ligações dissulfeto na tela
        if _pdb.__contains__('wt'):
            print_infos(message=f'{len(disulfides)} disulfide bond(s) found', type='info')
    
        # prepara input para o tleap
        with open(f'{prefix}.tleap1.in', 'w') as f:
            f.write("source leaprc.protein.ff14SBonlysc\n")
            f.write(f"protein = loadPdb {pdb}\n")
            if disulfides:
                for bond in disulfides:
                    f.write(f'{bond}\n')
            f.write(f"savePdb protein {pdb[:-4]}.pdb\n")
            f.write(f"quit")
        
        # executa o tleap
        subprocess.run(f"tleap -f {prefix}.tleap1.in > /dev/null 2>&1", stdout=subprocess.PIPE, shell=True)
        
        # otimização de estrutura com solvente implícito
        # ----------------------------------------------
        with open(f'{prefix}.tleap2.in', 'w') as f:
            if system_type == 'E:L':
                f.write( 'source leaprc.protein.ff14SBonlysc\n')
                f.write( 'source leaprc.water.tip3p\n')
                f.write( 'source leaprc.gaff2\n')
                f.write(f'REC = loadpdb {prefix}_02.pdb\n')
                f.write(f'LIG = loadmol2 {ligH_bcc}\n')
                f.write(f'loadamberparams {ligH_frcmod}\n')
                f.write(f'saveoff LIG data.lib\n')
                f.write(f'loadoff data.lib\n')
                f.write( 'COM = combine {REC LIG}\n')
                f.write( 'charge COM\n')
                if ionize is True:
                    f.write( 'addIons COM Na+ 0\n')
                    f.write( 'addIons COM Cl- 0\n')
                f.write(f'savepdb COM {prefix}_03.pdb\n')
                f.write(f'saveamberparm COM {prefix}_03.prmtop {prefix}_03.inpcrd\n')
                f.write( 'quit')
            if system_type == 'P:P':
                f.write( 'source leaprc.protein.ff14SBonlysc\n')
                f.write( 'source leaprc.water.tip3p\n')
                f.write( 'source leaprc.gaff2\n')
                f.write(f'REC = loadpdb {prefix}_02.pdb\n')
                f.write( 'charge REC\n')
                if ionize is True:
                    f.write( 'addIons REC Na+ 0\n')
                    f.write( 'addIons REC Cl- 0\n')
                f.write(f'savepdb REC {prefix}_03.pdb\n')
                f.write(f'saveamberparm REC {prefix}_03.prmtop {prefix}_03.inpcrd\n')
                f.write( 'quit')

        # executa o tleap
        subprocess.run(f'tleap -f {prefix}.tleap2.in', stdout=subprocess.PIPE, shell=True)
        
        # determina a carga total do complexo
        with open(f'{submit_dir}/leap.log', 'r') as f:
            lines = f.readlines()
            for line in lines:
                if line.__contains__('Total unperturbed charge:'):
                    charge = float(re.findall(r'-?\d+.\d+', line)[0])
                    break

        # prepara arquivo .conf para otimização de geometria no NAMD
        with open(f'{prefix}.namd-min.conf', 'w') as f:
            f.write( 'amber                  yes\n')
            f.write(f'parmfile               {prefix}_03.prmtop\n')
            f.write(f'ambercoor              {prefix}_03.inpcrd\n')
            f.write( 'switching              off\n')
            f.write( 'exclude                scaled1-4\n')
            f.write( '1-4scaling             0.833333\n')
            f.write( 'cutoff                 12\n')
            f.write( 'pairlistdist           14\n')
            f.write( 'timestep               2\n')
            f.write( 'rigidBonds             all\n')
            f.write( 'rigidTolerance         1e-08\n')
            f.write( 'nonbondedFreq          1\n')
            f.write( 'fullElectFrequency     2\n')
            f.write( 'stepspercycle          20\n')
            f.write(f'outputName             {os.path.basename(prefix)}_04\n')
            f.write(f'DCDFreq                250\n')
            f.write( 'xstFreq                100\n')
            f.write( 'outputEnergies         100\n')
            f.write( 'outputPressure         100\n')
            f.write( 'temperature            310\n')
            f.write( 'gbis                   on\n')
            f.write( 'alphacutoff            12.0\n')
            f.write( 'ionConcentration       0.15\n')
            f.write( 'minimize               10000\n')
        
        # otimização de geometria com o NAMD
        os.chdir(output_dir)
        
        # atualiza
        print_infos(message='running geometry optimization', type='info')
        
        if not os.path.isfile(f'{prefix}_04.dcd'):
            subprocess.run(f'namd3 {prefix}.namd-min.conf > /dev/null 2>&1', stdout=subprocess.PIPE, shell=True)        
        
        # extrai o último frame da trajetória de otimização de geometria
        with open(f'{prefix}.cpptraj.in', 'w') as f:
            f.write(f'trajout {prefix}_04_min_lastframe.pdb pdb onlyframes 40\nquit')
        
        # executa o cpptraj
        subprocess.run(f'cpptraj -p {prefix}_03.prmtop -y {prefix}_04.dcd -i {prefix}.cpptraj.in', stdout=subprocess.PIPE, shell=True)

        # ----------------------------------------------------------------------
        # reescreve a numeração dos resíduos de acordo com a estrutura wild_type
        #
        # GAMBIARRA ==> melhorar este bloco
        # ----------------------------------------------------------------------
        os.chdir(submit_dir)
        ref = {}
        with open(wt_resids, 'r') as f:
            atoms         = f.readlines()
            chain_id      = atoms[0].split()[1]
            seen_chains   = set(chain_id)
            ref[chain_id] = []
            for atom in atoms:
                atom      = atom.split()
                residue   = str(atom[0])
                cur_chain = str(atom[1])
                resnum    = int(atom[2])
                if cur_chain != chain_id:
                    if cur_chain not in seen_chains:
                        ref[cur_chain] = []
                ref[cur_chain].append([residue, resnum])
                chain_id = cur_chain
                seen_chains.add(chain_id)
        
        # ---
        pdb_parser = PDBParser(f'{prefix}_04_min_lastframe.pdb')
        pdb_parser.parse()
        atoms = pdb_parser.get_atoms()
        
        # ---
        seen_indices = set()
        for cur_chain in ref.keys():
            endat = len(ref[cur_chain])
            count = -1
            for atom in atoms:
                if atom['atom_index'] in seen_indices:
                    continue
                else:
                    if atom['atom_name'] == 'N':
                        count += 1
                    if atom['atom_type'] in ions:
                        atom['record_name'] = 'HETATM'
                        atom['chain_id'] = cur_chain
                        atom['residue_number'] = resnum + 1
                        resnum += 1
                        continue
                    if atom['residue_name'] == 'LIG':
                        count  += 1
                        atom['record_name'] = 'HETATM'
                        atom['chain_id'] = cur_chain
                        atom['residue_number'] = resnum + 1
                        continue
                    
                    if count <= endat - 1:
                        resnum = ref[cur_chain][count][1]
                        atom['chain_id'] = cur_chain
                        atom['residue_number'] = resnum
                        seen_indices.add(atom['atom_index'])
                    else:
                        break
        
        # ---
        _pdb = write_pdb(atoms, outfile=f'{output_dir}/{output_name}_05.pdb')
        pdbs.append(_pdb)
        pdb_parser = PDBParser(_pdb)        
        pdb_parser.parse()
        atoms = pdb_parser.get_atoms()

        # separa estruturas em "complex", "partner1" e "partner2"
        mylist = ['complex', 'partner1', 'partner2']
        com    = []
        p1     = []
        p2     = []
        for item in mylist:
            for atom in atoms:
                if item == 'complex':
                    if atom['chain_id'] in partner1 or atom['chain_id'] in partner2:
                        com.append(atom)
                if system_type == 'E:L':
                    if item == 'partner1':
                        if atom['chain_id'] in partner1 and not atom['residue_name'] == 'LIG':
                            p1.append(atom)
                    if item == 'partner2':
                        if atom['residue_name'] == 'LIG':
                            p2.append(atom)
                if system_type == 'P:P':
                    if item == 'partner1':
                        if atom['chain_id'] in partner1:
                            p1.append(atom)
                    if item == 'partner2':
                        if atom['chain_id'] in partner2:
                            p2.append(atom)

        # ---
        write_pdb(com, outfile=f'{prefix}_06.complex.pdb')
        write_pdb(p1,  outfile=f'{prefix}_06.partner1.pdb')
        write_pdb(p2,  outfile=f'{prefix}_06.partner2.pdb')

        # ---
        mopac_keywords = 'XYZ PDB PM7 GEO-OK ALLVEC LARGE VECTORS 1SCF MOZYME PL T=1D EPS=78.4 RSOLV=1.3 CUTOFF=9.0 LET DISP(1.0)'
        write_mop(com, mopac_keywords, outfile=f'{prefix}_07.complex.mop')
        write_mop(p1,  mopac_keywords, outfile=f'{prefix}_07.partner1.mop')
        write_mop(p2,  mopac_keywords, outfile=f'{prefix}_07.partner2.mop')

        # executa o mopac
        os.chdir(output_dir)
        
        print_infos(message='running QM calculation (complex)', type='info')
        runMOPAC(f'{prefix}_07.complex.mop', mopac_exec)

        print_infos(message='running QM calculation (complex)', type='info')
        runMOPAC(f'{prefix}_07.partner1.mop', mopac_exec)
        
        print_infos(message='running QM calculation (receptor)', type='info')
        runMOPAC(f'{prefix}_07.partner2.mop', mopac_exec)

        # ---    
        for file in files_to_remove:
            if os.path.isfile(file):
                os.remove(file)
    return pdbs

def runMOPAC(mop, mopac_exec):
    """
    """
    subprocess.run(f'{mopac_exec} {mop} > /dev/null 2>&1', shell=True, stdout=subprocess.PIPE)

def openbabel(partner2, lig_name):
    """
    """
    print_infos(message='add_hydrogens_to_lig is True', type='info')
    write_pdb(partner2, outfile=f'{lig_name}.pdb')
    subprocess.run(f'obabel -ipdb {lig_name}.pdb -omol2 -O {lig_name}_H.mol2 -p > /dev/null 2>&1', shell=True, stdout=subprocess.PIPE)
    return f'{lig_name}_H.mol2'

def antechamber(lig_file, lig_net_charge):
    """
    """
    ligH_bcc    = f'{lig_file[:-4]}_H.bcc.mol2'
    ligH_frcmod = f'{lig_file[:-4]}_H.bcc.frcmod'
    if not os.path.isfile(ligH_bcc):
        print_infos(message=f'calculating atom charges for the ligand', type='info')
        subprocess.run(f'antechamber \
        -fi mol2 \
        -fo mol2 \
        -i {lig_file} \
        -o {ligH_bcc} \
        -nc {lig_net_charge} \
        -at gaff2 \
        -c bcc \
        -rn LIG \
        -dr no \
        -pf yes > /dev/null 2>&1', stdout=subprocess.PIPE, shell=True)
    if os.path.isfile(ligH_bcc):
        subprocess.run(f'parmchk2 \
        -i {ligH_bcc} \
        -o {ligH_frcmod} \
        -f mol2 \
        -s gaff2 > /dev/null 2>&1', stdout=subprocess.PIPE, shell=True)
    else:
        print_infos(message='an error occurred while preparing the ligand', type='info')
        print_infos(message='try using the "--lig_net_charge" parameter', type='info')
        print_end()
    return ligH_bcc, ligH_frcmod

def aminoacids_1lettercode(resname):
    """
    """
    aa = {
        'ALA': 'A',
        'CYS': 'C',
        'ASP': 'D',
        'GLU': 'E', 'GLH': 'E',
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

def aminoacids_3lettercode(resname):
    """
    """
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

def write_resids(partners, outfile):
    """
    """
    lig_count = 0
    with open(outfile, 'w') as f:
        for atom in partners:
            resname = atom['residue_name']
            resnum  = atom['residue_number']
            chain   = atom['chain_id']
            condition = False
            if atom['record_name'] == 'ATOM' and atom['atom_name'] == 'CA':
                condition = True
            elif atom['residue_name'] in ligands:
                lig_count += 1
            elif atom['atom_type'] in ions:
                condition = False
            
            # ---
            if condition is True:
                f.write(f'{resname} {chain} {resnum}\n')
            if lig_count == 1:
                f.write(f'{resname} {chain} {resnum}\n')
    return outfile

def write_mop(partners, mopac_keywords, outfile):
    """
    """
    with open(outfile, 'w') as f:
        chain_id = partners[0]['chain_id']
        f.write(f"{mopac_keywords}\n\n\n")
        for atom in partners:
            cur_chain = atom['chain_id']
            if atom['residue_name'] in ions:
                f.write('TER\n')
            elif chain_id != cur_chain:
                f.write('TER\n')
                chain_id = cur_chain           
            if len(atom['atom_name']) == 3:
                atom_name = f"{atom['atom_name']:>4s}"
            else:
                atom_name = f"{atom['atom_name']:^4s}"
            f.write(
                f"{atom['record_name']:<6s}"
                f"{atom['atom_index']:>5d} "
                f"{atom_name} "
                f"{atom['residue_name']:<3s} "
                f"{atom['chain_id']:1s}"
                f"{atom['residue_number']:>4d}    "
                f"{atom['x']:8.3f}{atom['y']:8.3f}{atom['z']:8.3f}"
                f"{atom['occupancy']:6.2f}{atom['bfactor']:6.2f}           "
                f"{atom['atom_type']:>2s}\n")
    return outfile

def write_pdb(partners, outfile):
    """
    """
    with open(outfile, 'w') as f:
        chain_id = partners[0]['chain_id']
        for atom in partners:
            cur_chain = atom['chain_id']
            if chain_id != cur_chain:
                f.write('TER\n')
                chain_id = cur_chain
            f.write(
                f"{atom['record_name']:<6s}"
                f"{atom['atom_index']:>5d} "
                f"{atom['atom_name']:^4s} "
                f"{atom['residue_name']:<3s} "
                f"{atom['chain_id']:1s}"
                f"{atom['residue_number']:4d}    "
                f"{atom['x']:8.3f}{atom['y']:8.3f}{atom['z']:8.3f}"
                f"{atom['occupancy']:6.2f}{atom['bfactor']:6.2f}           "
                f"{atom['atom_type']:>2s}\n")
            if atom['residue_name'] in ions:
                f.write('TER\n')
    return outfile

def renumber_residues(pdb, outfile):
    """
    """
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
    """
    """
    cys_sg_atoms = []
    cys_residues = []
    j_coord_list = []
    ssbonds = []
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
    """
    """
    chains = set()
    with open(pdbfile, 'r') as file:
        for line in file:
            if line.startswith('ATOM'):
                chain_id = line[21]
                chains.add(chain_id)
    return chains

def find_ligands(pdbfile):
    """
    """
    ligands = set()
    with open(pdbfile, 'r') as file:
        for line in file:
            if line.startswith('HETATM'):
                hetatm = line[17:20]
                ligands.add(hetatm)
    return ligands

def find_atoms_closest_to_protein(atoms1, atoms2, dist_cutoff):
    """
    """
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
    """
    """
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
    """
    """
    i = np.array([atom1['x'], atom1['y'], atom1['z']])
    j = np.array([atom2['x'], atom2['y'], atom2['z']])
    return np.linalg.norm(i - j)

def istool(tool):
    """
    """
    return which(tool) is not None

def isdir(path):
    """
    """
    if os.path.isdir(path):
        return os.path.abspath(path)
    else:
        print_infos(message=f'path not found -> {path}', type='info')
        print_end()

def isfile(file):
    """
    """
    if os.path.isfile(file):
        return os.path.abspath(file)
    else:
        print_infos()

def print_end():
    """
    """
    exit(f'\n --- End process --- \n')

def print_infos(message, type):
    """
    """
    if type == 'info':
        print(f'      info: {message}')
    if type == 'warning':
        print(f'   warning: {message}')
    if type == 'structure':
        print(f' structure: {message}')
    if type == 'protocol':
        print(f'  protocol: {message}')

def processing_time(st):
    """
    """
    sum_x = 0
    for i in range(1000000):
        sum_x += i
    time.sleep(1)
    elapsed_time = time.time() - st
    print(' processing time:', time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))

def configure_requirements():
    """
    """
    mopac_exec   = '/opt/mopac/mopac-main/build/mopac'
    requirements = ['namd3', 'antechamber', 'cpptraj', 'parmchk2', 'obabel']
    for item in requirements:
        condition = istool(item)
        if condition is False:
            print(f' error: requirement not found -> {item}\n')
            print_end()
    if not os.path.isfile(mopac_exec):
        print(f' error: requirement not found -> {mopac_exec}\n')
        print_end()    
    return mopac_exec

def header():
    """
    """
    print('')
    print(' ================================================')
    print('  QRASy - Quantum Reactivity Alanine Scanning (y)')
    print()
    print(' Authors: Elton Chaves, Gerd Rocha, Igor Barden  ')
    print(' Contact: chavesejf@gmail.com')
    print('     DOI: N/A')
    print(' ================================================')
    print('')

if (__name__ == "__main__"):
    """
    """
    # define variável que armazena o diretório de submissão
    submit_dir = os.getcwd()

    # Define o tempo de início do script
    st = time.time()

    # imprime o cabeçalho na tela
    header()

    # ---
    mopac_exec = configure_requirements()

    # configura os argumentos do script
    # ---------------------------------
    parser    = argparse.ArgumentParser()
    mandatory = parser.add_argument_group('mandatory arguments')

    # obrigatórios
    mandatory.add_argument('--ipdb', nargs=1, type=isfile, required=True, metavar='',
    help='input file in the PDB format')
    mandatory.add_argument('--partner1', nargs=1, required=True, metavar='',
    help='chain ID of the partner1 (e.g.: receptor)')
    mandatory.add_argument('--partner2', nargs=1, required=True, metavar='',
    help='chain ID of the partner2 or ligand ID (e.g.: ligand)')
    mandatory.add_argument('--system_type', nargs=1, required=True, choices=['E:L', 'P:P'], default=['P:P'], metavar='',
    help='...')

    # opcionais
    parser.add_argument('--odir', nargs=1, type=isdir, default=[f'{submit_dir}'], metavar='',
    help=f'folder path to save the output files (default={submit_dir})')
    parser.add_argument('--mut_list', nargs=1, type=isfile, default=[None], metavar='',
    help='...')
    parser.add_argument('--int_dist_cutoff', nargs=1, type=float, default=[4.5], metavar='',
    help=f'...')
    parser.add_argument('--lig_net_charge', nargs=1, type=int, default=[0], metavar='',
    help=f'...')
    parser.add_argument('--ionize', action='store_true',
    help=f'...')
    
    # finaliza configuração dos argumentos e variáveis
    # ------------------------------------------------
    args                 = parser.parse_args()
    pdb_file             = args.ipdb[0]
    ligands              = find_ligands(pdb_file)
    partner1             = list(args.partner1)[0]
    partner2             = list(args.partner2)[0] if args.partner2[0] not in ligands else args.partner2[0]
    system_type          = args.system_type[0]
    mutant_list          = args.mut_list[0]
    int_dist_cutoff      = args.int_dist_cutoff[0]
    lig_net_charge       = args.lig_net_charge[0]
    ionize               = args.ionize
    gap_dist_cutoff      = 16
    backbone             = ['N', 'CA', 'C', 'O']
    ions                 = ['MG', 'CA', 'NA', 'CL', 'Cl-', 'FE', 'K', 'ZN', 'MN', 'BR', 'SE', 'AL', 'SB', 'BA', 'SI', 'BI', 'CU', 'AG', 'CO']
    output_name          = os.path.basename(pdb_file[:-4]).lower()
    output_dir           = f'{args.odir[0]}/outputs_qrasy/{output_name}'
    files_to_remove      = [
        f'{submit_dir}/data.lib',
        f'{submit_dir}/leap.log',   
        f'{submit_dir}/sqm.in',
        f'{submit_dir}/sqm.out',
        f'{submit_dir}/sqm.pdb']
    
    # =================
    # pré-processamento
    # =================
    wt, wt_resids, pdbs, ligH_bcc, ligH_frcmod = pre_processing(pdb_file, partner1, partner2, output_name, output_dir)
    
    # =================
    # pós-processamento
    # =================
    if os.path.isfile(wt) and len(pdbs) > 0:
        wt   = post_processing([wt], partner1, partner2, wt_resids, ligH_bcc, ligH_frcmod)
        pdbs = post_processing(pdbs, partner1, partner2, wt_resids, ligH_bcc, ligH_frcmod)
    else:
        print_infos(message='nothing to do', type='info')
    
    processing_time(st); print_end()

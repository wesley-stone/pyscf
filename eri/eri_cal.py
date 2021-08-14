import os
import sys
import re
import pickle
from types import FunctionType
from utils import create_if_not_exist
import numpy as np
import pandas as pd
from copy import deepcopy
from pyscf import gto, dft, lib
import random

ao_order = {
    1: np.array([
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1]
    ]),
    2: np.array([
        [2, 0, 0],
        [1, 1, 0],
        [1, 0, 1],
        [0, 2, 0],
        [0, 1, 1],
        [0, 0, 2]
    ])
}

def compose_cao(angl, env, atom_pos):
    #Note: follow the order from mol.ao_labels()
    assert(angl >= 0 and angl <=2)
    
    if angl == 0:
        rlts = [{
            'center': atom_pos,
            'coeff': env,
            'angl': np.zeros(3)
        }]
    else:
        rlts = [{
            'center': atom_pos,
            'coeff': env,
            'angl': a,
        } for a in ao_order[angl]]
    return rlts

def get_caos(atm_xyz, c_s=None, kernel=False, log_path=None):
    if log_path is None:
        mol = gto.Mole()
    else:
        mol = gto.Mole(log_path=log_path)
    mol.verbose=9
    mol.atom = atm_xyz
    if c_s is not None:
        mol.charge = int(c_s.split(' ')[0])
#    mol.spin = int(c_s.split(' ')[1])-1
    mol.cart=True
    mol.basis = '6-31G(d)'
    mol.build()

    binfo = deepcopy(mol._basis)
    bas_env = {k: gto.mole.make_bas_env(v) for k, v in binfo.items()}
    ainfo = deepcopy(mol._atom)

    # compose the configuration of each contracted ao
    contracted_ao = []
    for atm in ainfo:
        a_sym = atm[0]
        a_pos = atm[1:]
        _bas, _env = bas_env[a_sym]
        
        for b in _bas:
            angl = b[gto.mole.ANG_OF]
            coeff = _env.reshape(2, -1).T
            tmp_cao = compose_cao(angl, coeff, a_pos)
            contracted_ao.extend(tmp_cao)
    
    # check if the generated ao is consistent with pyscf
    assert(len(contracted_ao) == mol.nao)

    scf=dft.RKS(mol)
    scf.xc='M06-2X'
    if kernel:
        scf.kernel()
    
    return {
        'caos': contracted_ao,
        'mol': mol, 
        'scf': scf
    }

def random_point_in_box(range_xyz):
    xmin, xmax, ymin, ymax, zmin, zmax = range_xyz
    xyz = np.zeros(3)
    xyz[0] = random.uniform(xmin, xmax)
    xyz[1] = random.uniform(ymin, ymax)
    xyz[2] = random.uniform(zmin, zmax)
    return xyz

def L2distance(p1, p2):
    return np.sum((p1-p2)**2)

def generate_config(range_xyz, cnt, h2o_xyz):
    rlt = np.empty((cnt, 2, 3))
    pO = np.zeros(3)
    # min length from h to o
    h_to_o = L2distance(h2o_xyz[0], pO)
    # min length from h to h
    h_to_h = L2distance(h2o_xyz[0], h2o_xyz[1])
    for ci in range(cnt):
        p1 = random_point_in_box(range_xyz)
        while L2distance(p1, pO) < h_to_o:
            p1 = random_point_in_box(range_xyz)

        p2 = random_point_in_box(range_xyz)
        while L2distance(p2, p1) < h_to_h or L2distance(p2, pO) < h_to_o:
            p2 = random_point_in_box(range_xyz)
        rlt[ci, 0] = p1
        rlt[ci, 1] = p2
    return rlt

def generate_h2o_data(cfg, logger):
    cnt = cfg.data.cnt
    # atm position for h2o in a.u.
    h2o_xyz = np.array([
        [0,        0,  0.11779],
        [0,  0.75545, -0.47116],
        [0, -0.75545, -0.47116],
    ])
    # constrain H atom in the unit box.
    range_xyz = np.array([
        -1, 1, -1, 1, -1, 1
    ])
    atm_xyzs = generate_config(range_xyz, cnt, h2o_xyz)

    xyz_str = 'O 0 0 0;'
    atm_caos = []
    log_folder = 'scf/'
    create_if_not_exist(log_folder)
    logger.info('generating contracted ao...')
    for i, h2_xyz in enumerate(atm_xyzs):
        for h_xyz in h2_xyz:
            xyz_str = xyz_str + f'H {h_xyz[0]} {h_xyz[1]} {h_xyz[2]};'
        
        caos = get_caos(xyz_str, kernel=True, log_path=os.path.join(log_folder, f'{i}.dt'))['caos']
        atm_caos.append(caos)

    fns = [f for f in os.listdir(log_folder) 
           if f.endswith('.dt')]
    fns = sorted(fns, key=lambda x: int(x.split('.')[0]))

    data = []
    logger.info('extract dm and vj...')
    for caos, fn in zip(atm_caos, fns):
        abs_fn = os.path.join(log_folder, fn)
        with open(abs_fn, 'rb') as fp:
            mid_rlts = pickle.load(fp)
            for r in mid_rlts:
                vj = r['vj']
                dm = r['dm']
                data.append({
                    'caos': caos,
                    'dm': dm,
                    'vj': vj
                })
    
    with open(cfg.data.of_name, 'wb') as fp:
        pickle.dump(data, fp)
    return data


class JDatabase:
    def __init__(self, atm_xyzs) -> None:
        self.atm_xyzs = atm_xyzs
        atm_cfgs = [get_caos(xyz) for xyz in self.atm_xyzs]
        self.atm_obj = [cfg['mol'] for cfg in atm_cfgs]
        self.atm_scf = [cfg['scf'] for cfg in atm_cfgs]
        self.atm_caos = [cfg['caos'] for cfg in atm_cfgs]
        
    def getJ(self, i, dm):
        return self.atm_scf[i].get_j(self.atm_obj[i], dm, self.atm_caos[i])

def amino_eri(input_pth, output_pth):
    folder_pth = input_pth
    file_pths = os.listdir(folder_pth)

    cao_dict = {}
    
    for fname in file_pths:
        # folder structure:
        # xxx/amino_name/struc.xyz
        abs_file_pth = os.path.join(folder_pth, fname)
        if not os.path.isdir(abs_file_pth):
            continue
        txt_file_pth = os.path.join(abs_file_pth, 'struc.xyz')
        with open(txt_file_pth, 'r') as fp:
            txt = fp.readlines()
            atm_xyz = ''.join(txt[2:])
            print('+'*20, fname, '+'*20)
            
            cao = get_caos(atm_xyz, kernel=True)
            cao_dict[fname] = cao
    
    with open(output_pth, 'wb') as fp:
        pickle.dump(cao_dict, fp)

def ISOL24_eri(input_pth, output_pth):
    folder_pth = input_pth
    file_pths = os.listdir(folder_pth)

    cao_dict = {}
    
    c_s = pd.read_csv('/mnt/data/charge.csv', 'r')
    c_slist = list(c_s['cha'])

    for fname in file_pths:
        abs_file_pth = os.path.join(folder_pth, fname)
        molenum = int(re.findall(r"\d+\.?\d*",fname)[0]) - 1
        if not os.path.isdir(abs_file_pth):
            continue
        txt_file_pth = os.path.join(abs_file_pth, 'struc.xyz')
        with open(txt_file_pth, 'r') as fp:
            txt = fp.readlines()
            atm_xyz = ''.join(txt[2:])
            print('+'*20, fname, '+'*20)
            cao = get_caos(atm_xyz, c_slist[molenum], kernel=True)
            cao_dict[fname] = cao
        with open(output_pth, 'wb') as fp:
            pickle.dump(cao_dict[fname], fp)


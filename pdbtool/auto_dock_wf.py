import os
import argparse
import subprocess
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdMolTransforms import ComputeCentroid
from vina import Vina
v = Vina(sf_name='vina')


def getParser():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_pdb', help='input pdb file')
    parser.add_argument('--chain', default='A')
    parser.add_argument('--ligand', default=None)
    parser.add_argument('--output_dir', default='output')
    return parser

def check_dir(outputdir):
    if not os.path.isdir(outputdir):
        os.mkdir(outputdir)
    

def prepare_ligand_and_protein(inputfile, chain, ligand, outputdir):
    out_r_pdb = inputfile.replace('.pdb', f'_{chain}_receptor.pdb')
    out_r_pdb = os.path.join(outputdir, out_r_pdb)
    out_r_pdbqt = out_r_pdb.replace('.pdb', '.pdbqt')
    subprocess.call(f'pdb_selchain -{chain} {inputfile} | pdb_delhetatm | pdb_tidy > {out_r_pdb}', shell=True)
    subprocess.call(f'prepare_receptor -r {out_r_pdb} -o {out_r_pdbqt} -A checkhydrogens', shell=True)

    if ligand is not None:
        out_l_pdb = inputfile.replace('.pdb', f'_{chain}_ligand.pdb')
        out_l_pdb = os.path.join(outputdir, out_l_pdb)
        out_l_pdbqt = out_l_pdb.replace('.pdb', '.pdbqt')
        out_l_smi = out_l_pdb.replace('.pdb', '.smi')
        subprocess.call(f'pdb_selchain -{chain} {inputfile} | pdb_selresname -{ligand} > {out_l_pdb}', shell=True)
        subprocess.call(f'mk_prepare_ligand.py -i {out_l_pdb} -o {out_l_pdbqt} --add_hydrogen --pH 7.4', shell=True)
        #subprocess.call(f"obabel -ipdb {out_l_pdb} -osmi -O {out_l_smi}", shell=True)
        return {'receptor_pdbqt':out_r_pdbqt,
                'ligand_pdbqt':out_l_pdbqt,
                'ligand_pdb':out_l_pdb,
                #'ligand_smi':out_l_smi
                }
    return {'receptor_pdbqt':out_r_pdbqt}



if __name__=='__main__':
    parser = getParser()
    args = parser.parse_args()
    check_dir(args.output_dir)
    res = prepare_ligand_and_protein(args.input_pdb, args.chain, args.ligand, args.output_dir)
    if len(res) == 1:
        print('receptor PDBQT file is prepared')
        exit()

    ligand = Chem.MolFromPDBFile(res['ligand_pdb'])
    centroid = ComputeCentroid(ligand.GetConformer())
    # The all bond order of RDMol from the PDB is 1 but I think it doesn't problem to compute center of mass because it depends on location of atoms.
    # So I didnt call AssignBondOrdersFromTemplate
    print(f'Computed centroid {centroid.x:.03f}, {centroid.y:.03f}, {centroid.z:.03f}')
    v.set_receptor(res['receptor_pdbqt'])
    v.set_ligand_from_file(res['ligand_pdbqt'])
    v.compute_vina_maps(center=[centroid.x, centroid.y, centroid.z], box_size=[20, 20, 20])
    energy_minimized = v.optimize()
    print('Score after minimization : %.3f (kcal/mol)' % energy_minimized[0])
    v.write_pose(os.path.join(args.output_dir, 'receptor_ligand_minimized.pdbqt'), overwrite=True)

    # Dock the ligand
    v.dock(exhaustiveness=32, n_poses=20)
    v.write_poses(os.path.join(args.output_dir, 'dock_vina_out.pdbqt'), n_poses=10, overwrite=True)
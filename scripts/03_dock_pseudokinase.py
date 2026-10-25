# Note: single rigid receptor (AlphaFold model), Vina exhaustiveness 8 (16 for first 4 ligands); ensemble/MD refinement is future work
import sys
from rdkit import Chem; from rdkit.Chem import AllChem
from meeko import MoleculePreparation, PDBQTWriterLegacy
from vina import Vina
c=list(map(float,open('center.txt').read().split()))
done={l.split()[0] for l in open('scores.tsv')} if __import__('os').path.exists('scores.tsv') else set()
out=open('scores.tsv','a')
for line in open('ligs.tsv'):
    name,smi=line.strip().split('\t')
    if name in done: continue
    m=Chem.AddHs(Chem.MolFromSmiles(smi)); AllChem.EmbedMolecule(m,randomSeed=42); AllChem.MMFFOptimizeMolecule(m)
    s=PDBQTWriterLegacy.write_string(MoleculePreparation().prepare(m)[0])[0]
    v=Vina(sf_name='vina',seed=42,verbosity=0); v.set_receptor('kd.pdbqt'); v.set_ligand_from_string(s)
    v.compute_vina_maps(center=c,box_size=[24,24,24]); v.dock(exhaustiveness=8,n_poses=5)
    e=v.energies(n_poses=1)[0][0]; v.write_poses(f'pose_{name}.pdbqt',n_poses=1,overwrite=True)
    ha=m.GetNumHeavyAtoms(); out.write(f'{name}\t{e:.2f}\t{ha}\t{e/ha:.3f}\n'); out.flush(); print(name,e,flush=True)

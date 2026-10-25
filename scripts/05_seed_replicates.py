import sys
from rdkit import Chem; from rdkit.Chem import AllChem
from meeko import MoleculePreparation, PDBQTWriterLegacy
from vina import Vina
c=list(map(float,open('center.txt').read().split()))
L=dict(l.strip().split('\t') for l in open('ligs.tsv'))
done={tuple(l.split()[:2]) for l in open('reps.tsv')} if __import__('os').path.exists('reps.tsv') else set()
out=open('reps.tsv','a')
for name in sys.argv[1].split(','):
  for seed in (7,123):
    if (name,str(seed)) in done: continue
    m=Chem.AddHs(Chem.MolFromSmiles(L[name])); AllChem.EmbedMolecule(m,randomSeed=seed); AllChem.MMFFOptimizeMolecule(m)
    s=PDBQTWriterLegacy.write_string(MoleculePreparation().prepare(m)[0])[0]
    v=Vina(sf_name='vina',seed=seed,verbosity=0); v.set_receptor('kd.pdbqt'); v.set_ligand_from_string(s)
    v.compute_vina_maps(center=c,box_size=[24,24,24]); v.dock(exhaustiveness=8,n_poses=1)
    e=v.energies(n_poses=1)[0][0]; out.write(f'{name}\t{seed}\t{e:.2f}\n'); out.flush(); print(name,seed,e,flush=True)

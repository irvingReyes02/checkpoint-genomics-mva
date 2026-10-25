"""Structural context of BUBR1 p.Asn1002Lys on AlphaFold AF-O60566-F1 (v6).
Outputs pLDDT, 4.5 A contacts, relative SASA, distances to functional sites; trims kinase domain (766-1050)
and defines the docking box centred on the ATP-binding residues (772-780, K795, D882)."""
import gemmi, numpy as np, freesasa
s = gemmi.read_structure('data/AF-O60566-F1-model_v6.pdb'); s.remove_hydrogens()
ch = s[0][0]; R = {r.seqid.num: r for r in ch}
ca = lambda n: np.array(R[n]['CA'][0].pos.tolist())
ns = gemmi.NeighborSearch(s[0], s.cell, 6).populate()
nb = sorted({(m.to_cra(s[0]).residue.seqid.num, m.to_cra(s[0]).residue.name) for a in R[1002]
             for m in ns.find_atoms(a.pos, '\0', radius=4.5) if abs(m.to_cra(s[0]).residue.seqid.num - 1002) > 1})
print('N1002 pLDDT', R[1002][0].b_iso, '\ncontacts', nb)
for k in (775, 795, 882, 1008): print(k, R[k].name, 'CA-CA %.1f A' % np.linalg.norm(ca(k) - ca(1002)))
for i in range(len(ch) - 1, -1, -1):
    if not 766 <= ch[i].seqid.num <= 1050: del ch[i]
s.write_pdb('results/kd.pdb')
sasa = freesasa.calc(freesasa.Structure('results/kd.pdb')).residueAreas()['A']
print('relSASA N1002 %.2f' % sasa['1002'].relativeTotal)
c = np.mean([a.pos.tolist() for r in ch if r.seqid.num in list(range(772, 781)) + [795, 882] for a in r], 0)
open('results/center.txt', 'w').write(' '.join(map(str, c.round(3))))
# receptor: obabel results/kd.pdb -xr -h -p 7.4 -O results/kd.pdbqt

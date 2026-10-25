import gemmi,numpy as np,glob,os
rec=gemmi.read_structure('kd.pdb');ns=gemmi.NeighborSearch(rec[0],rec.cell,5).populate()
c=np.array(list(map(float,open('center.txt').read().split())))
KEY={772,773,774,775,776,777,778,779,780,795,882}
print('ligand\tcentroid_to_box_A\tn_contact_res\tkey_site_contacts')
for f in sorted(glob.glob('pose_*.pdbqt')):
    xyz=[];res=set()
    for l in open(f):
        if l.startswith(('ATOM','HETATM')) and l[77:79].strip()!='H' and l[77:79].strip()!='HD':
            p=np.array([float(l[30:38]),float(l[38:46]),float(l[46:54])]);xyz.append(p)
            for m in ns.find_atoms(gemmi.Position(*p),'\0',radius=4.0): res.add(m.to_cra(rec[0]).residue.seqid.num)
    d=np.linalg.norm(np.mean(xyz,0)-c)
    print(os.path.basename(f)[5:-6],f'{d:.1f}',len(res),','.join(map(str,sorted(res&KEY))),sep='\t')

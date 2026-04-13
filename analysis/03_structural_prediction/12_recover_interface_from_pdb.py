#!/usr/bin/env python3
"""12_recover_interface_from_pdb.py -- recover interface metrics from surviving PDBs.

Background:
  Post-AF2 cleanup extractor (09_extract_and_cleanup.py) was never invoked due to
  a path resolution bug, so result_model_*.pkl files were deleted without writing
  scores.json / pae_matrices.npz. Combined iptm+ptm from ranking_debug.json is the
  only pre-existing per-job metric. This script recovers what is possible from the
  surviving ranked_*.pdb files:

    - interface_plddt: mean pLDDT (B-factor of CA atoms) over chain-A residues
      whose Cb (or CA for glycine) is within CB_CUTOFF of any chain-B Cb.
    - interface_plddt_b: symmetric metric over chain B.
    - interface_plddt_joint: combined mean over both interface sides.
    - n_interface_residues_a / _b: counts of interface residues per chain.
    - mean_plddt_a / _b: whole-chain mean pLDDT per chain (sanity check).
    - per_model: list of the above per ranked_*.pdb (25 models).
    - summary_mean / summary_best: mean and best-of-25 across the above.

  What we CANNOT recover without pkl: raw ipTM (separate from pTM), PAE matrix,
  interchain_pae, ipSAE, LIS.

Usage:
  python3 12_recover_interface_from_pdb.py <results_root> [--out <tsv>] [--cb-cutoff 8.0]

  results_root: e.g. .../d1d3_fullscreen/results
  Writes one TSV row per target with per-target summary metrics.
"""
import argparse, os, glob, json, sys, math

CB_CUTOFF_DEFAULT = 8.0  # Angstroms, Cb-Cb (Cα for glycine)

def parse_pdb_atoms(path):
    # returns list of dicts: chain, resnum, resname, atom, x, y, z, bfac
    atoms=[]
    with open(path) as f:
        for ln in f:
            if not ln.startswith('ATOM'): continue
            atom=ln[12:16].strip()
            resname=ln[17:20].strip()
            chain=ln[21]
            try:
                resnum=int(ln[22:26])
                x=float(ln[30:38]); y=float(ln[38:46]); z=float(ln[46:54])
                bfac=float(ln[60:66])
            except ValueError: continue
            atoms.append((chain,resnum,resname,atom,x,y,z,bfac))
    return atoms

def per_residue_cb(atoms):
    # Returns {chain: {resnum: (x,y,z,plddt,resname)}} using CB, or CA for glycine.
    by_res={}
    for c,rn,resname,atom,x,y,z,bfac in atoms:
        key=(c,rn)
        if key not in by_res: by_res[key]=[resname,None,None]
        if atom=='CA':
            by_res[key][1]=(x,y,z,bfac,resname)
        if atom=='CB':
            by_res[key][2]=(x,y,z,bfac,resname)
    out={}
    for (c,rn),(resname,ca,cb) in by_res.items():
        coord = cb if cb is not None else ca
        if coord is None: continue
        out.setdefault(c,{})[rn]=coord
    return out

def interface_metrics(path, cutoff=CB_CUTOFF_DEFAULT):
    atoms=parse_pdb_atoms(path)
    if not atoms: return None
    res=per_residue_cb(atoms)
    chains=sorted(res.keys())
    if len(chains)<2: return None
    A,B=chains[0],chains[1]
    ra=res[A]; rb=res[B]
    # build coord lists
    la=list(ra.items()); lb=list(rb.items())
    cut2=cutoff*cutoff
    iface_a=set(); iface_b=set()
    for ia,(xa,ya,za,pa,na) in la:
        for ib,(xb,yb,zb,pb,nb) in lb:
            dx=xa-xb; dy=ya-yb; dz=za-zb
            if dx*dx+dy*dy+dz*dz <= cut2:
                iface_a.add(ia); iface_b.add(ib)
    def mean(x): return sum(x)/len(x) if x else 0.0
    plddt_a_all=[v[3] for v in ra.values()]
    plddt_b_all=[v[3] for v in rb.values()]
    plddt_a_if=[ra[i][3] for i in iface_a]
    plddt_b_if=[rb[i][3] for i in iface_b]
    return dict(
        chain_a=A, chain_b=B,
        n_res_a=len(ra), n_res_b=len(rb),
        mean_plddt_a=mean(plddt_a_all),
        mean_plddt_b=mean(plddt_b_all),
        n_interface_a=len(iface_a),
        n_interface_b=len(iface_b),
        interface_plddt_a=mean(plddt_a_if),
        interface_plddt_b=mean(plddt_b_if),
        interface_plddt_joint=mean(plddt_a_if+plddt_b_if),
        interface_residues_a=sorted(iface_a),
        interface_residues_b=sorted(iface_b),
    )

def summarize_target(subdir, cutoff):
    pdbs=sorted(glob.glob(os.path.join(subdir,'ranked_*.pdb')))
    if not pdbs: return None
    per=[]
    for p in pdbs:
        m=interface_metrics(p,cutoff)
        if m is None: continue
        m['model']=os.path.basename(p)
        per.append(m)
    if not per: return None
    def col(k): return [p[k] for p in per]
    def mean(x): return sum(x)/len(x) if x else 0.0
    def best(x): return max(x) if x else 0.0
    summary=dict(
        n_models=len(per),
        mean_plddt_a=mean(col('mean_plddt_a')),
        mean_plddt_b=mean(col('mean_plddt_b')),
        interface_plddt_joint_mean=mean(col('interface_plddt_joint')),
        interface_plddt_joint_best=best(col('interface_plddt_joint')),
        interface_plddt_a_mean=mean(col('interface_plddt_a')),
        interface_plddt_b_mean=mean(col('interface_plddt_b')),
        n_interface_a_mean=mean(col('n_interface_a')),
        n_interface_b_mean=mean(col('n_interface_b')),
    )
    return summary, per

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('results_root')
    ap.add_argument('--out',default=None)
    ap.add_argument('--cb-cutoff',type=float,default=CB_CUTOFF_DEFAULT)
    ap.add_argument('--json-dir',default=None,help='if set, write per-target JSON with per-model detail')
    a=ap.parse_args()
    root=a.results_root
    rows=[]
    for target in sorted(os.listdir(root)):
        tdir=os.path.join(root,target)
        if not os.path.isdir(tdir): continue
        # subdir is sflt1_vs_<target>
        subs=[d for d in glob.glob(os.path.join(tdir,'*')) if os.path.isdir(d)]
        if not subs: continue
        res=summarize_target(subs[0],a.cb_cutoff)
        if res is None: continue
        summary,per=res
        summary['target']=target
        rows.append(summary)
        if a.json_dir:
            os.makedirs(a.json_dir,exist_ok=True)
            with open(os.path.join(a.json_dir,f'{target}.json'),'w') as f:
                json.dump(dict(summary=summary,per_model=per),f,indent=1)
        print(f"[{len(rows):4d}] {target:<25} iPLDDT_joint_mean={summary['interface_plddt_joint_mean']:.1f} n_if_a={summary['n_interface_a_mean']:.1f}",file=sys.stderr)
    if not rows: return
    cols=['target','n_models','interface_plddt_joint_mean','interface_plddt_joint_best','interface_plddt_a_mean','interface_plddt_b_mean','n_interface_a_mean','n_interface_b_mean','mean_plddt_a','mean_plddt_b']
    out=a.out or os.path.join(os.path.dirname(root.rstrip('/')),'interface_recovery.tsv')
    with open(out,'w') as f:
        f.write('\t'.join(cols)+'\n')
        for r in rows:
            f.write('\t'.join(f"{r[c]:.3f}" if isinstance(r[c],float) else str(r[c]) for c in cols)+'\n')
    print(f'wrote {out} ({len(rows)} targets)',file=sys.stderr)

if __name__=='__main__': main()

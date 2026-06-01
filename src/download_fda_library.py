#!/usr/bin/env python3
"""Download FDA-approved drugs from ChEMBL — optimized with chunked writes."""
import csv, sys, time
from chembl_webresource_client.new_client import new_client

molecule = new_client.molecule
OUT = "data/fda_approved_full.csv"
CHUNK = 500  # ChEMBL pagination is fast for small ranges

print(f"Fetching FDA-approved drugs (max_phase=4)...", flush=True)

# Get total count first
total_raw = molecule.filter(max_phase=4).only(['molecule_chembl_id'])
total_count = len(total_raw[:20000])
print(f"Total candidates: {total_count}", flush=True)

fieldnames = ['chembl_id', 'name', 'smiles', 'mw', 'alogp', 'hba', 'hbd', 'psa',
              'first_approval', 'indication_class']
drugs = []
fetched = 0
written = 0
header_written = False

for offset in range(0, min(total_count, 20000), CHUNK):
    batch = molecule.filter(max_phase=4).only(
        ['molecule_chembl_id', 'pref_name', 'molecule_structures',
         'molecule_properties', 'first_approval', 'indication_class']
    )[offset:offset + CHUNK]
    
    fetched += len(batch)
    
    batch_drugs = []
    for mol in batch:
        smiles = None
        structs = mol.get('molecule_structures')
        if structs and 'canonical_smiles' in structs:
            smiles = structs['canonical_smiles']
        if not smiles:
            continue
        
        props = mol.get('molecule_properties') or {}
        try:
            mw = float(props.get('full_mol_weight') or props.get('mw_freebase') or 0)
        except (ValueError, TypeError):
            mw = 0
        
        if mw and (mw < 150 or mw > 800):
            continue
        
        batch_drugs.append({
            'chembl_id': mol['molecule_chembl_id'],
            'name': mol.get('pref_name', ''),
            'smiles': smiles,
            'mw': mw,
            'alogp': props.get('alogp', ''),
            'hba': props.get('hba', ''),
            'hbd': props.get('hbd', ''),
            'psa': props.get('psa', ''),
            'first_approval': mol.get('first_approval', ''),
            'indication_class': mol.get('indication_class', ''),
        })
    
    drugs.extend(batch_drugs)
    
    # Write incrementally
    mode = 'w' if not header_written else 'a'
    with open(OUT, mode, newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        if not header_written:
            writer.writeheader()
            header_written = True
        writer.writerows(batch_drugs)
    
    written += len(batch_drugs)
    print(f"  {fetched}/{total_count} fetched, {written} drug-like saved", flush=True)
    
    if len(batch) < CHUNK:
        break
    
    time.sleep(0.2)

print(f"\nFinal: {written} FDA-approved drug-like molecules saved to {OUT}", flush=True)

# SMILES-only version
smiles_out = "data/fda_approved_smiles.csv"
with open(smiles_out, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['Name', 'SMILES'])
    for d in drugs:
        name = d['name'] or d['chembl_id']
        writer.writerow([name, d['smiles']])
print(f"SMILES version: {smiles_out} ({len(drugs)} entries)", flush=True)
print("Done!", flush=True)

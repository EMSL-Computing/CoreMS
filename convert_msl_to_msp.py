#!/usr/bin/env python
"""Convert MetabRef JSON format to standard MSP format for GCMS libraries"""

import json
from pathlib import Path

def convert_metabref_json_to_msp(json_file, msp_file):
    """
    Convert MetabRef JSON format to standard MSP format with all metadata.
    
    Parameters
    ----------
    json_file : Path
        Input JSON file from MetabRef API
    msp_file : Path
        Output MSP file path
    """
    with open(json_file, 'r') as f:
        data = json.load(f)
    
    # Handle both direct list and {"GC-MS": [...]} wrappers
    if isinstance(data, dict) and "GC-MS" in data:
        compounds = data["GC-MS"]
    else:
        compounds = data
    
    with open(msp_file, 'w') as f_out:
        for compound in compounds:
            spectrum_data = compound.get("spectrum_data", {})
            
            # Write molecular metadata
            f_out.write(f"Name: {compound.get('molecule_name', 'Unknown')}\n")
            if compound.get('casno'):
                f_out.write(f"CAS: {compound['casno']}\n")
            if compound.get('formula'):
                f_out.write(f"Formula: {compound['formula']}\n")
            if compound.get('inchikey'):
                f_out.write(f"InChIKey: {compound['inchikey']}\n")
            if compound.get('inchi'):
                f_out.write(f"InChI: {compound['inchi']}\n")
            if compound.get('smiles'):
                f_out.write(f"SMILES: {compound['smiles']}\n")
            if compound.get('pubchem'):
                f_out.write(f"PubChem: {compound['pubchem']}\n")
            if compound.get('chebi'):
                f_out.write(f"ChEBI: {compound['chebi']}\n")
            if compound.get('kegg'):
                f_out.write(f"KEGG: {compound['kegg']}\n")
            if compound.get('refmet'):
                f_out.write(f"RefMet: {compound['refmet']}\n")
            if compound.get('iupac_name'):
                f_out.write(f"IUPAC_Name: {compound['iupac_name']}\n")
            
            # Write spectrum metadata
            if spectrum_data.get('ri'):
                f_out.write(f"RI: {spectrum_data['ri']}\n")
            if spectrum_data.get('rt'):
                f_out.write(f"RetentionTime: {spectrum_data['rt']}\n")
            if spectrum_data.get('comments'):
                f_out.write(f"Comment: {spectrum_data['comments']}\n")
            if spectrum_data.get('derivative'):
                f_out.write(f"Derivative: {spectrum_data['derivative']}\n")
            if spectrum_data.get('usi'):
                f_out.write(f"USI: {spectrum_data['usi']}\n")
            
            # Parse and write peaks
            mz_string = spectrum_data.get('mz', '')
            if mz_string:
                import re
                peaks = re.findall(r'\((\d+),(\d+)\)', mz_string)
                f_out.write(f"Num Peaks: {len(peaks)}\n")
                for mz, intensity in peaks:
                    f_out.write(f"{mz} {intensity}\n")
            
            f_out.write("\n")

if __name__ == "__main__":
    data_dir = Path("corems/molecular_id/data")
    
    # Convert GCMS library from JSON
    gcms_json = data_dir / "gcms_library_raw.json"
    gcms_msp = data_dir / "PNNLMetV20191015.msp"
    if gcms_json.exists():
        convert_metabref_json_to_msp(gcms_json, gcms_msp)
        print(f"Converted {gcms_json} to {gcms_msp}")
    
    # Convert FAMES library from JSON
    fames_json = data_dir / "fames_library_raw.json"
    fames_msp = data_dir / "FAMES_REF.msp"
    if fames_json.exists():
        convert_metabref_json_to_msp(fames_json, fames_msp)
        print(f"Converted {fames_json} to {fames_msp}")

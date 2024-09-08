import pandas as pd
import pymol
from pymol import cmd
from center_of_mass import com, get_com
import numpy as np
from scipy.spatial.distance import cdist

### HELPER FUNCTION ###
def get_partial_match(sequence, subsequence):
    """
    Searches for a partial match of the subsequence (where part of the subseqeunce appears
        in the sequence at the edge of the sequence)
    """
    subseq_padding = len(subsequence) // 2
    for size in range(len(subsequence) - 1, subseq_padding, -1):
        if sequence[:size] == subsequence[-size:]:
            return size - subseq_padding - 1
        if sequence[-size:] == subsequence[:size]:
            return len(sequence) - size + subseq_padding
    return -1

### PDB TO SEQ FILE
pdb_df = pd.read_csv("pdb_to_pymol_seq.csv")

# Maps the pdb to its sequence + offset
pdb_to_seq = {pdb: (seq, offset) for pdb, seq, offset in zip(pdb_df["pdb"], pdb_df["seq"], pdb_df["offset"])}

# First grab the information about pSTY binding
df = pd.read_csv("pSTY_binding.csv")

data = {"pdb": [], "protid": [], "site_true": [], "site_pdb": [], "start": [], "stop": [], "com_dist": [], "min_dist": [], "motif": [], "aa": [], "description": [], "sequence": [], "name": [], "pymol_aa": [], "cog": [], "bind_start_new": [], "bind_stop_new": []}
erroreous_data = {"pdb": [], "protid": [], "site_true": [], "site_pdb": [], "start": [], "stop": [], "motif": [], "aa": [], "description": [], "sequence": [], "name": []}

aa_to_name = {
    "Y": "TYR",
    "T": "THR",
    "S": "SER",
}

# Generated dictionary of 3 letter to 1 letter
three_to_one_letter = {k.upper(): v for k, v in {'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C', 'Gln': 'Q', 'Glu': 'E', 'Gly': 'G', 'His': 'H', 'Ile': 'I', 'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P', 'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V'}.items()}

not_found = []

for idx, row in df.iterrows():
    # Extract the necessary information
    pdb = row["PDB"]
    protid = row["ACC_ID"]
    site = row["Updated_Site"]
    site_true = row["Site"]
    start = row["begin"]
    stop = row["end"]
    motif = row["SITE_7_AA"].upper()
    sequence = row["sequence"]
    name = row["name"]
    description = row["description"]
    aa = motif[7]

    # Some of the motif alignment results are mislabelled...
    if aa not in ("Y", "S", "T"):
        print("MISLABELLED")

    # Reinit pymol and grab the pdb
    cmd.reinitialize()
    cmd.fetch(pdb)

    # Try to get the corresponding PDB
    try:
        seq, pdb_offset = pdb_to_seq[pdb]
    except:
        # we can't yet extract the seq from the pdb
        # NOTE: this should NEVER be the case since I've already filtered only the high quality PDBs
        not_found.append({"error": "pdb not found", "pdb": pdb, "uniprot": protid, "pymol_seq": "", "start": start, "stop": stop, "wanted_sequence": sequence, "name": name, "description": description, "motif": motif})

    # Try to find the motif. if cannot, run the partial match to check the ends. weird offsets for the 1-indexing
    EXP_SITE = seq.upper().find(motif.upper()) + 7 + pdb_offset

    if EXP_SITE == -1 + 7 + pdb_offset: 
        # Run partial match algorithm and quit if it's unsuccessful
        EXP_SITE = get_partial_match(seq, motif) + pdb_offset
        if EXP_SITE == -1 + pdb_offset:
            # Unfortunately not much to do here but can try to manually check
            not_found.append({"error": "motif not in pdb", "pdb": pdb, "uniprot": protid, "pymol_seq": seq, "start": start, "stop": stop, "wanted_sequence": sequence, "name": name, "description": description, "motif": motif})
            continue

    # Verify we have a STY
    if aa not in ("Y", "S", "T"):
        assert seq[EXP_SITE - pdb_offset] in ("Y", "S", "T")
        aa = seq[EXP_SITE - pdb_offset]
    else:
        assert seq[EXP_SITE - pdb_offset] == aa, f"{EXP_SITE}, {EXP_SITE - pdb_offset}, {pdb}, {aa}, {seq}, {seq[EXP_SITE - pdb_offset]} {site_true} {site}"

    aa_name = aa_to_name[aa]

    ans = set()
    chains = cmd.get_chains(pdb)
    cmd.select("site1", f"resi {EXP_SITE} and chain {chains[0]}")
    cmd.iterate("site1", "ans.add(resn)", space={"ans": ans})
    delta = EXP_SITE - site

    pymol_aa = list(ans)[0] if len(ans) > 0 else "ERROR"

    # get the site
    aa_to_oxy_name = {
        "Y": "OH",
        "S": "OG",
        "T": "OG1",
        "TYR": "OH",
        "SER": "OG",
        "THR": "OG1",
    }
    oxy_name = aa_to_oxy_name.get(aa, None) 

    binding_site_start = EXP_SITE + (start - site) - pdb_offset
    binding_site_stop = EXP_SITE + (stop - site) - pdb_offset

    # If the site isn't in the PDB (example: 6RCW, 462) 
    if binding_site_start < 0:
        not_found.append({"error": "site not in pdb", "pdb": pdb, "uniprot": protid, "pymol_seq": seq, "start": start, "stop": stop, "wanted_sequence": sequence, "name": name, "description": description, "motif": motif})
        continue
    pymol_seq = seq[binding_site_start:binding_site_stop + 1]
    if pymol_seq != sequence:
        not_found.append({"error": "protid and pymol seq don't match", "pdb": pdb, "uniprot": protid, "pymol_seq": seq, "start": start, "stop": stop, "wanted_sequence": sequence, "name": name, "description": description, "motif": motif})
        continue

    # Verify the sequence is correct
    assert pymol_seq == sequence, f"binding site {pymol_seq} and {sequence} do not match for pdb {pdb} {start=}, {stop=}, {binding_site_start}, {binding_site_stop}, {seq}, {EXP_SITE=}, {site=}, {pdb_offset=}"

    try:

        # calculate the COM (if the mass is non-zero) for the binding site and the site1
        center_of_mass = cmd.centerofmass(f"binding_site_{idx}")
        center_of_mass_site = cmd.centerofmass("site1") 

        # calculate the coordinates of the COM
        domain_coords = cmd.get_coords(f"binding_site_{idx}")

        # calculate the distance manually
        dist = ((center_of_mass_site[0] - center_of_mass[0])**2 + (center_of_mass_site[1] - center_of_mass[1])**2 + (center_of_mass_site[2] - center_of_mass[2])**2)**0.5

        # calculate distance from a set of points
        min_dist = np.min(cdist([center_of_mass_site], domain_coords))

        # calculate the cog
        cog_coord = domain_coords.mean(axis=0)
        cog_dist = ((center_of_mass_site[0] - cog_coord[0])**2 + (center_of_mass_site[1] - cog_coord[1])**2 + (center_of_mass_site[2] - cog_coord[2])**2)**0.5

        # update data with the info
        data["pdb"].append(pdb)
        data["protid"].append(protid)
        data["site_pdb"].append(EXP_SITE)
        data["site_true"].append(site_true)
        data["start"].append(start)
        data["stop"].append(stop)
        data["com_dist"].append(dist)
        data["min_dist"].append(min_dist)
        data["sequence"].append(sequence)
        data["name"].append(name)
        data["description"].append(description)
        data["motif"].append(motif)
        data["aa"].append(aa)
        data["pymol_aa"].append(pymol_aa)
        data["cog"].append(cog_dist)
        data["bind_start_new"].append(binding_site_start)
        data["bind_stop_new"].append(binding_site_stop)
    except Exception as e:
        # only error should be that the mass is zero
        print("ERROR::::", e)
        not_found.append({"error": e, "pdb": pdb, "uniprot": protid, "pymol_seq": seq, "start": start, "stop": stop, "wanted_sequence": sequence, "name": name, "description": description, "motif": motif})


# create the CSV
pd.DataFrame(data).to_csv(f"results_with_distance.csv", index=False)
pd.DataFrame(not_found).to_csv(f"pdb_without_distance.csv", index=False)

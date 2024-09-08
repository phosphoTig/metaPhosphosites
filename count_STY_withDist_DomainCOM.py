from pymol import cmd
import pickle
from scipy.spatial import cKDTree

def get_all_aa(aa, chain):
    """
    Selects a specific amino acid from a specific chain
    """
    name = f"{aa}_selection"
    cmd.select(name, f"resn {aa} and chain {chain}")
    return name

# Open an existing pickle file which maps PDBs to domains we're interested in
with open("pdb_to_domains.pkl", "rb") as f:
    pdb_to_domain = pickle.load(f)

pdb_domain_cnt_10 = {}
pdb_domain_cnt_25 = {}

for pdb, domain_arr in pdb_to_domain.items():
    cmd.reinitialize()
    cmd.fetch(pdb)

    chains = cmd.get_chains(pdb)

    # For now, grab all of the STY from chain 0
    tyr_select = get_all_aa("TYR", chains[0])
    thr_select = get_all_aa("THR", chains[0])
    ser_select = get_all_aa("SER", chains[0])

    # Add all of the coordinates
    tyr_coord, thr_coord, ser_coord = set(), set(), set()
    cmd.iterate_state(1, tyr_select, "ans.add((x, y, z))", space={"ans": tyr_coord})
    cmd.iterate_state(1, thr_select, "ans.add((x, y, z))", space={"ans": thr_coord})
    cmd.iterate_state(1, ser_select, "ans.add((x, y, z))", space={"ans": ser_coord})

    tyr_coord = list(tyr_coord)
    thr_coord = list(thr_coord)
    ser_coord = list(ser_coord)

    # Create a cKDTree for quick coordinate comparison
    tyr_kdtree = cKDTree((tyr_coord))
    thr_kdtree = cKDTree((thr_coord))
    ser_kdtree = cKDTree((ser_coord))

    # Map the coordinates to a residue
    coord_to_resn = {}
    cmd.iterate_state(1, tyr_select, "ans[(x, y, z)] = (resi)", space={"ans": coord_to_resn})
    cmd.iterate_state(1, thr_select, "ans[(x, y, z)] = (resi)", space={"ans": coord_to_resn})
    cmd.iterate_state(1, ser_select, "ans[(x, y, z)] = (resi)", space={"ans": coord_to_resn})

    # Iterate over the domains
    for begin, end in domain_arr:
        # Try to calculate the COM
        try:
            cmd.select(f"binding_site_{begin}_{end}", f"resi {begin}-{end} and chain {chains[0]}")
            coord = cmd.centerofmass(f"binding_site_{begin}_{end}")
        except Exception as e:
            # only error should be that the mass is zero
            print(e)
            continue
            

        # Get the number to STY 10 A away
        num_tyr = len({(coord_to_resn[tyr_coord[co]]) for co in tyr_kdtree.query_ball_point(coord, 10)})
        num_thr = len({(coord_to_resn[thr_coord[co]]) for co in thr_kdtree.query_ball_point(coord, 10)})
        num_ser = len({(coord_to_resn[ser_coord[co]]) for co in ser_kdtree.query_ball_point(coord, 10)})

        pdb_domain_cnt_10[f"{pdb}_{begin}_{end}"] = (num_ser, num_thr, num_tyr)

        # Get the number to STY 25 A away
        num_tyr = len({(coord_to_resn[tyr_coord[co]]) for co in tyr_kdtree.query_ball_point(coord, 25)})
        num_thr = len({(coord_to_resn[thr_coord[co]]) for co in thr_kdtree.query_ball_point(coord, 25)})
        num_ser = len({(coord_to_resn[ser_coord[co]]) for co in ser_kdtree.query_ball_point(coord, 25)})
        pdb_domain_cnt_25[f"{pdb}_{begin}_{end}"] = (num_ser, num_thr, num_tyr)

# Save the results

with open("25_cnt_w_chain.pkl", "wb") as f:
    pickle.dump(pdb_domain_cnt_25, f)

with open("12_cnt_w_chain.pkl", "wb") as f:
    pickle.dump(pdb_domain_cnt_10, f)

#Author = Shreya C., github.com/GenericP3rson 
#Reviewed/edits by = Tigist Tamir, github.com/phosphoTig

# Get binding domain information from uniprot for the csv
import requests
import pandas as pd
import sys

def get_binding_sites(uniprot_id):
    # UniProt API endpoint for retrieving protein information
    url = f'https://www.ebi.ac.uk/proteins/api/proteins/{uniprot_id}?sequence=true'

    try:
        # Make a GET request to the UniProt API
        response = requests.get(url, headers={'Accept': 'application/json'})
        response.raise_for_status()  # Check for HTTP request errors

        # Parse the JSON response
        protein_data = response.json()
        protein_id = protein_data["id"].replace("_HUMAN", "")

        # Check if there is information about features (including binding sites)
        if 'features' in protein_data and "sequence" in protein_data:
            sequence = protein_data['sequence']['sequence']
            binding_sites = []

            for feature in protein_data["features"]:
                if feature["type"] in ["BINDING", "ACT_SITE"]:
                    feature_type = feature["type"]
                    description = feature.get("description", "")
                    begin = int(feature["begin"])
                    end = int(feature["end"])
                    seq_segment = sequence[begin-1:end]
                    name = f'{protein_id}_{feature["ligand"]["name"]}' if feature_type == "BINDING" else f'{protein_id}_{feature_type}'

                    binding_sites.append({
                        "id": uniprot_id,
                        "description": description,
                        "name": name,
                        "begin": begin,
                        "end": end,
                        "sequence": seq_segment
                    })

            return binding_sites
        else:
            return None

    except requests.RequestException as e:
        print(f"Error: Unable to retrieve data for UniProt ID {uniprot_id}. Exception: {e}")
        return None

# Example usage
def main(start, stop):
    pathways = pd.read_csv("1_filtered_phosphosites.csv")
    uniprot_ids = list(set(pathways["ACC_ID"][start:stop]))
    
    all_binding_sites = []

    for uniprot_id in uniprot_ids:
        binding_sites = get_binding_sites(uniprot_id)
        if binding_sites:
            print(f"{uniprot_id}: {len(binding_sites)} binding sites found")
            all_binding_sites.extend(binding_sites)
    
    output_file = f"tmp_binding_site_info_{start}_{stop}.csv"
    pd.DataFrame(all_binding_sites).to_csv(output_file, index=False)
    print(f"Binding sites information saved to {output_file}")

if __name__ == "__main__":
    start = int(sys.argv[1])
    stop = int(sys.argv[2])
    main(start, stop)




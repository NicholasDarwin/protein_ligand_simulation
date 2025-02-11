import requests

def fetch_pdb(protein_name):
    url = f"https://data.rcsb.org/rest/v1/core/entry/{protein_name}"
    response = requests.get(url)
    if response.status_code == 200:
        pdb_id = response.json()["rcsb_id"]
        pdb_url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        pdb_response = requests.get(pdb_url)
        with open(f"data/{protein_name}.pdb", "wb") as file:
            file.write(pdb_response.content)
        return f"data/{protein_name}.pdb"
    else:
        return None

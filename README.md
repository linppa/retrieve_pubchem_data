# Retrive Pubchem Data

This script takes a list of drugs/compounds from the Comparative Toxicogenomics
Database & retrieves specific PubChem data
for each compound using the PUG REST API. The data is then saved in a .tsv file.

# Database Sources
- https://ctdbase.org/downloads/
- https://pubchem.ncbi.nlm.nih.gov/

# Usage

```bash
1. I used the CTD drug database as input, containing two columns, 'chemical_name' and 'cas_no'

2. Run 'retrieve_cid.py' to retrieve data for compounds. Update INPUT_FILE and OUTPUT_FILE in the script.

3. Run the output retrieved with previous compound script to check for substances that were not retrieved previously. Update INPUT_FILE and OUTPUT_FILE in the 'retrieve_sid.py' script.
```

# Output
Final output gives a file with compound/substance information from pubchem. The
columns are as follows: 

- chemical_name
- cas_noUNII	
- Pharmacology Classification	
- Pharmacology Description	
- Indication
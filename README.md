# Retrive Pubchem Data
- [Linda Quach](https://github.com/linppa)

This python script takes a list of drugs/compounds from the Comparative Toxicogenomics
Database (CTD) & retrieves specific PubChem json data for each compound/drug using the PUG
REST API. The data is then saved in a .tsv file.

# Database Sources
- https://ctdbase.org/downloads/
- https://pubchem.ncbi.nlm.nih.gov/

# Requirements
Please make sure to install all the necessary dependencies by running the following command:
```bash
pip install -r requirements.txt
```

# Usage
1. Download CTD drug database tsv as input, which contains two columns:
   `chemical_name` and `cas_no`. Alternatively, a cleaned version of the CTD is included in the `drug_input` folder.
   
2. Update `CTD_INPUT_FILE` in the script to the location of downloaded data.
   Update `OUTPUT_FILE` to the desired output location.

3. Run `retrieve_cid_chemical_name.py` to retrieve data for compounds. Update INPUT_FILE and OUTPUT_FILE in the script.

# Output
Final output gives a file with retrieved information from pubchem. The
retrived columns, if available on PubChem, are as follows: 

- `generic_name` -- record title or listing name of compound

- `url` -- URL to PubChem record

- `cas_no` -- A CAS Registry Number (also called CAS RN or CAS Number) is a proprietary registry number assigned by the Chemical Abstracts Service (CAS) division of the American Chemical Society (ACS).

- `unii_no` -- UNique Ingredient Identifier (UNII) code for this chemical. It is a non-proprietary registry number assigned by the U.S. Food and Drug Administration (FDA).

- `indications` -- A drug indication refers to any valid reason (e.g., signs/symptoms or diseases/disorders) to use a certain drug. 

- `pharm_class` -- Pharmacological action classes from the Medical Subject Headings (MeSH) thesaurus.

- `pharm_desc` -- Pharmacology classification detailed description

- `hazard_statements` -- GHS (Globally Harmonized System of Classification and Labelling of Chemicals) to identify hazardous chemicals and to inform users about these hazards.

# Alternative Scripts
- `retrieve_cid_generic_name.py` -- Retrieves PubChem data by generic name or CAS number.
- `retrieve_cid_dmf.py` -- Retrieves PubChem data by Drug Master File (DMF) name
  or CAS number.
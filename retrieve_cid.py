import requests
import csv
import pandas as pd
import time

CTD_INPUT_FILE = 'drug_input/ctd_drugs.tsv'
OUTPUT_FILE = 'drug_output/ctd_compounds_substances_.tsv'

# rate limit for pubchem
REQUESTS_PER_SECOND = 5

def get_cid_by_cas(cas_number):
    search_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{cas_number}/cids/JSON"
    response = requests.get(search_url)
    response_json = response.json()
    if 'IdentifierList' in response_json:
        return response_json['IdentifierList']['CID']
    else:
        raise ValueError(f"No CID for CAS {cas_number}")

def get_cid_by_name(drug_name):
    search_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{drug_name}/cids/JSON"
    response = requests.get(search_url)
    response_json = response.json()
    if 'IdentifierList' in response_json:
        return response_json['IdentifierList']['CID']
    else:
        raise ValueError(f"No CID for NAME {drug_name}")

def get_compound_info(cid):
    info_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{cid}/JSON/"
    response = requests.get(info_url)
    return response.json()

def extract_info(data):
    compound_info = {}
    compound_info['Drug Name'] = data['Record']['RecordTitle']
    sections = data['Record']['Section']

    for section in sections:
        # ======= 1. Names and Identifiers =========
        if section['TOCHeading'] == 'Names and Identifiers':
            for subsection in section['Section']:
                if subsection['TOCHeading'] == 'Other Identifiers':
                    for subsubsection in subsection['Section']:
                        if subsubsection['TOCHeading'] == 'CAS':
                            compound_info['CAS Number'] = subsubsection['Information'][0]['Value']['StringWithMarkup'][0]['String']
                        if subsubsection['TOCHeading'] == 'UNII':
                            compound_info['UNII'] = subsubsection['Information'][0]['Value']['StringWithMarkup'][0]['String']
    
        # ======= 7. Drug and Medication Information =========
        if section['TOCHeading'] == 'Drug and Medication Information':
            for subsection in section['Section']:
                if subsection['TOCHeading'] == 'Drug Indication':
                    compound_info['Indication'] = subsection['Information'][0]['Value']['StringWithMarkup'][0]['String']
        
        # ======= 10. Pharmacology and Biochemistry =========
        if section['TOCHeading'] == 'Pharmacology and Biochemistry':
            for subsection in section['Section']:
                if subsection['TOCHeading'] == 'MeSH Pharmacological Classification':
                    compound_info['Pharmacology Classification'] = subsection['Information'][0]['Name']
                    compound_info['Pharmacology Description'] = subsection['Information'][0]['Value']['StringWithMarkup'][0]['String']

    return compound_info

def main(input_filename, output_filename):
    drug_data = pd.read_csv(input_filename, sep='\t')
    results = []

    for index, row in drug_data.iterrows():
        cas_number = row['cas_no']
        chemical_name = row['chemical_name']
        try:
            if pd.notna(cas_number):
                cids = get_cid_by_cas(cas_number)
                if len(cids) > 1:
                    print(f"======================================================\n")
                    print(f"Multiple CIDs found.)\n")
                cid = cids[0]
            else:
                raise ValueError(f"Missing CAS for drug '{chemical_name}'\n")
            
            print(f"Processing {cas_number} | {chemical_name} | {cid}")
            data = get_compound_info(cid)
            compound_info = extract_info(data)
            combined_info = {**row.to_dict(), **compound_info}  # merge original row data with new data
            results.append(combined_info)
            time.sleep(1 / REQUESTS_PER_SECOND)                 # rate limiting
            print(f"======================================================\n")

        except Exception as e:
            print(f"======================================================\n")
            print(f"{e}; Checking drug name '{chemical_name}'\n")
            try:
                cids = get_cid_by_name(chemical_name)
                if len(cids) > 1:
                    print(f"Multiple CIDs found.)\n")
                cid = cids[0]

                print(f"Processing {cas_number} | {chemical_name} | {cid}")
                data = get_compound_info(cid)
                compound_info = extract_info(data)
                combined_info = {**row.to_dict(), **compound_info}  # merge original row data with new data
                results.append(combined_info)
                time.sleep(1 / REQUESTS_PER_SECOND)                 # rate limiting
                print(f"======================================================\n")

            except Exception as e2:
                print(f"{e2}")
                results.append(row.to_dict())

    # convert results to df & save
    results_df = pd.DataFrame(results)
    results_df.to_csv(output_filename, sep='\t', index=False)
    print(f"======================================================\n")
    print(f"Data saved to {output_filename}.\n")
    print(f"======================================================\n")

if __name__ == "__main__":
    main(CTD_INPUT_FILE, OUTPUT_FILE)
import requests
import pandas as pd
import time

CTD_INPUT_FILE = 'drug_output/ctd_compounds.tsv' # previous output from 'retrieve_cid.py'
OUTPUT_FILE = 'drug_output/ctd_compounds_substances.tsv'

# rate limit for pubchem
REQUESTS_PER_SECOND = 5

def get_sid_by_cas(cas_number):
    search_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/name/{cas_number}/sids/JSON"
    response = requests.get(search_url)
    response_json = response.json()
    if 'IdentifierList' in response_json:
        return response_json['IdentifierList']['SID']
    else:
        raise ValueError(f"Could not find SID for CAS number {cas_number}")

def get_sid_by_name(drug_name):
    search_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/name/{drug_name}/sids/JSON"
    response = requests.get(search_url)
    response_json = response.json()
    if 'IdentifierList' in response_json:
        return response_json['IdentifierList']['SID']
    else:
        raise ValueError(f"Could not find SID for drug name {drug_name}")

def get_substance_info(sid):
    info_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/substance/{sid}/JSON/"
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
        if pd.notna(row.get('Drug Name')):  # skip rows where 'Drug Name' is not null
            results.append(row.to_dict())   # append original row data
            continue

        cas_number = row['cas_no']
        chemical_name = row['chemical_name']
        try:
            if pd.notna(cas_number):
                sids = get_sid_by_cas(cas_number)
                if len(sids) > 1:
                    print(f"Multiple SIDs found for CAS number '{cas_number}', using the first one: {sids[0]}")
                sid = sids[0]
            else:
                raise ValueError(f"No CAS number provided, trying to search by drug name '{chemical_name}'")

            data = get_substance_info(sid)
            compound_info = extract_info(data)
            combined_info = {**row.to_dict(), **compound_info}  # merge original row data with new data
            results.append(combined_info)
            time.sleep(1 / REQUESTS_PER_SECOND)                 # rate limiting

        except Exception as e:
            print(f"Error processing CAS number {cas_number}: {e}, trying by drug name '{chemical_name}'")
            try:
                sids = get_sid_by_name(chemical_name)
                if len(sids) > 1:
                    print(f"Multiple SIDs found for drug name '{chemical_name}', using the first one: {sids[0]}")
                sid = sids[0]

                data = get_substance_info(sid)
                compound_info = extract_info(data)
                combined_info = {**row.to_dict(), **compound_info}  # merge original row data with new data
                results.append(combined_info)
                time.sleep(1 / REQUESTS_PER_SECOND)                 # rate limiting
            except Exception as e2:
                print(f"Error processing drug name '{chemical_name}': {e2}")
                results.append(row.to_dict())

    # convert results to df & save
    results_df = pd.DataFrame(results)
    results_df.to_csv(output_filename, sep='\t', index=False)
    print(f"Data saved to {output_filename}")

if __name__ == "__main__":
    main(CTD_INPUT_FILE, OUTPUT_FILE)

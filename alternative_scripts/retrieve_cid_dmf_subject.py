import requests
import pandas as pd
import time

# ================================================
# UPDATE THESE VARIABLES TO MATCH YOUR FILE PATHS
# ================================================
INPUT_FILE = 'drug_input/dmf_drugs.tsv'
OUTPUT_FILE = 'drug_output/dmf_drugs_output.tsv'

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

def extract_info(data, cid):
    compound_info = {}
    compound_info['generic_name'] = data['Record'].get('RecordTitle', None)
    compound_info['url'] = f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}"
    sections = data['Record'].get('Section', [])

    for section in sections:
        if section['TOCHeading'] == 'Names and Identifiers':
            for subsection in section.get('Section', []):
                if subsection['TOCHeading'] == 'Other Identifiers':
                    for subsubsection in subsection.get('Section', []):
                        if subsubsection['TOCHeading'] == 'CAS':
                            compound_info['cas_no'] = subsubsection['Information'][0]['Value']['StringWithMarkup'][0]['String']
                        if subsubsection['TOCHeading'] == 'UNII':
                            compound_info['unii_no'] = subsubsection['Information'][0]['Value']['StringWithMarkup'][0]['String']

        if section['TOCHeading'] == 'Drug and Medication Information':
            for subsection in section.get('Section', []):
                if subsection['TOCHeading'] == 'Drug Indication':
                    compound_info['indication'] = subsection['Information'][0]['Value']['StringWithMarkup'][0]['String']
        
        if section['TOCHeading'] == 'Pharmacology and Biochemistry':
            for subsection in section.get('Section', []):
                if subsection['TOCHeading'] == 'MeSH Pharmacological Classification':
                    compound_info['pharm_class'] = subsection['Information'][0]['Name']
                    compound_info['pharm_desc'] = subsection['Information'][0]['Value']['StringWithMarkup'][0]['String']

        if section['TOCHeading'] == 'Safety and Hazards':
            for subsection in section.get('Section', []):
                if subsection['TOCHeading'] == 'Hazards Identification':
                    for subsubsection in subsection.get('Section', []):
                        if subsubsection['TOCHeading'] == 'GHS Classification':
                            hazard_statements = []
                            for info in subsubsection.get('Information', []):
                                if info['Name'] == 'GHS Hazard Statements':
                                    for string_with_markup in info['Value']['StringWithMarkup']:
                                        hazard_statements.append(string_with_markup['String'])
                            compound_info['hazard_statements'] = '; '.join(hazard_statements)
    return compound_info

def main(input_filename, output_filename):
    drug_data = pd.read_csv(input_filename, sep='\t')
    results = []

    for index, row in drug_data.iterrows():
        cas_number = row.get('cas_no', None)
        dmf_subject = row.get('dmf_subject', None)
        try:
            if pd.notna(cas_number):
                cids = get_cid_by_cas(cas_number)
                if len(cids) > 1:
                    print(f"======================================================\n")
                    print(f"Multiple CIDs found for CAS {cas_number}.\n")
                cid = cids[0]
            elif pd.notna(dmf_subject):
                cids = get_cid_by_name(dmf_subject)
                if len(cids) > 1:
                    print(f"======================================================\n")
                    print(f"Multiple CIDs found for DMF subject '{dmf_subject}'.\n")
                cid = cids[0]
            else:
                raise ValueError(f"Missing CAS and DMF subject for drug in row {index}")

            print(f"Processing {cas_number} | {dmf_subject} | {cid}")
            data = get_compound_info(cid)
            compound_info = extract_info(data, cid)
            for key, value in compound_info.items():
                if key in row and pd.isna(row[key]):
                    row[key] = value
                elif key not in row:
                    row[key] = value
            time.sleep(1 / REQUESTS_PER_SECOND) # rate limiting
            print(f"Added data for {dmf_subject} or {cas_number}\n")
            print(f"======================================================\n")

        except Exception as e:
            print(f"Error processing {dmf_subject} or {cas_number}: {e}")
        results.append(row)

    # convert results to df & save
    results_df = pd.DataFrame(results)
    results_df.to_csv(output_filename, sep='\t', index=False)
    print(f"======================================================\n")
    print(f"Data saved to {output_filename}.\n")
    print(f"======================================================\n")

if __name__ == "__main__":
    main(INPUT_FILE, OUTPUT_FILE)
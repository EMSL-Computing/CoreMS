
from pathlib import Path

import pandas as pd

def merge_files(file_paths: list, output_filename:str, variable: str = 'Peak Height'):

        master_data_dict = []
        list_filenames = []
        for filepath in file_paths:
            
            filepath = Path(filepath)
            
            with filepath.open('r') as f:
                
                #data = json.loads(json.load(f))
                
                df = pd.read_csv(f)
                idx = df.groupby(['Molecular Formula'])['Confidence Score'].transform(max) == df['Confidence Score']

                df = df[idx]
                df.fillna(0, inplace=True)

                name_column = "{} ({})".format(variable, filepath.stem)
                
                df.rename({variable: name_column}, inplace=True, axis=1)

                list_filenames.append(name_column)
                master_data_dict.extend(df.to_dict('records'))

        formula_dict = {}
        for record in master_data_dict:
            molecular_formula = record.get('Molecular Formula')
            
            if molecular_formula in formula_dict.keys():
                formula_dict[molecular_formula].append(record)
            else:
                formula_dict[molecular_formula] = [record]
        
        def dict_mean(dict_list, average_keys):
            mean_dict = {}
            
            for key in average_keys:
                
                mean_dict[key] = sum(d[key] for d in dict_list) / len(dict_list)
            
            return mean_dict
            
        average_records = []
        
        average_keys = ['m/z', 'Calibrated m/z', 'Calculated m/z', 'Peak Area', 'Resolving Power', 'S/N', 'm/z Error (ppm)', 'm/z Error Score', 
                        'Isotopologue Similarity', 'Mono Isotopic Index', 'Confidence Score']
        average_keys.extend(list_filenames)

        for formula, records in formula_dict.items():
            
            #mean_dict = dict_mean(records, average_keys)
            mean_dict = {}
            for record in  records: 
                #get the selected variable
                for filename in list_filenames:
                    if filename in record.keys():
                        mean_dict[filename] = record[filename]
                
            for record in  records: 
                #than get the rest of the data
                for key in record.keys():
                    if key not in average_keys:
                        mean_dict[key] = record[key]
        
            average_records.append(mean_dict)                
        master_df = pd.DataFrame(average_records)
        
        master_df.set_index('Molecular Formula', inplace=True)
        print(master_df)

        master_df.to_csv('{}.csv'.format(output_filename))
        #grouped = master_df.groupby(["Molecular Formula", "Sample Name", "Peak Height"])

if __name__ == '__main__':
    
    file_paths = ['tests/tests_data/Auto_SRFA_QC.csv', 'tests/tests_data/Auto_SRFA_QC II.csv']
    output_path = 'test_aggregation'
    merge_files(file_paths, output_path)
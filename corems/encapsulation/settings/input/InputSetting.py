from corems.encapsulation.constant import Labels

class DataInputSetting:
        
        #add to this dict the VALUES to match your labels, THE ORDER WON"T MATTER
        #"column_translate" : {"m/z":"m/z", "Resolving Power":"Resolving Power", "Abundance":"Abundance" , "S/N":"S/N"}
        header_translate = {'m/z': Labels.mz, 
                            "Resolving Power":"Resolving Power",
                            "Res.":Labels.rp, 
                            'I':Labels.abundance,
                            "Abundance":"Abundance",
                            "Signal/Noise":"S/N",
                            "S/N":"S/N"}

def d_params(file_location): #pragma: no cover

        parameters = dict()

        parameters["Aterm"] = 0

        parameters["Bterm"] = 0

        parameters["Cterm"] = 0

        parameters["exc_high_freq"] = 0

        parameters["exc_low_freq"] = 0

        parameters["bandwidth"] = 0

        parameters['analyzer'] = 'Unknown'
        
        parameters['instrument_label'] = 'Unknown' 

        parameters['sample_name'] = 'Unknown'

        parameters["number_data_points"] = 0

        parameters["polarity"] = 'Unknown'

        parameters["filename_path"] = file_location

        """scan_number and rt will be need to lc ms"""

        parameters["mobility_scan"] = 0

        parameters["mobility_rt"] = 0

        parameters["scan_number"] = 0

        parameters["rt"] = 0

        return parameters
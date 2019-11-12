class DataInputSetting:
        
        #add to this dict the VALUES to match your labels, THE ORDER WON"T MATTER
        #"column_translate" : {"m/z":"m/z", "Resolving Power":"Resolving Power", "Abundance":"Abundance" , "S/N":"S/N"}
        header_translate = {"m/z":"m/z", 
                            "Resolving Power":"Resolving Power",
                            "Res.":"Resolving Power", 
                            'I':"Abundance",
                            "Abundance":"Abundance",
                            "Signal/Noise":"S/N",
                            "S/N":"S/N"}

def d_parms(file_location):

        parameters = dict()

        parameters["Aterm"] = None

        parameters["Bterm"] = None

        parameters["Cterm"] = None

        parameters["exc_high_freq"] = None

        parameters["exc_low_freq"] = None

        parameters["bandwidth"] = None

        parameters['analyzer'] = None
        
        parameters['instrument_label'] = None 

        parameters['sample_name'] = None

        parameters["number_data_points"] = 0

        parameters["polarity"] = None

        parameters["filename_path"] = file_location

        """scan_number and rt will be need to lc ms"""

        parameters["mobility_scan"] = 0

        parameters["mobility_rt"] = 0

        parameters["scan_number"] = 0

        parameters["rt"] = 0

        return parameters
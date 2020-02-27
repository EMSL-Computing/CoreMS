from corems.molecular_id.calc.ClusterFilter import ClusteringFilter

class MolecularFormulaSearchFilters:

    @staticmethod
    def filter_kendrick( ms_peak_indexes, mass_spectrum_obj):

        index_to_remove = []

        if mass_spectrum_obj.molecular_search_settings.use_runtime_kendrick_filter:
            
            index_to_remove = ClusteringFilter().filter_kendrick_by_index(ms_peak_indexes, mass_spectrum_obj)

            #for index in noise_indexes: self.mass_spectrum_obj[index].clear_molecular_formulas()
        for peak_index, mf_obj in index_to_remove:
                #print(peak_index, mf_obj)
                for iso_index, mf_iso in mf_obj.mspeak_mf_isotopologues_indexes:
                    mass_spectrum_obj[iso_index].remove_molecular_formula(mf_iso)    

                mass_spectrum_obj[peak_index].remove_molecular_formula(mf_obj)
                
                ms_peak_indexes.remove((peak_index, mf_obj))

        return ms_peak_indexes

    @staticmethod
    def check_min_peaks( ms_peak_indexes, mass_spectrum):
        
        if mass_spectrum.molecular_search_settings.use_min_peaks_filter:

            if not len(ms_peak_indexes) >= mass_spectrum.molecular_search_settings.min_peaks_per_class:
                
                for peak_index, mf_obj in ms_peak_indexes:
                
                    mass_spectrum[peak_index].remove_molecular_formula(mf_obj)

    @staticmethod
    def filter_isotopologue( ms_peak_indexes, mass_spectrum):
        
        isotopologue_count_threshold = mass_spectrum.molecular_search_settings.isotopologue_filter_threshold
        
        index_to_remove = []

        if mass_spectrum.molecular_search_settings.use_isotopologue_filter:

            for mspeak_index, mf_obj in ms_peak_indexes:
            
                if mf_obj.isotopologue_count_percentile < isotopologue_count_threshold:
                    
                    #need to make sure only the isotopologue is being removed before going forward,
                    #need to modify indexes storage to include molecular formula position inside mspeak
                                                    
                    #removes tuple obj from initial list to be used on next filter steps
                    ms_peak_indexes.remove((mspeak_index, mf_obj))
                    index_to_remove.append((mspeak_index, mf_obj))
                    index_to_remove.extend(mf_obj.mspeak_mf_isotopologues_indexes)

        #iterate over all indexes to be remove and remove the mf from the mspeak 
        
        for peak_index, mf_obj in index_to_remove:
                #print(peak_index, mf_obj)
                mass_spectrum[peak_index].remove_molecular_formula(mf_obj)
                
        return ms_peak_indexes 
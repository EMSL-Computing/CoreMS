import sys
import os

from corems.mass_spectrum.output.export import HighResMassSpecExport

def test_export_mass_spectrum(mass_spectrum_ftms):

    exportMS = HighResMassSpecExport('NEG_ESI_SRFA_CoreMS', mass_spectrum_ftms)

    exportMS._output_type = 'excel'
    exportMS.save()
    # Check that the file was created and remove it
    assert os.path.exists('NEG_ESI_SRFA_CoreMS.xlsx')
    os.remove('NEG_ESI_SRFA_CoreMS.xlsx')
    os.remove('NEG_ESI_SRFA_CoreMS.json')

    exportMS._output_type = 'csv'
    exportMS.save()
    assert os.path.exists('NEG_ESI_SRFA_CoreMS.csv')
    os.remove('NEG_ESI_SRFA_CoreMS.csv')
    os.remove('NEG_ESI_SRFA_CoreMS.json')

    exportMS._output_type = 'pandas'
    exportMS.save()
    assert os.path.exists('NEG_ESI_SRFA_CoreMS.pkl')
    os.remove('NEG_ESI_SRFA_CoreMS.pkl')
    os.remove('NEG_ESI_SRFA_CoreMS.json')

    exportMS._output_type = 'hdf5'
    exportMS.save()
    assert os.path.exists('NEG_ESI_SRFA_CoreMS.hdf5')
    os.remove('NEG_ESI_SRFA_CoreMS.hdf5')

    df = exportMS.get_pandas_df()    
    assert df.shape[0] > 10
    json_dump1 = exportMS.to_json() 

    mass_spectrum_ftms.to_excel('NEG_ESI_SRFA_CoreMS')
    assert os.path.exists('NEG_ESI_SRFA_CoreMS.xlsx')
    os.remove('NEG_ESI_SRFA_CoreMS.xlsx')
    os.remove('NEG_ESI_SRFA_CoreMS.json')
    df = mass_spectrum_ftms.to_dataframe()
    assert df.shape[0] > 10
    
    json_dump = mass_spectrum_ftms.to_json()
    assert len(json_dump) > 10
    assert json_dump1 == json_dump

    mass_spectrum_ftms.molecular_search_settings.output_score_method = "prob_score"
    mass_spectrum_ftms.to_csv('NEG_ESI_SRFA_CoreMS_prob_score')
    assert os.path.exists('NEG_ESI_SRFA_CoreMS_prob_score.csv')
    os.remove('NEG_ESI_SRFA_CoreMS_prob_score.csv')
    os.remove('NEG_ESI_SRFA_CoreMS_prob_score.json')

    mass_spectrum_ftms.to_pandas('NEG_ESI_SRFA_CoreMS')
    assert os.path.exists('NEG_ESI_SRFA_CoreMS.pkl')
    os.remove('NEG_ESI_SRFA_CoreMS.pkl')
    os.remove('NEG_ESI_SRFA_CoreMS.json')
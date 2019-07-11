
import os
import sys
sys.path.append(".")
from enviroms.emsl.yec.mass_spectrum.input.TextMassList import Read_MassList

__author__ = "Yuri E. Corilo"
__date__ = "Jul 02, 2019"


if __name__ == '__main__':

    directory = os.path.join(os.getcwd(), "res/")
    #filename_ESFA = "20190616_WK_ESFA_0pt2mgml_ESI_Neg_0pt8sFID_000001.d"
    #filename_SRFA = "20190205_WPK_SRFA_opt_000001.d"
    print(os.getcwd())
    file_name_ESFA_txt = directory + \
        "20190616_WK_ESFA_0pt2mgml_ESI_Neg_1pt4sFID_000001.ascii"

    full_file_path = file_name_ESFA_txt

    polariy = -1

    mass_spec = Read_MassList(full_file_path, polariy, delimiter="  ")

    mass_spec.plot_mz_domain_profile()

    print(mass_spec.number_average_molecular_weight())
    print(mass_spec.weight_average_molecular_weight())

    filtered_mass_peaks = mass_spec.filter_by_s2n(100)

    mass_spec_peaks_filtered = filtered_mass_peaks[0]

    '''
    print('Exp. Mass :',
          mass_spec_peaks_filtered.exp_mz, 
          '\nTheor. Mass :',
          mass_spec_peaks_filtered.molecular_formula.theoretical_mz,
          '\nMol. Formula :',
          mass_spec_peaks_filtered.molecular_formula.to_string(),
          '\nDBE :',
          mass_spec_peaks_filtered.molecular_formula.dbe, 
          '\nH/C :',
          mass_spec_peaks_filtered.molecular_formula.H_C,
          '\nClass :',
          mass_spec_peaks_filtered.molecular_formula.heteroatomic_class_label,
          '\nMass error :',
          mass_spec_peaks_filtered.molecular_formula.assigment_mass_error
          )
    '''

import os, sys
import pathlib

sys.path.append(".")
from enviroms.emsl.yec.mass_spectrum.input.TextMassList import Read_MassList

__author__ = "Yuri E. Corilo"
__date__ = "Jul 02, 2019"


if __name__ == "__main__":

    directory = os.path.join(os.getcwd(), "res/")
    # filename_ESFA = "20190616_WK_ESFA_0pt2mgml_ESI_Neg_0pt8sFID_000001.d"
    # filename_SRFA = "20190205_WPK_SRFA_opt_000001.d"

    file_name_ESFA_txt = (
        directory + "20190616_WK_ESFA_0pt2mgml_ESI_Neg_1pt4sFID_000001.ascii"
    )

    full_file_path = file_name_ESFA_txt

    polariy = -1

    mass_list_reader = Read_MassList(full_file_path, polariy, delimiter="  ")

    mass_spectrum = mass_list_reader.get_mass_spectrum(auto_process=True)

    mass_spectrum.plot_mz_domain_profile()

    print(
        "number_average_molecular_weight",
        mass_spectrum.number_average_molecular_weight(),
    )
    print(
        "weight_average_molecular_weight",
        mass_spectrum.weight_average_molecular_weight(),
    )

    filtered_mass_peaks = mass_spectrum.filter_by_s2n(100)

    filtered_mass_spec = filtered_mass_peaks[0]

    """
    after assigment
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
    """

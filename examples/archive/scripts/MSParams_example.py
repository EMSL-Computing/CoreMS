# This script demonstrates how to set global parameters and instantiate a mass spectrum object using them

from corems.encapsulation.factory.parameters import MSParameters
from corems.transient.input.brukerSolarix import ReadBrukerSolarix

# Set global parameters and instantiate a mass spectrum object using them
## Note that the default noise_threshold_method is 'log'
MSParameters.mass_spectrum.noise_threshold_method = 'relative_abundance'

parser = ReadBrukerSolarix("tests/tests_data/ftms/ESI_NEG_SRFA.d")
bruker_transient = parser.get_transient()
mass_spectrum_i = bruker_transient.get_mass_spectrum(plot_result=False, auto_process=False)
print(mass_spectrum_i.parameters.mass_spectrum.noise_threshold_method) # relative_abundance

# Create a new MSParameters instance with default parameters and assign it to the mass spectrum object
new_msparams = MSParameters(use_defaults=True)
mass_spectrum_i.parameters = new_msparams
print(mass_spectrum_i.parameters.mass_spectrum.noise_threshold_method) # log




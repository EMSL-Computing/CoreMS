import sys

from corems.ms_peak.factory.MSPeakClasses import ICRMassPeak


def test_mspeaks_fit(mass_spectrum_ftms):
    mass_spectrum_ftms[3].plot_simulation()
    mass_spectrum_ftms[3].plot_simulation(sim_type="gaussian", oversample_multiplier=10)
    mass_spectrum_ftms[3].plot_simulation()


def test_mspeak_calculations():
    kendrick_base = {"C": 1, "H": 2}
    ion_charge = +1
    mz_exp = 212.1234
    abundance = 200
    resolving_power = 1000000
    signal_to_noise = 200
    massspec_index = (300, 300, 300)
    index = 1
    mspeak = ICRMassPeak(
        ion_charge,
        mz_exp,
        abundance,
        resolving_power,
        signal_to_noise,
        massspec_index,
        index,
    )
    assert mspeak.resolving_power == 1000000

    mspeak.change_kendrick_base(kendrick_base)

    mspeak._calc_kmd(kendrick_base)
    mspeak.calc_area()

    assert round(mspeak.kendrick_mass, 3) == 211.887
    assert round(mspeak.kmd * 100, 0) == -89
    assert mspeak.knm == 211

    mspeak.set_calc_resolving_power(50, 3)
    assert round(mspeak.resolving_power, 0) == 9008907

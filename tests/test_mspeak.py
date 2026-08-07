import sys

from corems.ms_peak.factory.MSPeakClasses import ICRMassPeak


def test_mspeaks_fit(mass_spectrum_ftms):
    mass_spectrum_ftms[3].plot_simulation()
    mass_spectrum_ftms[3].plot_simulation(sim_type="gaussian", oversample_multiplier=10)
    mass_spectrum_ftms[3].plot_simulation()


def test_mspeak_calculations():
    kendrick_base = {"C": 1, "H": 2}
    polarity = +1
    mz_exp = 212.1234
    abundance = 200
    resolving_power = 1000000
    signal_to_noise = 200
    massspec_index = (300, 300, 300)
    index = 1
    mspeak = ICRMassPeak(
        polarity,
        mz_exp,
        abundance,
        resolving_power,
        signal_to_noise,
        massspec_index,
        index,
    )
    assert mspeak.resolving_power == 1000000
    assert mspeak.polarity == 1
    # Deprecated alias remains for backwards compatibility
    assert mspeak.ion_charge == mspeak.polarity

    mspeak.change_kendrick_base(kendrick_base)

    mspeak._calc_kmd(kendrick_base)
    mspeak.calc_area()

    assert round(mspeak.kendrick_mass, 3) == 211.887
    assert round(mspeak.kmd * 100, 0) == -89
    assert mspeak.knm == 211

    mspeak.set_calc_resolving_power(50, 3)
    assert round(mspeak.resolving_power, 0) == 9008907


def test_mspeak_polarity_alias_and_setter():
    """peak.ion_charge is a deprecated alias for peak.polarity."""
    mspeak = ICRMassPeak(+1, 100.0, 10.0, 1e5, 50.0, (0, 0, 0), 0)
    assert mspeak.polarity == 1
    assert mspeak.ion_charge == 1
    mspeak.ion_charge = -1
    assert mspeak.polarity == -1
    assert mspeak.ion_charge == -1


def test_add_mspeak_rejects_ion_charge_keyword(mass_spectrum_ftms):
    """add_mspeak must fail loudly if ion_charge= is passed."""
    import pytest

    with pytest.raises(TypeError, match="no longer accepts ion_charge"):
        mass_spectrum_ftms.add_mspeak(
            polarity=1,
            mz_exp=100.0,
            abundance=1.0,
            resolving_power=1e5,
            signal_to_noise=10.0,
            massspec_indexes=(0, 0, 0),
            ion_charge=1,
        )

    with pytest.raises(TypeError, match="no longer accepts ion_charge"):
        mass_spectrum_ftms.add_mspeak(
            1,
            100.0,
            1.0,
            1e5,
            10.0,
            (0, 0, 0),
            ion_charge=1,
        )

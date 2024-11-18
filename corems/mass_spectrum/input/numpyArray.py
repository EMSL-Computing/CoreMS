__author__ = "Yuri E. Corilo"
__date__ = "Oct 23, 2019"

from corems.encapsulation.constant import Labels
from corems.encapsulation.factory.parameters import default_parameters
from corems.mass_spectrum.factory.MassSpectrumClasses import (
    MassSpecCentroid,
    MassSpecProfile,
)


def ms_from_array_profile(
    mz,
    abundance,
    dataname: str,
    polarity: int = -1,
    auto_process: bool = True,
    data_type: str = Labels.simulated_profile,
):
    """Create a MassSpecProfile object from an array of m/z values and abundance values.

    Parameters
    ----------
    mz : numpy.ndarray
        Array of m/z values.
    abundance : numpy.ndarray
        Array of abundance values.
    dataname : str
        Name of the data.
    polarity : int, optional
        Polarity of the data. The default is -1.
    auto_process : bool, optional
        Flag to automatically process the data. The default is True.
    data_type : str, optional
        Type of the data. The default is Labels.simulated_profile.

    Returns
    -------
    MassSpecProfile
        The created MassSpecProfile object.
    """
    data_dict = {Labels.mz: mz, Labels.abundance: abundance}

    output_parameters = get_output_parameters(polarity, dataname)

    output_parameters[Labels.label] = data_type

    ms = MassSpecProfile(data_dict, output_parameters, auto_process=auto_process)

    return ms


def ms_from_array_centroid(
    mz,
    abundance,
    rp: list[float],
    s2n: list[float],
    dataname: str,
    polarity: int = -1,
    auto_process: bool = True,
):
    """Create a MassSpecCentroid object from an array of m/z values, abundance values, resolution power, and signal-to-noise ratio.

    Parameters
    ----------
    mz : numpy.ndarray
        Array of m/z values.
    abundance : numpy.ndarray
        Array of abundance values.
    rp : list(float)
        List of resolving power values.
    s2n : list(float)
        List of signal-to-noise ratio values.
    dataname : str
        Name of the data.
    polarity : int, optional
        Polarity of the data. The default is -1.
    auto_process : bool, optional

    Returns
    -------
    MassSpecCentroid
        The created MassSpecCentroid object.
    """
    data_dict = {
        Labels.mz: mz,
        Labels.abundance: abundance,
        Labels.s2n: s2n,
        Labels.rp: rp,
    }

    output_parameters = get_output_parameters(polarity, dataname)
    output_parameters[Labels.label] = Labels.corems_centroid

    return MassSpecCentroid(data_dict, output_parameters, auto_process)


def get_output_parameters(polarity: int, file_location: str):
    """Generate the output parameters for creating a MassSpecProfile or MassSpecCentroid object.

    Parameters
    ----------
    polarity : int
        Polarity of the data.
    file_location : str
        File location.

    Returns
    -------
    dict
        Output parameters.
    """
    d_params = default_parameters(file_location)

    d_params["analyzer"] = "Generic Simulated"

    d_params["instrument_label"] = "Generic Simulated"

    d_params["polarity"] = polarity

    d_params["filename_path"] = file_location

    d_params["mobility_scan"] = 0

    d_params["mobility_rt"] = 0

    d_params["scan_number"] = 0

    d_params["rt"] = 0

    d_params[Labels.label] = Labels.simulated_profile

    return d_params

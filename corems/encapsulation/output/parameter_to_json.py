import json

import toml
from pathlib import Path

from corems.encapsulation.output import parameter_to_dict
from corems.encapsulation.output.parameter_to_dict import get_dict_data_lcms


def dump_all_settings_json(filename="SettingsCoreMS.json", file_path=None):
    """
    Write JSON file into current directory with all the default settings for the CoreMS package.

    Parameters:
    ----------
    filename : str, optional
        The name of the JSON file to be created. Default is 'SettingsCoreMS.json'.
    file_path : str or Path, optional
        The path where the JSON file will be saved. If not provided, the file will be saved in the current working directory.
    """

    data_dict_all = parameter_to_dict.get_dict_all_default_data()

    if not file_path:
        file_path = Path.cwd() / filename

    with open(
        file_path,
        "w",
        encoding="utf8",
    ) as outfile:
        import re

        # pretty print
        output = json.dumps(
            data_dict_all, sort_keys=False, indent=4, separators=(",", ": ")
        )
        output = re.sub(r'",\s+', '", ', output)

        outfile.write(output)


def dump_ms_settings_json(filename="SettingsCoreMS.json", file_path=None):
    """
    Write JSON file into current directory with all the mass spectrum default settings for the CoreMS package.

    Parameters
    ----------
    filename : str, optional
        The name of the JSON file to be created. Default is 'SettingsCoreMS.json'.
    file_path : str or Path, optional
        The path where the JSON file will be saved. If not provided, the file will be saved in the current working directory.

    """
    data_dict = parameter_to_dict.get_dict_ms_default_data()
    if not file_path:
        file_path = Path.cwd() / filename

    with open(
        file_path,
        "w",
        encoding="utf8",
    ) as outfile:
        import re

        # pretty print
        output = json.dumps(
            data_dict, sort_keys=False, indent=4, separators=(",", ": ")
        )
        output = re.sub(r'",\s+', '", ', output)

        outfile.write(output)


def dump_gcms_settings_json(filename="SettingsCoreMS.json", file_path=None):
    """
    Write JSON file into current directory containing the default GCMS settings data.

    Parameters
    ----------
    filename : str, optional
        The name of the JSON file to be created. Default is 'SettingsCoreMS.json'.
    file_path : str or Path-like object, optional
        The path where the JSON file will be saved. If not provided, the file will be saved in the current working directory.
    """

    from pathlib import Path
    import json

    data_dict = parameter_to_dict.get_dict_gcms_default_data()

    if not file_path:
        file_path = Path.cwd() / filename

    with open(
        file_path,
        "w",
        encoding="utf8",
    ) as outfile:
        import re

        # pretty print
        output = json.dumps(
            data_dict, sort_keys=False, indent=4, separators=(",", ": ")
        )
        output = re.sub(r'",\s+', '", ', output)

        outfile.write(output)


def dump_all_settings_toml(filename="SettingsCoreMS.toml", file_path=None):
    """
    Write TOML file into the specified file path or the current directory with all the default settings for the CoreMS package.

    Parameters
    ----------
    filename : str, optional
        The name of the TOML file. Defaults to 'SettingsCoreMS.toml'.
    file_path : str or Path, optional
        The path where the TOML file will be saved. If not provided, the file will be saved in the current directory.

    """
    from pathlib import Path

    data_dict_all = parameter_to_dict.get_dict_all_default_data()

    if not file_path:
        file_path = Path.cwd() / filename

    with open(
        file_path,
        "w",
        encoding="utf8",
    ) as outfile:
        import re

        output = toml.dumps(data_dict_all)
        outfile.write(output)


def dump_ms_settings_toml(filename="SettingsCoreMS.toml", file_path=None):
    """
    Write TOML file into the current directory with all the mass spectrum default settings for the CoreMS package.

    Parameters
    ----------
    filename : str, optional
        The name of the TOML file to be created. Default is 'SettingsCoreMS.toml'.
    file_path : str or Path, optional
        The path where the TOML file should be saved. If not provided, the file will be saved in the current working directory.

    """
    data_dict = parameter_to_dict.get_dict_ms_default_data()

    if not file_path:
        file_path = Path.cwd() / filename

    with open(
        file_path,
        "w",
        encoding="utf8",
    ) as outfile:
        import re

        # pretty print
        output = toml.dumps(data_dict)
        outfile.write(output)


def dump_gcms_settings_toml(filename="SettingsCoreMS.toml", file_path=None):
    """
    Write TOML file into current directory containing the default GCMS settings data.

    Parameters
    ----------
    filename : str, optional
        The name of the TOML file. Defaults to 'SettingsCoreMS.toml'.
    file_path : str or Path, optional
        The path where the TOML file will be saved. If not provided, the file will be saved in the current working directory.

    """

    data_dict = parameter_to_dict.get_dict_gcms_default_data()

    if not file_path:
        file_path = Path.cwd() / filename

    with open(
        file_path,
        "w",
        encoding="utf8",
    ) as outfile:
        output = toml.dumps(data_dict)
        outfile.write(output)


def dump_lcms_settings_json(
    filename="SettingsCoreMS.json", file_path=None, lcms_obj=None
):
    """
    Write JSON file into current directory with all the LCMS settings data for the CoreMS package.

    Parameters
    ----------
    filename : str, optional
        The name of the JSON file. Defaults to 'SettingsCoreMS.json'.
    file_path : str or Path, optional
        The path where the JSON file will be saved. If not provided, the file will be saved in the current working directory.
    lcms_obj : object, optional
        The LCMS object containing the settings data. If not provided, the settings data will be retrieved from the default settings.

    """

    if lcms_obj is None:
        data_dict = parameter_to_dict.get_dict_lcms_default_data()
    else:
        data_dict = get_dict_data_lcms(lcms_obj)

    if not file_path:
        file_path = Path.cwd() / filename

    with open(
        file_path,
        "w",
        encoding="utf8",
    ) as outfile:
        outfile.write(json.dumps(data_dict, indent=4))


def dump_lcms_settings_toml(
    filename="SettingsCoreMS.toml", file_path=None, lcms_obj=None
):
    """
    Write TOML file into current directory with all the LCMS settings data for the CoreMS package.

    Parameters
    ----------
    filename : str, optional
        The name of the TOML file. Defaults to 'SettingsCoreMS.toml'.
    file_path : str or Path, optional
        The path where the TOML file will be saved. If not provided, the file will be saved in the current working directory.
    lcms_obj : object, optional
        The LCMS object containing the settings data. If not provided, the settings data will be retrieved from the default settings.

    """

    if lcms_obj is None:
        data_dict = parameter_to_dict.get_dict_lcms_default_data()
    else:
        data_dict = get_dict_data_lcms(lcms_obj)

    if not file_path:
        file_path = Path.cwd() / filename

    with open(
        file_path,
        "w",
        encoding="utf8",
    ) as outfile:
        output = toml.dumps(data_dict)
        outfile.write(output)

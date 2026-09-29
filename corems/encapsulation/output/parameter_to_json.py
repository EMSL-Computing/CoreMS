import json

import tomlkit
from pathlib import Path

from corems.encapsulation.output import parameter_to_dict
from corems.encapsulation.output.parameter_to_dict import get_dict_data_lcms


def _toml_array_item(item):
    """Prepare one array element for TOML.

    The previous encoder wrote a null array element as the string ``None``.
    """
    if item is None:
        return "None"
    return _toml_value(item)


def _toml_value(value):
    """Prepare ``value`` for TOML.

    TOML has no null. Table entries whose value is ``None`` are omitted, which
    matches the encoder previously used for settings files.
    """
    if isinstance(value, dict):
        return {
            key: _toml_value(item)
            for key, item in value.items()
            if item is not None
        }
    if isinstance(value, tuple):
        return tuple(_toml_array_item(item) for item in value)
    if isinstance(value, list):
        return [_toml_array_item(item) for item in value]
    return value


def _dumps_toml(data):
    """Serialize ``data`` to a TOML string."""
    return tomlkit.dumps(_toml_value(data))


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

        output = _dumps_toml(data_dict_all)
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
        output = _dumps_toml(data_dict)
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
        output = _dumps_toml(data_dict)
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
        output = _dumps_toml(data_dict)
        outfile.write(output)


def dump_lcms_collection_settings_json(
    filename="SettingsCoreMS.json", file_path=None, lcms_collection=None
):
    """Write JSON file with LCMS collection settings data.

    Parameters
    ----------
    filename : str, optional
        The name of the JSON file. Defaults to 'SettingsCoreMS.json'.
    file_path : str or Path, optional
        The path where the JSON file will be saved. If not provided, the file will be saved in the current working directory.
    lcms_collection : LCMSCollection, optional
        The LCMS collection object containing the settings data. If not provided, the settings data will be retrieved from the default settings.
    """
    from corems.encapsulation.output.parameter_to_dict import (
        get_dict_data_lcms_collection,
        get_dict_lcms_collection_default_data,
    )

    if lcms_collection is None:
        data_dict = get_dict_lcms_collection_default_data()
    else:
        data_dict = get_dict_data_lcms_collection(lcms_collection)

    if not file_path:
        file_path = Path.cwd() / filename

    with open(
        file_path,
        "w",
        encoding="utf8",
    ) as outfile:
        outfile.write(json.dumps(data_dict, indent=4))


def dump_lcms_collection_settings_toml(
    filename="SettingsCoreMS.toml", file_path=None, lcms_collection=None
):
    """Write TOML file with LCMS collection settings data.

    Parameters
    ----------
    filename : str, optional
        The name of the TOML file. Defaults to 'SettingsCoreMS.toml'.
    file_path : str or Path, optional
        The path where the TOML file will be saved. If not provided, the file will be saved in the current working directory.
    lcms_collection : LCMSCollection, optional
        The LCMS collection object containing the settings data. If not provided, the settings data will be retrieved from the default settings.
    """
    from corems.encapsulation.output.parameter_to_dict import (
        get_dict_data_lcms_collection,
        get_dict_lcms_collection_default_data,
    )

    if lcms_collection is None:
        data_dict = get_dict_lcms_collection_default_data()
    else:
        data_dict = get_dict_data_lcms_collection(lcms_collection)

    if not file_path:
        file_path = Path.cwd() / filename

    with open(
        file_path,
        "w",
        encoding="utf8",
    ) as outfile:
        output = _dumps_toml(data_dict)
        outfile.write(output)

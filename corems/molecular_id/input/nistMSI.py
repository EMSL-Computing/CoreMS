__author__ = "Yuri E. Corilo"
__date__ = "Feb 12, 2020"

from threading import Thread
from pathlib import Path

from corems.molecular_id.factory.EI_SQL import EI_LowRes_SQLite


class ReadNistMSI(Thread):
    """A class for reading NIST MSI files and storing the data in a SQLite database.

    Parameters
    ----------
    file_path : str
        The path to the NIST MSI file.
    url : str, optional
        The URL for the SQLite database. Default is 'sqlite://'.

    Raises
    ------
    FileExistsError
        If the specified file does not exist.

    Attributes
    ----------
    file_path : str
        The path to the NIST MSI file.
    url : str
        The URL for the SQLite database.
    sqlLite_obj : EI_LowRes_SQLite
        The SQLite object for storing the compound data.

    Methods
    -------
    * run().
        Runs the thread and initializes the SQLite object.
    * get_sqlLite_obj().
        Returns the SQLite object.
    * get_compound_data_dict_list().
        Parses the NIST MSI file and returns a list of compound data dictionaries.
    """

    def __init__(self, file_path, url="sqlite://"):
        Thread.__init__(self)
        file_path = Path(file_path)

        if not file_path.exists():
            raise FileExistsError("File does not exist: " + file_path)

        self.file_path = file_path

        self.url = url

    def run(self):
        """Runs the thread and initializes the SQLite object."""
        self.sqlLite_obj = self.get_sqlLite_obj()

    def get_sqlLite_obj(self):
        """Returns the SQLite object.

        Returns
        -------
        EI_LowRes_SQLite
            The SQLite object for storing the compound data.
        """
        compound_data_dict_list = self.get_compound_data_dict_list()

        sqlLite_obj = EI_LowRes_SQLite(url=self.url)

        for data_dict in compound_data_dict_list:
            if not data_dict.get("NUM PEAKS"):
                data_dict["NUM PEAKS"] = len(data_dict.get("mz"))
            if not data_dict.get("CASNO"):
                data_dict["CASNO"] = data_dict.get("CAS")
                if not data_dict["CASNO"]:
                    data_dict["CASNO"] = 0
            # print(data_dict)
            try:
                sqlLite_obj.add_compound(data_dict)
            except:
                print(data_dict.get("NAME"))

        return sqlLite_obj

    def get_compound_data_dict_list(self):
        """Parses the NIST MSI file and returns a list of compound data dictionaries.

        Returns
        -------
        list
            A list of compound data dictionaries.
        """
        list_dict_data = []

        with open(self.file_path) as msifile:
            content = msifile.readlines()

            i = 0

            dict_data = dict()
            dict_data["mz"] = list()
            dict_data["abundance"] = list()
            # for line in content:
            #   print(line, line=="\n" )

            while i < len(content):
                split_line = content[i].split(":")

                # empty line
                if len(content[i]) == 1:
                    i += 1
                    if dict_data.get("NAME"):
                        list_dict_data.append(dict_data)

                    # print(dict_data)
                    dict_data = dict()
                    dict_data["mz"] = list()
                    dict_data["abundance"] = list()

                # metadata, name, ri, rt etc
                elif len(split_line) >= 2:
                    label = split_line[0]
                    value = ":".join(split_line[1:]).strip("\n").strip("")
                    dict_data[label] = value
                    i += 1

                # mz and abundance pairs
                elif len(split_line) == 1:
                    for s in content[i].strip("\n").strip("").split("(")[1:]:
                        values = s.split(" ")

                        if values[0] == "":
                            mz = values[1]
                        else:
                            mz = values[0]

                        abun = values[-2].strip(")")

                        dict_data["mz"].append(mz)
                        dict_data["abundance"].append(abun)

                    i += 1
                # something else
                else:
                    i += 1

        return list_dict_data

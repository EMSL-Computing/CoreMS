"""Utility functions for Bruker data processing."""

from pathlib import Path
from s3path import S3Path


def get_scan_attributes(scan_attr, imaging_info_attr) -> dict:
    """
    Get the scan attributes from the scan.xml or ImagingInfo.xml file.
    
    Parameters
    ----------
    d_directory_location : str, Path, or S3Path
        Directory containing the XML files
        
    Returns
    -------
    dict
        Dictionary containing the scan number as key and a tuple of retention time, TIC, 
        and optionally maxpeak and spotname as values.


    TODO: We need to reformat the dictionary to actually include keys and values so it is self-descriptive.
    TODO: This will break the code, so a new version is needed. 
    TODO: Will need to make sure theres tests which capture this change.
    
    """
    from bs4 import BeautifulSoup
    
    scan_xml_exists = scan_attr.exists()
    imaging_info_exists = imaging_info_attr.exists()

    if scan_xml_exists:
        try:
            soup = BeautifulSoup(scan_attr.open(), "xml")
            list_rt = [float(rt.text) for rt in soup.find_all("minutes")]
            list_tic = [float(tic.text) for tic in soup.find_all("tic")]
            list_scan = [int(scan.text) for scan in soup.find_all("count")]
            
            # Check if maxpeak exists (more comprehensive version)
            # TODO: Enable this, but it could break code so a new version is needed
            enable_maxpeak = False
            if enable_maxpeak:
                maxpeak_elements = soup.find_all("maxpeak")
                if maxpeak_elements:
                    list_maxpeak = [float(maxpeak.text) for maxpeak in maxpeak_elements]
                    dict_scan_rt_tic = dict(zip(list_scan, zip(list_rt, list_tic, list_maxpeak)))
                else:
                    dict_scan_rt_tic = dict(zip(list_scan, zip(list_rt, list_tic)))

            dict_scan_rt_tic = dict(zip(list_scan, zip(list_rt, list_tic)))

            return dict_scan_rt_tic
        except Exception as e:
            raise FileNotFoundError(f"Error reading scan.xml: {e}")
    elif imaging_info_exists:
        try:
            soup = BeautifulSoup(imaging_info_attr.open(), "xml")
            list_rt = [float(rt.text) for rt in soup.find_all("minutes")]
            list_tic = [float(tic.text) for tic in soup.find_all("tic")]
            list_maxpeak = [float(maxpeak.text) for maxpeak in soup.find_all("maxpeak")]
            list_scan = [int(scan.find("count").text) for scan in soup.find_all("scan")]
            list_spotname = [
                scan.find("spotName").text for scan in soup.find_all("scan")
            ]
            dict_scan_rt_tic = dict(zip(list_scan, zip(list_rt, list_tic, list_maxpeak, list_spotname)))
            return dict_scan_rt_tic
        except Exception as e:
            raise FileNotFoundError(f"Error reading ImagingInfo.xml: {e}")
    else:
        raise FileNotFoundError(
            "Dataset does not contain a 'scan.xml' or 'ImagingInfo.xml' file."
        )
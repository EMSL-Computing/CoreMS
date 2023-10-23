import requests
import os
import json 
from pathlib import Path


username = os.environ.get('DMS_USER_NAME') or 'XXXX'
password = os.environ.get('DMS_USER_PASSWORD') or 'XXXXX'

s = requests.Session()
s.auth = (username, password)

url_request = "https://dms2.pnl.gov/data/ax/json/list_report/dataset/{}".format("4-AA_L-AsparticAcid") 
r = s.get(url_request)

data_list = r.json()
for data in data_list:
    print(data)


filein = Path("data/ftms_nom_data_products.json")
outfile = Path("data/ftms_nom_nmdc_data.csv")

with outfile.open("w") as csv_data:

    csv_data.write('{},{} \n'.format('DMS Dataset ID', "Dataset"))

    with filein.open() as json_data:

        list_data = json.load(json_data)
        for data in list_data:
            filename = (data.get("name").replace(".csv", ""))
            if filename[0:6] == 'Brodie':

                url_request = "https://dms2.pnl.gov/data/ax/json/list_report/dataset/{}".format(filename) 
                r = s.get(url_request)

                data_list = r.json()
                for data in data_list:
                    (data.get('ID'))

                csv_data.write('{},{}\n'.format(data.get('ID'), filename))


import hashlib
import json
from pathlib import Path

registration_path = Path("db/gcms_metabolomics_data_products.json")
metadata_path = Path("db/gcms_metabolomics_metadata_products.json")

metabolomics_data_path = Path(registration_path)
data_dir = Path("db/nmdc_data")

if metabolomics_data_path.exists():

    url_list = []
    with metabolomics_data_path.open('r') as metabolomics_data_products:
        products = json.load(metabolomics_data_products)
        for product in products:

            name = product.get("name")
            out_file_path = data_dir / name
            md5_checksum = hashlib.md5(out_file_path.open('rb').read()).hexdigest()
            product["md5_checksum"] = md5_checksum
            product["id"] = "{}:{}".format("nmdc", md5_checksum)
            product["file_size_bytes"] = out_file_path.stat().st_size
            url_list.append(product.get("url").replace('.csv', ".json").replace('results', 'metadata'))

    with registration_path.open('w') as metabolomics_data_products:
        metabolomics_data_products.write(json.dumps(products, indent=1))

    with metadata_path.open('w') as metabolomics_data_products:
        metabolomics_data_products.write(json.dumps(url_list, indent=1))

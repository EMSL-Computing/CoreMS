# TODO KRH: This is a work in progress. It is not yet functional but serves as the development testbed for the
# TODO KRH: LCMS metabolomics workflow.

from corems.molecular_id.search.database_interfaces import MSPInterface

msp_file_path = "/Users/heal742/LOCAL/05_NMDC/02_MetaMS/metams/data/databases/20250407_gnps_curated.msp"

my_msp = MSPInterface(file_path=msp_file_path)

my_msp_FE = my_msp._to_flashentropy()
print("finished")


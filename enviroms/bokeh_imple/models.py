from bokeh.core.properties import List, String, Dict, Int
from bokeh.models import LayoutDOM

class FileInput(LayoutDOM):
__implementation__ = 'static/js/extensions_file_input.coffee'
__javascript__ = './input_widget/static/js/papaparse.js'

value = String(help="""
Selected input file.
""")

file_name = String(help="""
Name of the input file.
""")

accept = String(help="""
Character string of accepted file types for the input. This should be
written like normal html.
""")

data = List(Dict(keys_type=String, values_type=Int), default=[], help="""
List of dictionary containing the inputed data. This the output of the parser.
""")
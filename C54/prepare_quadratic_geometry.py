#!/usr/bin/env python3
from pathlib import Path
import json
P=Path(__file__).resolve().parent/'data';D=json.loads((P/'quadratic_prepass_function.json').read_text());assert D['complete']
(P/'quadratic_geometry_data.m').write_text('QuadraticFieldPolynomialData:='+str(D['field_polynomial'])+';\nQuadraticFunctionData:='+str(D['filter_function'])+';\nQuadraticCuspSignature:='+str(D['base_cusp_signature'])+';\nQuadraticCuspValue:='+str(D['base_cusp_value'])+';\n')

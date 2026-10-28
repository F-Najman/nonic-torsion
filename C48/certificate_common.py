from pathlib import Path
import ast
P=Path(__file__).resolve().parent

def literal(path):
 s=path.read_text().replace('[*','[').replace('*]',']');out={}
 for n in ast.parse(s).body:
  assert isinstance(n,ast.Assign) and len(n.targets)==1 and isinstance(n.targets[0],ast.Name)
  out[n.targets[0].id]=ast.literal_eval(n.value)
 return out

def magma(path):
 out={}
 for stmt in path.read_text().split(';'):
  if ':=' not in stmt:continue
  key,value=stmt.split(':=',1);value=value.strip().replace('[*','[').replace('*]',']');value={'true':'True','false':'False'}.get(value,value)
  out[key.strip()]=ast.literal_eval(value)
 return out

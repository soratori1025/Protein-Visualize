import urllib.request, json
with urllib.request.urlopen('http://localhost:8000/api/secondary-structure/uniprot/P31645') as r:
  data = json.loads(r.read())
  for r in data.get('regions', []):
      print(f"{r.get('type')}: {r.get('description')} ({r.get('start')}-{r.get('end')})")

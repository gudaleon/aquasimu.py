# aquasim_py v1 (skeleton)
Port plan of AQUASIM 2.1g aqcore to Python. See aquasim_py/__init__.py for the module map.
Test: `python3 -c "from aquasim_py import aqufile as a; t=a.load('model.aqu'); print(a.find_section(t,'COMPSYS')[0])"`

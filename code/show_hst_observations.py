from astroquery.mast.missions import MastMissions

missions = MastMissions(mission='hst')

target = 'NGC-1851'
radius = '3'
radius_units = 'arcminutes'
select_cols = [
]

# Use query_object method to resolve the object name into coordinates
results = missions.query_object(
    target,
    radius=radius, 
    radius_units=radius_units, 
    select_cols=select_cols,
    sci_pep_id='12311, 13297',
    sci_spec_1234='F275W,F275W;*,*;F275W,*;F275W;*,F336W,F336W;*,*;F336W,*;F336W;*',
    sci_aec='S')

print(results)


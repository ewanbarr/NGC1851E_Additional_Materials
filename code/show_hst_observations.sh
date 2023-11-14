curl -d "{\"target\":[\"NGC-1851\"],\"radius\":3,\"radius_units\":\"arcminutes\",\"conditions\":[{\"sci_pep_id\":\"12311, 13297\"},{\"sci_spec_1234\":\"F275W,F275W;*,*;F275W,*;F275W;*,F336W,F336W;*,*;F336W,*;F336W;*\"},{\"sci_aec\":\"S\"}],\"limit\":5000}" -H "Content-Type: application/json" -X POST "https://mast.stsci.edu/search/hst/api/v0.1/search"

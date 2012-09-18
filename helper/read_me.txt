iAF1260 information contained in the SBML file format

reactions (format: ‘R_<reaction abbreviation >’) 
- reaction name, reversibility, reaction stoichiometry, gene-protein-reaction (GPR) association, subsystem, E.C. number
metabolites (format: ‘M_<metabolite abbreviation>_<compartment abbreviation>’) 
- metabolite name, compartment, charge, formula (appended to the end of the name, <metabolite name>_FORMULA)
a flux distribution associated with a steady-state modeling simulation 
- lower bound, upper bound, objective coefficient, flux value, reduced cost


SBML File Properties	
file name	Ec-iAF1260-flux1.xml
organism	E. coli K-12 MG1655
model	iAF1260
Biomass Objective Function (BOF)	Ec_biomass_iAF1260_core_59p81M (E. coli biomass objective function (iAF1260) - core - with 59.81 GAM estimate)
flux balance analysis objective	maximize BOF
Growth Associated Maintenance (GAM)	59.81 mmol ATP gDW-1
Non-Growth Associated Maintenance (NGAM)	8.39 mmol ATP gDW-1 hr-1
media conditions	computational minimal media
carbon source	8 mmol glucose gDw-1 hr-1
aerobic or anaerobic	18.5 mmol O2 gDw-1 hr-1
additional constraints	
reactions constrained to zero	CAT, SPODM, SPODMpp, FHL
flux split between reaction pairs	none
	
	
SBML File Properties	
file name	Ec-iAF1260-flux2.xml
organism	E. coli K-12 MG1655
model	iAF1260
Biomass Objective Function (BOF)	Ec_biomass_iAF1260_core_59p81M (E. coli biomass objective function (iAF1260) - core - with 59.81 GAM estimate)
flux balance analysis objective	maximize BOF
Growth Associated Maintenance (GAM)	59.81 mmol ATP gDW-1
Non-Growth Associated Maintenance (NGAM)	8.39 mmol ATP gDW-1 hr-1
media conditions	computational minimal media
carbon source	11.0 mmol glucose gDw-1 hr-1
aerobic or anaerobic	18.2 mmol O2 gDw-1 hr-1
additional constraints	
reactions constrained to zero	152 reactions identified to be unavailable to the cell under glucose aerobic conditions
flux split between reaction pairs	NDH-1:NDH-2 is 1:1 (3 pairs of reactions - different quinone usage)  NADH10:NADH17pp, NADH5:NADH16pp, NADH9:NADH18pp

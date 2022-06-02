

"""# Libraries import"""

from __future__ import print_function
import cobra 
from cobra import Model, Reaction, Metabolite
model = Model('R.Opacus DSM43205')
#import escher
#from escher import Builder
import pandas as pd
import pprint
import csv

model.add_metabolites([
Metabolite('rb15bp', name= 'D-Ribulose 1,5-bisphosphate' , formula='C5H8O11P2' , compartment='c', charge = -4),
Metabolite('CO2', name= 'CO2', formula = 'C1O2', compartment= 'c', charge = 0),
Metabolite('PGA', name = '3-Phospho-D-glycerate', formula = 'C3H4O7P', compartment='c', charge=-3), #acido 3pg
Metabolite('H2O' , name='H2O', formula='H2O1', charge = 0, compartment= 'c'),
Metabolite('H', name = 'Hydron', formula='H', charge= 1),
Metabolite('dpg13' , name='3-Phospho-D-glyceroyl phosphate',formula='C3H4O10P2', compartment='c', charge = -4 ),
Metabolite('ATP' , name='Adenosine 5-triphosphate',formula='C10H12N5O13P3', compartment='c', charge = -4 ),
Metabolite('ADP' , name='Adenosine 5-diphosphate',formula='C10H12N5O10P2', compartment='c', charge = -3 ),
Metabolite('g3p', name='Glyceraldehyde 3-phosphate',formula='C3H5O6P', compartment='c', charge = -2 ),
Metabolite('nad', name='Nicotinamide adenine dinucleotide',formula='C21H26N7O14P2', compartment='c', charge = -1),
Metabolite('nadh', name='Nicotinamide adenine dinucleotide - reduced',formula='C21H27N7O14P2', compartment='c', charge = -2 ),
#pi = Metabolite('pi', name='Phosphate',formula='H1O4P', compartment='c', charge = 0),
Metabolite('fdp', name='D-Fructose 1,6-bisphosphate',formula='C6H10O12P2', compartment='c', charge = -4),
Metabolite('dhap', name='Glycerone phosphate',formula='C3H5O6P', compartment='c', charge = -2),
Metabolite('f6p', name='D-Fructose 6-phosphate',formula='C6H11O9P', compartment='c', charge = -2),
Metabolite('pi_2', name='Phosphate',formula='H1O4P', compartment='c', charge = -2), #Ver si esto es legal. 
Metabolite('e4p', name='D-Erythrose 4-phosphate',formula='C4H7O7P', compartment='c', charge = -2),
Metabolite('xu5p_D', name='D-Xylulose 5-phosphate',formula='C5H9O8P', compartment='c', charge = -2),
Metabolite('s17bp', name='Sedoheptulose 1,7-bisphosphate',formula='C7H12O13P2', compartment='c', charge = -4),
Metabolite('s7p',name='Sedoheptulose 7-phosphate',formula='C7H13O10P', compartment='c', charge = -2 ),
Metabolite('r5p', name='Alpha-D-Ribose 5-phosphate',formula='C5H9O8P', compartment='c', charge = -2),
Metabolite('ru5p_D', name='D-Ribulose 5-phosphate',formula='C5H9O8P', compartment='c', charge = -2),
Metabolite('actp', name='Acetyl phosphate',formula='C2H3O5P', compartment='c', charge = -2),
Metabolite('Glc_aD', name='Alpha-D-glucose',formula='C6H12O6', compartment='c', charge = 0),
Metabolite('g6p_A', name='Alpha-D-glucose 6 phosphate',formula='C6H11O9P', compartment='c', charge = -2),
Metabolite('f6p_B', name='Beta-D-Fructose 6-phosphate',formula='C6H11O9P', compartment='c', charge = -2),
Metabolite('fdp_B', name='Beta-D-Fructose 1,6-bisphosphate',formula='C6H10O12P2', compartment='c', charge = -4),
Metabolite('pg2', name='D-Glycerate 2-phosphate',formula='C3H4O7P', compartment='c', charge = -3),
Metabolite('pep', name='Phosphoenolpyruvate',formula='C3H2O6P', compartment='c', charge = -3),
Metabolite('pyr', name='Pyruvate',formula='C3H3O3', compartment='c', charge = -1),
Metabolite('g6p_B', name='beta-D-Glucose 6-phosphate',formula='C6H11O9P', compartment='c', charge = -2),
Metabolite('pgl6', name='6-phospho-D-glucono-1,5-lactone',formula='C6H9O9P', compartment='c', charge = -2),
Metabolite('pgc6', name='6-Phospho-D-gluconate',formula='C6H10O10P', compartment='c', charge = -3),
Metabolite('ddg6p_2', name='2-Dehydro-3-deoxy-D-gluconate 6-phosphate',formula='C6H8O9P', compartment='c', charge = -3),
#Metabolite('glyc3p', name=' Glycerol 3-phosphate',formula='C3H7O6P', compartment='c', charge = -2), #Formulas sacadas de CHEBI Y METANETX
#Metabolite('acoa', name=' Acyl-CoA',formula='	C22H31N7O17P3RS', compartment='c', charge = -4),
#Metabolite('1-Acylglycerol-3P', name=' 1-acyl-sn-glycero-3-phosphate',formula='C4H6O7PR', compartment='c', charge = -2), #CHEBI:57970
Metabolite('coa', name=' Coenzyme A',formula='C21H32N7O16P3S', compartment='c', charge = -4), #-4 
#Metabolite('Phosphidate', name=' Phosphatidate',formula='C5H5O8PR2', compartment='c', charge = -2), #
#Metabolite('DAG', name='1,2-Diacyl-sn-glycerol',formula='C5H6O5R2', compartment='c', charge = 0), #CHEBI:17815
#Metabolite('TAG', name='Triacylglycerol',formula='C6H5O6R3', compartment='c', charge = 0),
Metabolite('meoh', name='Methanol',formula='CH4O1', compartment='c', charge = 0),
#Metabolite('fald', name='Formaldehyde',formula='CH2O', compartment='c', charge = 0),
#Metabolite('for', name='Formate',formula='CH1O2', compartment='c', charge = -1),
Metabolite('accoa', name='Acetyl-CoA',formula='C23H34N7O17P3S', compartment='c', charge = -4),
Metabolite('cit', name='citrate',formula='C6H5O7', compartment='c', charge = -3),
Metabolite('oaa', name='oxaloacetate',formula='C4H2O5', compartment='c', charge = -2),
Metabolite('acon_C', name='Cis-Aconitate',formula='C6H3O6', compartment='c', charge = -3),
Metabolite('icit', name='Isocitrate',formula='C6H5O7', compartment='c', charge = -3),
Metabolite('osuc', name='Oxalosuccinate',formula='C6H3O7', compartment='c', charge = -3),
Metabolite('nadp', name='Nicotinamide adenine dinucleotide phosphate',formula='C21H25N7O17P3', compartment='c', charge = -3),
Metabolite('nadph', name='Nicotinamide adenine dinucleotide phosphate - reduced',formula='C21H26N7O17P3', compartment='c', charge = -4),
Metabolite('akg', name='2-Oxoglutarate',formula='C5H4O5', compartment='c', charge = -2),
Metabolite('mal_L', name='L-Malate',formula='C4H4O5', compartment='c', charge = -2),
Metabolite('fad', name='Flavin adenine dinucleotide oxidized',formula='C27H31N9O15P2', compartment='c', charge = -2), #-3 opción
Metabolite('fadh2', name='Flavin adenine dinucleotide reduced',formula='C27H33N9O15P2', compartment='c', charge = -2),
Metabolite('fum', name='Fumarate',formula='C4H2O4', compartment='c', charge = -2),
Metabolite('succ_c', name='Succinate',formula='C4H4O4', compartment='c', charge = -2),
Metabolite('succoa', name='Succinyl-CoA',formula='C25H35N7O19P3S', compartment='c', charge = -5),
Metabolite('q8_c', name='Ubiquinone-8',formula='C49H74O4', compartment='c', charge = 0),
Metabolite('q8h2_c' , name = 'Ubiquinol-8', formula = 'C49H76O4' , compartment = 'c' , charge = 0),
Metabolite('hco3_c' , name = 'Bicarbonate', formula = 'HCO3' , compartment = 'c' , charge = -1),
Metabolite('gdp_c' , name = 'gdp_c', formula = 'C10H12N5O11P2' , compartment = 'c' , charge = -3),
Metabolite('gtp_c' , name = 'gtp', formula = 'C10H12N5O14P3' , compartment = 'c' , charge = -4),       
Metabolite('glx' , name = 'Glyoxylate', formula = 'C2H1O3' , compartment = 'c' , charge = -1),
Metabolite('mql8' , name='Menaquinol 8', formula='C51H74O2', charge = 0, compartment= 'c'),
Metabolite('mqn8' , name='Menaquinone 8', formula='C51H72O2', charge = 0, compartment= 'c'),
Metabolite('H2S' , name='Hydrogen sulfide', formula='H2S', charge = -1, compartment= 'c'),
Metabolite('SO3' , name='Sulfite', formula='SO3', charge = -2, compartment= 'c'),
Metabolite('fdxrd' , name='Reduced ferredoxin', formula='Fe2S2', charge = -1, compartment= 'c'), #consultar
Metabolite('fdxo_2_2' , name='Oxidized ferredoxin', formula='Fe2S2', charge = 0, compartment= 'c'), #consultar
Metabolite('no2' , name='Nitrite', formula='NO2', charge = -1, compartment= 'c'),
Metabolite('no3' , name='Nitrate', formula='NO3', charge = -1, compartment= 'c'),
Metabolite('h_e' , name='Hydron', formula='H1', charge = 1, compartment= 'e'),
Metabolite('o2', name= 'O2', formula = 'O2', compartment= 'c', charge = 0),
Metabolite('gln__L', name= 'L-Glutamine', formula = 'C5H10N2O3', compartment= 'c', charge = 0),
Metabolite('glu__L', name= 'L-Glutamate', formula = 'C5H8NO4', compartment= 'c', charge = -1),
#Metabolite('g6p', name= 'D-Glucose 6-phosphate', formula = 'C6H11O9P', compartment= 'c', charge = -2),
Metabolite('ficytC', name= 'Ferricytochrome c', formula = 'C42H54FeN8O6S2', compartment= 'c', charge = 3),
Metabolite('focytC', name= 'Ferrocytochrome C', formula = 'C42H54FeN8O6S2', compartment= 'c', charge = 2),
Metabolite('o2e', name= 'Oxygen', formula = 'O2', compartment= 'e', charge = 0),
Metabolite('CO2e', name= 'CO2', formula = 'C1O2', compartment= 'e', charge = 0),
Metabolite('H2Oe', name= 'H2Oe', formula = 'H2O1', compartment= 'e', charge = 0),
Metabolite('glc__D_e', name= 'beta-D-Glucose', formula = 'C6H12O6', compartment= 'e', charge = 0),
Metabolite('meoh_e', name= 'Methanol', formula = 'C1H4O1', compartment= 'e', charge = 0),
Metabolite('pyr_e', name='Pyruvate extracellular',formula='C3H3O3', compartment='e', charge = -1),
Metabolite('oaa_e', name = 'oxaloacetate extracellular', formula ='C4H2O5', compartment = 'e', charge = -2),
Metabolite('nh4_c', name = 'Ammonium', formula = 'NH4' , compartment = 'c', charge = 1),
Metabolite('nh4_e',name = 'Ammonium', formula = 'NH4' , compartment = 'e', charge = 1 ),
Metabolite('succ_e',name = 'Succinate', formula = 'C4H4O4', compartment='e', charge = -2),
Metabolite('akg_e', name='2-Oxoglutarate',formula='C5H4O5', compartment='e', charge = -2),
Metabolite('gln__L_e', name= 'L-Glutamine', formula = 'C5H10N2O3', compartment= 'e', charge = 0),
Metabolite('glu__L_e', name= 'L-Glutamate', formula = 'C5H8NO4', compartment= 'e', charge = -1),
Metabolite('amp_c', name= 'AMP', formula = 'C10H12N5O7P', compartment= 'c', charge = -2),


])

number_of_metabolites=range(len(model.metabolites))
idx=[model.metabolites[i].id for i in number_of_metabolites]

for i in number_of_metabolites: 
  dummy=model.metabolites.get_by_id(idx[i])
  globals()[idx[i]]=dummy

"""# Reactions, exchanges, sinks and demands

## Intracellular reactions
"""

#1 reaccion 
reaction = Reaction('RBPC') #BiGG
reaction.name = '3-phospho-D-glycerate carboxy-lyase (dimerizing; D-ribulose-1,5-bisphosphate-forming)'
reaction.subsystem = 'Carbon Fixation - Calvin-Benson-Bassham'
reaction.lower_bound = 0.  
reaction.upper_bound = 1000.  #Irreversible
reaction.add_metabolites({
    rb15bp: 1.0,
    CO2: 1.0,
    H2O: 1.0, #Carga variada 
    PGA: -2.0,
    H: -2.0})#Carga variada 

reaction.reaction
reaction.gene_reaction_rule = '(A8G16_RS25910)' 
reaction.genes
model.add_reactions([reaction])

#2 reaccion 
reaction = Reaction('PGK')
reaction.name = 'Phosphoglycerate kinase'
reaction.subsystem = 'Carbon Fixation - Calvin-Benson-Bassham'
reaction.lower_bound = -1000. 
reaction.upper_bound = 1000.  #SEED
reaction.add_metabolites({
	PGA: -1.0,
	ATP: -1.0,
	dpg13:  1.0,
	ADP: 1.0,
	})
reaction.reaction
reaction.gene_reaction_rule = '(A8G16_RS25850 or A8G16_RS33840)' 
reaction.genes
model.add_reactions([reaction])

#3 reaccion 
reaction = Reaction('GAPDHh')
reaction.name = 'Glyceraldehyde 3-phosphate dehydrogenase (NAD)'
reaction.subsystem = 'Carbon Fixation - Calvin-Benson-Bassham'
reaction.lower_bound = -1000. #Reversible SEED
reaction.upper_bound = 1000.  
reaction.add_metabolites({
	dpg13: -1.0,
	H: -1.0,
	nadh: -1.0,
	g3p: 1.0,
	nad: 1.0,
	pi_2: 1.0,
	})
reaction.reaction
reaction.gene_reaction_rule = '(A8G16_RS25845 or A8G16_RS02100 or A8G16_RS33845)' 
reaction.genes
model.add_reactions([reaction])

reaction = Reaction('FBA')
reaction.name = 'Fructose-bisphosphate aldolase'
reaction.subsystem = 'Carbon Fixation - Calvin-Benson-Bassham'
reaction.lower_bound = -1000. #reversible 4.1.2.13
reaction.upper_bound = 1000.  
reaction.add_metabolites({
	dhap: -1.0,
	g3p: -1.0,
	fdp: 1.0,
	})
reaction.reaction
reaction.gene_reaction_rule = '(A8G16_RS25905 or A8G16_RS30370 or A8G16_RS41170)'
reaction.genes
model.add_reactions([reaction])

reaction = Reaction('FBP')
reaction.name = 'Fructose-bisphosphatase'
reaction.subsystem = 'Carbon Fixation - Calvin-Benson-Bassham'
reaction.lower_bound = -1000. #reversible 3.1.3.11
reaction.upper_bound = 1000.  
reaction.add_metabolites({
	f6p: 1.0,
	H2O: -1.0,
	fdp: -1.0,
	pi_2: 1.0
	})
reaction.reaction
reaction.gene_reaction_rule = '(A8G16_RS00015 or A8G16_RS03405 or A8G16_RS25835)' 
model.add_reactions([reaction])

reaction = Reaction('TKT2')
reaction.name = 'Transketolase'
reaction.subsystem = 'Carbon Fixation - Calvin-Benson-Bassham'
reaction.lower_bound = -1000. #reversible
reaction.upper_bound = 1000.  
reaction.add_metabolites({
	f6p: -1.0,
	g3p: -1.0,
	e4p: 1.0,
	xu5p_D: 1.0
	})
reaction.reaction
reaction.gene_reaction_rule = '(A8G16_RS25895 or A8G16_RS00390 or A8G16_RS33795)'
model.add_reactions([reaction])


reaction = Reaction('FBA3')
reaction.name = 'Sedoheptulose 1,7-bisphosphate D-glyceraldehyde-3-phosphate-lyase'
reaction.subsystem = 'Carbon Fixation - Calvin-Benson-Bassham'
reaction.lower_bound = -1000. #reversible 4.1.2.13
reaction.upper_bound = 1000.  
reaction.add_metabolites({
	s17bp: 1.0,
	dhap: -1.0,
	e4p: -1.0,
	})
reaction.reaction
reaction.gene_reaction_rule = '(A8G16_RS25905 or A8G16_RS30370 or A8G16_RS41170)' 
model.add_reactions([reaction])

reaction = Reaction('SBP_1')
reaction.name = 'Sedoheptulose-bisphosphatase'
reaction.subsystem = 'Carbon Fixation - Calvin-Benson-Bassham'
reaction.lower_bound = -1000. #3.1.3.11
reaction.upper_bound = 1000.  
reaction.add_metabolites({
	s17bp: -1.0,
	H2O: -1.0,
	pi_2: 1.0,
	s7p: 1.0,

	})
reaction.reaction
reaction.gene_reaction_rule = '(A8G16_RS00015 or A8G16_RS03405 or A8G16_RS25835)' #Buscar Gen
model.add_reactions([reaction])

reaction = Reaction('TKT1')
reaction.name = 'Transketolase'
reaction.subsystem = 'Carbon Fixation - Calvin-Benson-Bassham'
reaction.lower_bound = -1000. 
reaction.upper_bound = 1000.  
reaction.add_metabolites({
	xu5p_D: 1.0,
	r5p: 1.0,
	g3p: -1.0,
	s7p: -1.0,

	})
reaction.reaction
reaction.gene_reaction_rule = '(A8G16_RS25895 or A8G16_RS00390 or A8G16_RS33795)' 
model.add_reactions([reaction])


reaction = Reaction('RPIh')
reaction.name = 'Ribose-5-phosphate isomerase'
reaction.subsystem = 'Carbon Fixation - Calvin-Benson-Bassham'
reaction.lower_bound = -1000. #reversible
reaction.upper_bound = 1000.  
reaction.add_metabolites({
	ru5p_D: 1.0,
	r5p: -1.0,
	})
reaction.reaction
reaction.gene_reaction_rule = '(A8G16_RS23885)' 
model.add_reactions([reaction])


reaction = Reaction('PRUK')
reaction.name = 'Phosphoribulokinase'
reaction.subsystem = 'Carbon Fixation - Calvin-Benson-Bassham'
reaction.lower_bound = 0. #irreversible
reaction.upper_bound = 1000.  
reaction.add_metabolites({
	ru5p_D: -1.0,
	ATP: -1.0,
	rb15bp: 1.0,
	ADP: 1.0,
	H: 1.0,
	})
reaction.reaction
reaction.gene_reaction_rule = '(A8G16_RS25840)' 
model.add_reactions([reaction])

reaction = Reaction('RPE')
reaction.name = 'Ribulose 5-phosphate 3-epimerase'
reaction.subsystem = 'Carbon Fixation - Calvin-Benson-Bassham'
reaction.lower_bound = -1000. #reversible
reaction.upper_bound = 1000.  
reaction.add_metabolites({
	xu5p_D: -1.0,
	ru5p_D: 1.0,
	})
reaction.reaction
reaction.gene_reaction_rule = '(A8G16_RS25890 or A8G16_RS15870)' 
model.add_reactions([reaction])


reaction = Reaction('TPI')
reaction.name = 'Triose-phosphate isomerase'
reaction.subsystem = 'Carbon Fixation - Calvin-Benson-Bassham'
reaction.lower_bound = -1000. #reversible
reaction.upper_bound = 1000.  
reaction.add_metabolites({
	g3p: -1.0,
	dhap: 1.0,
	})
reaction.reaction
reaction.gene_reaction_rule = '(A8G16_RS33835)' 
model.add_reactions([reaction])


reaction = Reaction('PKETX')
reaction.name = 'Phosphoketolase (xylulose-5-phosphate utilizing)'
reaction.subsystem = 'Carbon Fixation - Calvin-Benson-Bassham'
reaction.lower_bound = -1000. #reversible confirmar
reaction.upper_bound = 1000.  
reaction.add_metabolites({
	g3p: 1.0,
	H2O: 1.0,
	actp: 1.0,
	pi_2: -1.0,
	xu5p_D: -1.0,
	})
reaction.reaction
reaction.gene_reaction_rule = '(A8G16_RS31985)'
model.add_reactions([reaction])

#EMP

reaction = Reaction('GLUKA')
reaction.name = 'Glucokinase/hexokinase (glc-A)'
reaction.subsystem = 'Embden-Meyerhof-Parnas'
reaction.lower_bound = 0. #irreversible
reaction.upper_bound = 1000.  
reaction.add_metabolites({
	ADP: 1.0,
	H: 1.0,
	g6p_A: 1.0,
	ATP: -1.0,
	Glc_aD: -1.0,
	})
reaction.reaction
reaction.gene_reaction_rule = '(GLUKAgen)' #Buscar Gen
model.add_reactions([reaction])


reaction = Reaction('PGI')
reaction.name = 'glucose-6-phosphate isomerase'
reaction.subsystem = 'Embden-Meyerhof-Parnas'
reaction.lower_bound = -1000. #reversible 5.3.1.9
reaction.upper_bound = 1000.  
reaction.add_metabolites({
	f6p_B: 1.0,
	g6p_B: -1.0,
	})
reaction.reaction
reaction.gene_reaction_rule = '(PGIAGen_)' #Buscar Gen
model.add_reactions([reaction])

reaction = Reaction('PFKh')
reaction.name = '6-Phosphofructokinase'
reaction.subsystem = 'Embden-Meyerhof-Parnas'
reaction.lower_bound = 0. #ireversible
reaction.upper_bound = 1000.  
reaction.add_metabolites({
	fdp_B: 1.0,
	ADP: 1.0,
	H: 1.0,
	ATP: -1.0,
	f6p_B: -1.0,
		})
reaction.reaction
reaction.gene_reaction_rule = '(PFKhgen)' #Buscar Gen
model.add_reactions([reaction])

reaction = Reaction('BFBP')
reaction.name = 'Beta-D-Fructose 1,6-bisphosphate 1-phosphohydrolase'
reaction.subsystem = 'Embden-Meyerhof-Parnas'
reaction.lower_bound = 0. #ireversible 3.1.3.11/ FBP/SBP_1
reaction.upper_bound = 1000.  
reaction.add_metabolites({
	fdp_B: -1.0,
	H2O: -1.0,
	pi_2: 1.0,
	f6p_B: 1.0,
		})
reaction.reaction
reaction.gene_reaction_rule = '(A8G16_RS00015 or A8G16_RS03405 or A8G16_RS25835)' #CHECK
model.add_reactions([reaction])

reaction = Reaction('FBAh')
reaction.name = 'fructose-bisphosphate aldolase'
reaction.subsystem = 'Embden-Meyerhof-Parnas'
reaction.lower_bound = -1000. #reversible 4.1.2.13
reaction.upper_bound = 1000.  
reaction.add_metabolites({
	fdp_B: -1.0,
	dhap: 1.0,
	g3p: 1.0,
		})
reaction.reaction
reaction.gene_reaction_rule = '(A8G16_RS25905 or A8G16_RS30370 or A8G16_RS41170)' #Buscar Gen
model.add_reactions([reaction])

reaction = Reaction('GAPDHh-EMP')
reaction.name = 'Glyceraldehyde 3-phosphate dehydrogenase (NAD)'
reaction.subsystem = 'Embden-Meyerhof-Parnas' #1.2.1.12
reaction.lower_bound = -1000. 
reaction.upper_bound = 1000.  
reaction.add_metabolites({
	dpg13: -1.0,
	H: -1.0,
	nadh: -1.0,
	g3p: 1.0,
	nad: 1.0,
	pi_2: 1.0,
	})
reaction.reaction
reaction.gene_reaction_rule = '(A8G16_RS25845 or A8G16_RS02100 or A8G16_RS33845)' 
reaction.genes
model.add_reactions([reaction])

reaction = Reaction('PGK-EMP')
reaction.name = 'Phosphoglycerate kinase'
reaction.subsystem = 'Embden-Meyerhof-Parnas' #2.7.2.3
reaction.lower_bound = -1000. 
reaction.upper_bound = 1000.  
reaction.add_metabolites({
	PGA: -1.0,
	ATP: -1.0,
	dpg13:  1.0,
	ADP: 1.0,
	})
reaction.reaction
reaction.gene_reaction_rule = '(A8G16_RS25850 or A8G16_RS33840)' 
reaction.genes
model.add_reactions([reaction])

reaction = Reaction('TPI-EMP')
reaction.name = 'Triose-phosphate isomerase'
reaction.subsystem = 'Embden-Meyerhof-Parnas'
reaction.lower_bound = -1000. #reversible 	5.3.1.1
reaction.upper_bound = 1000.  
reaction.add_metabolites({
	g3p: -1.0,
	dhap: 1.0,
	})
reaction.reaction
reaction.gene_reaction_rule = '(A8G16_RS33835)' 
model.add_reactions([reaction])


reaction = Reaction('PGM')
reaction.name = 'Phosphoglycerate mutase'
reaction.subsystem = 'Embden-Meyerhof-Parnas' #5.4.2.11/5.4.2.12 DEP/IND respc. Se ponen genes de ambas enzimas. Consultar 
reaction.lower_bound = -1000. 
reaction.upper_bound = 1000.  
reaction.add_metabolites({
	PGA: 1.0,
	pg2: -1.0,
	})
reaction.reaction
reaction.gene_reaction_rule = '(A8G16_RS01505 or A8G16_RS24150 or A8G16_RS17480)' #Buscar Gen
reaction.genes
model.add_reactions([reaction])

reaction = Reaction('ENO')
reaction.name = 'Enolase'
reaction.subsystem = 'Embden-Meyerhof-Parnas' #4.2.1.11 Reversible
reaction.lower_bound = -1000. 
reaction.upper_bound = 1000.  
reaction.add_metabolites({
	H2O: 1.0,
	pep: 1.0,
	pg2: -1.0,
	})
reaction.reaction
reaction.gene_reaction_rule = '(A8G16_RS14105 or A8G16_RS24170)' #Buscar Gen
reaction.genes
model.add_reactions([reaction])

reaction = Reaction('PYK')
reaction.name = 'Pyruvate kinase' #2.7.1.40
reaction.subsystem = 'Embden-Meyerhof-Parnas'
reaction.lower_bound = 0. #Irreversible
reaction.upper_bound = 1000.  
reaction.add_metabolites({
	ADP: -1.0,
	H: -1.0,
	pep: -1.0,
	ATP: 1.0,
	pyr: 1.0,
	})
reaction.reaction
reaction.gene_reaction_rule = '(A8G16_RS27205 or A8G16_RS39215)' #Buscar Gen
reaction.genes
model.add_reactions([reaction])

#ED

reaction = Reaction('G6PBDH')
reaction.name = 'Glucose 6-phosphate dehydrogenase'
reaction.subsystem = 'Entner–Doudoroff pathway'
reaction.lower_bound = 0. #irreversible EC.1.1.1.49
reaction.upper_bound = 1000.  
reaction.add_metabolites({
	g6p_B: -1.0,
	nadp: -1.0,
	pgl6: 1.0,
	H: 1.0,
	nadph:1.0,
	})
reaction.reaction
reaction.gene_reaction_rule = '(A8G16_RS08400 or A8G16_RS13485 or A8G16_RS33805)' 
model.add_reactions([reaction])

reaction = Reaction('PGL')
reaction.name = '6-phosphogluconolactonase'
reaction.subsystem = 'Entner–Doudoroff pathway'
reaction.lower_bound = -1000. #Irreversible EC Number: 3.1.1.31
reaction.upper_bound = 1000.  
reaction.add_metabolites({
	pgl6: -1.0,
	H2O: -1.0,
	pgc6: 1.0,
	H: 1.0,
	})
reaction.reaction
reaction.gene_reaction_rule = '(A8G16_RS00100 or A8G16_RS33815)' #Buscar Gen
model.add_reactions([reaction])

reaction = Reaction('EDD')
reaction.name = '6-phosphogluconate dehydratase'
reaction.subsystem = 'Entner–Doudoroff pathway'
reaction.lower_bound = 0. #Irreversible EC.4.2.1.12
reaction.upper_bound = 1000.  
reaction.add_metabolites({
	pgc6: -1.0,
	H2O: 1.0,
	ddg6p_2: 1.0,
	})

reaction.reaction
reaction.gene_reaction_rule = '(A8G16_RS08405)' #Buscar Gen
model.add_reactions([reaction])

reaction = Reaction('EDA')
reaction.name = '2-dehydro-3-deoxy-phosphogluconate aldolase'
reaction.subsystem = 'Entner–Doudoroff pathway'
reaction.lower_bound = 0. #Irreversible EC. 4.1.2.14
reaction.upper_bound = 1000.  
reaction.add_metabolites({
	pyr: 1.0,
	g3p: 1.0,
	ddg6p_2: -1.0,
	})

reaction.reaction
reaction.gene_reaction_rule = '(A8G16_RS08410)' #Buscar Gen
model.add_reactions([reaction])

#Ruta de Kennedy
"""
reaction = Reaction('GPAT1_2')
reaction.name = 'glycerol-3-phosphate O-acyltransferase 1/2'
reaction.subsystem = 'Kennedy Pathway'
reaction.lower_bound = -1000. #reversible
reaction.upper_bound = 1000.  
reaction.add_metabolites({
	glyc3p: -1.0,
	AcylCoA: -1.0,
	coa: 1.0,
	Acylglycerol13P: 1.0,
	})
reaction.reaction
reaction.gene_reaction_rule = '(GPAT1_2gen)' #Buscar Gen
model.add_reactions([reaction])

reaction = Reaction('plsC/LPA')
reaction.name = ''
reaction.subsystem = 'Kennedy Pathway'
reaction.lower_bound = -1000. #reversible
reaction.upper_bound = 1000.  
reaction.add_metabolites({
	coa: -1.0,
	AcylCoA: 1.0,
	phosphidate: -1.0,
	Acylglycerol13P: 1.0,
	})
reaction.reaction
reaction.gene_reaction_rule = '(plsCgen)' #Buscar Gen
model.add_reactions([reaction])

reaction = Reaction('PAP')
reaction.name = 'phosphatidate phosphatase'
reaction.subsystem = 'Kennedy Pathway'
reaction.lower_bound = -1000. #reversible
reaction.upper_bound = 1000.  
reaction.add_metabolites({
	diacyl_glycerol: 1.0,
	pi_2: 1.0,
	phosphidate: -1.0,
	H2O: -1.0,
	})
reaction.reaction
reaction.gene_reaction_rule = '(papgen)' #Buscar Gen
model.add_reactions([reaction])

reaction = Reaction('DGAT1')
reaction.name = 'diacylglycerol O-acyltransferase 1'
reaction.subsystem = 'Kennedy Pathway'
reaction.lower_bound = -1000. #reversible
reaction.upper_bound = 1000.  
reaction.add_metabolites({
	diacyl_glycerol: 1.0,
	AcylCoA: 1.0,
	coa: -1.0,
	TAG: -1.0,
	})
reaction.reaction
reaction.gene_reaction_rule = '(DGAT1gen)' #Buscar Gen
model.add_reactions([reaction])
"""
#PP  

reaction = Reaction('G6PI_1')
reaction.name = 'glucose-6-phosphate isomerase'
reaction.subsystem = 'Pentose Phosphate pathway'
reaction.lower_bound = -1000. #reversible EC. 5.3.1.9
reaction.upper_bound = 1000.  
reaction.add_metabolites({
	g6p_B: 1.0,
	g6p_A: -1.0,
	})
reaction.reaction
reaction.gene_reaction_rule = '(G6PI_1gen)' #Buscar Gen
model.add_reactions([reaction])

reaction = Reaction('PGDHh')
reaction.name = '6-phosphogluconate dehydrogenase'
reaction.subsystem = 'Pentose Phosphate pathway'
reaction.lower_bound = 0. #Irreversible Ec. 1.1.1.44/1.1.1.343 Se utilizan los dos genes 
reaction.upper_bound = 1000.  
reaction.add_metabolites({
	pgc6: -1.0,
	nadp: -1.0,
	nadph: 1.0,
	ru5p_D: 1.0,
	CO2: 1.0,
	})
reaction.reaction
reaction.gene_reaction_rule = '(A8G16_RS27605 or A8G16_RS00020)' #Buscar Gen
model.add_reactions([reaction])

reaction = Reaction('RPE-PPP')
reaction.name = 'Ribulose 5-phosphate 3-epimerase'
reaction.subsystem = 'Pentose Phosphate pathway'
reaction.lower_bound = -1000. #reversible 5.1.3.1
reaction.upper_bound = 1000.  
reaction.add_metabolites({
	xu5p_D: -1.0,
	ru5p_D: 1.0,
	})
reaction.reaction
reaction.gene_reaction_rule = '(A8G16_RS25890 or A8G16_RS15870)' 
model.add_reactions([reaction])

reaction = Reaction('RPIh-PPP')
reaction.name = 'Ribose-5-phosphate isomerase'
reaction.subsystem = 'Pentose Phosphate pathway'
reaction.lower_bound = -1000. #reversible 5.3.1.6
reaction.upper_bound = 1000.  
reaction.add_metabolites({
	ru5p_D: 1.0,
	r5p: -1.0,
	})
reaction.reaction
reaction.gene_reaction_rule = '(A8G16_RS23885)' 
model.add_reactions([reaction])

reaction = Reaction('TAh')
reaction.name = 'transaldolase'
reaction.subsystem = 'Pentose Phosphate pathway'
reaction.lower_bound = -1000. #reversible EC. 2.2.1.2
reaction.upper_bound = 1000.  
reaction.add_metabolites({
	g3p: -1.0,
	s7p: -1.0,
	e4p: 1.0,
	f6p_B: 1.0,
	})
reaction.reaction
reaction.gene_reaction_rule = '(A8G16_RS00035 or A8G16_RS21695 or A8G16_RS33800)' #Buscar Gen
model.add_reactions([reaction])

reaction = Reaction('TKT2-PPP')
reaction.name = 'Transketolase'
reaction.subsystem = 'Pentose Phosphate pathway'
reaction.lower_bound = -1000. #reversible
reaction.upper_bound = 1000.  
reaction.add_metabolites({
	f6p_B: -1.0,
	g3p: -1.0,
	e4p: 1.0,
	xu5p_D: 1.0
	})
reaction.reaction
reaction.gene_reaction_rule = '(A8G16_RS25895 or A8G16_RS00390 or A8G16_RS33795)'
model.add_reactions([reaction])

reaction = Reaction('TKT1-PPP')
reaction.name = 'Transketolase'
reaction.subsystem = 'Pentose Phosphate pathway'
reaction.lower_bound = -1000. 
reaction.upper_bound = 1000.  
reaction.add_metabolites({
	xu5p_D: 1.0,
	r5p: 1.0,
	g3p: -1.0,
	s7p: -1.0,

	})
reaction.reaction
reaction.gene_reaction_rule = '(A8G16_RS25895 or A8G16_RS00390 or A8G16_RS33795)' 
model.add_reactions([reaction])

"""
#Methanol degradation

reaction = Reaction('ALCD1')
reaction.name = 'methanol dehydrogenase'
reaction.subsystem = ' Methanol degradation'
reaction.lower_bound = -1000. #reversible
reaction.upper_bound = 1000.  
reaction.add_metabolites({
	meoh: -1.0,
	nad: -1.0,
	fald: 1.0,
	nadh: 1.0,
	H: 1.0,

	})
reaction.reaction
reaction.gene_reaction_rule = '(ALCD1Gen_)' #Buscar Gen
model.add_reactions([reaction])

reaction = Reaction('ALDD1')
reaction.name = 'formaldehyde dehydrogenase'
reaction.subsystem = ' Methanol degradation'
reaction.lower_bound = -1000. #reversible
reaction.upper_bound = 1000.  
reaction.add_metabolites({
	H2O: -1.0,
	fald: -1.0,
	nad: -1.0,
	nadh: 1.0,
	H: 2.0,
	for_:1.0,
	})
reaction.reaction
reaction.gene_reaction_rule = '(ALDD1Gen_)' #Buscar Gen
model.add_reactions([reaction])

reaction = Reaction('FDH')
reaction.name = 'formate dehydrogenase'
reaction.lower_bound = -1000. #reversible
reaction.upper_bound = 1000.  
reaction.add_metabolites({
	CO2: 1.0,
	nad: -1.0,
	nadh: 1.0,
	for_:-1.0,
	})
reaction.reaction
reaction.gene_reaction_rule = '(fdhgen)' #Buscar Gen
model.add_reactions([reaction])
"""

#TCA Cycle

reaction = Reaction('CS')
reaction.name = 'citrate synthase'
reaction.subsystem = ' TCA Cycle'
reaction.lower_bound = 0. #irreversible 2.3.3.16
reaction.upper_bound = 1000.  
reaction.add_metabolites({
	accoa: -1.0,
	H2O: -1.0,
	oaa: -1.0,
	cit:1.0,
	coa: 1.0,
	H:1.0,
	})
reaction.reaction
reaction.gene_reaction_rule = '(A8G16_RS23015 or A8G16_RS23035)' 
model.add_reactions([reaction])

reaction = Reaction('ACONTa')
reaction.name = 'Aconitase (Aconitase (half-reaction A, Citrate hydro-lyase)'
reaction.subsystem = ' TCA Cycle'
reaction.lower_bound = -1000. #reversible  4.2.1.3
reaction.upper_bound = 1000.  
reaction.add_metabolites({
	cit: -1.0,
	acon_C: 1.0,
	H2O: 1.0,
	})
reaction.reaction
reaction.gene_reaction_rule = '(A8G16_RS21735 or A8G16_RS39530)' 
model.add_reactions([reaction])

reaction = Reaction('ACONTb')
reaction.name = 'Aconitase (half-reaction B, Isocitrate hydro-lyase)'
reaction.subsystem = ' TCA Cycle'
reaction.lower_bound = -1000. #reversible 4.2.1.3
reaction.upper_bound = 1000.  
reaction.add_metabolites({
	icit: 1.0,
	acon_C: -1.0,
	H2O: -1.0,
	})
reaction.reaction
reaction.gene_reaction_rule = '(ACONTb_gen)' 
model.add_reactions([reaction])

reaction = Reaction('ICITRED')
reaction.name = 'Isocitrate:NADP+ oxidoreductase (decarboxylating)'
reaction.subsystem = ' TCA Cycle'
reaction.lower_bound = -1000. #reversible EC.1.1.1.42
reaction.upper_bound = 1000.  
reaction.add_metabolites({
	icit: -1.0,
	nadp: -1.0,
	nadph: 1.0,
	osuc: 1.0,
	H: 1.0,
	})
reaction.reaction
reaction.gene_reaction_rule = '(A8G16_RS25530 or A8G16_RS28305)' 
model.add_reactions([reaction])

#Half reaction       

reaction = Reaction('OSUCCL')
reaction.name = 'Oxalosuccinate carboxy-lyase (2-oxoglutarate-forming)'
reaction.subsystem = ' TCA Cycle'
reaction.lower_bound = -1000. #reversible EC.1.1.1.42
reaction.upper_bound = 1000.  
reaction.add_metabolites({
	osuc: -1.0,
	H: -1.0,
	akg: 1.0,
	CO2: 1.0,
	})
reaction.reaction
reaction.gene_reaction_rule = '(A8G16_RS25530 or A8G16_RS28305)' 
model.add_reactions([reaction])

reaction = Reaction('MDH')
reaction.name = 'Malate dehydrogenase '
reaction.subsystem = ' TCA Cycle'
reaction.lower_bound = -1000. #reversible EC. 1.1.1.37
reaction.upper_bound = 1000.  
reaction.add_metabolites({
	nad: -1.0,
	mal_L: -1.0,
	oaa: 1.0,
	nadh: 1.0,
	H: 1.0,
	})
reaction.reaction
reaction.gene_reaction_rule = '(A8G16_RS11075)' 
model.add_reactions([reaction])

reaction = Reaction('FUM')
reaction.name = 'Fumarase'
reaction.subsystem = ' TCA Cycle'
reaction.lower_bound = -1000. #reversible EC.4.2.1.2
reaction.upper_bound = 1000.  
reaction.add_metabolites({
	fum: -1.0,
	H2O: -1.0,
	mal_L: 1.0,
	})
reaction.reaction
reaction.gene_reaction_rule = '(A8G16_RS03235)' 
model.add_reactions([reaction])

reaction = Reaction('SUCOAS')
reaction.name = 'Succinyl-CoA synthetase (ADP-forming)'
reaction.subsystem = ' TCA Cycle'
reaction.lower_bound = -1000. #reversible EC. 6.2.1.5* /2.8.3.18 
reaction.upper_bound = 1000.  
reaction.add_metabolites({
	ATP: -1.0,
	coa: -1.0,
	succ_c: -1.0,
	ADP: 1.0, 
	pi_2: 1.0, 
	succoa: 1.0,
	})
reaction.reaction
reaction.gene_reaction_rule = '(A8G16_RS32760 or A8G16_RS32765)' 
model.add_reactions([reaction])

reaction = Reaction('SUCDi')
reaction.name = 'Succinate dehydrogenase (irreversible)' 
reaction.subsystem = ' TCA Cycle'
reaction.lower_bound = 0. #reversible #1.3.5.1
reaction.upper_bound = 1000.  
reaction.add_metabolites({
	succ_c: -1.0,
	q8_c: -1.0,
	fum: 1.0,
	q8h2_c: 1.0,
	})
reaction.reaction
reaction.gene_reaction_rule = '(A8G16_RS11065 or A8G16_RS11070 or A8G16_RS27415)' 
model.add_reactions([reaction])

reaction = Reaction('AKGDH')
reaction.name = '2-Oxogluterate dehydrogenase'
reaction.subsystem = ' TCA Cycle'
reaction.lower_bound = 0. #reversible #1.8.1.4
reaction.upper_bound = 1000.  
reaction.add_metabolites({
	akg: -1.0,
	coa: -1.0,
	nad: -1.0,
	CO2: 1.0,
	nadh: 1.0,
	succoa: 1.0,
	})
reaction.reaction
reaction.gene_reaction_rule = '(A8G16_RS01180 or A8G16_RS14355 or A8G16_RS21205 or A8G16_RS23340 or A8G16_RS28185)' 
model.add_reactions([reaction])

#Anaplerotic Pathway

reaction = Reaction('PC')
reaction.name = 'Pyruvate carboxylase'
reaction.subsystem = ' Anaplerotic Pathway'
reaction.lower_bound = 0. #reversible
reaction.upper_bound = 1000.  
reaction.add_metabolites({
	ATP: -1.0,
	pyr: -1.0,
	hco3_c: -1.0,
	ADP: 1.0,
	pi_2: 1.0,
	oaa: 1.0,
	H: 1.0,
	})
reaction.reaction
reaction.gene_reaction_rule = '(A8G16_RS22465)' 
model.add_reactions([reaction])

reaction = Reaction('PPC')
reaction.name = 'Phosphoenolpyruvate carboxylase'
reaction.subsystem = ' Anaplerotic Pathway'
reaction.lower_bound = 0. #Irreversible
reaction.upper_bound = 1000.  
reaction.add_metabolites({
	pi_2: 1.0,
	oaa: 1.0,
	H: 1.0,
	H2O: -1.0,
	CO2: -1.0,
	pep: -1.0,
	})
reaction.reaction
reaction.gene_reaction_rule = '(A8G16_RS22465)' 
model.add_reactions([reaction])


reaction = Reaction('PEPCK_re')
reaction.name = 'Phosphoenolpyruvate carboxykinase (GTP)'
reaction.subsystem = ' Anaplerotic Pathway'
reaction.lower_bound = 0. #Irreversible Ec. 4.1.1.32 Prokka indica que es GTP dependiente NO PPCK
reaction.upper_bound = 1000.  
reaction.add_metabolites({
	oaa: -1.0,
	gtp_c: -1.0,
	CO2: 1.0,
	gdp_c: 1.0,
	pep:1.0,
	})
reaction.reaction
reaction.gene_reaction_rule = '(A8G16_RS11825 or A8G16_RS11840 or A8G16_RS15295 or A8G16_RS15300 or A8G16_RS29775 or A8G16_RS30415 or A8G16_RS31685)' 
model.add_reactions([reaction])

reaction = Reaction('ME1')
reaction.name = 'Malic enzyme (NAD)'
reaction.subsystem = ' Anaplerotic Pathway'
reaction.lower_bound = 0. #reversible EC.1.1.1.38
reaction.upper_bound = 1000.  
reaction.add_metabolites({
	nad: -1.0,
	mal_L: -1.0,
	nadh: 1.0,
	CO2: 1.0,
	pyr: 1.0,
	})
reaction.reaction
reaction.gene_reaction_rule = '(A8G16_RS14495 or A8G16_RS16130 or A8G16_RS24755 or A8G16_RS24755 or A8G16_RS29645)' 
model.add_reactions([reaction])

#Glyoxylate pathway

reaction = Reaction('ICL')
reaction.name = 'Isocitrate lyase'
reaction.subsystem = ' Glyoxylate Pathway'
reaction.lower_bound = 0. #Irreversible EC. 4.1.3.1
reaction.upper_bound = 1000.   
reaction.add_metabolites({
	icit: -1.0,
	succ_c: 1.0,
	glx: 1.0,
	})
reaction.reaction
reaction.gene_reaction_rule = '(A8G16_RS01270)' 
model.add_reactions([reaction])


reaction = Reaction('MALS')
reaction.name = 'Malate synthase' 
reaction.subsystem = ' Glyoxylate Pathway'
reaction.lower_bound = 0. #irreversible EC 2.3.3.9
reaction.upper_bound = 1000.  
reaction.add_metabolites({
	accoa: -1.0,
	H2O: -1.0,
	glx: -1.0,
	mal_L: 1.0,
	H: 1.0,
	coa:1.0,
	})
reaction.reaction
reaction.gene_reaction_rule = '(A8G16_RS28620 or A8G16_RS38450 or A8G16_RS39500)' 
model.add_reactions([reaction])


#ETC

reaction = Reaction('NADH6')
reaction.name = 'NADH dehydrogenase (ubiquinone-8 & 3.5 protons)' 
reaction.subsystem = ' Electron Transport Chain'
reaction.lower_bound = 0. #irreversible EC 1.6.5.2
reaction.upper_bound = 1000.  
reaction.add_metabolites({
	H: -4.5,
	q8_c: -1.0,
	nadh: -1.0,
	q8h2_c: 1.0,
	nad: 1.0,
	h_e:3.5,
	})
reaction.reaction
reaction.gene_reaction_rule = '(A8G16_RS01195 or A8G16_RS10975 or A8G16_RS11565 or A8G16_RS20225 or A8G16_RS29330)' 
model.add_reactions([reaction])

reaction = Reaction('NADH10')
reaction.name = 'NADH dehydrogenase (menaquinone-8 & 0 protons)' 
reaction.subsystem = ' Electron Transport Chain'
reaction.lower_bound = 0. #buscar EC
reaction.upper_bound = 1000.  
reaction.add_metabolites({
	H: -1.0,
	mqn8: -1.0,
	nadh: -1.0,
	mql8: 1.0,
	nad: 1.0,
	})
reaction.reaction
reaction.gene_reaction_rule = '(NADH10gen)' 
model.add_reactions([reaction])

reaction = Reaction('sulfite-reductase-(NADPH)')
reaction.name = 'hydrogen-sulfide:NADP+ oxidoreductase' 
reaction.subsystem = ' Electron Transport Chain'
reaction.lower_bound = -1000. #reversible EC 1.8.1.2
reaction.upper_bound = 1000.  
reaction.add_metabolites({
	nadp: -3.0,
	H2O: -3.0,
	H2S: -1.0,
	H: 4.0,
	SO3: 1.0,
	nadph:3.0,
	})
reaction.reaction
reaction.gene_reaction_rule = '(sulfgen)' 
model.add_reactions([reaction])


reaction = Reaction('FNOR')
reaction.name = 'Ferredoxin---NADP+ reductase' 
reaction.subsystem = 'Electron Transport Chain'
reaction.lower_bound = -1000. #reversible 1.18.1.2
reaction.upper_bound = 1000.  
reaction.add_metabolites({
	H: -1.0,
	nadp: -1.0,
	fdxrd: -2.0,
	nadph: 1.0,
	fdxo_2_2: 2.0,
	})
reaction.reaction
reaction.gene_reaction_rule = '(A8G16_RS31920)' 
model.add_reactions([reaction])

reaction= Reaction('SUCD1')
reaction.name = 'Succinate dehydrogenase' 
reaction.subsystem = 'Electron Transport Chain'
reaction.lower_bound = 0. #reversible 1.3.5.1
reaction.upper_bound = 1000.  
reaction.add_metabolites({
	fadh2: -1.0,
	q8_c: -1.0,
	fad: 1.0,
	q8h2_c: 1.0,
	})
reaction.reaction
reaction.gene_reaction_rule = '(A8G16_RS11065 or A8G16_RS11070 or A8G16_RS27415)' 
model.add_reactions([reaction])

reaction = Reaction('SUCDi-ETC')
reaction.name = 'Succinate dehydrogenase (irreversible)' 
reaction.subsystem = 'Electron Transport Chain'
reaction.lower_bound = 0. #irreversible #1.3.5.1
reaction.upper_bound = 1000.  
reaction.add_metabolites({
	succ_c: -1.0,
	q8_c: -1.0,
	fum: 1.0,
	q8h2_c: 1.0,
	})
reaction.reaction
reaction.gene_reaction_rule = '(A8G16_RS11065 or A8G16_RS11070 or A8G16_RS27415)' 
model.add_reactions([reaction])

reaction = Reaction('CYTBD')
reaction.name = 'Cytochrome oxidase bd (ubiquinol-8: 2 protons)' 
reaction.subsystem = 'Electron Transport Chain'
reaction.lower_bound = 0. #irreversible #7.1.1.7 Cytochrome bd-I ubiquinol oxidase subunit 2
reaction.upper_bound = 1000.  
reaction.add_metabolites({
	H: -2.0,
	o2: -0.5,
	q8h2_c: -1.0,
	q8_c: 1.0,
	H2O: 1.0,
	h_e: 2.0,
	})
reaction.reaction
reaction.gene_reaction_rule = '(A8G16_RS11065 or A8G16_RS11070 or A8G16_RS27415)' 
model.add_reactions([reaction])

reaction = Reaction('CYTBD2pp')
reaction.name = 'Cytochrome oxidase bd (menaquinol-8: 2 protons) (periplasm)' 
reaction.subsystem = 'Electron Transport Chain'
reaction.lower_bound = -1000. 
reaction.upper_bound = 1000.  
reaction.add_metabolites({
	H: -2.0,
	mql8: -1.0,
	o2: - 0.5,
	H2O: 1.0,
	mqn8: 1.0,
	h_e: 2.0,
	})
reaction.reaction
reaction.gene_reaction_rule = '(CYTBD2ppGen)' 
model.add_reactions([reaction])

reaction = Reaction('CYOO2pp')
reaction.name = 'Cytochrome-c oxidase (2 protons translocated) periplasm' 
reaction.subsystem = 'Electron Transport Chain' # EC 1.9.3.1
reaction.lower_bound = -1000. 
reaction.upper_bound = 1000.  
reaction.add_metabolites({
	H: -4.0,
	focytC: -2.0,
	o2: - 0.5,
	H2O: 1.0,
	ficytC: 2.0,
	h_e: 2.0,
	})
reaction.reaction
reaction.gene_reaction_rule = '(CY002ppGen)' 
model.add_reactions([reaction])


reaction = Reaction('NO3R1')
reaction.name = 'Nitrate reductase (Ubiquinol-8)' 
reaction.subsystem = 'Electron Transport Chain'
reaction.lower_bound = -1000. #Reversible # 1.7.99.4 se remplaza por ec:1.9.6.1 periplasmic ?
reaction.upper_bound = 1000.  
reaction.add_metabolites({
	H: -2.0,
	no3: -1.0,
	q8h2_c: -1.0,
	H2O: 1.0,
	no2:1.0,
	q8_c: 1.0,
	h_e: 2.0
	})
reaction.reaction
reaction.gene_reaction_rule = '(A8G16_RS23470 or A8G16_RS26740)' 
model.add_reactions([reaction])

reaction = Reaction('NO3R2')
reaction.name = 'Nitrate reductase (Menaquinol-8)' 
reaction.subsystem = 'Electron Transport Chain'
reaction.lower_bound = -1000. #
reaction.upper_bound = 1000.  
reaction.add_metabolites({
	H: -2.0,
	no3: -1.0,
	mql8: -1.0,
	H2O: 1.0,
	no2:1.0,
	mqn8: 1.0,
	h_e: 2.0
	})
reaction.reaction
reaction.gene_reaction_rule = '(A8G16_RS23470 or A8G16_RS26740)' 
model.add_reactions([reaction])


reaction = Reaction('FRD2')
reaction.name = 'Fumarate reductase' 
reaction.subsystem = 'Electron Transport Chain'
reaction.lower_bound = -1000. #Reversible EC 1.3.99.1 se remplaza por 1.3.5.1 CHECK
reaction.upper_bound = 1000.  
reaction.add_metabolites({
	fum: -1.0,
	mql8: -1.0,
	mqn8: 1.0,
	succ_c: 1.0,
	})
reaction.reaction
reaction.gene_reaction_rule = '(A8G16_RS11065 or A8G16_RS11070 or A8G16_RS27415)' 
model.add_reactions([reaction])

reaction = Reaction('ATPS4r')
reaction.name = 'ATP synthase (four protons for one ATP)' 
reaction.subsystem = 'Electron Transport Chain'
reaction.lower_bound = 0. #Irreversible 
reaction.upper_bound = 1000.  
reaction.add_metabolites({
	ADP: -1.0,
	pi_2: -1.0,
	h_e: -4.0,
	ATP: 1.0,
    H:3.0,
    H2O:1.0,
	})
reaction.reaction
reaction.gene_reaction_rule = '(XX or YY or ZZ)' 
model.add_reactions([reaction])

#intermediate reactions

reaction = Reaction('PDH')
reaction.name = 'Pyruvate dehydrogenase' 
reaction.subsystem = 'Intermediate reactions'
reaction.lower_bound = 0. #Irreversible EC  1.4.1.3
reaction.upper_bound = 1000. 
reaction.add_metabolites({
	nad: -1.0,
	coa: -1.0,
	pyr: -1.0,
	nadh: 1.0,
  CO2:1.0,
  accoa:1.0,
	})
reaction.reaction
reaction.gene_reaction_rule = '(A8G16_RS04155 or A8G16_RS35685 or A8G16_RS39220)' 
model.add_reactions([reaction])

reaction = Reaction('GLUDy')
reaction.name = 'Glutamate dehydrogenase (NADP)' 
reaction.subsystem = 'Intermediate reactions'
reaction.lower_bound = -1000. #Irreversible EC  1.2.4.1
reaction.upper_bound = 1000. 
reaction.add_metabolites({
	nadp: -1.0,
	H2O: -1.0,
	glu__L: -1.0,
	nadph: 1.0,
    akg:1.0,
    H:1.0,
	nh4_c: 1.0,
	})
reaction.reaction
reaction.gene_reaction_rule = '(A8G16_RS06845)' 
model.add_reactions([reaction])

reaction = Reaction('GLUN')
reaction.name = 'Glutaminase' 
reaction.subsystem = 'Intermediate reactions'
reaction.lower_bound = 0. #Irreversible EC  3.5.1.2
reaction.upper_bound = 1000. 
reaction.add_metabolites({
    gln__L: -1.0,
    H2O: -1.0,
    glu__L: 1.0,
    nh4_c: 1.0,
    })
reaction.reaction
reaction.gene_reaction_rule = '(A8G16_RS38985)' 
model.add_reactions([reaction])

reaction = Reaction('GLUSy')
reaction.name = 'Glutamate synthase (NADPH)' 
reaction.subsystem = 'Intermediate reactions'
reaction.lower_bound = 0. #Irreversible EC  1.4.1.13
reaction.upper_bound = 1000. 
reaction.add_metabolites({
    akg: -1.0,
    gln__L: -1.0,
    H: -1.0,
    nadph: -1.0,
    glu__L: 2.0,
    nadp: 1.0,
    })
reaction.reaction
reaction.gene_reaction_rule = '(A8G16_RS14890 or A8G16_RS14895 or A8G16_RS27210 or A8G16_RS27215)' 
model.add_reactions([reaction])

reaction = Reaction('GLNS')
reaction.name = 'Glutamine synthetase' 
reaction.subsystem = 'Intermediate reactions'
reaction.lower_bound = 0. #Irreversible EC  6.3.1.2
reaction.upper_bound = 1000. 
reaction.add_metabolites({
    ATP: -1.0,
    glu__L: -1.0,
    nh4_c: -1.0,
    gln__L: 1.0,
    ADP: 1.0,
    H: 1.0,
    pi_2: 1.0,
    })
reaction.reaction
reaction.gene_reaction_rule = '(A8G16_RS12745 or A8G16_RS22745 or A8G16_RS38440)' 
model.add_reactions([reaction])

reaction = Reaction('PPS')
reaction.name = 'Phosphoenolpyruvate synthase' 
reaction.subsystem = 'Intermediate reactions'
reaction.lower_bound = -1000. #reversible EC  2.7.9.2
reaction.upper_bound = 1000. 
reaction.add_metabolites({
    ATP: -1.0,
    H2O: -1.0,
    pyr: -1.0,
    pep: 1.0,
    ADP: 1.0,
    amp_c: 1.0,
    pi_2: 1.0,
    })
reaction.reaction
reaction.gene_reaction_rule = '(A8G16_RS12745 or A8G16_RS22745 or A8G16_RS38440)' 
model.add_reactions([reaction])


# Biomass

reaction = Reaction('Biomass')
reaction.name = 'Biomass_R_Opacus_43205_core_w_GAM' 
reaction.subsystem = 'Biomass'
reaction.lower_bound = 0. #Irreversible
reaction.upper_bound = 1000. 
reaction.add_metabolites({
	PGA: -1.21,
	accoa: -3.02,
	ATP:  -45.40,
	e4p: - 0.29,
	f6p: - 0.06,
	g3p: - 0.1,
	g6p_B: - 0.17,
	gln__L: 0,#-0.21,
	glu__L: 0,#- 5.01,
	H2O: -22.21,
	nad: -3.54,
	nadph: -19.65,
	oaa: -1.44,
	pep: -0.42,
	pyr: -2.29,
	r5p: -0.72,
	ADP: 45.40,
	akg: 3.32,
	coa:3.02,
	H: 29.29,
	nadh: 3.54,
	nadp: 19.65,
	pi_2: 45.40,
	})
reaction.reaction
reaction.gene_reaction_rule = '(GenBiomass)' 
model.add_reactions([reaction])

model.reactions

"""## Exchanges, demands and sinks"""

model.add_boundary(model.metabolites.get_by_id("CO2e"), type="exchange")
model.add_boundary(model.metabolites.get_by_id("o2e") , type = "exchange")
model.add_boundary(model.metabolites.get_by_id("h_e") , type = "exchange")
model.add_boundary(model.metabolites.get_by_id("H2Oe") , type = "exchange")
model.add_boundary(model.metabolites.get_by_id("glc__D_e") , type = "exchange")
model.add_boundary(model.metabolites.get_by_id("meoh_e") , type = "exchange")
model.add_boundary(model.metabolites.get_by_id("pyr_e") , type = "exchange")
model.add_boundary(model.metabolites.get_by_id("oaa_e") , type = "exchange")
model.add_boundary(model.metabolites.get_by_id("nh4_e") , type = "exchange")
model.add_boundary(model.metabolites.get_by_id("succ_e") , type = "exchange")
model.add_boundary(model.metabolites.get_by_id("akg_e") , type = "exchange")
model.add_boundary(model.metabolites.get_by_id("glu__L_e") , type = "exchange")
model.add_boundary(model.metabolites.get_by_id("gln__L_e") , type = "exchange")
#model.add_boundary(model.metabolites.get_by_id(""))
#model.add_boundary(model.metabolites.get_by_id("PHB_c") , type = "sink")
#model.add_boundary(model.metabolites.get_by_id("TAG_c") , type = "sink")

#Reacciones de intercambio

reaction = Reaction('O2t')
reaction.name = 'O2 transport diffusion' 
reaction.subsystem = 'exchange'
reaction.lower_bound = -1000. #Reversible 
reaction.upper_bound = 1000.  
reaction.add_metabolites({
	model.metabolites.get_by_id('o2e'): -1.0,
	model.metabolites.get_by_id('o2'): 1.0,
	})
reaction.reaction
reaction.gene_reaction_rule = '(Spontaneous)' 
model.add_reactions([reaction])

# CO2 transport reaction
reaction = Reaction('CO2t')
reaction.name = 'CO2 transport diffusion' 
reaction.subsystem = 'exchange'
reaction.lower_bound = -1000. #Reversible 
reaction.upper_bound = 1000.  
reaction.add_metabolites({
	model.metabolites.get_by_id('CO2e'): -1.0,
	model.metabolites.get_by_id('CO2'): 1.0,
	})
reaction.reaction
reaction.gene_reaction_rule = '(Spontaneous)' 
model.add_reactions([reaction])

# H2O transport reaction
reaction = Reaction('H2Ot')
reaction.name = 'H2O transport diffusion' 
reaction.subsystem = 'exchange'
reaction.lower_bound = -1000. #Reversible 
reaction.upper_bound = 1000.  
reaction.add_metabolites({
	model.metabolites.get_by_id('H2Oe'): -1.0,
	model.metabolites.get_by_id('H2O'): 1.0,
	})
reaction.reaction
reaction.gene_reaction_rule = '(Spontaneous)' 
model.add_reactions([reaction])


# Glucose transport reaction by GLC PTS
reaction = Reaction('GLCpts')
reaction.name = 'D-glucose transport via PEP-Pyr-PTS' 
reaction.subsystem = 'exchange'
reaction.lower_bound = 0. #Reversible 
reaction.upper_bound = 1000.  
reaction.add_metabolites({
	model.metabolites.get_by_id('pep'): -1.0,
	model.metabolites.get_by_id('glc__D_e'): -1.0,
    model.metabolites.get_by_id('g6p_B'): 1,
    model.metabolites.get_by_id('pyr'): 1,
	})
reaction.reaction
reaction.gene_reaction_rule = '(Arreglar)' 
model.add_reactions([reaction])

# Methanol transport by diffusion
reaction = Reaction('MeOHt')
reaction.name = 'Methanol reversible transport via diffusion' 
reaction.subsystem = 'exchange'
reaction.lower_bound = 0. #Reversible 
reaction.upper_bound = 1000.  
reaction.add_metabolites({
	model.metabolites.get_by_id('meoh_e'): -1.0,
	model.metabolites.get_by_id('meoh'): 1.0,
	})
reaction.reaction
reaction.gene_reaction_rule = '(Spontaneous)' 
model.add_reactions([reaction])

# Pyruvate transport by simport
reaction = Reaction('PYRt2')
reaction.name = 'Pyruvate transport in via proton symport' 
reaction.subsystem = 'exchange'
reaction.lower_bound = -1000. #Reversible 
reaction.upper_bound = 1000.  
reaction.add_metabolites({
	model.metabolites.get_by_id('h_e'): -1.0,
	model.metabolites.get_by_id('pyr_e'): -1.0,
  model.metabolites.get_by_id('H'): 1.0,
  model.metabolites.get_by_id('pyr'): 1.0
	})
reaction.reaction
reaction.gene_reaction_rule = '(Spontaneous)' 
model.add_reactions([reaction])

reaction = Reaction('OAAt')
reaction.name = 'Oxaloacetate transport' 
reaction.subsystem = 'exchange'
reaction.lower_bound = -1000. #Reversible 
reaction.upper_bound = 1000.  
reaction.add_metabolites({
	model.metabolites.get_by_id('oaa'): -1.0,
	model.metabolites.get_by_id('oaa_e'): 1.0,
	})
reaction.reaction
reaction.gene_reaction_rule = '(Spontaneous)' 
model.add_reactions([reaction])

reaction = Reaction('NH4t')
reaction.name = 'Ammonia reversible transport' 
reaction.subsystem = 'exchange'
reaction.lower_bound = -1000. #Reversible 
reaction.upper_bound = 1000.  
reaction.add_metabolites({
	model.metabolites.get_by_id('nh4_e'): -1.0,
	model.metabolites.get_by_id('nh4_c'): 1.0,
	})
reaction.reaction
reaction.gene_reaction_rule = '(Spontaneous)' 
model.add_reactions([reaction])

reaction = Reaction('SUCCt3')
reaction.name = 'Succinate transport out via proton antiport' 
reaction.subsystem = 'exchange'
reaction.lower_bound = 0. #Reversible 
reaction.upper_bound = 1000.  
reaction.add_metabolites({
    model.metabolites.get_by_id('succ_e'): 1.0,
        model.metabolites.get_by_id('h_e'): 1.0,
    model.metabolites.get_by_id('H'): -1.0,
    model.metabolites.get_by_id('succ_c'): -1.0,
    })
reaction.reaction
reaction.gene_reaction_rule = '(Spontaneous)' 
model.add_reactions([reaction])

reaction = Reaction('SUCCt2_2')
reaction.name = 'Succinate transport via proton symport (2 H)' 
reaction.subsystem = 'exchange'
reaction.lower_bound = 0. #Reversible 
reaction.upper_bound = 1000.  
reaction.add_metabolites({
    model.metabolites.get_by_id('succ_e'): -1.0,
    model.metabolites.get_by_id('h_e'): -2.0,
    model.metabolites.get_by_id('H'): 2.0,
    model.metabolites.get_by_id('succ_c'): 1.0,
    })
reaction.reaction
reaction.gene_reaction_rule = '(Spontaneous)' 
model.add_reactions([reaction])

reaction = Reaction('AKGt2r')
reaction.name = '2 oxoglutarate reversible transport via symport' 
reaction.subsystem = 'exchange'
reaction.lower_bound = -1000. #Reversible 
reaction.upper_bound = 1000.  
reaction.add_metabolites({
    model.metabolites.get_by_id('akg_e'): -1.0,
    model.metabolites.get_by_id('h_e'): -1.0,
    model.metabolites.get_by_id('H'): -1.0,
    model.metabolites.get_by_id('akg'): 1.0,
    })
reaction.reaction
reaction.gene_reaction_rule = '(Spontaneous)' 
model.add_reactions([reaction])

reaction = Reaction('GLUt2r')
reaction.name = 'L glutamate transport via proton symport reversible' 
reaction.subsystem = 'exchange'
reaction.lower_bound = -1000. #Reversible 
reaction.upper_bound = 1000.  
reaction.add_metabolites({
    model.metabolites.get_by_id('glu__L_e'): -1.0,
    model.metabolites.get_by_id('h_e'): -1.0,
    model.metabolites.get_by_id('H'): -1.0,
    model.metabolites.get_by_id('glu__L'): 1.0,
    })
reaction.reaction
reaction.gene_reaction_rule = '(Spontaneous)' 
model.add_reactions([reaction])

reaction = Reaction('GLNabc')
reaction.name = 'L-glutamine transport via ABC system' 
reaction.subsystem = 'exchange'
reaction.lower_bound = 0. #Reversible 
reaction.upper_bound = 1000.  
reaction.add_metabolites({
    model.metabolites.get_by_id('gln__L_e'): -1.0,
    model.metabolites.get_by_id('ATP'): -1.0,
    model.metabolites.get_by_id('H2O'): -1.0,
    model.metabolites.get_by_id('gln__L'): 1.0,
    model.metabolites.get_by_id('ADP'): 1.0,
    model.metabolites.get_by_id('H'): 1.0,
    model.metabolites.get_by_id('pi_2'): 1.0,
    })
reaction.reaction
reaction.gene_reaction_rule = '(Spontaneous)' 
model.add_reactions([reaction])

#Conviertiendo en Excel

writer = pd.ExcelWriter('Tesis_R.Opacus_FromPython.xlsx')

lista_Metabolitos_id = []
lista_Metabolitos_name = []
lista_Metabolitos_formula = []
lista_Metabolitos_carga = []

for x in model.metabolites:
    lista_Metabolitos_id.append(x.id)
    lista_Metabolitos_name.append(x.name)
    lista_Metabolitos_formula.append(x.formula)
    lista_Metabolitos_carga.append(x.charge)

metabolitos_id = pd.DataFrame(lista_Metabolitos_id, columns  = ["id_BiGG"])
metabolitos_name = pd.DataFrame(lista_Metabolitos_name, columns  = ["Nombre"])
metabolitos_formula = pd.DataFrame(lista_Metabolitos_formula, columns  = ["formula"])
metabolitos_carga = pd.DataFrame(lista_Metabolitos_carga, columns  = ["carga"])
metabolitos = pd.concat([metabolitos_id,metabolitos_name,metabolitos_formula,metabolitos_carga], axis =1)
metabolitos.to_excel( writer, index = False, sheet_name = "Metabolitos")


lista_reacciones_id = []
lista_reacciones_nombre = []
lista_reaccion = []
lista_reacciones_gen = []
lista_subsistema = []

for x in model.reactions:
	lista_reacciones_id.append(x.id)
	lista_reacciones_nombre.append(x.name)
	lista_reaccion.append(x.reaction)
	lista_subsistema.append(x.subsystem)

for i in model.reactions:
	lista_reacciones_gen.append(i.gene_reaction_rule)
 
reaccion_id = pd.DataFrame(lista_reacciones_id, columns  = ["Identificador BIGG Reaccion"])
reaccion_nombre = pd.DataFrame(lista_reacciones_nombre, columns  = ["Nombre Reaccion"])
reaccion_formula = pd.DataFrame(lista_reaccion, columns  = ["Reaccion"])
reaccion_subsistema = pd.DataFrame(lista_subsistema, columns =["Subsistema"])
reaccion_gen = pd.DataFrame(lista_reacciones_gen, columns  = ["Gen Asociado"])
reacciones = pd.concat([reaccion_id,reaccion_nombre,reaccion_formula,reaccion_subsistema,reaccion_gen], axis =1)
reacciones.to_excel(writer, index = False, sheet_name = "Reacciones" )

writer.save()
writer.close()



# Iterate through the the objects in the model
print("Reactions")
print("---------")
for x in model.reactions:
    print("%s : %s" % (x.id, x.reaction))

#print("")
#print("Metabolites")
#print("-----------")
#for x in model.metabolites:
#	print('%s : %s' % (x.id, x.formula))

#print("")
#print("Genes")
#print("-----")

for x in model.genes:
   associated_ids = (i.id for i in x.reactions)
   print("%s is associated with reactions: %s" %
         (x.id, "{" + ", ".join(associated_ids) + "}"))


for x in model.reactions:
	balance = x.check_mass_balance()
	if len(balance) == 0:
		print ("La reacción",x.id,"esta bien balanceada", balance)
	else:
		print("La reacción",x.id,"esta mal balanceada", balance, "revisar cargas y coeficientes")

print ("Exchanges reactions are:")
for x in model.exchanges:
	print(x.id)


print ("Demands reactions are:")
for x in model.demands:
	print(x.id)

"""
print ("Sinks reactions are:")
for x in model.sinks:
	print(x.id)
"""

cobra.io.save_json_model(model, "Rhodococcus Opacus 43205.json")
"""
from cobra.util.solver import linear_reaction_coefficientss
print(linear_reaction_coefficients(model))
"""

"""
model.objetive = 'Biomass_R_Opacus_43205_core_w_GAM'
model.reactions.get_by_id("Biomass_R_Opacus_43205_core_w_GAM").upper_bound = 1000.
linear_reaction_coefficients(model)
print(model.optimize().objective_value)
"""

model

"""# Media configuration and model solution

### Set objective function
"""

model.reactions.get_by_id( "Biomass").objective_coefficient = 0
model.reactions.get_by_id( "GLCpts").objective_coefficient = 0
model.reactions.get_by_id( "PYRt2").objective_coefficient = 0
model.reactions.get_by_id( "OAAt").objective_coefficient = 0
model.reactions.get_by_id( "EX_pyr_e").objective_coefficient = 0

"""### Set uptake bounds"""

model.reactions.get_by_id("EX_glc__D_e").lower_bound= -10
model.reactions.get_by_id("EX_glc__D_e").upper_bound= -5

"""### Solve"""

FBA_sol=model.optimize()

print(FBA_sol)

"""## Solution summary"""

print(model.summary(FBA_sol))

"""## Solution report in human language"""

rxn_ids = [r.id for r in model.reactions]
sol = FBA_sol.fluxes[ (abs(FBA_sol.fluxes) > 1e-6) ]
sol_nice = sol.to_frame()
sol_nice["reaction"] = sol.index
sol_nice["LB"] = [model.reactions.get_by_id(rid).lower_bound for rid in sol.index]
sol_nice["UB"] = [model.reactions.get_by_id(rid).upper_bound for rid in sol.index]
sol_nice ["subsystem"] = [model.reactions.get_by_id(rid).subsystem for rid in sol.index]
print(sol_nice.sort_values(by="fluxes"))


import simulation_class as sc
import pytest

## Do real tests

def test_MesoRDsimulation():
	new_sim = sc.MesoRDsimulation('simulation_example',0.01)
	
def test_MesoRDsimulation_exception():
	with pytest.raises(NotADirectoryError):#, match='must be 0 or None'):
		new_sim2 = sc.MesoRDsimulation('sitin_amle',0.01)
	

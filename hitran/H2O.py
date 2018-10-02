
from hapi import *

db_begin('data')

for i in range (0,100):
	a = i*20 + 300
	if a < 100:
		b = "CH4-data/CH4_0%dK_1atm.txt" % (a)	
	else:
		b = "CH4-data/CH4_%dK_1atm.txt" % (a)
	print "doing temperature: %d" % (a)
	
	nu,coef = absorptionCoefficient_Lorentz(SourceTables='CH4',Environment={'T':a,'p':1},GammaL='gamma_self',HITRAN_units=False,File=b)

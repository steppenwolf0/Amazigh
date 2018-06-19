import os
import math
import random

init = 30
stop = 40
for x in range (init,stop):
	var='gnome-terminal -x sh -c "'+' ./bayesOpt '+str(x)+'"'
	print(var)
	os.system(var)

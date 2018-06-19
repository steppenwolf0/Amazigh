import os
import math
import random

generation=0
for x in range (0,10):
	var='oarsub \' ./bayesOpt '+str(x)+'\''
	#var='sh gen2.sh 100 '+str(w1)+' '+str(w2)+' '+str(w3)+' '+str(w4)+' '+str(h1)+' '+str(h2)+' '+str(h3)+' '+str(wd1)+' '+str(wd2)+' '+str(wd3)+' '+str(x)+' '+str(generation)
	print(var)
	os.system(var)

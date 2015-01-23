from pylab import *

expression = loadtxt("normalised_gene_expression")
expression = expression
#print expression.mean(axis=1)[:,newaxis]
#print std(expression, axis=1)
#print [expression.mean(axis=1)[:,newaxis]]
#print [std(expression, axis=1)]

savetxt('mean_std', vstack(([std(expression, axis=1)], [mean(expression, axis=1)])).T)




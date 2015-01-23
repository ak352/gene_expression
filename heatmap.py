import scipy
import scipy.cluster.hierarchy as hier
from os.path import expanduser
from pylab import *
import sys

num_genes = int(sys.argv[1])

def red_black_green():
    cdict = {
       'green': ((0.0, 0.0, 0.0),
               (0.5, 0.0, 0.0),
               (1.0, 1.0, 1.0)),
       'blue': ((0.0, 0.0, 0.0),
                (1.0, 0.0, 0.0)),
       'red': ((0.0, 0.0, 1.0),
                 (0.5, 0.0, 0.0),
                 (1.0, 0.0, 0.0))
       }

    my_cmap = matplotlib.colors.LinearSegmentedColormap(
        'my_colormap', cdict, 100)

    return my_cmap




def fcluster( pts, ncluster, method="complete", criterion="maxclust" ):
    """ -> (pts, Y pdist, Z linkage, T fcluster, clusterlists)   
        ncluster = n1 + n2 + ... (including n1 singletons)
        av cluster size = len(pts) / ncluster
    """
    pts = np.asarray(pts)
    Y = scipy.spatial.distance.pdist( pts , 'euclidean')  # ~ N^2 / 2
    #Y = scipy.spatial.distance.pdist( pts , 'correlation')  # ~ N^2 / 2
    Z = hier.linkage( Y, method )  # N-1                         
    T = hier.fcluster( Z, ncluster, criterion=criterion )
        # clusters = clusterlists(T)
    return (pts, Y, Z, T)

#infile = "%s/Downloads/E-MTAB-264/CHB_p3_expression.txt" % expanduser("~")
infile = "GDS3268.soft"
expressions  = []
probes = 0
num_probes_excluded = 0
with open(infile) as f:
    next(f)
    next(f)
    for line in f:
        probes += 1
        line = line[:-1].split("\t")
        try:
            measured = [float(k) for k in line[1:] if k !='null']
            if not measured:
                num_probes_excluded += 1
                continue
            avg_exp = median(measured)
            #print "average expression = ", avg_exp
            expressions.append([float(k) if k!='null' else avg_exp for k in line[1:]])
        except ValueError:
            print line
            
            sys.exit()
        if probes == num_genes:
            break
expr_array = array(expressions)

#Normalise gene expression
#row_sum = expr_array.mean(axis=1)
#expr_array = expr_array - row_sum[:, newaxis]

savetxt('normalised_gene_expression', expr_array)

    
print infile
print "Number of probes = %d" % len(expressions)
print "Number of samples = %d" % len(expressions[0])
print "Number of probes excluded = %d" % num_probes_excluded

#plot_data = ma.masked_equal(expr_array, 0)
#imshow(expr_array) #, cmap=cm.get_cmap("Reds"), interpolation="nearest")
#pcolor(expr_array)
#colorbar(im)
#show()

pts,Y,Z,T = fcluster(expr_array, 5)
R = hier.dendrogram( Z , orientation='right')
idxs = [int(k) for k in R["ivl"]]
#print i
a2 = expr_array[idxs,:]

figure()
subplot(1,2,1)
X,Y = meshgrid(range(expr_array.shape[1]+1), range(expr_array.shape[0]+1))
im2 = pcolormesh(X,Y,expr_array, cmap='hot') #cm.get_cmap("Reds"))
colorbar(im2)
title('Unsorted')

#figure()
subplot(1,2,2)
X,Y = meshgrid(range(a2.shape[1]+1), range(a2.shape[0]+1))
im = pcolormesh(X,Y,a2, cmap='hot') #cm.get_cmap("Reds"))
colorbar(im)
title('Sorted')


## Clustering by sample too.
a2 = a2.T
pts1,Y1,Z1,T1 = fcluster(a2, 2)
figure()
R1 = hier.dendrogram( Z1 , orientation='right')
idxs1 = [int(k) for k in R1["ivl"]]
a3 = a2[idxs1, :]
a3 = a3.T

fig = figure(figsize=(8,8))
axd = fig.add_axes([0.09,0.1,0.2,0.6])
R = hier.dendrogram( Z , orientation='right')
axd.set_xticks([])
axd.set_yticks([])
axf = fig.add_axes([0.29, 0.71, 0.6, 0.2])
R1 = hier.dendrogram(Z1)
axf.set_xticks([])
axf.set_yticks([])

axm = fig.add_axes([0.29,0.1,0.6,0.6]) # x-pos, y-pos, width, height
X,Y = meshgrid(range(a3.shape[1]+1), range(a3.shape[0]+1))
im = pcolormesh(X,Y,a3, cmap=red_black_green()) #'hot') #cm.get_cmap("Reds"))
axis('tight')
#fig.colorbar(im)
axm.get_yaxis().set_visible(False)
#title('Sorted')
#axm.set_xticklabels([''] + list(df_rowclust.columns))
#axm.set_yticklabels([''] + list(df_rowclust.index))

show()


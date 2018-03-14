from sklearn.cluster import DBSCAN
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from sklearn import metrics
from sklearn.datasets.samples_generator import make_blobs
from sklearn.preprocessing import StandardScaler
import numpy as np
import plotly
from plotly.graph_objs import  Layout,Scatter3d
import plotly.graph_objs as go
import seaborn as sns; sns.set()
def pairCorrelationFunction_3D(x, y, z, S, rMax, dr):
    """Compute the three-dimensional pair correlation function for a set of
    spherical particles contained in a cube with side length S.  This simple
    function finds reference particles such that a sphere of radius rMax drawn
    around the particle will fit entirely within the cube, eliminating the need
    to compensate for edge effects.  If no such particles exist, an error is
    returned.  Try a smaller rMax...or write some code to handle edge effects! ;)
    Arguments:
        x               an array of x positions of centers of particles
        y               an array of y positions of centers of particles
        z               an array of z positions of centers of particles
        S               length of each side of the cube in space
        rMax            outer diameter of largest spherical shell
        dr              increment for increasing radius of spherical shell
    Returns a tuple: (g, radii, interior_indices)
        g(r)            a numpy array containing the correlation function g(r)
        radii           a numpy array containing the radii of the
                        spherical shells used to compute g(r)
        reference_indices   indices of reference particles
    """
    from numpy import zeros, sqrt, where, pi, mean, arange, histogram

    # Find particles which are close enough to the cube center that a sphere of radius
    # rMax will not cross any face of the cube
    #print('Checks: ')
    #print('    x   |   y   |   z   |  S  | rMax |  dr')
    #print(x[0] , y[0] , z[0] ,S, rMax,dr)
    bools1 = abs(x) > rMax
    bools2 = abs(x) < (S - rMax)
    bools3 = abs(y) > rMax
    bools4 = abs(y) < (S - rMax)
    bools5 = abs(z) > rMax
    bools6 = abs(z) < (S - rMax)

    interior_indices, = where(bools1 * bools2 * bools3 * bools4 * bools5 * bools6)
    num_interior_particles = len(interior_indices)

    if num_interior_particles < 1:
        raise  RuntimeError ("No particles found for which a sphere of radius rMax\
                will lie entirely within a cube of side length S.  Decrease rMax\
                or increase the size of the cube.")

    edges = arange(0., rMax + 1.1 * dr, dr)
    num_increments = len(edges) - 1
    g = zeros([num_interior_particles, num_increments])
    radii = zeros(num_increments)
    numberDensity = len(x) / S**3

    # Compute pairwise correlation for each interior particle
    for p in range(num_interior_particles):
        index = interior_indices[p]
        d = sqrt((x[index] - x)**2 + (y[index] - y)**2 + (z[index] - z)**2)
        d[index] = 2 * rMax

        (result, bins) = histogram(d, bins=edges, normed=False)
        g[p,:] = result / numberDensity

    # Average g(r) for all interior particles and compute radii
    g_average = zeros(num_increments)
    for i in range(num_increments):
        radii[i] = (edges[i] + edges[i+1]) / 2.
        rOuter = edges[i + 1]
        rInner = edges[i]
        g_average[i] = mean(g[:, i]) / (4.0 / 3.0 * pi * (rOuter**3 - rInner**3))

    return (g_average, radii, interior_indices)
    # Number of particles in shell/total number of particles/volume of shell/number density
    # shell volume = 4/3*pi(r_outer**3-r_inner**3)


def set_mol_info():
    pass
    return

def dbscan_w_plot(all_pos):
    X = all_pos
    #= StandardScaler().fit_transform(all_pos)

    # #############################################################################
    # Compute DBSCAN

    db = DBSCAN( eps=0.8 ,min_samples = 4 ).fit(X)
    core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
    core_samples_mask[db.core_sample_indices_] = True
    labels = db.labels_

    # Number of clusters in labels, ignoring noise if present.
    n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)

    print('Estimated number of clusters: %d' % n_clusters_)
    print("Silhouette Coefficient: %0.3f" %
          metrics.silhouette_score(X, labels))

    # #############################################################################
    # Plot result

    new_xyz = [] #np.zeros((X.shape[0],3))
    new_cols = [] #np.zeros((X.shape[0], 4))
    new_labels = []
    new_op = []

    fig = plt.figure()
    ax = Axes3D(fig)
    # Black removed and is used for noise instead.
    plotly_data = []
    unique_labels = set(labels)
    colors = [plt.cm.Spectral(each) for each in np.linspace(0, 1, len(unique_labels))]

    #print("Color variety: ", len(colors))
    #print(colors)

    for k, col in zip(unique_labels, colors):
        if k == -1:
            # Black used for noise.
            col = [0, 0, 0, 0.3]

        class_member_mask = (labels == k)

        xy = X[class_member_mask & core_samples_mask]

        for tmp_xy in xy:
            new_xyz.append(tmp_xy)
            #tmp_color = list(col)
            #tmp_color = [tmp_color[0],tmp_color[1],tmp_color[2],0.3]
            new_cols.append(tuple(col))
            if k == -1:
                new_labels.append('Unclustered')
                new_op.append(float(0.3))
            else:
                new_labels.append('Cluster: '+ str(k))
                new_op.append(float(1.0))


        ax.scatter(xy[:, 0], xy[:, 1],xy[:, 2], 'o', c=tuple(col),
                   edgecolors='k', s=14)


        xy = X[class_member_mask & ~core_samples_mask]

        for tmp_xy in xy:
            new_xyz.append(tmp_xy)
            #col = tuple(col[0],col[1],col[2],0.3)
            new_cols.append(tuple(col))
            if k == -1:
                new_labels.append('Unclustered')
                new_op.append(float(0.3))
            else:
                new_labels.append('Cluster: '+ str(k))
                new_op.append(float(1.0))

        ax.scatter(xy[:, 0], xy[:, 1], xy[:, 2], 'o', c=tuple(col),
                   edgecolors='k', s=6)

    new_xyz = np.asarray(new_xyz)
    print(len(new_op))

    #names = set_mol_info(X)
    remove_idx = []
    for lbl_tmp, opa_tmp, col_tmp, xyz_tmp, idx_tmp in zip(new_labels,new_op,new_cols,new_xyz, range(len(new_labels))):
        if lbl_tmp != 'Unclustered':
            remove_idx.append(idx_tmp)

    new_labels = new_labels.pop(for r in remove_idx)
    new_xyz =
    new_cols =
    new_op =

    plotly.offline.plot({
        "data": [Scatter3d(x=new_xyz[:, 0], y=new_xyz[:, 1], z=new_xyz[:, 2],
                           marker=dict(color=new_cols,line = dict(width=1)),
                           mode='markers',text =new_labels , opacity = dict(opacity = new_op) )],  #
        #"data" : plotly_data ,
        "layout": Layout(title="Active protein water molecules heatmap")
    })

    plt.title('Estimated number of clusters: %d' % n_clusters_)


    ########################################################################################


def distance_matrix(molecule_pop_lst_active, all_pos):
    one_prot = np.zeros(int(molecule_pop_lst_active[0]))
    for i in range(int(molecule_pop_lst_active[0])):
        for j in range(int(molecule_pop_lst_active[0])):
            if (i < j):
                a = all_pos[i, :]
                b = all_pos[j, :]
                one_prot[j][i] = np.linalg.norm(a - b)
from xlrd import open_workbook
from scipy.interpolate import spline
import plotly
from plotly.graph_objs import  Layout,Scatter3d
from mpl_toolkits.mplot3d import proj3d
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from sklearn.cluster import DBSCAN
from sklearn import metrics
from sklearn.datasets.samples_generator import make_blobs
from sklearn.preprocessing import StandardScaler
#from mayavi import mlab
from sklearn.neighbors import KernelDensity
from sklearn.decomposition import PCA
from sklearn.cross_decomposition import PLSRegression
import scipy as sp
from scipy import stats
from matplotlib.widgets import CheckButtons
import matplotlib.pyplot as plt
#import seaborn as sns; sns.set()
import matplotlib as mpl
#mpl.style.use('classic')

#from sheatheer_jones import *

from rdf_pdf import * #pairCorrelationFunction_3D

wb = open_workbook('watermap_statistical_xyz_analysis.xlsx')



def create_rdf(molecule_pop_lst, all_pos, domain_size, dr, particle_radius, rMax, dr_set, line_smoothing):
    g_r_lst = []
    r_lst = []
    for i in range(len(molecule_pop_lst)):
        if i == 0:
            g_r, r, reference_indices = pairCorrelationFunction_3D(all_pos[:int(molecule_pop_lst[i] - 1), 0],
                                                                   all_pos[:int(molecule_pop_lst[i] - 1), 1],
                                                                   all_pos[:int(molecule_pop_lst[i] - 1), 2],
                                                                   domain_size, rMax, dr)
            if line_smoothing:
                xnew = np.linspace(r.min(), r.max(), 300)  #
                g_r_lst.append(spline(r, g_r, xnew))
                r_lst.append(xnew)
            else:
                g_r_lst.append(g_r)
                r_lst.append(r)

        else:
            # print(int(molecule_pop_lst[i - 1]),' -+- ',int(molecule_pop_lst[i] - 1))
            g_r, r, reference_indices = pairCorrelationFunction_3D(
                all_pos[int(molecule_pop_lst[i - 1]):int(molecule_pop_lst[i] + molecule_pop_lst[i - 1] - 1), 0],
                all_pos[int(molecule_pop_lst[i - 1]):int(molecule_pop_lst[i] + molecule_pop_lst[i - 1] - 1), 1],
                all_pos[int(molecule_pop_lst[i - 1]):int(molecule_pop_lst[i] + molecule_pop_lst[i - 1] - 1), 2],
                domain_size, rMax, dr_set)
            if line_smoothing:
                xnew = np.linspace(r.min(), r.max(), 300)  #
                g_r_lst.append(spline(r, g_r, xnew))
                r_lst.append(xnew)
            else:
                g_r_lst.append(g_r)
                r_lst.append(r)

    return g_r_lst, r_lst

def find_average(g_r_lst, r_lst, names_lst , minimax = True):
    # Discard min &max

    factor = len(names_lst)
    avg = [0 for u in range(len(r_lst[0]))]
    #print('RLIST_RANGE: ', len(r_lst[0]))
    for i in range(len(r_lst[0])):
        # for j in range(factor):
        tmp = []
        for row in g_r_lst:
            #print(row[i])
            tmp.append(float(row[i]))
        # print(tmp)
        if(minimax):
            tmp.remove(max(tmp))
            tmp.remove(min(tmp))
        #print(i)
        avg[i] = sum(tmp) / float(len(tmp))

    return r_lst[0], avg

class w_mol:

    def __init__(self, x,y,z,prot_name, _id , dh ,tds , dg):
        self.x = x
        self.y = y
        self.z = z
        self._id = _id
        self.prot_name = prot_name
        self.dh = dh
        self.tds = tds
        self.dg = dg


colum_offset = 12
row_offset = 0 #  mol_name
mol_num_pos = 1 #  mol_total

sheet_idx = 0
names_lst_active = []
molecule_pop_lst_active = []
mol_pos_active = []
mol_obj_active = []

names_lst_inactive = []
molecule_pop_lst_inactive = []
mol_pos_inactive = []
mol_obj_inactive = []

for sheet in wb.sheets():
    print(wb.sheet_names())
    if sheet_idx == 0:
        number_of_rows = sheet.nrows
        number_of_columns = sheet.ncols


        items = []

        rows = []
        tmp_buff = []
        mol_idx = - 1
        for col in range(number_of_columns):
            tmp_arr = np.zeros(3)
            if (col % 13 == 0):
                molecule_name = (sheet.cell(0, col).value)
                molecule_pop = (sheet.cell(1, col).value)
                print(molecule_name,'   -    ', molecule_pop)
                names_lst_active.append(molecule_name)
                molecule_pop_lst_active.append(molecule_pop)
                mol_pos_active.append([])
                mol_idx += 1


            # xyz - coordinate
            if (col % 13 == 1):
                for row in range(1, int(molecule_pop_lst_active[-1] + 1)):

                    x_value = (sheet.cell(row, col).value)
                    y_value = (sheet.cell(row, col + 1).value)
                    z_value = (sheet.cell(row, col + 2).value)


                    #print('x:',x_value,' y:',y_value,' z:',z_value)
                    tmp_arr[0] = float(x_value)
                    tmp_arr[1] = float(y_value)
                    tmp_arr[2] = float(z_value)
                    dh_val = (sheet.cell(row, col + 5).value)
                    tds_val = (sheet.cell(row, col + 6).value)
                    dg_val = (sheet.cell(row, col + 7).value)

                    mol_obj_active.append(w_mol(x_value,y_value,z_value,names_lst_active[mol_idx],names_lst_active[mol_idx] +'_'+ str(row),\
                                                float(dh_val),float(tds_val),float(dg_val)))
                    #print(tmp_arr)
                    tmp_buff.append(tmp_arr)
                    mol_pos_active[mol_idx].append(tmp_arr)

        sheet_idx += 1

    else:
        number_of_rows = sheet.nrows
        number_of_columns = sheet.ncols

        items = []

        rows = []
        tmp_buff = []
        mol_idx = - 1
        for col in range(number_of_columns):
            tmp_arr = np.zeros(3)
            if (col % 13 == 0):
                molecule_name = (sheet.cell(0, col).value)
                molecule_pop = (sheet.cell(1, col).value)
                print(molecule_name, '   -    ', molecule_pop)
                names_lst_inactive.append(molecule_name)
                molecule_pop_lst_inactive.append(molecule_pop)
                mol_pos_inactive.append([])
                mol_idx += 1

            # xyz - coordinate
            if (col % 13 == 1):
                for row in range(1, int(molecule_pop_lst_inactive[-1] + 1)):
                    x_value = (sheet.cell(row, col).value)
                    y_value = (sheet.cell(row, col + 1).value)
                    z_value = (sheet.cell(row, col + 2).value)

                    # print('x:',x_value,' y:',y_value,' z:',z_value)
                    tmp_arr[0] = float(x_value)
                    tmp_arr[1] = float(y_value)
                    tmp_arr[2] = float(z_value)

                    dh_val = (sheet.cell(row, col + 5).value)
                    tds_val = (sheet.cell(row, col + 6).value)
                    dg_val = (sheet.cell(row, col + 7).value)

                    mol_obj_inactive.append(w_mol(x_value, y_value, z_value, names_lst_inactive[mol_idx],names_lst_inactive[mol_idx] + '_' + str(row), \
                                                float(dh_val), float(tds_val), float(dg_val)))

                    tmp_buff.append(tmp_arr)
                    mol_pos_inactive[mol_idx].append(tmp_arr)


all_pos_active = np.zeros([int(len(mol_obj_active)),3])
np_names_active = np.zeros([int(len(mol_obj_active)),1])
idx = 0
for moll in mol_obj_active:

    #print (moll.x, moll.y ,moll.z)
    all_pos_active[idx][0]=  moll.x
    all_pos_active[idx][1] = moll.y
    all_pos_active[idx][2] = moll.z

    names_lst_active.append(moll._id)

    idx += 1


all_pos_inactive = np.zeros([int(len(mol_obj_inactive)),3])
np_names_inactive = np.zeros([int(len(mol_obj_inactive)),1])
idx = 0
for moll in mol_obj_inactive:

    #print (moll.x, moll.y ,moll.z)
    all_pos_inactive[idx][0] = moll.x
    all_pos_inactive[idx][1] = moll.y
    all_pos_inactive[idx][2] = moll.z

    names_lst_inactive.append(moll._id)
    idx += 1

all_pos = np.vstack((all_pos_active, all_pos_inactive))

print(np.shape(all_pos_active))
print(np.shape(all_pos_inactive))
print(np.shape(all_pos))

center = np.mean(all_pos,axis=0)

print ('Center: ', center)


#for i in range(np.shape(all_pos)[1]):






def plot_KDE_2D(mol_obj, label):
    xyz = center
    X = []

    for mol in mol_obj:
        b = np.array([mol.x, mol.y, mol.z])
        # print(b)
        dist = np.linalg.norm(xyz - b)
        X.append(dist)

    # Silverman's Rule of Thumb
    BW = 1.06 * np.power(np.std(np.asarray(X),axis=0),-0.2) / 2
    print(BW)
    for kernel in ['epanechnikov']:#, 'tophat', 'epanechnikov']:
        kde = KernelDensity(kernel=kernel, bandwidth=BW).fit(np.asarray(X).reshape(-1,1))
        X_plot = np.linspace(0, 25, 10000)[:, np.newaxis]
        log_dens = kde.score_samples(X_plot)

        plt.plot(X_plot[:, 0], np.exp(log_dens), '-',
                 label="{1} BW = '{0:.2f}'".format(BW ,label))

    plt.legend(loc='upper left')

def plot_KDE_3D(mol_obj, label , f, axarr ):




    X = []

    for mol in mol_obj:
        b = np.array([mol.dh, mol.dg, mol.tds])
        X.append(b)

    # Silverman's Rule of Thumb
    BW = 1.06 * np.power(np.std(np.asarray(X)[:,0],axis=0),-0.2) / 100


    X = np.asarray(X)
    xyz_min = np.min(X , axis= 0 )
    xyz_max = np.max(X , axis= 0 )

   # BW = hsj(X[:,0])
    print('MINIMAX Shape;', np.shape(xyz_max) , xyz_min ,xyz_max)
    print(BW)

    for kernel in ['gaussian']:#, 'tophat', 'epanechnikov']:
        kde = KernelDensity(kernel=kernel, bandwidth=BW).fit(np.asarray(X).reshape(-1,1))
        X_plot_x = np.linspace(xyz_min[0], xyz_max[0], 10000)[:, np.newaxis]
        X_plot_y = np.linspace(xyz_min[1], xyz_max[1], 10000)[:, np.newaxis]
        X_plot_z = np.linspace(xyz_min[2], xyz_max[2], 10000)[:, np.newaxis]
        X_plot = np.transpose(np.array([X_plot_x,X_plot_y,X_plot_z]))
        #print('Shape;', np.shape(X_plot))
        log_dens = kde.score_samples(X_plot.reshape(-1,1))


        #print(np.shape(log_dens))
        #plt.plot(X_plot_x, np.exp(log_dens[0:1000]), '-',
         #        label="{1} BW = '{0:.2f}'".format(BW ,label))

        #ax.scatter(X[:,0],X[:,1],X[:,2], c=log_dens[0:len(X[:,0])],
               #    label=label)

    #axarr[0].plot(X_plot_x, np.exp(log_dens[:len(X_plot_x)]))
    #axarr[1].plot(X_plot_y, np.exp(log_dens[len(X_plot_x):2*len(X_plot_x)]))
    #axarr[2].plot(X_plot_z, np.exp(log_dens[2*len(X_plot_x):3*len(X_plot_x)]))

    axarr[0].hist(X[:,0], 50, normed=1,  alpha=0.75)
    axarr[1].hist(X[:,1], 50, normed=1,  alpha=0.75)
    axarr[2].hist(X[:,2], 50, normed=1,  alpha=0.75)

    plt.legend(loc='upper left')


def PLS(mol_a, mol_i):
    X = []
    Y = []
    for mol in mol_a:
        b = np.array([mol.x - center[0], mol.y - center[1], mol.z - center[2], mol.dg, mol.dh, mol.tds])

        X.append(b)
        Y.append(np.array([1,1]))

    #Y = []
    for mol in mol_i:
        b = np.array([mol.x - center[0], mol.y - center[1], mol.z - center[2], mol.dg, mol.dh, mol.tds])

        X.append(b)
        Y.append(np.array([0, 0]))

    pls2 = PLSRegression(n_components=2)
    x_scores , y_scores = pls2.fit_transform(X,Y)

    plt.figure(10, figsize=(5, 5))
    plt.scatter(x_scores[:,0],x_scores[:,1])
    plt.scatter(y_scores[:, 0], y_scores[:, 1])
    plt.show()
    print(np.shape(x_scores))
    print(np.shape(y_scores))



#plot_KDE_2D(mol_obj_active, 'Active')
#plot_KDE_2D(mol_obj_inactive, 'Inactive')
PLS(mol_obj_active,mol_obj_inactive)
#fig = plt.figure()
plt.show()
#ax = fig.add_subplot(111, projection='3d')

f, axarr = plt.subplots(3, sharex=True)

#plot_KDE_3D(mol_obj_active, 'Active' ,f, axarr )
#plot_KDE_3D(mol_obj_inactive, 'Inactive',f, axarr )

plt.legend(loc='upper left')
plt.ylabel('Density')
plt.xlabel(r'$\AA$')
plt.title('KDE of H2O distances from center')
plt.show()

##########


X = []
for mol in mol_obj_active:
    b = np.array([mol.x - center[0], mol.y- center[1], mol.z- center[2], mol.dg, mol.dh, mol.tds])
    X.append(b)

pca = PCA(n_components=2)

X_r = pca.fit(X).transform(X)

plt.scatter(X_r[:, 0], X_r[:, 1], s= 2,c = 'b',
                 label="Active")


X = []
for mol in mol_obj_inactive:
    b = np.array([mol.x - center[0], mol.y- center[1], mol.z- center[2], mol.dg, mol.dh, mol.tds])

    X.append(b)


pca = PCA(n_components=2)

X_r = pca.fit(X).transform(X)

plt.scatter(X_r[:, 0], X_r[:, 1], s= 2,c = 'r',
                 label="Inactive")
plt.legend(loc='upper left')

plt.title('PCA w/ 2 components from 6d')
plt.show()

#################
X = []
for mol in mol_obj_inactive:
    b = np.array([ mol.dg, mol.dh, mol.tds])
    X.append(b)

X= np.asarray(X)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(X[:, 0], X[:, 1], X[:, 2])
X = []
for mol in mol_obj_active:
    b = np.array([ mol.dg, mol.dh, mol.tds])
    X.append(b)

X= np.asarray(X)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(X[:, 0], X[:, 1], X[:, 2])
plt.show()



pca = PCA(n_components=2)

X_r = pca.fit(X).transform(X)


xc, yc = np.mean(X_r[:, 0]), np.mean(X_r[:, 1])

print('Mean: ' +str(xc) +','+ str(yc))

plt.scatter(X_r[:, 0], X_r[:, 1], s= 2,#c = 'b',
                 label="Inactive")

X = []
for mol in mol_obj_inactive:
    b = np.array([ mol.dg, mol.dh, mol.tds])
    X.append(b)


pca = PCA(n_components=2)

X_r = pca.fit(X).transform(X)

plt.scatter(X_r[:, 0], X_r[:, 1], s= 2,#c = 'r',
                 label="Inactive")
xc, yc = np.mean(X_r[:, 0]), np.mean(X_r[:, 1])

print('Mean: ' +str(xc) +','+ str(yc))

plt.title('PCA w/ 2 components from 3d (energies)')
plt.legend(loc='upper left')

plt.show()





plt.legend(loc='upper left')
plt.ylabel('Density')
plt.xlabel(r'$\AA$')
plt.title('KDE of H2O distances from center')
plt.show()
"""
from Bio import *
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
structure = PDBParser().get_structure('5bmy', '5bmy.pdb')


model = structure[0]
chain = model['A']

for i in chain.get_residues():
    print(i.get_resnames())
"""
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
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
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

class molecule_set:

    active_molecules = []
    inactive_molecules = []

    def load_dataset(self):

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

                            molecule_set.active_molecules.append(w_mol(x_value,y_value,z_value,names_lst_active[mol_idx],names_lst_active[mol_idx] +'_'+ str(row),\
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

                            molecule_set.inactive_molecules.append(w_mol(x_value, y_value, z_value, names_lst_inactive[mol_idx],names_lst_inactive[mol_idx] + '_' + str(row), \
                                                        float(dh_val), float(tds_val), float(dg_val)))

                            tmp_buff.append(tmp_arr)
                            mol_pos_inactive[mol_idx].append(tmp_arr)

    def get_xyz(self, act):

        all_pos_active = np.zeros([int(len(molecule_set.active_molecules)), 3])
        np_names_active = np.zeros([int(len(molecule_set.active_molecules)), 1])
        idx = 0
        for moll in molecule_set.active_molecules:
            # print (moll.x, moll.y ,moll.z)
            all_pos_active[idx][0] = moll.x
            all_pos_active[idx][1] = moll.y
            all_pos_active[idx][2] = moll.z
            idx += 1

        all_pos_inactive = np.zeros([int(len(molecule_set.inactive_molecules)),3])
        np_names_inactive = np.zeros([int(len(molecule_set.inactive_molecules)),1])
        idx = 0
        for moll in molecule_set.inactive_molecules:

            #print (moll.x, moll.y ,moll.z)
            all_pos_inactive[idx][0] = moll.x
            all_pos_inactive[idx][1] = moll.y
            all_pos_inactive[idx][2] = moll.z
            idx += 1


        if act == 'active':
            return all_pos_active
        elif act == 'inactive':
            return all_pos_inactive
        else:
            all_pos = np.vstack((all_pos_active, all_pos_inactive))
            return all_pos

    def get_com(self):
        return np.mean((self.get_xyz('all')),axis=0)

    def plot_KDE_2D_per_mol(self, mol_obj, label):
        xyz = self.get_com()
        X = []

        for mol in mol_obj:
            b = np.array([mol.x, mol.y, mol.z])
            # print(b)
            dist = np.linalg.norm(xyz - b)
            X.append(dist)

        # Silverman's Rule of Thumb
        BW = 1.06 * np.power(np.std(np.asarray(X), axis=0), -0.2) / 2
        print(BW)
        for kernel in ['epanechnikov']:  # , 'tophat', 'epanechnikov']:
            kde = KernelDensity(kernel=kernel, bandwidth=BW).fit(np.asarray(X).reshape(-1, 1))
            X_plot = np.linspace(0, 25, 10000)[:, np.newaxis]
            log_dens = kde.score_samples(X_plot)

            plt.plot(X_plot[:, 0], np.exp(log_dens), '-',
                     label="{1} BW = '{0:.2f}'".format(BW, label))

        plt.legend(loc='upper left')

    def plot_KDE_2D(self):
        self.plot_KDE_2D_per_mol(molecule_set.active_molecules, 'Active')
        self.plot_KDE_2D_per_mol(molecule_set.inactive_molecules, 'Active')

        plt.legend(loc='upper left')
        plt.ylabel('Density')
        plt.xlabel(r'$\AA$')
        plt.title('KDE of H2O distances from center')
        plt.show()

    def plot_PCA_6D(self, active, inactive):
        X = []
        for mol in active:
            b = np.array([mol.x - center[0], mol.y - center[1], mol.z - center[2], mol.dg, mol.dh, mol.tds])
            X.append(b)

        pca = PCA(n_components=2)

        X_r = pca.fit(X).transform(X)

        plt.scatter(X_r[:, 0], X_r[:, 1], s=2, c='b',
                    label="Active ")#+str(pca.get_precision()))

        Y = []
        for mol in inactive:
            b = np.array([mol.x - center[0], mol.y - center[1], mol.z - center[2], mol.dg, mol.dh, mol.tds])

            Y.append(b)

        pca = PCA(n_components=2)

        Y_r = pca.fit(Y).transform(Y)

        plt.scatter(Y_r[:, 0], Y_r[:, 1], s=2, c='r',
                    label="Inactive ")#+str(pca.get_precision()))
        plt.legend(loc='upper left')

        plt.title('PCA w/ 2 components from 6d')
        print(pca.get_precision())
        plt.show()

        return X_r , Y_r

    def num_PCA_6D(self):
        X = []
        Y = []
        for mol in molecule_set.active_molecules:
            b = np.array([mol.x - center[0], mol.y - center[1], mol.z - center[2], mol.dg, mol.dh, mol.tds])
            X.append(b)

        pca = PCA(n_components=2)
        X_r = pca.fit(X).transform(X)

        for mol in molecule_set.inactive_molecules:
            b = np.array([mol.x - center[0], mol.y - center[1], mol.z - center[2], mol.dg, mol.dh, mol.tds])
            Y.append(b)

        pca = PCA(n_components=2)
        Y_r = pca.fit(Y).transform(Y)

        return X_r, Y_r

    def find_active_by_index(self, idx):
        return molecule_set.active_molecules[idx]

    def find_inactive_by_index(self, idx):
        return molecule_set.inactive_molecules[idx]

    def find_by_pos(self, xyz):
        for mol in [molecule_set.inactive_molecules,molecule_set.active_molecules]:
            b = np.array([mol.x, mol.y, mol.z])
            # print(b)
            dist = np.linalg.norm(xyz - b)
            if dist < 0.001:
                return mol._id, str(mol.dh), str(mol.tds), str(mol.dg)

    def get_molecules_from_PCA_area(self , x_low, x_high , y_low, y_high):

        a_pos, in_pos = self.num_PCA_6D()
        a_mols, ia_mols = [] ,[]
        idx = 0
        for mm in a_pos:
            if mm[0] < x_high and mm[0] > x_low and mm[1] < y_high and mm[1] > y_low :
                a_mols.append(self.find_active_by_index(idx))
            idx += 1
        print('Index A: ', idx)
        idx = 0
        for mm in a_pos:
            if mm[0] < x_high and mm[0] > x_low and mm[1] < y_high and mm[1] > y_low:
                ia_mols.append(self.find_inactive_by_index(idx))
            idx += 1
        print('Index I: ', idx)
        return a_mols, ia_mols

    def get_xyz_ext(self, external):

        all_pos_active = np.zeros([int(len(external)), 3])
        idx = 0
        for moll in external:
            # print (moll.x, moll.y ,moll.z)
            all_pos_active[idx][0] = moll.x
            all_pos_active[idx][1] = moll.y
            all_pos_active[idx][2] = moll.z
            idx += 1

        return all_pos_active

    def get_all_ext(self, external):

        all_pos_active = np.zeros([int(len(external)), 6])
        idx = 0
        for moll in external:
            # print (moll.x, moll.y ,moll.z)
            all_pos_active[idx][0] = moll.x
            all_pos_active[idx][1] = moll.y
            all_pos_active[idx][2] = moll.z
            all_pos_active[idx][3] = moll.dh
            all_pos_active[idx][4] = moll.tds
            all_pos_active[idx][5] = moll.dg
            idx += 1

        return all_pos_active

    def mols_distance_from_com(self,mol_obj):
        xyz = self.get_com()
        X = []

        for mol in mol_obj:
            b = np.array([mol.x, mol.y, mol.z])
            # print(b)
            dist = np.linalg.norm(xyz - b)
            if dist < 10.0:
                X.append(mol)

        return X

    def plot_LDA(self, active, inactive):
        X = []
        Y = []
        for mol in active:
            b = np.array([mol.x - center[0], mol.y - center[1], mol.z - center[2], mol.dg, mol.dh])#, mol.tds])
            X.append(b)
            Y.append(np.array([1]))

        for mol in inactive:
            b = np.array([mol.x - center[0], mol.y - center[1], mol.z - center[2], mol.dg, mol.dh])#, mol.tds])
            X.append(b)
            Y.append(np.array([2]))

        plt.figure(23)
        #Y = np.asarray(Y)
        lda = LDA(n_components=2)

        X_r2 = lda.fit(X,Y).transform(X)

        #plt.scatter(X_r2[Y == 1, 0],X_r2[Y == 1, 1], s=2, c='b',
         #           label="Active ")  # +str(pca.get_precision()))
        #plt.scatter(X_r2[Y == 2, 0], X_r2[Y == 2, 1], s=2, c='r',
           #         label="Inactive ")  # +str(pca.get_precision()))
        print('-------------------------------')
        print(np.shape(X_r2))
        plt.scatter(X_r2[:, 0], X_r2[:, 1], s=2, c='r', label="Inactive ")

        plt.legend(loc='upper left')

        plt.title('LDA w/ 2 components from 6d')
        plt.show()

        #return X_r, Y_r




wm = molecule_set()
wm.load_dataset()

center = np.mean(wm.get_xyz('all'),axis=0)

#print ('Center: ', center)

print(wm.get_com())

wm.plot_KDE_2D()

wm.plot_PCA_6D(wm.active_molecules,wm.inactive_molecules)


wm.plot_PCA_6D(wm.mols_distance_from_com(wm.active_molecules),wm.mols_distance_from_com(wm.inactive_molecules))

wm.plot_LDA(wm.active_molecules,wm.inactive_molecules)

a_pca, ia_pca = wm.get_molecules_from_PCA_area(-20,-3,-3,15)








np_a_pca = wm.get_all_ext(a_pca)
np_ia_pca = wm.get_all_ext(ia_pca)
print(np.shape(np_a_pca))
print(np.shape(np_ia_pca))
f, axarr = plt.subplots(3, sharex=True)

axarr[0].hist(np_a_pca[:,0], 50, normed=1,  alpha=0.75)
axarr[1].hist(np_a_pca[:,1], 50, normed=1,  alpha=0.75)
axarr[2].hist(np_a_pca[:,2], 50, normed=1,  alpha=0.75)
axarr[0].hist(np_ia_pca[:,0], 50, normed=1,  alpha=0.75)
axarr[1].hist(np_ia_pca[:,1], 50, normed=1,  alpha=0.75)
axarr[2].hist(np_ia_pca[:,2], 50, normed=1,  alpha=0.75)


plt.show()

f, axarr = plt.subplots(3, sharex=True)

axarr[0].hist(np_a_pca[:,3], 50, normed=1,  alpha=0.75)
axarr[1].hist(np_a_pca[:,4], 50, normed=1,  alpha=0.75)
axarr[2].hist(np_a_pca[:,5], 50, normed=1,  alpha=0.75)
axarr[0].hist(np_ia_pca[:,3], 50, normed=1,  alpha=0.75)
axarr[1].hist(np_ia_pca[:,4], 50, normed=1,  alpha=0.75)
axarr[2].hist(np_ia_pca[:,5], 50, normed=1,  alpha=0.75)
plt.show()
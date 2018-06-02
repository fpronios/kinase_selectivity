from xlrd import open_workbook
from scipy.interpolate import spline
import plotly

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
import seaborn as sns; sns.set()
import matplotlib as mpl
from sklearn.model_selection import cross_validate, KFold
from sklearn import tree
import graphviz

#from sheatheer_jones import *
#from scipy.spatial import Delaunay
#mpl.style.use('classic')

from rdf_pdf import * #pairCorrelationFunction_3D

wb = open_workbook('watermap_statistical_xyz_analysis.xlsx')




def kde2D(x, y, bandwidth,min_max= [0,1,0,1] , xbins=100j, ybins=100j,  **kwargs):
    """Build 2D kernel density estimate (KDE)."""

    # create grid of sample locations (default: 100x100)
    xx, yy = np.mgrid[min_max[0]:min_max[1]:xbins,
                        min_max[2]:min_max[3]:ybins]

    xy_sample = np.vstack([yy.ravel(), xx.ravel()]).T
    xy_train  = np.vstack([y, x]).T

    kde_skl = KernelDensity(bandwidth=bandwidth, **kwargs)
    kde_skl.fit(xy_train)

    # score_samples() returns the log-likelihood of the samples
    z = np.exp(kde_skl.score_samples(xy_sample))
    return xx, yy, np.reshape(z, xx.shape)


def get_axis_KDE( i_mol_obj, a_mol_obj, axis ,z_low= 1, z_high = 1 , com= [0,0,0]):



    z_min = 500
    z_max = -500
    x_min = 500
    x_max = -500
    y_min = 500
    y_max = -500

    for mol in a_mol_obj:
        if z_min > mol.z :
            z_min = mol.z
        if z_max < mol.z :
            z_max = mol.z

    for mol in a_mol_obj:
        if x_min > mol.x :
            x_min = mol.x
        if x_max < mol.x :
            x_max = mol.x

    for mol in a_mol_obj:
        if y_min > mol.y :
            y_min = mol.y
        if y_max < mol.y :
            y_max = mol.y

    for zz in range(int(z_max+abs(z_min))):
        Xa, Ya = [], []
        Xi, Yi = [], []

        for mol in a_mol_obj:
            if mol.z < z_min + zz +1 and mol.z > z_min + zz :
                if axis == 'x':
                    Xa.append(mol.x)
                #elif axis == 'y':
                    Ya.append(mol.y)
                else:
                    Xa.append(mol.z)

        for mol in i_mol_obj:
            if mol.z < z_min + zz +1 and mol.z > z_min + zz :
                if axis == 'x':
                    Xi.append(mol.x)
                #elif axis == 'y':
                    Yi.append(mol.y)
                else:
                    Xi.append(mol.z)

        #fig, (ax1, ax2) = plt.subplots(1, 2)
        xa, ya = np.asarray(Xa), np.asarray(Ya)
        xi, yi = np.asarray(Xi), np.asarray(Yi)

        # Set up the figure
        if len(Xi) > 1 and len(Xa) > 1:
            f, ax = plt.subplots(figsize=(8, 8))
            ax.set_aspect("equal")

            # Draw the two density plots
            ax = sns.kdeplot(xi, yi,
                             cmap="Reds", shade=True, shade_lowest=False, alpha = 0.65)
            ax = plt.scatter(xi, yi, marker="+", c = 'r', label = "Inactive")
            ax = sns.kdeplot(xa, ya,
                             cmap="Blues", shade=True, shade_lowest=False , alpha = 0.65)
            ax = plt.scatter(xa, ya, marker="+", c='b',  label = "Active")
            ax = plt.title('$Z$ between [%0.2f, %0.2f] $\AA$, {mean = %0.2f}' % (z_min+zz ,z_min+zz + 1,com[2],))
            ax = plt.ylabel("$Y$: $\AA$")
            ax = plt.xlabel("$X$: $\AA$")

            ax = plt.scatter(com[0], com[1], s=80,marker="+", facecolor='green' , label = 'Center of Mass')

            plt.legend(loc='upper right')
            plt.savefig("saved_figs/combined/comb_z_%d_%0.2f_%0.2f_.png" % (100+zz,z_min + zz, z_min + zz + 1,))
            #plt.show()

        sns.set(style="white")

        if len(Xi) > 1 and len(Xa) > 1000:
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 5))

            len_xa, len_xi, len_ya, len_yi = len(Xa), len(Xi), len(Ya), len(Yi)
            plt.suptitle('$Z$ between [%0.2f, %0.2f] $\AA$, {mean = %0.2f} \nNumber of Active : %d, Inactive : %d' % (z_min + zz, z_min + zz + 1,com[2],len_xa,len_xi,))

            ax1 = sns.kdeplot(xa, color='b', label='Active', ax=ax1)
            ax1 = sns.rugplot(xa, color='b', ax=ax1)
            ax1 = sns.kdeplot(xi, color='r', label='Inactive', ax=ax1)
            ax1 = sns.rugplot(xi, color='r', ax=ax1)
            ax1 = sns.rugplot([com[0]], color='g', height=0.1, ax=ax1, label = 'Center of Mass')

            #ax = plt.title('$Z$ in (%0.2f, %0.2f)' % (z_min + zz, z_min + zz + 1,))
            ax2 = sns.kdeplot(ya, color='b', label='Active', ax=ax2)
            ax2 = sns.rugplot(ya, color='b', ax=ax2)
            ax2 = sns.kdeplot(yi, color='r', label='Inactive', ax=ax2)
            ax2 = sns.rugplot(yi, color='r', ax=ax2)
            ax2 = sns.rugplot([com[1]], color='g', height=0.1, ax=ax2 , label = 'Center of Mass')

            ax1.set_title("$X$ axis")
            ax1.set_ylabel("Density")
            ax1.set_xlabel("$\AA$")
            ax2.set_title("$Y$ axis")
            ax2.set_xlabel("$\AA$")
            ax1 = plt.legend(loc = 'upper right')
            ax2 = plt.legend(loc = 'upper right')
            #fig.
            plt.subplots_adjust(left=0.05, bottom=None, right=0.95, top=0.8,
                            wspace=None, hspace=None)

            plt.savefig("saved_figs/xy_hist/comb_z_%d_%0.2f_%0.2f_.png" % (100+zz,z_min + zz, z_min + zz + 1,))
            #plt.show()

        if len(Xi) > 1 and len(Xa) > 1000:


            g = sns.jointplot(xa, ya, kind="kde", size=8, space=0 )
            g.plot_joint(plt.scatter, c="w", s=35, linewidth=1, marker="+")
            g.ax_joint.collections[0].set_alpha(0)
            g.fig.suptitle('$Z$ in (%0.2f, %0.2f)' % (z_min+zz ,z_min+zz + 1,))
            g.set_axis_labels("$X$", "$Y$")
            g.savefig("saved_figs/iso_levels/active_z_%d_%0.2f_%0.2f_.png" % (100+zz,z_min + zz, z_min + zz + 1,))
            #plt.show()

            g = sns.jointplot(xi, yi, kind="kde", size=8, space=0)
            g.plot_joint(plt.scatter, c="w", s=35, linewidth=1, marker="+")
            g.ax_joint.collections[0].set_alpha(0)
            g.fig.suptitle('$Z$ in (%0.2f, %0.2f)' % (z_min + zz, z_min + zz + 1,))
            g.set_axis_labels("$X$", "$Y$")
            g.savefig("saved_figs/iso_levels/inactive_z_%d_%0.2f_%0.2f_.png" % (100+zz,z_min + zz, z_min + zz + 1,))

    # sns.jointplot(x,y, kind="kde", size=7, space=0, ax=g , color='r')


    # Silverman's Rule of Thumb
    #BW = 1.06 * np.power(np.std(np.asarray(Xa), axis=0), -0.2) / 2

    #xx, yy, zz = kde2D(x, y, BW *5, [x_min,x_max,y_min,y_max])
    #fig, (ax1, ax2) = plt.subplots(1, 2)
    #ax1.pcolormesh(xx, yy, zz)
    #ax1.scatter(x, y, s=2, facecolor='white')
    #ax1.scatter(com[0], com[1], s=4, facecolor='red')
    #ax1.set_title("Active, z: ( %.2f ,%.2f )" %(z_low, z_high, ))
    #ax1.set_xlim(x_min, x_max)
    #ax1.set_ylim(y_min, y_max)
    #plt.show()

    #x, y = np.asarray(Xi), np.asarray(Yi)
    #xx, yy, zz = kde2D(x, y, BW * 5, [x_min, x_max, y_min, y_max])
    #ax2.pcolormesh(xx, yy, zz)
    #ax2.set_title("Inactive, z: ( %.2f ,%.2f )" %(z_low, z_high, ))
    #ax2.scatter(x, y, s=2, facecolor='white')
    #ax2.scatter(com[0], com[1], s=4, facecolor='red')
    #ax2.set_xlim(x_min,x_max)
    #ax2.set_ylim(y_min, y_max)



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
    prot_names_inact = []
    prot_names_act = []
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

    def get_dgs(self, act):

        all_pos_active = np.zeros([int(len(molecule_set.active_molecules)), 3])
        np_names_active = np.zeros([int(len(molecule_set.active_molecules)), 1])
        idx = 0
        for moll in molecule_set.active_molecules:
            # print (moll.x, moll.y ,moll.z)
            all_pos_active[idx][0] = moll.dg
            all_pos_active[idx][1] = moll.tds
            all_pos_active[idx][2] = moll.dh
            idx += 1

        all_pos_inactive = np.zeros([int(len(molecule_set.inactive_molecules)),3])
        np_names_inactive = np.zeros([int(len(molecule_set.inactive_molecules)),1])
        idx = 0
        for moll in molecule_set.inactive_molecules:

            #print (moll.x, moll.y ,moll.z)
            all_pos_inactive[idx][0] = moll.dg
            all_pos_inactive[idx][1] = moll.tds
            all_pos_inactive[idx][2] = moll.dh
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

    def get_mean_enrg(self):
        return  np.mean((self.get_dgs('all')), axis=0)

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

        return X

    def plot_KDE_2D(self):
        self.plot_KDE_2D_per_mol(molecule_set.active_molecules, 'Active')
        self.plot_KDE_2D_per_mol(molecule_set.inactive_molecules, 'Inactive')

        plt.legend(loc='upper left')
        plt.ylabel('Density')
        plt.xlabel(r'$\AA$')
        plt.title('KDE of H2O distances from center')
        plt.show()

    def plot_PCF_by_prot(self):
        molecule_set.prot_names_inact = list(set([prot.prot_name for prot in molecule_set.inactive_molecules]))
        molecule_set.prot_names_act = list(set([prot.prot_name for prot in molecule_set.active_molecules]))
        #self.plot_KDE_2D_per_mol(molecule_set.active_molecules, 'Active')
        #self.plot_KDE_2D_per_mol(molecule_set.inactive_molecules, 'Inactive')

        self.get_PCF_ltt_by_prot(molecule_set.active_molecules,molecule_set.prot_names_act, 'Active')
        self.get_PCF_ltt_by_prot(molecule_set.inactive_molecuvtkles, molecule_set.prot_names_inact, 'Inactive')

        plt.legend(loc='upper left')
        plt.ylabel('Density')
        plt.xlabel(r'$\AA$')
        plt.title('KDE of H2O distances from center')
        plt.show()

    def get_PCF_ltt_by_prot(self,mol_obj, prot_names,label):
        xyz = self.get_com()
        for prot in prot_names:
            X = []
            IDX = []
            for idx, mol in enumerate(mol_obj):
                b = np.array([mol.x, mol.y, mol.z])
                # print(b)
                dist = np.linalg.norm(xyz - b)
                if dist < 10:
                    X.append(dist)
                    IDX.append(idx)

            mols_ltt = len(X)
            mol_pcf = [[] for i in range(mols_ltt)]
            for idx, mol in enumerate(mol_obj):
                if idx in IDX:
                    for idx_inner, mol_inner in enumerate(mol_obj):
                        if idx_inner in IDX and idx_inner != idx and mol_inner.prot_name == mol.prot_name and mol.prot_name == prot:
                            a = np.array([mol.x, mol.y, mol.z])
                            b = np.array([mol_inner.x, mol_inner.y, mol_inner.z])
                            # print(b)
                            dist = np.linalg.norm(a - b)

                            mol_pcf[IDX.index(idx)].append(dist)
            to_remove = []
            for idx,  mol_empty_test in enumerate(mol_pcf):
                if len(mol_empty_test)==0:
                    #print(idx)
                    to_remove.append(idx)

            mol_pcf = [i for j, i in enumerate(mol_pcf) if j not in to_remove]

            X = np.array(mol_pcf)
            print("mol_pcf")
            X = X.ravel()
            print(len(X))
            BW = 1.06 * np.power(np.std(np.asarray(X), axis=0), -0.2) / 2
            print(BW)
            for kernel in ['epanechnikov']:  # , 'tophat', 'epanechnikov']:
                kde = KernelDensity(kernel=kernel, bandwidth=BW).fit(np.asarray(X).reshape(-1, 1))
                X_plot = np.linspace(0, 25, 10000)[:, np.newaxis]
                log_dens = kde.score_samples(X_plot)
                if label == 'Active':
                    l_type = '-'
                    c_color = 'r'
                else:
                    l_type = '*-'
                    c_color = 'b'
                plt.plot(X_plot[:, 0], np.exp(log_dens), c =c_color,
                         label="{1} BW = '{0:.2f}'".format(BW,  label + ": " + prot))

            plt.legend(loc='upper right', prop={'size': 2})

        return X

    def get_PCF_ltt(self, mol_obj, label):
        xyz = self.get_com()
        X = []
        IDX = []
        for idx , mol in enumerate(mol_obj):
            b = np.array([mol.x, mol.y, mol.z])
            # print(b)
            dist = np.linalg.norm(xyz - b)
            if dist < 10:
                X.append(dist)
                IDX.append(idx)

        mols_ltt = len(X)
        mol_pcf = [[] for i in range (mols_ltt)]
        for idx, mol in enumerate(mol_obj):
            if idx in IDX:
                for idx_inner, mol_inner in enumerate(mol_obj):
                    if idx_inner in IDX and idx_inner != idx:
                        a = np.array([mol.x, mol.y, mol.z])
                        b = np.array([mol_inner.x, mol_inner.y, mol_inner.z])
                        # print(b)
                        dist = np.linalg.norm(a - b)

                        mol_pcf[IDX.index(idx)].append(dist)

        X = np.array(mol_pcf)
        print("mol_pcf")
        X = X.ravel()

        BW = 1.06 * np.power(np.std(np.asarray(X), axis=0), -0.2) / 2
        print(BW)
        for kernel in ['epanechnikov']:  # , 'tophat', 'epanechnikov']:
            kde = KernelDensity(kernel=kernel, bandwidth=BW).fit(np.asarray(X).reshape(-1, 1))
            X_plot = np.linspace(0, 25, 100000)[:, np.newaxis]
            log_dens = kde.score_samples(X_plot)

            plt.plot(X_plot[:, 0], np.exp(log_dens), '-',
                     label="{1} BW = '{0:.2f}'".format(BW, label))

        plt.legend(loc='upper left')

        return  X

    def plot_KDE_2D_per_mol_enrg(self, mol_obj, label):
        xyz = self.get_mean_enrg()
        X = []

        for mol in mol_obj:
            b = np.array([mol.dh, mol.tds, mol.dg])
            # print(b)
            dist = np.linalg.norm(xyz - b)
            if np.linalg.norm(xyz) > np.linalg.norm(b) :
                dist = - dist
            X.append(mol.dg)

        # Silverman's Rule of Thumb
        BW = 1.06 * np.power(np.std(np.asarray(X), axis=0), -0.2) / 2
        print(BW)
        for kernel in ['epanechnikov']:  # , 'tophat', 'epanechnikov']:
            kde = KernelDensity(kernel=kernel, bandwidth=BW).fit(np.asarray(X).reshape(-1, 1))
            X_plot = np.linspace(-25, 25, 10000)[:, np.newaxis]
            log_dens = kde.score_samples(X_plot)

            plt.plot(X_plot[:, 0], np.exp(log_dens), '-',
                     label="{1} BW = '{0:.2f}'".format(BW, label))

        plt.legend(loc='upper right')

        return X

    def plot_PCF(self):
        self.get_PCF_ltt(molecule_set.active_molecules, 'Active')
        self.get_PCF_ltt(molecule_set.inactive_molecules, 'Inactive')

        plt.legend(loc='upper right')
        plt.ylabel('Density')
        plt.xlabel('Energy Delta')
        plt.title('KDE of H2O distances from mean energies')
        plt.show()

    def plot_KDE_2D_enrg(self):
        self.plot_KDE_2D_per_mol_enrg(molecule_set.active_molecules, 'Active')
        self.plot_KDE_2D_per_mol_enrg(molecule_set.inactive_molecules, 'Inactive')

        plt.legend(loc='upper right')
        plt.ylabel('Density')
        plt.xlabel('Energy Delta')
        plt.title('KDE of H2O distances from mean energies')
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

    def mols_distance_from_enrg(self,mol_obj):
        xyz = self.get_mean_enrg()
        X = []

        for mol in mol_obj:
            b = np.array([mol.dh, mol.tds, mol.dg])
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

    def get_per_protein_ctc(self):
        molecule_set.prot_names_inact = list(set([prot.prot_name for prot in molecule_set.inactive_molecules]))
        molecule_set.prot_names_act = list(set([prot.prot_name for prot in molecule_set.active_molecules]))
        print(molecule_set.prot_names_inact)
        print(molecule_set.prot_names_act)

        usfl_act = [[] for i in range(len(molecule_set.prot_names_act))]
        usfl_inact = [[] for i in range(len(molecule_set.prot_names_inact))]

        for prot, i in enumerate(molecule_set.prot_names_act):
            usfl_act[i].append([mol for mol in molecule_set.active_molecules if (mol.protein_name == prot)])

    def get_n_closest_by_prot(self,mol_obj, act = True,n = 25):
        molecule_set.prot_names_act = list(set([prot.prot_name for prot in molecule_set.active_molecules]))
        molecule_set.prot_names_inact = list(set([prot.prot_name for prot in molecule_set.inactive_molecules]))
        if act:
            prot_names = molecule_set.prot_names_act
        else:
            prot_names = molecule_set.prot_names_inact

        xyz = self.get_com()
        n_by_prot = []
        for prot in prot_names:
            X = []
            IDX = []
            for idx, mol in enumerate(mol_obj):
                b = np.array([mol.x, mol.y, mol.z])
                # print(b)
                dist = np.linalg.norm(xyz - b)
                X.append(dist)
                IDX.append(idx)

            #x_np = np.asarray(X)
            #idx_np = np.asarray(IDX)

            Z = [x for _, x in sorted(zip(X, IDX))]
            IDX = Z

            mols_ltt = len(X)
            mol_pcf = [[] for i in range(mols_ltt)]
            n_test = 0
            for idx, mol in enumerate(mol_obj):
                if idx in IDX:
                    for idx_inner, mol_inner in enumerate(mol_obj):
                        if n_test < n and idx_inner in IDX and idx_inner != idx and mol_inner.prot_name == mol.prot_name and mol.prot_name == prot:
                            a = np.array([mol_inner.x, mol_inner.y, mol_inner.z])
                            mol_pcf[IDX.index(idx)].append(a)
                            n_test += 1
            to_remove = []
            for idx,  mol_empty_test in enumerate(mol_pcf):
                if len(mol_empty_test)==0:
                    #print(idx)
                    to_remove.append(idx)

            mol_pcf = [i for j, i in enumerate(mol_pcf) if j not in to_remove]

            X = np.array(mol_pcf)
            #print("N closest mols %s" % prot)
            #print(X)
            n_by_prot.append(X)

        return n_by_prot


wm = molecule_set()
wm.load_dataset()

center = np.mean(wm.get_xyz('all'),axis=0)

#print ('Center: ', center)


def mayavi_3d_density(wm):
    xi, yi, zi = wm.get_xyz_ext(wm.mols_distance_from_com(
        wm.inactive_molecules)).T  # tot_mols[:inactive_rows, 0],tot_mols[:inactive_rows, 1],tot_mols[:inactive_rows, 2]
    xa, ya, za = wm.get_xyz_ext(wm.mols_distance_from_com(
        wm.active_molecules)).T  # tot_mols[inactive_rows:, 0],tot_mols[inactive_rows:, 1],tot_mols[inactive_rows:, 2]

    import hdbscan

    from scipy import stats
    from mayavi import mlab
    kde = stats.gaussian_kde(wm.get_xyz_ext(wm.mols_distance_from_com(wm.active_molecules)).T)

    density = kde(wm.get_xyz_ext(wm.mols_distance_from_com(wm.active_molecules)).T)

    # Plot scatter with mayavi
    figure = mlab.figure('DensityPlot')
    pts = mlab.points3d(xa, ya, za, density, scale_mode='none', scale_factor=1)
    mlab.axes()
    mlab.show()

    density = kde(wm.get_xyz_ext(wm.mols_distance_from_com(wm.inactive_molecules)).T)
    # Evaluate kde on a grid
    xmin, ymin, zmin = xi.min(), yi.min(), zi.min()
    xmax, ymax, zmax = xi.max(), yi.max(), zi.max()
    xi, yi, zi = np.mgrid[xmin:xmax:30j, ymin:ymax:30j, zmin:zmax:30j]
    coords = np.vstack([item.ravel() for item in [xi, yi, zi]])
    density = kde(coords).reshape(xi.shape)

    # Plot scatter with mayavi
    figure = mlab.figure('DensityPlot')

    grid = mlab.pipeline.scalar_field(xi, yi, zi, density)
    min = density.min()
    max = density.max()
    mlab.pipeline.volume(grid, vmin=min, vmax=min + .5 * (max - min))

    mlab.axes()
    mlab.show()

inactive_rows = wm.get_all_ext(wm.inactive_molecules).shape[0]
inactive_cols = wm.get_all_ext(wm.inactive_molecules).shape[1]
active_rows = wm.get_all_ext(wm.active_molecules).shape[0]

tot_mols = np.ones((inactive_rows + active_rows,inactive_cols+1))
tot_mols[:inactive_rows,:-1] =  wm.get_all_ext(wm.inactive_molecules)
tot_mols[:inactive_rows,-1] = 0

tot_mols[inactive_rows:,:-1] = wm.get_all_ext(wm.active_molecules)
tot_mols[inactive_rows:,-1] = 1

print([wm.get_com()[0],wm.get_com()[1],wm.get_com()[2]])

#get_axis_KDE( wm.inactive_molecules, axis = 'x')
from plotly.graph_objs import  Layout,Scatter3d, Marker,Mesh3d, Figure,Data




def chain_grouping(n_cl):
    n = len(n_cl[0][0])

    #for prot in n_cl:
    #print(len(n_cl[1][0]))
    #    print("&&&&&&&&&&&&&&&&&&&&&    ")
     #   print(prot[0])
    distance_list = [[0 for j in range(len(n_cl))] for i in range(n)]
    xyz_list = [[0 for j in range(len(n_cl))] for i in range(n)]
    for mol in range(n):
        for en, prot in enumerate(n_cl):
            compare_to = prot[0][mol]
            min = 1111
            for mol_search in range(n):
                b = prot[0][mol_search]
                #print("To test:", b)
                #print("Compare to: ", compare_to)
                dist = np.linalg.norm(compare_to - b)
                if dist < min and dist != 0.0:
                    min = dist
                    distance_list[mol][en] = dist
                    xyz_list[mol][en] = b


    #print(len(distance_list[0][0]))
    pre_mean_list = []
    for i in range(n):
        #print(distance_list[i])
        #print(xyz_list[i])
        pre_mean_list.append(np.asarray(xyz_list[i]))

    mean_lst = []
    for i in range(n):

        mean_lst.append(np.mean(pre_mean_list[i], axis = 0))

    #print(mean_lst[0])
    #print(mean_lst[1])
    mean_np = np.asarray(mean_lst)
    #print(mean_np)
    return mean_np#group_com

#######
n_closest = wm.get_n_closest_by_prot(wm.active_molecules,act= True, n = 25)
mean_act = chain_grouping(n_closest)

n_closest = wm.get_n_closest_by_prot(wm.inactive_molecules,act= False, n = 25)
mean_inact = chain_grouping(n_closest)

#mean_act = np.load("mean_act.npy")
#mean_inact = np.load("mean_inact.npy")
np.save("mean_act",mean_act)
np.save("mean_inact",mean_inact)
print(np.shape(mean_act))
print(np.shape(mean_inact))

xa, ya, za = mean_act[:, 0],mean_act[:, 1],mean_act[:, 2]
xi, yi, zi = mean_inact[:, 0],mean_inact[:, 1],mean_inact[:, 2]

print()

alpha = 4
global_opacity = 0.6
#print(wm.get_xyz_ext(wm.mols_distance_from_com(wm.inactive_molecules,"NONE")))


input = ("DONE")
pointsi=Scatter3d(mode = 'markers',
                 name = 'active_points',
                 x =xi,
                 y= yi,
                 z= zi,
                 marker = Marker( size=4, color='#458B00' )
)
simplexesi = Mesh3d(alphahull = alpha,
                   name = 'active_simplexes',
                   x =xi,
                   y= yi,
                   z= zi,
                   color='rgba(25,100,255, 0.5)', #set the color of simplexes in alpha shape
                   opacity=global_opacity
)
pointsa=Scatter3d(mode = 'markers',
                 name = 'inactive_points',
                 x =xa,
                 y= ya,
                 z= za,
                 marker = Marker( size=4, color='#FF8B00' )
)
simplexesa = Mesh3d(alphahull = alpha,
                   name = 'inactive_simplexes',
                   x =xa,
                   y= ya,
                   z= za,
                   color='rgba(255,0,0, 0.5)', #set the color of simplexes in alpha shape
                   opacity=global_opacity
)
#x_style = dict( zeroline=False, range=[-2.85, 4.25], tickvals=np.linspace(-2.85, 4.25, 5)[1:].round(1))
#y_style = dict( zeroline=False, range=[-2.65, 1.32], tickvals=np.linspace(-2.65, 1.32, 4)[1:].round(1))
#z_style = dict( zeroline=False, range=[-3.67,1.4], tickvals=np.linspace(-3.67, 1.4, 5).round(1))
layout=Layout(title='Alpha shape of Active vs Inactive proteins alpha = %d' % alpha#,
              #width=1200,
             #height=1000
              #scene = Scene(
              #xaxis = x_style,
              #yaxis = y_style,
              #zaxis = z_style
            # )
             )
fig1=Figure(data=Data([pointsi, simplexesi,pointsa, simplexesa]), layout=layout)
#fig1 = go.Figure(data=data1, layout=layout)
plotly.offline.plot(fig1, filename='plotly_out/alpha_shape_a_%d_op_%0.2f_n_close.html' % (alpha,global_opacity,))

print(wm.get_com())




"""

get_axis_KDE( wm.inactive_molecules, wm.active_molecules, 'x' , wm.get_com()[2] - 0.25 ,wm.get_com()[2] + 0.25,[wm.get_com()[0],wm.get_com()[1],wm.get_com()[2]])

wm.get_per_protein_ctc()

X_a = wm.plot_KDE_2D_per_mol(wm.active_molecules,'Active')
Y_a = wm.plot_KDE_2D_per_mol_enrg(wm.active_molecules,'Active')
X_i = wm.plot_KDE_2D_per_mol(wm.inactive_molecules,'Active')
Y_i = wm.plot_KDE_2D_per_mol_enrg(wm.inactive_molecules,'Active')


dist_a = []
enrg_delta_a = []
for x, y in zip(X_a,Y_a):
    if x < 10.0:
        dist_a.append(x)
        enrg_delta_a.append(y)


dist_i = []
enrg_delta_i  = []
for x, y in zip(X_i,Y_i):
    if x < 10.0:
        dist_i.append(x)
        enrg_delta_i.append(y)


fig = plt.figure(4)
plt.clf()
plt.cla()
plt.scatter(dist_a,enrg_delta_a, label = 'Active' , s = 4)

plt.scatter(dist_i,enrg_delta_i, label = 'Inactive', s = 4)
plt.xlabel(r'$\AA$')
plt.ylabel("Energy Delta")
plt.title('Energy vs distance')
plt.show()



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

"""



#wm.get_PCF_ltt(wm.active_molecules, "Test")
#wm.plot_PCF_by_prot()
#wm.plot_PCF()

#wm.get_per_protein_ctc()
#prot_names_inact = list(set([prot.prot_name for prot in wm.inactive_molecules]))
#prot_names_act = list(set([prot.prot_name for prot in wm.active_molecules]))

#wm.prot_names_act = prot_names_act
#wm.prot_names_inact = prot_names_inact

#tri = Delaunay(tot_mols[inactive_rows:,:-4], 3)

#print(tri.simplices)

#layout = dict(title = 'Active RDF')
#fig = dict(data = (points[:,0], points[:,1], tri.simplices.copy()), layout = layout)
#plotly.offline.plot(fig, filename='delauney_3D.html')

#wm.plot_KDE_2D_enrg()


#wm.plot_KDE_2D()

#wm.plot_PCA_6D(wm.active_molecules,wm.inactive_molecules)


#wm.plot_PCA_6D(wm.mols_distance_from_com(wm.active_molecules),wm.mols_distance_from_com(wm.inactive_molecules))

#wm.plot_LDA(wm.active_molecules,wm.inactive_molecules)

#a_pca, ia_pca = wm.get_molecules_from_PCA_area(-20,-3,-3,15)


#plotly.offline.plot({
#    "data": [Scatter3d(x=x, y=y,z=z, marker = dict(color = density), mode = 'markers')],
#
#    "layout": Layout(title="Active protein water molecules heatmap")
#})

#plt.show()
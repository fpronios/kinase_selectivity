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

import scipy as sp
from scipy import stats
from matplotlib.widgets import CheckButtons
import matplotlib.pyplot as plt
import seaborn as sns; sns.set()


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

    def __init__(self, x,y,z,prot_name, _id):
        self.x = x
        self.y = y
        self.z = z
        self._id = _id
        self.prot_name = prot_name


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

                    mol_obj_active.append(w_mol(x_value,y_value,z_value,names_lst_active[mol_idx],names_lst_active[mol_idx] + str(row)))
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

                    mol_obj_inactive.append(w_mol(x_value, y_value, z_value, names_lst_inactive[mol_idx], names_lst_inactive[mol_idx] + str(row)))
                    # print(tmp_arr)
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



# input pos matrix
dbscan_w_plot(all_pos_active)

#X = all_pos#= StandardScaler().fit_transform(all_pos)
# Particle setup
domain_size = 60.0
# Calculation setup
dr = 0.1
### Random arrangement of particles ###
particle_radius = 0.01
rMax = 10
#Compute pair correlation
dr_set = 0.1
line_smoothing = False


g_r_lst_active , r_lst_active = create_rdf(molecule_pop_lst_active, all_pos_active, domain_size , dr , particle_radius ,rMax,dr_set,line_smoothing)

g_r_lst_inactive , r_lst_inactive = create_rdf(molecule_pop_lst_inactive, all_pos_inactive, domain_size , dr , particle_radius ,rMax,dr_set,line_smoothing)

# Visualize

fig, ax = plt.subplots()
line_lst = []

for r_iter, g_iter in zip(r_lst_inactive, g_r_lst_inactive):
    line_iter, = ax.plot(r_iter , g_iter)
    line_lst.append(line_iter)

plt.legend(tuple(names_lst_active))
plt.xlabel('r')
plt.ylabel('g(r)')
plt.xlim( (0, rMax) )
plt.ylim( (0, 1.6 * domain_size))


#plt.plot(r, g_r)#, color='black')
plt.subplots_adjust(left=0.2)

leg = ax.legend(loc='top right', fancybox=True, shadow=True)
#leg.get_frame().set_alpha(0.4)

rax = plt.axes([0.05, 0.3, 0.1, 0.35])
binary_def = [True for i in range (len(names_lst_active))]


check = CheckButtons(rax, tuple(names_lst_active), tuple(binary_def))


def func(label):
    for lbl, idx in zip(names_lst_active, range (len(names_lst_active))):
        if label == lbl:
            line_lst[idx].set_visible(not line_lst[idx].get_visible())

    plt.draw()

check.on_clicked(func)





x_a, y_a,  = find_average(g_r_lst_active,r_lst_active,names_lst_active)
x_i, y_i,  = find_average(g_r_lst_inactive,r_lst_inactive,names_lst_inactive)
plt.figure(3)

plt.plot(x_a, y_a)
plt.plot(x_i, y_i)
plt.legend(('Active','Inactive'))
plt.xlabel('r')
plt.ylabel('g(r)')
plt.title('mean RDF of Active vs Inactive')

x_a, y_a, = find_average(g_r_lst_active,r_lst_active,names_lst_active, False)
x_i, y_i,  = find_average(g_r_lst_inactive,r_lst_inactive,names_lst_inactive, False)
plt.figure(4)

plt.plot(x_a, y_a)
plt.plot(x_i, y_i)
plt.legend(('Active','Inactive'))
plt.xlabel('r')
plt.ylabel('g(r)')
plt.title('mean RDF of Active vs Inactive')


x = all_pos_active[:,0]
y = all_pos_active[:,1]
z = all_pos_active[:,2]
xyz = np.vstack([x,y,z])

density = stats.gaussian_kde(xyz)(xyz)

idx = density.argsort()
x, y, z, density = x[idx], y[idx], z[idx], density[idx]

#kde = stats.gaussian_kde(np.asanyarray(g_r_lst_active).ravel())
#density = kde(all_pos_active)

#kde = KernelDensity(kernel='gaussian', bandwidth=0.5).fit([x_a,y_a])
#log_dens = kde.score_samples(x_a.reshape(1, -1))
##ax.plot(np.exp(log_dens), '-',
 #       label="kernel = '{0}'".format('gaussian'))
#
fig6 = plt.figure(6)
ax = fig6.add_subplot(111, projection='3d')
ax.scatter(x, y, z, c=density)

ax.grid(True)




#plotly.offline.plot({
#    "data": [Scatter3d(x=x, y=y,z=z, marker = dict(color = density), mode = 'markers')],
#
#    "layout": Layout(title="Active protein water molecules heatmap")
#})

#plt.show()
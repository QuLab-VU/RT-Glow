import csv
import numpy as np


#=====static data=======
static_cells5000 = []
static_cells2500 = []
static_cells1250 = []
static_cells625 = []
static_cells312 = []
static_cells156 = []
static_cells78 = []
static_cells39 = []
static_cells0 = []

ct=0
with open('data/H841_Static_Lum_5k.csv', newline='', encoding='utf-8-sig') as csvfile_s:
    full_data_s = csv.reader(csvfile_s, delimiter=',')
    for row in full_data_s:
        if ct == 0:
            static_cells = row
        elif ct > 0:
            static_cells5000.append(float(row[0]))
            static_cells2500.append(float(row[1]))
            static_cells1250.append(float(row[2]))
            static_cells625.append(float(row[3]))
            static_cells312.append(float(row[4]))
            static_cells156.append(float(row[5]))
            static_cells78.append(float(row[6]))
            static_cells39.append(float(row[7]))
            static_cells0.append(float(row[8]))
        ct = ct+1


#====dynamic data====

#tdata_d = []
tdata = []
exp_data_DMSO = [[], [], []]
exp_data_10uM = [[], [], []]
exp_data_2_5uM = [[], [], []]
exp_data_625nM = [[], [], []]
exp_data_156_25nM = [[], [], []]
exp_data_39nM = [[], [], []]
exp_data_9_76nM = [[], [], []]
exp_data_2_44nM = [[], [], []]
exp_data_0_61nM = [[], [], []]
exp_data_0_15nM = [[], [], []]

ct=0
with open('data/H841_RAW_TAK901.csv', newline='') as csvfile_d:
    full_data_d = csv.reader(csvfile_d, delimiter=',')
    for row in full_data_d:
        #print(row)
        if ct > 0:
            tdata.append(float(row[1]))
            exp_data_10uM[0].append(float(row[2]))
            exp_data_10uM[1].append(float(row[3]))
            exp_data_10uM[2].append(float(row[4]))
            exp_data_2_5uM[0].append(float(row[5]))
            exp_data_2_5uM[1].append(float(row[6]))
            exp_data_2_5uM[2].append(float(row[7]))
            exp_data_625nM[0].append(float(row[8]))
            exp_data_625nM[1].append(float(row[9]))
            exp_data_625nM[2].append(float(row[10]))
            exp_data_156_25nM[0].append(float(row[11]))
            exp_data_156_25nM[1].append(float(row[12]))
            exp_data_156_25nM[2].append(float(row[13]))
            exp_data_39nM[0].append(float(row[14]))
            exp_data_39nM[1].append(float(row[15]))
            exp_data_39nM[2].append(float(row[16]))
            exp_data_9_76nM[0].append(float(row[17]))
            exp_data_9_76nM[1].append(float(row[18]))
            exp_data_9_76nM[2].append(float(row[19]))
            exp_data_2_44nM[0].append(float(row[20]))
            exp_data_2_44nM[1].append(float(row[21]))
            exp_data_2_44nM[2].append(float(row[22]))
            exp_data_0_61nM[0].append(float(row[23]))
            exp_data_0_61nM[1].append(float(row[24]))
            exp_data_0_61nM[2].append(float(row[25]))
            exp_data_0_15nM[0].append(float(row[26]))
            exp_data_0_15nM[1].append(float(row[27]))
            exp_data_0_15nM[2].append(float(row[28]))
            exp_data_DMSO[0].append(float(row[29]))
            exp_data_DMSO[1].append(float(row[30]))
            exp_data_DMSO[2].append(float(row[31]))
        ct = ct+1



#=======================
#=====prepare data======
#=======================

def prepare_static_data(data_array, background_luminescence):
    datapoints = len(data_array)
    cleanup = [0] * datapoints
    for j in range(datapoints):
        cleanup[j] = data_array[j] - background_luminescence
        if cleanup[j] < 0:
            cleanup[j] = 0

    return cleanup


def prepare_dynamic_data(data_array, background_luminescence):
    datapoints = len(data_array[0])
    duplicates = len(data_array)
    mean_init_lum = 0

    for i in range(duplicates):
        mean_init_lum = mean_init_lum + data_array[i][0]
    mean_init_lum = mean_init_lum/duplicates

    #background_luminescence = mean_init_lum-initial_luminescence_diff

    clean_mean = [0] * datapoints
    clean_exp1 = [0] * datapoints
    clean_exp2 = [0] * datapoints
    clean_exp3 = [0] * datapoints
    for j in range(datapoints):
        for i in range(duplicates):
            clean_mean[j] = clean_mean[j] + data_array[i][j]
        print(clean_mean[j] / duplicates, background_luminescence, clean_mean[j] / duplicates - background_luminescence)
        clean_mean[j] = clean_mean[j] / duplicates - background_luminescence
        #clean_exp1[j] = data_array[0][j] - background_luminescence
        #clean_exp2[j] = data_array[1][j] - background_luminescence
        #clean_exp3[j] = data_array[2][j] - background_luminescence
        if clean_mean[j] < 0:
            clean_mean[j] = 0
        #if clean_exp1[j] < 0:
        #    clean_exp1[j] = 0
        #if clean_exp2[j] < 0:
        #    clean_exp2[j] = 0
        #if clean_exp3[j] < 0:
        #    clean_exp3[j] = 0

    return [clean_mean, clean_exp1, clean_exp2, clean_exp3]




#=====static data=======

static_initial_ct = [float(x) for x in static_cells]
rev = [x for x in reversed(static_initial_ct)]
rev = np.array(rev)

static_cells_mean = [static_cells0[0], static_cells39[0], static_cells78[0], static_cells156[0],
                     static_cells312[0], static_cells625[0], static_cells1250[0], static_cells2500[0], static_cells5000[0]]
static_cells_exp1 = [static_cells0[1], static_cells39[1], static_cells78[1], static_cells156[1],
                     static_cells312[1], static_cells625[1], static_cells1250[1], static_cells2500[1], static_cells5000[1]]
static_cells_exp2 = [static_cells0[2], static_cells39[2], static_cells78[2], static_cells156[2],
                     static_cells312[2], static_cells625[2], static_cells1250[2], static_cells2500[2], static_cells5000[2]]
static_cells_exp3 = [static_cells0[3], static_cells39[3], static_cells78[3], static_cells156[3],
                     static_cells312[3], static_cells625[3], static_cells1250[3], static_cells2500[3], static_cells5000[3]]

background_luminescence = static_cells_mean[0]

cleaned_static = np.array(prepare_static_data(static_cells_mean, background_luminescence))
cleaned_static_exp1 = np.array(prepare_static_data(static_cells_exp1, background_luminescence))
cleaned_static_exp2 = np.array(prepare_static_data(static_cells_exp2, background_luminescence))
cleaned_static_exp3 = np.array(prepare_static_data(static_cells_exp3, background_luminescence))
#print(cleaned_static)
#print(cleaned_static_exp1)
#print(cleaned_static_exp2)
#print(cleaned_static_exp3)
#print(rev)

f=open('data/cleaned_static_data.txt','w')
#f.write('init_cell_count cleaned_static_data cleaned_experimental_replicate1 cleaned_experimental_replicate2 cleaned_experimental_replicate3\n')
f.write('init_cell_count cleaned_static_data\n')
for i in range(len(rev)):
    #f.write(str(rev[i])+' '+str(cleaned_static[i])+' '+str(cleaned_static_exp1[i])+' '+str(cleaned_static_exp2[i])+' '+str(cleaned_static_exp3[i])+'\n')
    f.write(str(rev[i])+' '+str(cleaned_static[i])+'\n')
f.close()



#====dynamic data=======

tdata = [x/24 for x in tdata]
# initial luminescence is cleaned up static luminescence at a cell count of 300:
initial_luminescence = cleaned_static[4]
#print(initial_luminescence)
#exit()

cleaned_DMSO= np.array(prepare_dynamic_data(exp_data_DMSO, background_luminescence))
cleaned_10uM = np.array(prepare_dynamic_data(exp_data_10uM, background_luminescence))
cleaned_2_5uM = np.array(prepare_dynamic_data(exp_data_2_5uM, background_luminescence))
cleaned_625nM = np.array(prepare_dynamic_data(exp_data_625nM, background_luminescence))
cleaned_156_25nM = np.array(prepare_dynamic_data(exp_data_156_25nM, background_luminescence))
cleaned_39nM = np.array(prepare_dynamic_data(exp_data_39nM, background_luminescence))
cleaned_9_76nM = np.array(prepare_dynamic_data(exp_data_9_76nM, background_luminescence))
cleaned_2_44nM = np.array(prepare_dynamic_data(exp_data_2_44nM, background_luminescence))
cleaned_0_61nM = np.array(prepare_dynamic_data(exp_data_0_61nM, background_luminescence))
cleaned_0_15nM = np.array(prepare_dynamic_data(exp_data_0_15nM, background_luminescence))

print(cleaned_DMSO[0])
print(cleaned_10uM[0])
print(cleaned_2_5uM[0])
print(cleaned_625nM[0])
print(cleaned_156_25nM[0])
print(cleaned_39nM[0])
print(cleaned_9_76nM[0])
print(cleaned_2_44nM[0])
print(cleaned_0_61nM[0])
print(cleaned_0_15nM[0])

####DATA EXACTLY AS IS THAT WE ARE FITTING TO (ONLY THE AVERAGE FILE - create file folder with the input files)
####subtract 118 of everything

f2=open('data/cleaned_dynamic_data.txt','w')
f2.write('tdata DMSO_mean '
         '10uM_mean '
         'cleaned_2.5uM_mean '
         'cleaned_625nM_mean '
         'cleaned_156.25nM_mean '
         'cleaned_39nM_mean '
         'cleaned_9.76nM_mean '
         'cleaned_2.44nM_mean '
         'cleaned_0.61nM_mean '
         'cleaned_0.15nM_mean\n')
'''
f2.write('tdata cleaned_DMSO_mean cleaned_DMSO_exp1 cleaned_DMSO_exp2 cleaned_DMSO_exp3 '
         'cleaned_10uM_mean cleaned_10uM_exp1 cleaned_10uM_exp2 cleaned_10uM_exp3 '
         'cleaned_2.5uM_mean cleaned_2.5uM_exp1 cleaned_2.5uM_exp2 cleaned_2.5uM_exp3 '
         'cleaned_625nM_mean cleaned_625nM_exp1 cleaned_625nM_exp2 cleaned_625nM_exp3 '
         'cleaned_156.25nM_mean cleaned_156.25nM_exp1 cleaned_156.25nM_exp2 cleaned_156.25nM_exp3 '
         'cleaned_39nM_mean cleaned_39nM_exp1 cleaned_39nM_exp2 cleaned_39nM_exp3 '
         'cleaned_9.76nM_mean cleaned_9.76nM_exp1 cleaned_9.76nM_exp2 cleaned_9.76nM_exp3 '
         'cleaned_2.44nM_mean cleaned_2.44nM_exp1 cleaned_2.44nM_exp2 cleaned_2.44nM_exp3 '
         'cleaned_0.61nM_mean cleaned_0.61nM_exp1 cleaned_0.61nM_exp2 cleaned_0.61nM_exp3 '
         'cleaned_0.15nM_mean cleaned_0.15nM_exp1 cleaned_0.15nM_exp2 cleaned_0.15nM_exp3\n') '''
for i in range(len(tdata)):
    #print(str(cleaned_10uM[i]), str(cleaned_2_5uM[i]))

    f2.write(str(tdata[i])+' '+str(cleaned_DMSO[0][i])+' '
             +str(cleaned_10uM[0][i])+' '
             +str(cleaned_2_5uM[0][i])+' '
             +str(cleaned_625nM[0][i])+' '
             +str(cleaned_156_25nM[0][i])+' '
             +str(cleaned_39nM[0][i])+' '
             +str(cleaned_9_76nM[0][i])+' '
             +str(cleaned_2_44nM[0][i])+' '
             +str(cleaned_0_61nM[0][i])+' '
             +str(cleaned_0_15nM[0][i])+'\n')
    '''
    f2.write(str(tdata[i])+' '+str(cleaned_DMSO[0][i])+' '+str(cleaned_DMSO[1][i])+' '+str(cleaned_DMSO[2][i])+' '+str(cleaned_DMSO[3][i])+' '
             +str(cleaned_10uM[0][i])+' '+str(cleaned_10uM[1][i])+' '+str(cleaned_10uM[2][i])+' '+str(cleaned_10uM[3][i])+' '
             +str(cleaned_2_5uM[0][i])+' '+str(cleaned_2_5uM[1][i])+' '+str(cleaned_2_5uM[2][i])+' '+str(cleaned_2_5uM[3][i])+' '
             +str(cleaned_625nM[0][i])+' '+str(cleaned_625nM[1][i])+' '+str(cleaned_625nM[2][i])+' '+str(cleaned_625nM[3][i])+' '
             +str(cleaned_156_25nM[0][i])+' '+str(cleaned_156_25nM[1][i])+' '+str(cleaned_156_25nM[2][i])+' '+str(cleaned_156_25nM[3][i])+' '
             +str(cleaned_39nM[0][i])+' '+str(cleaned_39nM[1][i])+' '+str(cleaned_39nM[2][i])+' '+str(cleaned_39nM[3][i])+' '
             +str(cleaned_9_76nM[0][i])+' '+str(cleaned_9_76nM[1][i])+' '+str(cleaned_9_76nM[2][i])+' '+str(cleaned_9_76nM[3][i])+' '
             +str(cleaned_2_44nM[0][i])+' '+str(cleaned_2_44nM[1][i])+' '+str(cleaned_2_44nM[2][i])+' '+str(cleaned_2_44nM[3][i])+' '
             +str(cleaned_0_61nM[0][i])+' '+str(cleaned_0_61nM[1][i])+' '+str(cleaned_0_61nM[2][i])+' '+str(cleaned_0_61nM[3][i])+' '
             +str(cleaned_0_15nM[0][i])+' '+str(cleaned_0_15nM[1][i])+' '+str(cleaned_0_15nM[2][i])+' '+str(cleaned_0_15nM[3][i])+'\n') '''
f2.close()
###############################################################
# Computational analysis of FRET rates between Chlorophylls
# @author : Han Zhang, Jianping GUO 
# use : python3 Forster.py H_matrix_file spec_tmppath
#       eg: python Forster.py 216_H_matrix.txt  ./spec_tmp/
###############################################################

# ------------------------------------------------------------#
# H_matrix_file contains the overall effective Hamiltonian of the PSI system;
# The diagonal element is the site energy;
# The non-diagonal element is the coupling value between pigment i and j;
# 
# filelist.txt contains a list of pigments of PSI, in the following format:
# resname-chain-reid
# eg : CLA-A-601
# ------------------------------------------------------------#

import math
import numpy as np
import sys
    	
class Forster():
    def __init__(self, S = 160, fwhm = 240):
        self.S = float(S)
        self.fwhm = float(fwhm)

    def cal_gau_spec(self, start, end, step, energy, out_name, type):
        start = float(start)
        end = float(end)
        step = int(step)
        step = float((end - start) / step)
        delta = self.fwhm 
        strengths = []          
        x_data = [] 
        if type == 'abs': 
            for i in np.arange(start, end, step):
                flag = math.exp(-(energy - i) ** 2.0 / (2.0 * (delta ** 2.0))) / (math.sqrt(2.0 * math.pi) * delta)
                strengths.append(flag)
                x_data.append(i)
        elif type == 'flu':
            for i in np.arange(start, end, step):
                flag = math.exp(-(energy - self.S - i) ** 2.0 / (2.0 * (delta ** 2.0))) / (
                            math.sqrt(2.0 * math.pi) * delta)
                strengths.append(flag)
                x_data.append(i)
        else:
            print('abs or flu')
        strengths = np.array(strengths).reshape(-1, 1)
        x_data = np.array(x_data).reshape(-1, 1)
        final_data = np.hstack((x_data, strengths))
        np.savetxt(out_name, final_data, fmt="%15.10f")

    def all_gau_spec(self, energy, path, start, end, step):
        file1 = path + str(energy)+ '_abs.txt'
        file2 = path + str(energy)+ '_flu.txt'
        self.cal_gau_spec(start, end, step, energy, file1, 'abs')
        self.cal_gau_spec(start, end, step, energy, file2, 'flu')

    def cal_J(self, fi, fj):
        spectra1 = np.loadtxt(fi)
        size = spectra1.size
        spectra2 = np.loadtxt(fj)
        x = []
        for i in range(spectra1.shape[0]):
            a = spectra1[i, 1] * spectra2[i, 1]
            x.append(a)    # 强度
        x = np.array(x)
        f1 = np.hstack((spectra1[:, 0], x)).reshape(2, int(size / 2)) 
        J_ij = abs(np.trapz(f1[1, :], f1[0, :])) # 梯形法近似定积分
        return J_ij

    def cal_site_F_k(self, V, J_ab):
        h_bar = 1.05457266e-34  # unit: J*s
        cm_J = 1 / 8065.541 * 1.6021766208e-19  # convert cm-1 into J
        k = (2 * math.pi / h_bar) * ((V * cm_J) ** 2) * (J_ab / 1.9864465642168332e-23)  # unit: s-1
        k_ps = k / (10 ** 12)  # unit ps-1
        return k_ps


if __name__ == '__main__':
    file = sys.argv[1]
    matrix = np.loadtxt(file)
    shape = matrix.shape[0]
    F = Forster()
    tmp_path = sys.argv[2]
    
    for i in range(shape):
        F.all_gau_spec(matrix[i][i], tmp_path, 10000, 20000, 10000)
    print('spec_cal done.')
    
 
    J = np.zeros([shape, shape])
    for i in range(shape):
        for j in range(shape):
            if i != j:
                flufile = tmp_path + str(matrix[i][i]) + '_flu.txt'
                absfile = tmp_path + str(matrix[j][j]) + '_abs.txt'
                print(flufile, absfile)
                J[i][j] = F.cal_J(flufile, absfile)
                print("I J is ", i, j, J[i][j])
            else:
                continue
    print('J_cal done.')

    f = open("filelist.txt", 'r')
    filelist = f.readlines()
    f.close()
    k = np.zeros([shape, shape])
    for i in range(shape):
        for j in range(shape):
            if i != j:
                if filelist[i][:3] == filelist[j][:3]:
                    k[i][j] = F.cal_site_F_k(matrix[i][j], J[i][j])
                    print('i j is ', i, j, k[i][j])
                else:
                    print(filelist[i][:3],filelist[j][:3])
                    J_ref = 0.0009789061621090092
                    k[i][j] = F.cal_site_F_k(matrix[i][j], J_ref)
                    print('i j is ', i, j, k[i][j])
            else:
                continue
    print('k_cal done.')

    J_out_name = 'J.txt'  # spectral_overlap_file
    k_out_name = 'k.txt'  # transfer_rates_file
    np.savetxt(J_out_name, J, fmt="%20.10f")
    np.savetxt(k_out_name, k, fmt="%20.10f")

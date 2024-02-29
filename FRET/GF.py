############################################################
# The generalized Forster theory is used to calculate the energy transfer rate 
# and time constant between color clusters
# @author : Han Zhang, Jianping GUO 
# use : python  GF.py H_matrixfile Donor_atart Doner_end Accept_start Accept_end spectmp outdir
############################################################

# ------------------------------------------------------------#
# H_matrix_file contains the overall effective Hamiltonian of the PSI system;
# The diagonal element is the site energy;
# The non-diagonal element is the coupling value between pigment i and j;
# ------------------------------------------------------------#

from Forster import Forster
import numpy as np
import sys
import math
from scipy import constants

def cal_GF_k(filename, start1, end1, start2, end2, start, end, step, spec_tmp, outdir):
    F = Forster() 
    a = np.loadtxt(filename)
    shape = np.shape(a)[0]
    D = a[start1:end1, start1:end1]
    D_d, D_v = np.linalg.eigh(D)
    A = a[start2:end2, start2:end2]
    A_d, A_v = np.linalg.eigh(A)

    # cal coupling
    C = a[start1:end1, start2:end2]
    V_m = np.zeros(D_v.shape[0] * A_v.shape[0]).reshape(D_v.shape[0], A_v.shape[0])
    for i in range(D_v.shape[0]):
        D_epi_1 = D_v[:, i].reshape(D_v.shape[0], 1)
        for j in range(A_v.shape[0]):
            A_epi_1 = A_v[:, j].reshape(1, A_v.shape[0])
            V1 = np.dot(D_epi_1, A_epi_1)
            V_m[i][j] = np.sum(V1 * C)
            print("V_m i j is ", V_m[i][j], i, j)

    # cal Boltzmann distribution
    Z = 0.0
    D_them = []
    # Partition function
    for num in D_d:
        tmp = np.exp(- num / (constants.k * 300 / (0.1602177 * 1e-18) * 8065.541))
        D_them.append(tmp)
        Z += tmp
    print("Partition function is: ", Z)
    n = 0.0
    # Boltzmann distribution
    D_boltzmann = []
    for num in D_them:
        D_boltzmann.append(num / Z)
        n += num / Z

    # cal spectral overlap
    for num in D_d:
        file1 = spec_tmp + str(num)+ '_flu.txt'
        F.cal_gau_spec(start, end, step, float(num), file1, 'flu' ) 
    for num in A_d:
        file2 = spec_tmp + str(num)+ '_abs.txt'
        F.cal_gau_spec(start, end, step, float(num), file2, 'abs' )
    J = np.zeros(D_v.shape[0] * A_v.shape[0]).reshape(D_v.shape[0], A_v.shape[0])
    for i in range(D_v.shape[0]):
        for j in range(A_v.shape[0]):
            J[i][j] = F.cal_J(spec_tmp + str(D_d[i])+ '_flu.txt', spec_tmp + str(A_d[j])+ '_abs.txt')
    print("Spectral overlap: \n", J)

    # cal GF kate(unit: ps^-1)
    K = np.zeros(D_v.shape[0] * A_v.shape[0]).reshape(D_v.shape[0], A_v.shape[0])
    h_bar = 1.05457266e-34  # unit: J*s
    cm_J = 1 / 8065.541 * 1.6021766208e-19  # convert cm-1 into J
    # 
    f = open("filelist-gf", 'r')
    filelist = f.readlines()
    f.close()
    Dfilelist = filelist[start1:end1]
    Afilelist = filelist[start2:end2]
    for i in range(D.shape[0]):
        for j in range(A.shape[0]):
            if Dfilelist[i][:3] == Afilelist[j][:3]:
                K[i][j] = (2 * math.pi / h_bar) * D_boltzmann[i] * ((V_m[i][j] * cm_J) ** 2) \
                           * (J[i][j] / 1.9864465642168332e-23) / (10 ** 12)
            else:
                print(Dfilelist[i][:3],Afilelist[j][:3])
                J_ref = 0.0009789061621090092  # Liu data 
                K[i][j] = (2 * math.pi / h_bar) * D_boltzmann[i] * ((V_m[i][j] * cm_J) ** 2) \
                          * (J_ref / 1.9864465642168332e-23) / (10 ** 12)

    # exciton inf file
    np.savetxt(outdir + "tmp_d.txt", D, fmt="%20.10f")
    np.savetxt(outdir + "tmp_a.txt", A, fmt="%20.10f")
    np.savetxt(outdir + "tmp_c.txt", C, fmt="%15.10f")
    np.savetxt(outdir + "tmp_d_d.txt", D_d, fmt="%15.10f")
    np.savetxt(outdir + "tmp_d_v.txt", D_v.T, fmt="%15.10f")
    np.savetxt(outdir + "tmp_a_d.txt", A_d, fmt="%15.10f")
    np.savetxt(outdir + "tmp_a_v.txt", A_v.T, fmt="%15.10f")
    np.savetxt(outdir + "tmp_coup.txt", V_m, fmt="%15.10f")
    np.savetxt(outdir + "tmp_j.txt", J, fmt="%15.10f")
    np.savetxt(outdir + "tmp_k.txt", K, fmt="%15.10f")
    file8 = open(outdir + "tmp_d.txt", "r").readlines()
    inf8 = ""
    for i in file8:
        inf8 += i

    file1 = open(outdir + "tmp_d_d.txt", "r").readlines()
    file2 = open(outdir + "tmp_d_v.txt", "r").readlines()
    inf1 = ""
    for i in range(D_d.shape[0]):
        inf1 += '   ' + file1[i].strip() + '   ' + str(round(D_boltzmann[i], 10)) + '     ' + file2[i]
    file9 = open(outdir + "tmp_a.txt", "r").readlines()
    inf9 = ""
    for i in file9:
        inf9 += i

    file3 = open(outdir + "tmp_a_d.txt", "r").readlines()
    file4 = open(outdir + "tmp_a_v.txt", "r").readlines()
    inf2 = ""
    for i in range(A_d.shape[0]):
        inf2 += '   ' + file3[i].strip() + '       ' + file4[i]

    file10 = open(outdir + "tmp_c.txt", "r").readlines()
    inf10 = ""
    for i in file10:
        inf10 += i
    file5 = open(outdir + "tmp_coup.txt", "r").readlines()
    inf5 = ""
    for i in file5:
        inf5 += i
    file6 = open(outdir + "tmp_j.txt", "r").readlines()
    inf6 = ""
    for i in file6:
        inf6 += i
    file7 = open(outdir + "tmp_k.txt", "r").readlines()
    inf7 = ""
    for i in file7:
        inf7 += i

    with open(outdir + "exciton_inf.txt", "w") as f:
        f.write("".join("Donor matrix(cm^-1):\n" + inf8 + "\n"))
        f.write("".join("Donor exciton states(cm^-1):\n" + inf1 + "\n"))
        f.write("".join("Acceptor matrix(cm^-1):\n" + inf9 + "\n"))
        f.write("".join("Acceptor exciton states(cm^-1):\n" + inf2 + "\n"))
        f.write("".join("Donor-Acceptor coups(cm^-1):\n" + inf10 + "\n"))
        f.write("".join("Exciton coups(cm^-1):\n" + inf5 + "\n"))
        f.write("".join("Spectral overlaps:\n" + inf6 + "\n"))
        f.write("".join("Transfer rates(ps^-1):\n" + inf7 + "\n"))
        f.write("Kate all(ps^-1): " + str(np.sum(K)) + '\n')
        f.write("Time constants(ps): " + str(1 / np.sum(K)))
    f.close()

    return K, np.sum(K)

###########################################################################
start, end, step = 10000, 20000, 10000
print(sys.argv)
K, all_K = cal_GF_k(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]),start, end, step, sys.argv[6], sys.argv[7])
print("Rate all(ps^-1): ", np.sum(K))
print("Time constants(ps): ", 1 / np.sum(K))



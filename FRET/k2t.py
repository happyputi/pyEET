#######################################################
# Trans k matrix to t matrix
# t matrix is time matrix
# @author : Jianping GUO 
# use: python3 k2t.py k_filename t_filename filelistname
#######################################################


import numpy as np
import sys

def tarns_km_tm(k_filename, t_filename):
    k_matrix = np.loadtxt(k_filename)
    x_shape = k_matrix.shape[0]
    y_shape = k_matrix.shape[1]
    t_matrix = np.zeros((x_shape,y_shape))
    for i in range(x_shape):
        for j in range(y_shape):
            if k_matrix[i][j]==0.0:
                t_matrix[i][j]=100.0
            else:
                t_matrix[i][j]=1/k_matrix[i][j]
    np.savetxt(t_filename, t_matrix, fmt="%20.10f")
    print(t_filename, 'is done!\n')

def getcontact(t_filename, filelistname):
    f = open(filelistname, 'r')
    filelist = f.readlines()
    f.close()
    t_matrix = np.loadtxt(t_filename)
    x_shape = t_matrix.shape[0]
    y_shape = t_matrix.shape[1]
    f_less1 = open("t_less1.dat", 'w')
    f_1t10 = open("t_b1t10.dat",'w')
    f_10t20 = open("t_b10t20.dat",'w')
    for i in range(x_shape):
        for j in range(i+1, y_shape):
            if t_matrix[i][j] < t_matrix[j][i]:
                T = t_matrix[i][j]
                if T <= 1.0:
                    line = ' '.join(filelist[i][:9].split('-')) + ' '+ ' '.join(filelist[j][:9].split('-')) + ' ' + str(T) + '\n' 
                    f_less1.write(line)
                elif 1.0 < T <= 10.0:
                    line = ' '.join(filelist[i][:9].split('-')) + ' '+ ' '.join(filelist[j][:9].split('-')) + ' ' + str(T) + '\n'
                    f_1t10.write(line)
                elif 10.0 < T <= 20.0:
                    line = ' '.join(filelist[i][:9].split('-')) + ' '+ ' '.join(filelist[j][:9].split('-')) + ' ' + str(T) + '\n'
                    f_10t20.write(line)
                else:
                    pass
            else:
                T = t_matrix[j][i]
                if T <= 1.0:
                    line = ' '.join(filelist[j][:9].split('-')) + ' '+ ' '.join(filelist[i][:9].split('-')) + ' ' + str(T) + '\n'
                    f_less1.write(line)
                elif 1.0 < T <= 10.0:
                    line = ' '.join(filelist[j][:9].split('-')) + ' '+ ' '.join(filelist[i][:9].split('-')) + ' ' + str(T) + '\n'
                    f_1t10.write(line)
                elif 10.0 < T <= 20.0:
                    line = ' '.join(filelist[j][:9].split('-')) + ' '+ ' '.join(filelist[i][:9].split('-')) + ' ' + str(T) + '\n'
                    f_10t20.write(line)
                else:
                    pass

    f_less1.close() 
    f_1t10.close()
    f_10t20.close()
    print("contactfile is make!\n")

if __name__ == '__main__':
    print("trans k matrix to t matrix\n") 
    tarns_km_tm(sys.argv[1], sys.argv[2])
    print("get contactfile\n")
    getcontact(sys.argv[2], sys.argv[3])
    
    

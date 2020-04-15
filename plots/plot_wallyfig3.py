import matplotlib.pyplot as plt
import numpy as np

full_off_con = []
full_off_lin = []

full_noO_con = []
full_noO_lin = []

_050_off_con = []
_050_off_lin = []
_050_noO_con = []
_050_noO_lin = []

_100_off_con = []
_100_off_lin = []
_100_noO_con = []
_100_noO_lin = []

_200_off_con = []
_200_off_lin = []
_200_noO_con = []
_200_noO_lin = []

_400_off_con = []
_400_off_lin = []
_400_noO_con = []
_400_noO_lin = []

_500_off_con = []
_500_off_lin = []
_500_noO_con = []
_500_noO_lin = []

_600_off_con = []
_600_off_lin = []
_600_noO_con = []
_600_noO_lin = []

_750_off_con = []
_750_off_lin = []
_750_noO_con = []
_750_noO_lin = []

f2p = []
f2nf2p_con = []
f2nf2p_lin = []

xs = []; F2ds_con = []; F2ds_lin = []
with open("deuterium-data/deut-conv-gt.txt","r") as f:
    for line in f:
        parse = line.strip().split(" ")
        xs.append(      float(parse[0]) )

        full_off_con.append(    float(parse[1]  ))
        full_off_lin.append(    float(parse[2]  ))

        full_noO_con.append(    float(parse[3]  ))
        full_noO_lin.append(    float(parse[4]  ))

        _050_off_con.append(    float(parse[5]  ))
        _050_off_lin.append(    float(parse[6]  ))
        _050_noO_con.append(    float(parse[7]  ))
        _050_noO_lin.append(    float(parse[8]  ))

        _100_off_con.append(    float(parse[9]  ))
        _100_off_lin.append(    float(parse[10]  ))
        _100_noO_con.append(    float(parse[11]  ))
        _100_noO_lin.append(    float(parse[12]  ))

        _200_off_con.append(    float(parse[13]  ))
        _200_off_lin.append(    float(parse[14]  ))
        _200_noO_con.append(    float(parse[15]  ))
        _200_noO_lin.append(    float(parse[16]  ))

        _400_off_con.append(    float(parse[17]  ))
        _400_off_lin.append(    float(parse[18]  ))
        _400_noO_con.append(    float(parse[19]  ))
        _400_noO_lin.append(    float(parse[20]  ))

        _500_off_con.append(    float(parse[21]  ))
        _500_off_lin.append(    float(parse[22]  ))
        _500_noO_con.append(    float(parse[23]  ))
        _500_noO_lin.append(    float(parse[24]  ))

        _600_off_con.append(    float(parse[25]  ))
        _600_off_lin.append(    float(parse[26]  ))
        _600_noO_con.append(    float(parse[27]  ))
        _600_noO_lin.append(    float(parse[28]  ))

        _750_off_con.append(    float(parse[29]  ))
        _750_off_lin.append(    float(parse[30]  ))
        _750_noO_con.append(    float(parse[31]  ))
        _750_noO_lin.append(    float(parse[32]  ))

        f2p.append( float(parse[33]) )
        f2nf2p_con.append(  float(parse[34]) )
        f2nf2p_lin.append(  float(parse[35]) )
 

full_off_con = np.asarray(full_off_con)
full_off_lin = np.asarray(full_off_lin)
                          
full_noO_con = np.asarray(full_noO_con)
full_noO_lin = np.asarray(full_noO_lin)
                          
_050_off_con = np.asarray(_050_off_con)
_050_off_lin = np.asarray(_050_off_lin)
_050_noO_con = np.asarray(_050_noO_con)
_050_noO_lin = np.asarray(_050_noO_lin)
                          
_100_off_con = np.asarray(_100_off_con)
_100_off_lin = np.asarray(_100_off_lin)
_100_noO_con = np.asarray(_100_noO_con)
_100_noO_lin = np.asarray(_100_noO_lin)
                          
_200_off_con = np.asarray(_200_off_con)
_200_off_lin = np.asarray(_200_off_lin)
_200_noO_con = np.asarray(_200_noO_con)
_200_noO_lin = np.asarray(_200_noO_lin)
                          
_400_off_con = np.asarray(_400_off_con)
_400_off_lin = np.asarray(_400_off_lin)
_400_noO_con = np.asarray(_400_noO_con)
_400_noO_lin = np.asarray(_400_noO_lin)
                          
_500_off_con = np.asarray(_500_off_con)
_500_off_lin = np.asarray(_500_off_lin)
_500_noO_con = np.asarray(_500_noO_con)
_500_noO_lin = np.asarray(_500_noO_lin)
                          
_600_off_con = np.asarray(_600_off_con)
_600_off_lin = np.asarray(_600_off_lin)
_600_noO_con = np.asarray(_600_noO_con)
_600_noO_lin = np.asarray(_600_noO_lin)
                          
_750_off_con = np.asarray(_750_off_con)
_750_off_lin = np.asarray(_750_off_lin)
_750_noO_con = np.asarray(_750_noO_con)
_750_noO_lin = np.asarray(_750_noO_lin)

f2p = np.asarray(f2p)
f2nf2np_con = np.asarray(f2nf2p_con)
f2nf2np_lin = np.asarray(f2nf2p_lin)



# Make (F2p + F2n)
f2pPf2n = f2p * ( np.ones(len(f2p)) +f2nf2p_con )
# Correct the returned F2d/(F2p+F2n) to F2d
f2d = full_off_con * f2pPf2n


plt.plot( xs, (f2d - f2pPf2n)/f2d ,color='red')
plt.plot( xs, ( _200_off_con*f2pPf2n - 0.038*f2pPf2n) / f2d  ,color='green')
plt.plot( xs, ( ( full_off_con - _200_off_con ) * f2pPf2n - (1-0.038)*f2pPf2n ) / f2d, color='blue')
plt.axhline(y=0,linestyle='--',linewidth=3,color='black')

plt.xlim([0.2,0.85])
plt.ylim([-0.06,0.06])
plt.show()

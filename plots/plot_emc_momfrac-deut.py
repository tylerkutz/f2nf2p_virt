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

xs = []; F2ds_con = []; F2ds_lin = []
with open("deut-conv-lt.txt","r") as f:
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



plt.figure(1)
plt.plot(xs,full_off_lin,color='orange',linewidth=3,linestyle='-',label='Conv. Linear')
plt.plot(xs,full_off_con,color='green',linewidth=3,linestyle='-',label='Conv. Const.')
plt.plot(xs,full_noO_con,color='purple',linewidth=3,linestyle='--',label='Conv. No Off.')
plt.plot(xs,_200_noO_con,color='purple',linewidth=3,linestyle=':',label='Conv. No Off. < 200MeV')
plt.plot(xs,_100_noO_con,color='blue',linewidth=3,linestyle=':',label='Conv. No Off. < 100MeV')
plt.legend(numpoints=1,loc=2,fontsize=15)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.ylim([0.5,1.2])
plt.xlim([0.2,0.95])
plt.xlabel('x',fontsize=16)
plt.ylabel('F2d/(F2p+F2n)',fontsize=16)
plt.axhline(y=1,linestyle='--',color='black',linewidth=3)
plt.tight_layout()
plt.savefig('plot-f2d-conv.pdf')

plt.figure(5)
plt.plot(xs,    _200_off_con    ,linestyle='--'    ,color='blue',   label='W/ Off < 200MeV')
plt.plot(xs,    _400_off_con    ,linestyle='--'    ,color='green',   label='W/ Off < 400MeV')
plt.plot(xs,    _750_off_con    ,linestyle='--'    ,color='purple',   label='W/ Off < 750MeV')
plt.plot(xs,    full_off_con    ,linestyle='-'     ,color='red',   label='W/ Off (All p)')

plt.legend(numpoints=1,loc=2,fontsize=15)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.ylim([0.9,1.2])
plt.xlim([0.2,0.95])
plt.xlabel('x',fontsize=16)
plt.ylabel('F2d/(F2p+F2n)',fontsize=16)
plt.axhline(y=1,linestyle='--',color='black',linewidth=3)
plt.tight_layout()

plt.figure(6)
plt.plot(xs,    _200_noO_con    ,linestyle='--'    ,color='blue',   label='W/o Off < 200MeV')
plt.plot(xs,    _400_noO_con    ,linestyle='--'    ,color='green',   label='W/o Off < 400MeV')
plt.plot(xs,    _750_noO_con    ,linestyle='--'    ,color='purple',   label='W/o Off < 750MeV')
plt.plot(xs,    full_noO_con    ,linestyle='-'     ,color='red',   label='W/o Off (All p)')

plt.legend(numpoints=1,loc=2,fontsize=15)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.ylim([0.9,1.2])
plt.xlim([0.2,0.95])
plt.xlabel('x',fontsize=16)
plt.ylabel('F2d/(F2p+F2n)',fontsize=16)
plt.axhline(y=1,linestyle='--',color='black',linewidth=3)
plt.tight_layout()


plt.figure(4)
plt.plot(xs,    full_off_con - full_noO_con ,linestyle='-',color='red', label='W/ - W/o Off (All p)')
plt.plot(xs,    _200_off_con - _200_noO_con ,linestyle='--',color='blue', label='W/ - W/o Off < 200')
plt.plot(xs,    _400_off_con - _400_noO_con ,linestyle='--',color='green', label='W/ - W/o Off < 400')
plt.plot(xs,    _750_off_con - _750_noO_con ,linestyle='--',color='purple', label='W/ - W/o Off < 750')
plt.legend(numpoints=1,loc=3,fontsize=15)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
#plt.ylim([0.9,1.2])
plt.xlim([0.2,0.95])
plt.xlabel('x',fontsize=16)
#plt.ylabel('F2d/(F2p+F2n)',fontsize=16)
#plt.axhline(y=1,linestyle='--',color='black',linewidth=3)
plt.tight_layout()

plt.figure(3)
plt.plot(xs,    (_200_off_con - _200_noO_con)/(full_off_con - full_noO_con ) ,linestyle='--',color='blue', label='(W/-W/o Off < 200)/(W-W/o Off All)')
plt.plot(xs,    (_400_off_con - _400_noO_con)/(full_off_con - full_noO_con ) ,linestyle='--',color='green', label='<400')
plt.plot(xs,    (_750_off_con - _750_noO_con)/(full_off_con - full_noO_con ) ,linestyle='--',color='purple', label='<750')
plt.legend(numpoints=1,loc=3,fontsize=15)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.ylim([0.,1.])
plt.xlim([0.2,0.95])
plt.xlabel('x',fontsize=16)
#plt.ylabel('F2d/(F2p+F2n)',fontsize=16)
#plt.axhline(y=1,linestyle='--',color='black',linewidth=3)
plt.tight_layout()




plt.show()

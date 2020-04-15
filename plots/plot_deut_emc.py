import matplotlib.pyplot as plt
import numpy as np
def readTheory(x,F2d,filename):
    with open(filename,"r") as f:
        for line in f:
            parse = line.strip().split(" ")
            x.append(       float(parse[0]) )
            F2d.append(     float(parse[1]) )

xs = []; F2ds_con = []; F2ds_lin = []
with open("deut-conv.txt","r") as f:
    for line in f:
        parse = line.strip().split(" ")
        xs.append(      float(parse[0]) )
        F2ds_con.append(     float(parse[1]) )
        F2ds_lin.append(    float(parse[2]) )
x_prl = []; F2ds_L = []; F2ds_H = [];
with open("deut-prl.txt","r") as f:
    for line in f:
        parse = line.strip().split(" ")
        x_prl.append(      float(parse[0]) )
        F2ds_L.append(  1./float(parse[1]) )
        F2ds_H.append(  1./float(parse[2]) )

with open("deut-bonus.txt","r") as f:
    data = f.readlines()
x_bon   = np.asarray(data[0].strip().split(" "),dtype=float)
f2d_bon = np.asarray(data[2].strip().split(" "),dtype=float)
e1_bon  = np.asarray(data[3].strip().split(" "),dtype=float)
e2_bon  = np.asarray(data[4].strip().split(" "),dtype=float)
e3_bon  = np.asarray(data[5].strip().split(" "),dtype=float)
e_bon = np.sqrt( np.power(e1_bon,2) + np.power(e2_bon,2) + np.power(e3_bon,2) )

x_cj = [];  F2d_cj = []
x_mst = []; F2d_mst = []
x_kp = [];  F2d_kp = []
readTheory(x_cj,F2d_cj,"deut-cj15.txt")
readTheory(x_kp,F2d_kp,"deut-kp.txt")
readTheory(x_mst,F2d_mst,"deut-mst.txt")

plt.errorbar(x_bon,f2d_bon,yerr=e_bon,linestyle='',color='black',marker='o')

plt.plot(x_cj,F2d_cj,color='green',linewidth=3,linestyle='--',label='CJ15')
plt.plot(x_mst,F2d_mst,color='red',linewidth=3,linestyle='--',label='MST')
plt.plot(x_kp,F2d_kp,color='blue',linewidth=3,linestyle='--',label='KP')
plt.plot(xs,F2ds_lin,color='orange',linewidth=3,linestyle='-',label='Conv. Linear')
plt.plot(xs,F2ds_con,color='purple',linewidth=3,linestyle='-',label='Conv. Const.')
plt.fill_between(x_prl,F2ds_L,F2ds_H,color='cornflowerblue',label='PRL Work')

plt.legend(numpoints=1,loc=2,fontsize=15)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.ylim([0.9,1.2])
plt.xlim([0.2,0.95])
plt.xlabel('x',fontsize=16)
plt.ylabel('F2d/(F2p+F2n)',fontsize=16)
plt.axhline(y=1,linestyle='--',color='black',linewidth=3)
plt.tight_layout()
plt.savefig('plot-f2d.pdf',bbox_inches='tight')


plt.figure(9)

plt.show()

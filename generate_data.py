from cmath import sqrt
from matplotlib.pyplot import close
import numpy as np

with open("double_star_system.dat",'w') as dss:
    # 设定双星系统w=1，质量m1=1，m2=3,R=4
    timeline=np.linspace(0,7.1)
    x1=[]
    y1=[]
    z1=[]
    x2=[]
    y2=[]
    z2=[]
    for t in timeline:
        x1.append(str(round(3*np.cos(t),5)))
        y1.append(str(round(3*np.sin(t),5)))
        z1.append("0")
        x2.append(str(round(-np.cos(t),5)))
        y2.append(str(round(-np.sin(t),5)))
        z2.append("0")

    for i in range(len(timeline)):
        dss.write(x1[i]+"\t\t"+y1[i]+"\t\t"+z1[i]+"\t\t")
        dss.write(x2[i]+"\t\t"+y2[i]+"\t\t"+z2[i]+"\t\t")
        dss.write(str(timeline[i])+"\n")


pi = 3.1415926

with open("3_body_system_equal_mass.dat",'w') as tbs:
    #一种三体运动的特解
    timeline=np.linspace(0,7.1)
    x1=[]
    y1=[]
    x2=[]
    y2=[]
    x3=[]
    y3=[]
    for t in timeline:
        r=2
        x1.append(str(round(r*np.cos(t),5)))
        y1.append(str(round(r*np.sin(t),5)))
        x2.append(str(round(r*np.cos(t+2*pi/3),5)))
        y2.append(str(round(r*np.sin(t+2*pi/3),5)))
        x3.append(str(round(r*np.cos(t-2*pi/3),5)))
        y3.append(str(round(r*np.sin(t-2*pi/3),5)))

    for i in range(len(timeline)):
        tbs.write(x1[i]+"\t\t"+y1[i]+"\t\t"+"0"+"\t\t")
        tbs.write(x2[i]+"\t\t"+y2[i]+"\t\t"+"0"+"\t\t")
        tbs.write(x3[i]+"\t\t"+y3[i]+"\t\t"+"0"+"\t\t")
        tbs.write(str(timeline[i])+"\n")



with open("3_body_system.dat",'w') as tbs:
    #一种三体运动的特解
    timeline=np.linspace(0,7.1,1000)
    x1=[]
    y1=[]
    x2=[]
    y2=[]
    x3=[]
    y3=[]
    for t in timeline:
        r1=1
        r2=np.sqrt(7/3)


        x1.append(str(round(r1*np.cos(t),5)))
        y1.append(str(round(r1*np.sin(t),5)))
        x2.append(str(round(r2*np.cos(t+2.2845207),5)))
        y2.append(str(round(r2*np.sin(t+2.2845207),5)))
        x3.append(str(round(r2*np.cos(t-2.2845207),5)))
        y3.append(str(round(r2*np.sin(t-2.2845207),5)))


    for i in range(len(timeline)):
        tbs.write(x1[i]+"\t\t"+y1[i]+"\t\t"+"0"+"\t\t")
        tbs.write(x2[i]+"\t\t"+y2[i]+"\t\t"+"0"+"\t\t")
        tbs.write(x3[i]+"\t\t"+y3[i]+"\t\t"+"0"+"\t\t")
        tbs.write(str(timeline[i])+"\n")


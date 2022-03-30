import numpy as np
from numpy.linalg import *
from scipy import interpolate
import matplotlib.pyplot as plt


class Star:
    star_num=0

    # x     array       观测坐标点
    # y     array
    # z     array
    # t     array       观测时间
    # mass  float       星体质量
    # tck   tuple array 样条插值的特征量
    # order float       星体编号


    def __init__(self) -> None:
        self.x=[]
        self.y=[]
        self.z=[]
        self.t=[]
        self.mass=1 #default value
        self.tck=[(),(),()]
        self.order=Star.star_num
        Star.star_num += 1

    def append_data(self,x,y,z,t):
        self.x.append(x)
        self.y.append(y)
        self.z.append(z)
        self.t.append(t)
    
    def give_mass(self,mass):
        self.mass=mass

    def get_data(self,flag):
        if flag=="x":
            return self.x
        elif flag=="y":
            return self.y
        elif flag=="z":
            return self.z
        elif flag=="t":
            return self.t
        elif flag=="order":
            return self.order
        elif flag=="tck":
            return self.tck
        

class Observer:
    # star_list         array   存放Star的实例化对象
    # star_num          int     星星数量
    # timeline_step     float   绘图时的步长
    # x                 array   绘图时的x坐标点(一次索引的下标是星星编号，二次索引的下标是时间)
    # y                 array   绘图时的y坐标点
    # z                 array   绘图时的z坐标点

    def __init__(self,star_num) -> None:
        self.star_list=[]
        self.star_num=star_num
        for i in range(star_num):
            self.star_list.append(Star())
    
    def read_data_from_file(self,filename='data.dat'):
        with open(filename,'r') as file:
            for lines in file:
                word=lines.split()
                
                if len(word) != 3*self.star_num + 1:
                    assert 0==1
                
                for i in range(0,self.star_num):
                    self.star_list[i].append_data(float(word[3*i]),float(word[3*i+1]),float(word[3*i+2]),float(word[3*self.star_num]))
    
    def interpolate_BSpline(self):
        for i in range (self.star_num):
            tck1=interpolate.splrep(self.star_list[i].get_data("t"),self.star_list[i].get_data("x"),k=5,s=0.2)
            tck2=interpolate.splrep(self.star_list[i].get_data("t"),self.star_list[i].get_data("y"),k=5,s=0.2)
            tck3=interpolate.splrep(self.star_list[i].get_data("t"),self.star_list[i].get_data("z"),k=5,s=0.2)
            self.star_list[i].tck=[tck1,tck2,tck3]

    def center_div(self,func_values,step):
        center_div=[]
        for j in range(1,len(func_values)-1):
            center_div.append((func_values[j+1]-func_values[j-1]/(2*step)))
        return center_div

# 求解齐次线性方程组，原理是求解特征向量后得到一个非零解。调整affordable_miss来得到模糊解
    def solve_null_space(self,associate_matrix, affordable_miss):
        _, miss, character = np.linalg.svd(associate_matrix)
        null_space = np.compress(miss <= affordable_miss, character, axis=0)
        return null_space.T

#————————以下两个方法可以被 precision_analysis_2D 代替————————#        
    def precision_analysis_2Body(self,affordable_miss=1e-2):
        # 计算中心差分
        div_x=[]
        div_y=[]
        
        for order in range(len(self.x)):
            div_x.append(self.center_div(self.x[order],self.timeline_step))
            div_y.append(self.center_div(self.y[order],self.timeline_step))
        
        # 计算偏微分=0的方程的各项系数
        # [V11 V12][m1]-[0]
        # [V21 V22][m2]-[0] 

        V11=[]
        V12=[]
        V22=[]
        for t in range(len(div_x[0])):
            V11.append(div_x[0][t]**2+div_y[0][t]**2)
            V12.append(div_x[0][t]*div_x[1][t]+div_y[0][t]*div_y[1][t])
            V22.append(div_x[1][t]**2+div_y[1][t]**2)

        final_solution=np.array([[0.],[0.]])

        for t in range(len(V11)):
            associate_matrix=np.array([     [ V11[t], V12[t] ],[ V12[t], V22[t] ]     ])
            local_solution=self.solve_null_space(associate_matrix,affordable_miss)
            final_solution+=local_solution
         
        final_solution=final_solution/final_solution[0][0]
        print(final_solution)
        return final_solution

    def precision_analysis_3Body(self,affordable_miss=1e-2):
        # 计算中心差分
        div_x=[]
        div_y=[]
        
        for order in range(len(self.x)):
            div_x.append(self.center_div(self.x[order],self.timeline_step))
            div_y.append(self.center_div(self.y[order],self.timeline_step))
        
        # 计算偏微分=0的方程的各项系数
        # [V11 V12 V13][m1] [0]
        # [V21 V22 V23][m2]=[0]
        # [V31 V32 V33][m3] [0]

        V11=[]
        V12=[]
        V13=[]
        V22=[]
        V23=[]
        V33=[]
        for t in range(len(div_x[0])):
            V11.append(div_x[0][t]**2+div_y[0][t]**2)
            V12.append(div_x[0][t]*div_x[1][t]+div_y[0][t]*div_y[1][t])
            V13.append(div_x[0][t]*div_x[2][t]+div_y[0][t]*div_y[2][t])
            V22.append(div_x[1][t]**2+div_y[1][t]**2)
            V23.append(div_x[1][t]*div_x[2][t]+div_y[1][t]*div_y[2][t])
            V33.append(div_x[2][t]**2+div_y[2][t]**2)

        final_solution=np.array([[0.],[0.],[0.]])

        for t in range(len(V11)):
            associate_matrix=np.array([[ V11[t], V12[t],V13[t]],[ V12[t], V22[t] ,V23[t]]],[V13[t],V23[t],V33[t]])
            local_solution=self.solve_null_space(associate_matrix,affordable_miss)
            final_solution+=local_solution
         
        final_solution=final_solution/final_solution[0][0]
        print(final_solution)
        return final_solution
#--------以上两个方法可以被 precision_analysis_2D 代替--------#

    def precision_analysis_2D(self,affordable_miss=1e-2):
        # 计算中心差分
        div_x=[]
        div_y=[]
        
        for order in range(len(self.x)):
            div_x.append(self.center_div(self.x[order],self.timeline_step))
            div_y.append(self.center_div(self.y[order],self.timeline_step))
        
        final_solution=np.zeros((self.star_num,1))
        associate_matrix=np.zeros((self.star_num,self.star_num))
        for t in range(len(div_x[0])):
            for i in range(self.star_num):
                for j in range(self.star_num):
                    associate_matrix[i][j]=(div_x[i][t]*div_x[j][t]+div_y[i][t]*div_y[j][t])
            
            local_solution=self.solve_null_space(associate_matrix,affordable_miss)
            if local_solution.size!=self.star_num:
                continue
            final_solution+=local_solution

        final_solution=final_solution/final_solution[0][0]
        print(final_solution)
        return final_solution

    def plot2_BSpline(self,predict_scale,mode):
        t_max=max(self.star_list[0].get_data("t"))
        t_min=min(self.star_list[0].get_data("t"))

        if mode=="predict":
            timeline,step=np.linspace(t_max,t_max+predict_scale*(t_max-t_min),retstep=True)
        elif mode=="transition":
            timeline,step=np.linspace(t_max-predict_scale*(t_max-t_min),t_max+predict_scale*(t_max-t_min),retstep=True)
        else:
            timeline,step=np.linspace(t_min-predict_scale*(t_max-t_min),t_max+predict_scale*(t_max-t_min),retstep=True)

        self.timeline_step=step #存储步长，用于计算差分。

        x=[]
        y=[] #用于存储生成的各星体在对应时刻的横纵坐标列表（列表的列表）。
        for order in range(self.star_num):
            x.append(interpolate.splev(timeline,self.star_list[order].tck[0]))
            y.append(interpolate.splev(timeline,self.star_list[order].tck[1]))
            plt.plot(x[order],y[order]) #绘图
        
        #将数组存储起来，便于调用
        self.x=x
        self.y=y

        plt.plot(np.cos(timeline),np.sin(timeline)) # 描画圆周运动对比图
        # plt.show()

    def plot3_BSpline(self,predict_scale,mode):
        t_max=max(self.star_list[0].get_data("t"))
        t_min=min(self.star_list[0].get_data("t"))

        if mode=="predict":
            timeline,step=np.linspace(t_max,t_max+predict_scale*(t_max-t_min),retstep=True)
        elif mode=="transition":
            timeline,step=np.linspace(t_max-predict_scale*(t_max-t_min),t_max+predict_scale*(t_max-t_min),retstep=True)
        else:
            timeline,step=np.linspace(t_min-predict_scale*(t_max-t_min),t_max+predict_scale*(t_max-t_min),retstep=True)

        self.timeline_step=step #存储步长，用于计算差分。

        x=[]
        y=[] 
        z=[]#用于存储生成的各星体在对应时刻的横纵坐标列表（列表的列表）。
        for order in range(self.star_num):
            x.append(interpolate.splev(timeline,self.star_list[order].tck[0]))
            y.append(interpolate.splev(timeline,self.star_list[order].tck[1]))
            z.append(interpolate.splev(timeline,self.star_list[order].tck[2]))
            plt.plot(x[order],y[order],z[order])

        self.x=x
        self.y=y
        self.z=z    
        # plt.show()  

# 圆周运动图形模拟
    def run2D(self,filename,predict_scale=0,mode="display"):
        self.read_data_from_file(filename=filename)
        self.interpolate_BSpline()
        self.plot2_BSpline(predict_scale=predict_scale,mode=mode)
        final_solution=self.precision_analysis_2D(affordable_miss=0.01)
        plt.xlabel("x")
        plt.ylabel("y")
        plt.text(0,0,str(final_solution.T))
        plt.show()

#--------以下为圆周运动模拟--------#
# observer=Observer(2)
# observer.run2D("circle.dat",predict_scale=-0.05,mode="display")

#--------以下为双星系统模拟--------#
# observer=Observer(2)
# observer.run2D("double_star_system.dat",predict_scale=-0.05,mode="display")

#--------以下为三星系统模拟--------#
observer=Observer(3)
observer.run2D("3_body_system.dat",predict_scale=-0.05,mode="display")



# 三体运动的（非动力学）拟合以及质量比估计

**摘要**：本项目使用**样条插值法**对星体位置的理想天文观测数据进行插值拟合，得到多星体运动的非动力学拟合结果，并使用**数值微分法**由动量守恒方法估计各个星体之间的质量比值，使得星体运动轨迹的拟合、预测以及质量比估计不依赖与星体之间作用力的形式。并且可将实际星体的质量与估计质量进行比较，来评估预测的精度。



## 问题概述

本项目系南京大学2020级计算物理导论课程作业——第一次作业，经助教允许可以将作业代码公开发布。

本项目尝试绕过行星动力学，用数值方法对无解析解的三体问题进行一些拟合、预测与估计。

显然，对于万有引力下的多体问题有动力学方程：
$$
m_i\ddot{\boldsymbol{r}}_i=G\sum_j\frac{m_im_j}{|\boldsymbol r_i-\boldsymbol r_j|^3}(\boldsymbol r_i-\boldsymbol r_j)
$$
由这些方程和初始条件，我们可以用数值方法解微分方程来解决三体问题。

如果我们并不了解我们需要研究的星体之间的动力学呢？（比如星体之间的作用力形式未知）

那么我们只能本着一件事：天文观测数据（实验）来进行预测。



### 非动力学拟合

广泛地说，如果我们有一系列天文观测者观测得到的N星坐标以及观测时间构成的四元组，这就标志着当前N星的一个状态，也就是：
$$
\psi(t)=\{(x_i,y_i,z_i),t\}
$$
其中各坐标显然都是时间的函数。

这里我们要做出一个假设：时空状态的变化是连续的，因为星星们（按照观测）总是这样连续地改变位置，时间也是均匀流动的。

在这种情况下，我们把N星的坐标和时间数据作插值处理，得到3N个样条插值函数。

这样我们就得到了星体的运动方程。我们将拟合函数作图，就可以得到N体运动的轨迹。

运用简单的中心差分方法，我们就能得到星体系统的运行速度等各项运动参数，而不依赖于N星系统当地的力的形式。



### 质量比估计和精确度判断

另外，我们还担心这样的拟合是否合适？我们可以适当地引入一些经典动力学的广泛标准，比如一些守恒律作为评判标准。

由于对于太空中不受外力的物体，不论三个星体之间力的形式如何，他们总要遵守动量守恒定律，即：
$$
m_1\boldsymbol v_1+m_2\boldsymbol v_2+m_3\boldsymbol v_3=\boldsymbol P
$$
一般情况下，我们并不知道星体的质量如何。因此我们可以使用观测数据对星体的质量之比进行估计。

回到动量守恒方程，其三个分量应该分别守恒，模平方也自然守恒。

我们尝试对拟合函数作差分得到速度，将动量模平方对质量求各偏导数，可以得到使总守恒量变化最小的质量方程，即：
$$
V_{ij}\boldsymbol m_j=0\\
where\ v_{ij}=\dot x_i\dot x_j+\dot y_i\dot y_j+\dot z_i\dot z_j
$$
问题归结为求上述系数矩阵 V 的零空间。我们用特征向量，特征值方法可以求的一个模糊解。

（由于计算精度，V不一定精确满足不满秩的条件，而是粗略满足，因此通过放大误差条件来得到解）

不妨设定第一个星体的质量为1，可以解得其他星体与第一个星体的质量之比。多次取点求平均值，就能得到各个星体质量比的估计值。

如果我们知道各个星体的质量的话，那我们就可以用这个预测值和真实值进行比较，以评估我们的拟合精确度。



## 程序概述

程序编写上，我们编写一个N星系统的模拟程序，然后基于计算简洁性，对双星系统进行物理模拟，然后与我们样条模拟的结果进行对比。

然后我们从已有的三体问题的特殊解上采集数据点，放到我们的程序中进行计算处理，来看看是否在一段时间内符合三体问题的解，以及质量估计是否正确

另外，由于星体数小于等于3的问题可以是平面问题，我们就选取平面解（即将Z置为0）来方便展示。

程序结构如下：

```
3_body system
	|___N_body_system.py
	|___generate_data.py
	|___circle.dat
	|___double_star_system.dat
	|___3_body_system.dat
	|___graph
		|___circle.png
		|___double_star_system.png
		|___3_body_system.png
```

N_body_system.py 是主模拟程序，从根目录下读取观测数据，输出拟合图像。

generate_data.py 是数据生成程序，除 circle.dat 是用计算器计算以外，其他都是由此程序生成。

下面3个 .dat 文件是拟合用的数据。



## 测试样例

测试样例的数据文件位于根目录中，拟合图像为 graph 文件夹中的同名文件。

### 圆周运动模拟

为观察简便起见，我们先研究一个二体问题，中心天体不动，将行星运动设为：
$$
(x,y)=(\cos t,\sin t)
$$
并输入约30个观测数据。

取 predict_scale =-0.05 mode="display" 绘制得到如下图案：

其中黄色是预测轨迹，绿色是标准圆周运动的轨迹。

可以看到对路径的拟合效果还是不错的，最大部分的水平偏差也未超过轨道半径的二十分之一。

程序也正确地输出了质量估计值：【1，0】，这符合我们中心天体不动（质量无穷大）的假设。



### 双星系统

由动力学我们可以解出双星系统的运行轨迹，由于质心系也是惯性系，为绘图方便不妨取质心为原点:	
$$
(x_1,y_1)=\frac{m_2R}{m_1+m_2}(\cos wt,\sin wt),\ (x_2,y_2)=-\frac{m_1R}{m_1+m_2}(\cos wt,\sin wt)
$$
代入数据时，假设m1=1，m2=3，R=4，w=1，由程序 generate_data.py 生成了 50 组数据。

程序正确给出质量估计值，接近【1：3】。



### 三体问题的特殊解模拟

查阅到一种三体问题的特殊解，运动方程如下：
$$
r_1=\frac{\sqrt{3}}{4}a\cos(wt)\\
r_2=\frac{\sqrt{7}}{4}a\cos(wt+\arccos(-\sqrt{\frac{3}{7}}))\\
r_3=\frac{\sqrt{7}}{4}a\cos(wt-\arccos(-\sqrt{\frac{3}{7}}))\\
w=2\sqrt{\frac{Gm}{a^3}}
$$
此解为两个质量为 m 的星体和一个质量为 2m 的星体构成等边三角形运动，系安徽省2015年物理高考题。

由于长度系数只与绘图比例有关，我们不妨取w=1，a=4/$\sqrt{3}$，仍然调用 generate_data.py 生成1000组观测数据。

可以看到，质量很接近我们设定的【2:1:1】，说明我们的结果拟合得很好，而且质量比估计也很准确。



## 思考

程序可以改进的地方：

1. 对于标准三体问题以至于更多体问题的没有进行样例测试（这部分可以后面使用微分方程法解决三体问题数值解之后生成数据代入）

2. 样条插值函数对于拟合情况效果很好，但是对于未来轨迹的预测则只能在极小一段时间内生效。可能需要输入更多个周期的数据预测曲线才会变得好看一些。

    

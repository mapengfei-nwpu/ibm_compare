## 1.引言

本文讲述了有限元框架下的IBM及其实现，包括了固体和流体都采用欧拉描述的Level set方法和传统的固体采用拉格朗日描述，流体采用欧拉描述的方法。

早在上世纪70年代，Peskin就提出了IBM，并用于求解流固耦合问题[1,2]。流体区域和固体区域求解的都是流体方程，固体对流体的反馈力通过源项作用给流体方程[3]。这样的模型是把固体看成是没有体积的纤维，每一根纤维都对流体施加反馈力。显然这和实际的物理模型是有出入的，不是我们最终想要的。体积不为零的超弹性物体首次在[6]中解决了。通过将方程转换成弱形式，然后再返回到强形式，推导出了更加通用的IBM。但是这里的要求是固体的粘性和流体的粘性相同。随后，更进一步，[7]提出了IBM的完全变分形式，可以在有限元框架下求解更广泛的问题，比如固体和流体的密度可以不相同，固体的粘性和流体的粘性可以不相同，甚至没有粘性，是纯弹性材料。甚至是固体可压。此外，采用完全变分形式的IBM不需要delta函数插值，并且通过了大量算例[8]的测试。本文主要采用的是[6]中的推导出的IBM。

IBM除了在建模方面有了长足的发展，流体和固体的交互也有了很大的进步。最初的IBM都是通过delta函数插值的方式进行流体和固体之间的信息传递。在数值计算时，delta函数需要通过连续函数近似，我们不仅需要寻找合适的delta函数，还要寻找合适的插值方式。[9]提出了一种利用高斯积分进行插值，使得固体网格的密度可以比流体网格小。也有人在探索不需要delta插值的IBM。[10]推到出了完全欧拉格式的IBM，固体和流体可以公用同一套网格，因此不需要delta函数插值。并且[9]利用了level set的思想追踪固体的区域，弥补了欧拉描述不能跟踪固体点的缺点。[7]利用了变分形式的IBM可以直接建立起流体和固体之间的联系，不需要引入delta函数。但是这种方法在实际计算的时候涉及到了两套不同网格上的有限元空间，无法直接通过现有的有限元框架求解。

[9]和[10]描述的方法都是基于有限差分的，这两种方法目前还有没有限元形式的实现，本文尝试着将它们通过有限元离散求解。

## 1.有限元形式下的IBM

$\Omega$为固定的控制单元，$B_t$是浸没在流体中的固体，$\Omega \backslash B_{t}$为流体区域。固体的运动可以描述为$\zeta: B \rightarrow B_{t}, \boldsymbol{x}=\zeta(\boldsymbol{s}, t)$，其中$B$是固体的参考状态。在控制单元内部，$B_t$的运动和流体的运动均可通过相同的连续方程和动量方程描述：

$$
\frac{\partial \rho}{\partial t}+\nabla \cdot(\rho \boldsymbol{u})=0 \quad
$$
$$
\rho\left[\frac{\partial \boldsymbol{u}}{\partial t}
+(\boldsymbol{u}\cdot \nabla) \boldsymbol{u} \right]= \nabla \cdot \boldsymbol{\sigma}+\rho \boldsymbol{b}
$$

其中$\rho(\boldsymbol{x}, t)$是密度，$\boldsymbol{u}(\boldsymbol{x}, t)$是速度，$\sigma(\boldsymbol{x}, t)$是柯西应力张量，$\boldsymbol{b}(\boldsymbol{x}, t)$是体积力。方程(1)和(2)对流体和固体都成立，区别在于固体和流体的本构方程不同，也就是柯西应力张量的不同。在讨论IBM的时候，不失一般性，我们可以假设固体和流体的密度恒为单位1，且没有外源力，那么原来的方程可以改成(3)和(4)。IBM的重点在于对柯西应力张量的描述。
$$
\nabla \cdot \boldsymbol{u}=0 \quad
$$
$$
\frac{\partial \boldsymbol{u}}{\partial t}
+(\boldsymbol{u}\cdot \nabla) \boldsymbol{u}= \nabla \cdot \boldsymbol{\sigma}
$$

**流体的本构方程** 在不可压的条件下，流体的应力张量可以写成
$$
\sigma_f=-p |+\sigma_{\mathrm{f}}^{v} \quad
$$
$$
\quad \sigma_{\mathrm{f}}^{v}=\mu_{\mathrm{f}}\left(\nabla \boldsymbol{u}+\nabla \boldsymbol{u}^{\mathrm{T}}\right)
$$
其中，$\sigma_f$表示流体区域的柯西应力张量，$\sigma^v_f$为粘性项。

**固体的本构方程** 同样的，我们只考虑不可压情况下的柯西应力张量，相对于流体，固体多出了因形变产生的应力$\sigma_s^e$。
$$
\sigma_s=-p |+\sigma_{\mathrm{s}}^{e}+\sigma_{\mathrm{s}}^{v}
$$

其中$\sigma_{\mathrm{s}}^{e}=J^{-1} \mathrm{P}_{\mathrm{s}}^{e} \mathrm{F}^{\mathrm{T}}$，$\mathrm{P}_{\mathrm{s}}^{e}=\frac{\partial W_{\mathrm{s}}^{e}(\mathrm{F})}{\partial \mathrm{F}}$， $\sigma_{\mathrm{s}}^{v}=\mu_{\mathrm{s}}\left(\nabla \boldsymbol{u}+\nabla \boldsymbol{u}^{\mathrm{T}}\right)$。

流体和固体的应力张量只是在各自内部适用，这里并没有考虑交界面上的应力张量，即假定$\sigma n|_{\partial B}=0$。[4]具体描述了这种方法，[9]描述了如何将固体内部的力和固体表面的力合并在一起。

**变分形式** 确定了柯西应力张量以后，原问题就可以转换成变分形式的问题。给定空间

$$
\begin{array}{c}
\mathbf{V}=H_{0}^{1}(\Omega)^{d}:=\left\{\mathbf{v} \in L^{2}(\Omega)^{d} \text { s.t. } \nabla \mathbf{v} \in L^{2}(\Omega)^{d \times d},\left.\mathbf{v}\right|_{\partial \Omega}=\mathbf{0}\right\} \\
Q=L_{0}^{2}(\Omega):=\left\{q \in L^{2}(\Omega) \text { s.t. } \int_{\Omega} q=0\right\}
\end{array}
$$

原问题可以写成：

寻找$u\in V$和$p\in Q$使得$\forall \boldsymbol v\in V$和$\forall q \in Q$都有
$$
\int_\Omega \frac{\partial \boldsymbol u}{\partial t} \boldsymbol v dx+\int_\Omega \boldsymbol u \cdot \nabla \boldsymbol u \boldsymbol v dx
=\int_\Omega \nabla\cdot\boldsymbol\sigma \boldsymbol vdx \\
\int_\Omega \nabla \cdot \boldsymbol u q dx=0
$$
观察动量方程变分形式的右端，我们可以推导

$$
\int_\Omega \nabla\cdot\sigma vdx=\int_{\Omega\backslash B_t} \nabla \cdot \sigma_f v dx+\int_{B_t} \nabla \cdot \sigma_s v dx\\
=\int_\Omega \nabla \cdot\sigma_f v dx+\int_{B_t}\nabla\cdot(\sigma_s-\sigma_f)vdx
$$
假设固体的粘性和流体的粘性相同，即$\mu_s=\mu_f$，那么固体和流体的应力张量可以消去相同的部分，得到
$$
\int_{B_t}\nabla\cdot(\sigma_s-\sigma_f)vdx=\int_{B_t}\nabla\cdot\sigma_s^e vdx=\int _{B_t}\nabla \cdot (J^{-1} \mathrm{P}_{\mathrm{s}}^{e} \mathrm{F}^{\mathrm{T}})vdx=\int_B\nabla \cdot\mathrm{P}_{\mathrm{s}}^{e} vdX
$$
最后一个等式是由于固体不可压，所以$J=1$，且有$F^{-T}=dX/dx$。对于这一项，参考状态下的积分和实时状态下的积分结果一样。随后，原来的动量方程可以写成
$$
\int_{\Omega}(\frac{\partial \boldsymbol{u}}{\partial t}
+(\boldsymbol{u}\cdot \nabla) \boldsymbol{u})vdx=\int_{\Omega} \nabla \cdot \boldsymbol{\sigma}_fvdx+\int_{B_t}\nabla \cdot ( \mathrm{P}_{\mathrm{s}}^{e} \mathrm{F}^{\mathrm{T}})vdx
$$
$$
\int_{\Omega}(\frac{\partial \boldsymbol{u}}{\partial t}
+(\boldsymbol{u}\cdot \nabla) \boldsymbol{u})vdx=\int_{\Omega} \nabla \cdot \boldsymbol{\sigma}_fvdx+\int_{B}\nabla \cdot \mathrm{P}_{\mathrm{s}}^{e} vdX
$$
这两个式子分别对应了欧拉法和拉格朗日法计算弹性方程。传统的IBM采用的是下面的这个公式，Heltai等人介绍的完全变分格式IBM和Maitre等人提出的level set格式的IBM采用的是上面的这个公式。~~有的论文中提及了固体和流体交界处的应力可能存在跳跃条件，还没看懂。~~

现在IBM就很明了了，在变分形式下，流固耦合方程只是比纯粹的流体方程的求解多出了一项$\int_B\nabla \cdot P_s^e\boldsymbol v dX$或者是$\int_{B_t}\nabla \cdot ( \mathrm{P}_{\mathrm{s}}^{e} \mathrm{F}^{\mathrm{T}})vdx$。如果在有限元框架下求解这个方程是不需要delta插值的，直接积分即可，但是这样实际操作比较困难，测试函数和$P_s^eF^T$所在的有限元空间处于不同的网格中，因此对这一项的积分无法通过现有的实现方式实现，目前只有Heltai等人实现了这种方法。此外，有限差分框架下的IBM有两种方式可以求解这一项，一种是最经典的delta函数插值的方式，固体采用拉格朗日描述，另一种是LS-IBM方法，固体采用欧拉描述。

## 2.固体应力的计算

通常IBM可以分解成以下三部分的求解：

1. 求解流体的速度
2. 求解固体的形变梯度和应力
3. 将应力通过插值等方式写进流体的控制方程

这三个部分可以整体求解，或者完全分开求解，或者2，3一起求解，1单独求解，每一种方法都能找到对应的文章。对于FE-IBM，本文采用的是完全分开求解的方式。对于FE-LS-IBM，本文采用的是2，3一起求解，1单独求解。

##### Delta 函数插值

任意一个定义在B空间上的函数都可以通过delta插值转移到$\Omega$上。
$$
\mathbf{V}(\mathbf{s}, t)=\mathbf{v}(\mathbf{X}(\mathbf{s}, t))=\int_{\Omega} \mathbf{v}(\mathbf{x}) \boldsymbol{\delta}(\mathbf{x}-\mathbf{X}(\mathbf{s}, t)) \mathrm{d} \mathbf{x} \quad \forall \mathbf{s} \in \mathcal{B}
$$
通过这个公式，在B上的积分项可以推导成
$$
\int_B(\nabla\cdot P)\cdot V(X,t)dX\\
=\int_B(\nabla\cdot P)\cdot \int _\Omega v\delta(x-X)dx dX\\
=\int_\Omega\int_B(\nabla\cdot P)\cdot  \delta(x-X)dXvdx
$$
其中，令$G=\nabla\cdot P$，我们可以通过B上的有限元空间计算出$G$。（前面假设了固体边界上没有外法向的力的作用）
$$
\int_BGVdX=\int_B(\nabla\cdot P)VdX=\int_BP:\nabla VdX
$$

##### 输运方程和水平集函数

IBLS方法的核心是推导出固体的形变梯度随时间变化的函数，它借助两个方程。第一个是输运方程，流场中任意一点的拉格朗日坐标是个常量，因此随体导数为零，根据这个公式可以导出追踪欧拉坐标的方程，对坐标梯度求逆就可以推导出流场中任意点的形变梯度。第二个是水平集方程，它用于追踪固体区域，将非固体区域的形变梯度归零。有了这两个公式以后，我们就可以知道任意时刻任意位置的形变梯度，进而求出任意位置的力。

有两种方法计算形变梯度，一种是用输运方程追踪坐标轨迹，然后求逆就可以得到形变梯度。另一种方法是直接求形变梯度的输运方程，由于算例的边界速度不为零，我们采用第二种方法。令F就是我们之前设定的形变梯度，我们可以得到下面的方程[10]

$$
\int_\Omega F_{t}vdx+\int_\Omega\mathrm{u} \cdot \nabla(F)vdx=\int_\Omega\nabla u^{T} F v dx
$$

令$I$为单位矩阵，$O$为零矩阵，在t=0时刻，初始的形变梯度为

$$
{F}_{0}=\left\{\begin{array}{ll}
{O} & \text { in } \Omega \backslash {B}_{t} \\
{I} & \text { in } {B}_{t}
\end{array}\right.
$$

因为形变梯度方程是双曲型的方程，传统有限元无法求解。有限元框架下求解双曲型方程有两种方法，一种是GLS，一种是间断有限元。我们选取了SUPG方法，它属于GLS的一种，时间离散使用了隐格式。

因为形变梯度是欧拉描述的形变梯度，我们可以直接计算应力项
$$
\int_{B_t}\nabla \cdot ( \mathrm{P}_{\mathrm{s}}^{e} \mathrm{F}^{\mathrm{T}})vdx=\int_{\Omega}\nabla \cdot ( \mathrm{P}_{\mathrm{s}}^{e} \mathrm{F}^{\mathrm{T}})vdx
$$


~~但是三维的形变梯度追踪方程有9个变量，每个变量有M个自由度，最差的情况，可能要做(2(9m-1)9m)^9m次乘法与加法，照这样来算，计算量约为流体求解的81倍。从计算效率方面考虑，我不认为可以推广到三维。有限差分之所以可以推广到三维是因为可以显式求解，不用求解庞大的线性方程组。如果采用DG，则可以显式求解，逐个单元求解。~~

### 3.时间离散

Chorin的投影方法和SAV格式

#### 数值模拟

###### 1. 方腔流

需要列出一个速度的表格，或者画图。相同时刻，取x=0.5和y=0.5画两幅图，进行比较。

投影+IBLS，投影+IBFE，SAV+IBLS，SAV+IBFE，



###### 2. 圆盘震动

同上



[1] Charles S. Peskin. Numerical analysis of blood ow in the heart

[2] C.S. Peskin. Flow patterns around heart valves: A numerical method

[3] C.S. Peskin, The immersed boundary method

[4] On the hyper-elastic formulation of the immersed boundary method

[5] A finite element approach for the immersed boundary method

[6] On the hyper-elastic formulation of the immersed boundary method

[7] Variational Implementation of Immersed Finite Element Methods

[8] Benchmarking the immersed finite element method for fluid–structure interaction problems

[9] Hybrid finite difference/finite element immersed boundary method

[10] Eulerian formulation and level set models for incompressible fluid-structure interaction
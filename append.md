## 1.有限元形式下的IBM

$\Omega$为固定的控制单元，$B_t$是浸没在流体中的固体，$\Omega \backslash B_{t}$为流体区域。固体的运动可以描述为$\zeta: B \rightarrow B_{t}, \boldsymbol{x}=\zeta(\boldsymbol{s}, t)$，其中$B$是固体的参考状态。在控制单元内部，$B_t$的运动和流体的运动均可通过相同的连续方程和动量方程描述：

$$
\frac{\partial \rho}{\partial t}+\nabla \cdot(\rho \boldsymbol{u})=0 \quad
$$
$$
\rho\left[\frac{\partial \boldsymbol{u}}{\partial t}
+(\boldsymbol{u}\cdot \nabla) \boldsymbol{u} \right]= \nabla \cdot \boldsymbol{\sigma}+\rho \boldsymbol{b}
$$

其中$\rho(\boldsymbol{x}, t)$是密度，$\boldsymbol{u}(\boldsymbol{x}, t)$是速度，$\sigma(\boldsymbol{x}, t)$是柯西应力张量，$\boldsymbol{b}(\boldsymbol{x}, t)$是体积力。方程(1)和(2)对流体和固体都成立，区别在于固体和流体的本构方程不同，也就是柯西应力张量的不同。在讨论IBM的时候，不失一般性，我们可以假设固体和流体的密度恒为单位1，且没有外源力，那么原来的方程可以改成(3)和(4)。
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

（流体和固体的应力张量只是在各自内部适用，这里并没有考虑交界面上的应力张量。IBFD/FE的文章中讲述了这种方法）

**动量方程的弱形式** 在动量方程左右两端乘上任意函数$v$，再在整个区域上积分，我们就可以得到
$$
\int_\Omega \nabla\cdot\sigma vdx=\int_{\Omega\backslash B_t} \nabla \cdot \sigma_f v dx+\int_{B_t} \nabla \cdot \sigma_s v dx\\
=\int_\Omega \nabla \cdot\sigma_f v dx+\int_{B_t}\nabla\cdot(\sigma_s-\sigma_f)vdx
$$
在固体的粘性和流体的粘性相同的情况下，可以得到
$$
\int_{B_t}\nabla\cdot(\sigma_s-\sigma_f)vdx=\int_{B_t}\nabla\cdot\sigma_s^e vdx=\int _{B_t}\nabla \cdot (J^{-1} \mathrm{P}_{\mathrm{s}}^{e} \mathrm{F}^{\mathrm{T}})vdx=\int_BdX
$$
因为固体不可压，在参考状态下积分和在实时状态下的积分结果一样。

（此处提及固体和流体交界处的应力可能存在跳跃条件，那样的话就不能分部积分了）

现在IBM就很明了了，只是比流体方程的求解多出了一项$\int_B\nabla \cdot P_s^eF^Tdx$，如果在有限元框架下求解这个方程是不需要delta插值的，直接积分即可。但是在实际操作过程中，测试函数和$P_s^eF^T$所在的有限元空间处于不同的网格中，对这一项的积分要么在全域上进行，要么就寻找测试函数非零的区域。传统的IB方法中有两种方法可以解这一项，一是IBLS方法，另一个是delta函数插值。

**再返回强形式** 





主要的区别

IBLS是利用拉格朗日法描述固体的运动，IBFE是用欧拉法描述固体的运动，利用levelset方程捕捉流体和固体之间的界面。



#### Delta 函数插值

传统的IBFE的方法是在固体的拉格朗日坐标系中求出力，然后在根据固体的欧拉坐标和流体的拉格朗日坐标，利用delta函数，将力施加到流体方程上。

这里面很关键的一步是固体是如何获取速度，如何计算出形变梯度和力，如何将力反馈给流体。

固体的应力求解需要在拉格朗日框架下进行，相乘的基函数对应的是拉格朗日意义上的基函数，在这样的坐标系下计算出

#### 输运方程和水平集函数

IBLS方法的核心是推导出固体的形变梯度随时间变化的函数，它借助两个方程。第一个是输运方程，流场中任意一点的拉格朗日坐标是个常量，因此随体导数为零，根据这个公式可以导出追踪欧拉坐标的方程，对坐标梯度求逆就可以推导出流场中任意点的形变梯度。第二个是水平集方程，它用于追踪固体区域，将非固体区域的形变梯度归零。有了这两个公式以后，我们就可以知道任意时刻任意位置的形变梯度，进而求出任意位置的力。

有两种方法计算形变梯度，一种是用输运方程追踪坐标轨迹，然后求逆就可以得到形变梯度。另一种方法是直接求形变梯度的输运方程。

因为形变梯度方程是双曲型的方程，传统有限元无法求解。有限元框架下求解双曲型方程有两种方法，一种是GLS，一种是间断有限元。我们选取了SUPG方法，它属于GLS的一种，时间离散使用了向后欧拉或者CN这样的隐格式。

但是三维的形变梯度追踪方程有9个变量，每个变量有M个自由度，最差的情况，可能要做(2(9m-1)9m)^9m次乘法与加法，照这样来算，计算量约为流体求解的81倍。从计算效率方面考虑，我不认为可以推广到三维。有限差分之所以可以推广到三维是因为可以显式求解，不用求解庞大的线性方程组。如果采用DG，则可以显式求解，逐个单元求解。

### 流体的求解格式

Chorin的投影方法

SAV方法



#### 数值模拟

###### 1. 方腔流

需要列出一个速度的表格，或者画图。相同时刻，取x=0.5和y=0.5画两幅图，进行比较。

SAV+IBLS, SAV+IBFE, 投影+IBLS, 投影+IBFE



###### 2. 圆盘震动

同上
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
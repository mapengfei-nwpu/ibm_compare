from dolfin import *

# 零Dirichlet边界条件
# 区域都是[0,1]x[0,1], T=0.1, delta=0.1

# 1.1
solution_1 = {}
ux = "-exp(t) * x^2 * (x - 1)^2 * y * (y - 1) * (2 * y - 1) / 256" 
uy = "exp(t) * x * (x - 1) * (2 * x - 1) * y^2 * (y - 1)^2/256"
p = "exp(t) * (x^3 - 1/4)" 

# 1.2
ux = "exp(t) * sin(pi * x)^2 * sin(2 * pi * y)"
uy = "-exp(t) * sin(2 * pi * x) * sin(pi * y)^2"
p = "exp(t) * (sin(pi * y) - 2 / pi)"

# 1.3
ux = "20 * x^2 * (x - 1)^2 * y * (y - 1) * (2 * y - 1) * t"
uy = "-20 * x * (x - 1) * (2 * x - 1) * y^2 * (y - 1)^2 * t"
p = "10 * (2 * x - 1) * (2 * y - 1)"

# 1.4
# 压强空间收敛阶不稳定
ux = "2 * pi * sin(pi * x)^2 * sin(pi * y) * cos(pi * y) * cos(t)"
uy = "-2 * pi * sin(pi * x) * cos(pi * x) * sin(pi * y)^2 * cos(t)"
p = "cos(pi * x) * cos(pi * y)"

# 非零Dirichlet边界条件
# 2.1 [0,2]x[-1,1]
ux = "2 * cos(pi * y) * sin(pi * x) * sin(t)"
uy = "-2 * sin(pi * y) * cos(pi * x) * sin(t)" 
p = "2 * sin(pi * y) * sin(pi * x) * cos(t)"

# 2.2 [0,1]x[-0.25,0]
ux = "(x^2 * y^2 + exp(-y)) * cos(2 * pi * t)"
uy = "(-2 * x * y^3/3 + 2 - pi * sin(pi * x)) * cos(2 * pi * t)" 
p = "-(2 - pi * sin(pi * x)) * cos(2 * pi * y) * cos(2 * pi * t)"

# 1.1

fx = "3*x^2*exp(t) + (x^2*exp(t)*(x - 1)^2*(y - 1))/6400 + (x^2*exp(t)*(2*y - 1)*(x - 1)^2)/12800 + (x^2*y*exp(t)*(x - 1)^2)/6400 + (x^2*y*exp(t)*(2*y - 1)*(y - 1))/12800 + (y*exp(t)*(2*y - 1)*(x - 1)^2*(y - 1))/12800 + (x*y*exp(t)*(2*x - 2)*(2*y - 1)*(y - 1))/6400 - (x^2*y*exp(t)*(2*y - 1)*(x - 1)^2*(y - 1))/256 - (x^3*y^2*exp(2*t)*(2*x - 1)*(x - 1)^3*(y - 1)^2*(6*y^2 - 6*y + 1))/65536 + (x^3*y^2*exp(2*t)*(2*y - 1)*(x - 1)^2*(y - 1)*(2*x^2 - 3*x + 1)*(2*y^2 - 3*y + 1))/32768"

fy = "(x*y^2*exp(t)*(2*x - 1)*(x - 1)*(y - 1)^2)/256 - (y^2*exp(t)*(2*x - 1)*(y - 1)^2)/12800 - (x*y^2*exp(t)*(y - 1)^2)/6400 - (x*y^2*exp(t)*(2*x - 1)*(x - 1))/12800 - (x*exp(t)*(2*x - 1)*(x - 1)*(y - 1)^2)/12800 - (x*y*exp(t)*(2*x - 1)*(2*y - 2)*(x - 1))/6400 - (y^2*exp(t)*(x - 1)*(y - 1)^2)/6400 - (x^2*y^3*exp(2*t)*(2*y - 1)*(x - 1)^2*(y - 1)^3*(6*x^2 - 6*x + 1))/65536 + (x^2*y^3*exp(2*t)*(2*x - 1)*(x - 1)*(y - 1)^2*(2*x^2 - 3*x + 1)*(2*y^2 - 3*y + 1))/32768"

# 1.2

fx = "(exp(t)*(50*sin(pi*x)^2*sin(2*pi*y) - pi^2*cos(pi*x)^2*sin(2*pi*y) + 3*pi^2*sin(pi*x)^2*sin(2*pi*y) + 100*pi*exp(t)*cos(pi*x)*sin(pi*x)^3*sin(2*pi*y)^2 - 100*pi*exp(t)*cos(2*pi*y)*sin(pi*x)^2*sin(2*pi*x)*sin(pi*y)^2))/50"
 
fy = "(exp(t)*(50*pi*cos(pi*y) - 50*sin(2*pi*x)*sin(pi*y)^2 + pi^2*cos(pi*y)^2*sin(2*pi*x) - 3*pi^2*sin(2*pi*x)*sin(pi*y)^2 + 100*pi*exp(t)*cos(pi*y)*sin(2*pi*x)^2*sin(pi*y)^3 - 100*pi*exp(t)*cos(2*pi*x)*sin(pi*x)^2*sin(pi*y)^2*sin(2*pi*y)))/50"

# 1.3

fx = "40*y - (2*t*x^2*(2*y - 1)*(x - 1)^2)/5 - (4*t*x^2*y*(x - 1)^2)/5 - (4*t*x^2*(x - 1)^2*(y - 1))/5 - (2*t*x^2*y*(2*y - 1)*(y - 1))/5 - (2*t*y*(2*y - 1)*(x - 1)^2*(y - 1))/5 + 20*x^2*y*(2*y - 1)*(x - 1)^2*(y - 1) - (4*t*x*y*(2*x - 2)*(2*y - 1)*(y - 1))/5 - 400*t^2*x^3*y^2*(2*x - 1)*(x - 1)^3*(y - 1)^2*(6*y^2 - 6*y + 1) + 800*t^2*x^3*y^2*(2*y - 1)*(x - 1)^2*(y - 1)*(2*x^2 - 3*x + 1)*(2*y^2 - 3*y + 1) - 20"
 
fy = "40*x + (2*t*y^2*(2*x - 1)*(y - 1)^2)/5 + (4*t*x*y^2*(y - 1)^2)/5 + (4*t*y^2*(x - 1)*(y - 1)^2)/5 + (2*t*x*y^2*(2*x - 1)*(x - 1))/5 + (2*t*x*(2*x - 1)*(x - 1)*(y - 1)^2)/5 - 20*x*y^2*(2*x - 1)*(x - 1)*(y - 1)^2 + (4*t*x*y*(2*x - 1)*(2*y - 2)*(x - 1))/5 - 400*t^2*x^2*y^3*(2*y - 1)*(x - 1)^2*(y - 1)^3*(6*x^2 - 6*x + 1) + 800*t^2*x^2*y^3*(2*x - 1)*(x - 1)*(y - 1)^2*(2*x^2 - 3*x + 1)*(2*y^2 - 3*y + 1) - 20"

# 1.4 

fx = "(3*pi^3*cos(pi*y)*sin(pi*x)^2*sin(pi*y)*cos(t))/25 - pi*cos(pi*y)*sin(pi*x) - 2*pi*cos(pi*y)*sin(pi*x)^2*sin(pi*y)*sin(t) + 4*pi^3*cos(pi*x)*sin(pi*x)^3*sin(pi*y)^4*cos(t)^2 - (pi^3*cos(pi*x)^2*cos(pi*y)*sin(pi*y)*cos(t))/25 + 4*pi^3*cos(pi*x)*cos(pi*y)^2*sin(pi*x)^3*sin(pi*y)^2*cos(t)^2"

fy = "2*pi*cos(pi*x)*sin(pi*x)*sin(pi*y)^2*sin(t) - (3*pi^3*cos(pi*x)*sin(pi*x)*sin(pi*y)^2*cos(t))/25 - pi*cos(pi*x)*sin(pi*y) + 4*pi^3*cos(pi*y)*sin(pi*x)^4*sin(pi*y)^3*cos(t)^2 + (pi^3*cos(pi*x)*cos(pi*y)^2*sin(pi*x)*cos(t))/25 + 4*pi^3*cos(pi*x)^2*cos(pi*y)*sin(pi*x)^2*sin(pi*y)^3*cos(t)^2"
 
# 2.1

fx = "2*cos(pi*y)*sin(pi*x)*cos(t) + 4*pi*cos(pi*x)*sin(pi*x) - 4*pi*cos(pi*x)*sin(pi*x)*cos(t)^2 + (pi^2*cos(pi*y)*sin(pi*x)*sin(t))/25 + 2*pi*cos(pi*x)*sin(pi*y)*cos(t)"
 
fy = "4*pi*cos(pi*y)*sin(pi*y) - 2*cos(pi*x)*sin(pi*y)*cos(t) - 4*pi*cos(pi*y)*sin(pi*y)*cos(t)^2 - (pi^2*cos(pi*x)*sin(pi*y)*sin(t))/25 + 2*pi*cos(pi*y)*sin(pi*x)*cos(t)"
# 2.2

fx = "cos(2*pi*t)^2*(- 2*y*x^2 + exp(-y))*((2*x*y^3)/3 + pi*sin(x*pi) - 2) - (y^2*cos(2*pi*t))/50 - (cos(2*pi*t)*(2*x^2 + exp(-y)))/100 - 2*pi*sin(2*pi*t)*(exp(-y) + x^2*y^2) + 2*x*y^2*cos(2*pi*t)^2*(exp(-y) + x^2*y^2) + pi^2*cos(2*pi*t)*cos(pi*x)*cos(2*pi*y)"
 
fy = "2*pi*sin(2*pi*t)*((2*x*y^3)/3 + pi*sin(x*pi) - 2) - (pi^3*cos(2*pi*t)*sin(pi*x))/100 - cos(2*pi*t)^2*((2*y^3)/3 + pi^2*cos(x*pi))*(exp(-y) + x^2*y^2) + (x*y*cos(2*pi*t))/25 - 2*pi*cos(2*pi*t)*sin(2*pi*y)*(pi*sin(pi*x) - 2) + 2*x*y^2*cos(2*pi*t)^2*((2*x*y^3)/3 + pi*sin(x*pi) - 2)"
 

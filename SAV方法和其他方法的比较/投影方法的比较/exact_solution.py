# 1.1
solution_1 = {}
ux = "-exp(t) * x[0]*x[0] * (x[0] - 1) * (x[0] - 1) * x[1] * (x[1] - 1) * (2 * x[1] - 1) / 256"
uy = "exp(t) * x[0] * (x[0] - 1) * (2 * x[0] - 1) * x[1]* x[1] * (x[1] - 1)* (x[1] - 1)/256"
p = "exp(t) * (x[0]*x[0]*x[0] - 1/4)"
fx = "3*x[0]*x[0]*exp(t) + (x[0]*x[0]*exp(t)*(x[0] - 1)*(x[0] - 1)*(x[1] - 1))/6400 + (x[0]*x[0]*exp(t)*(2*x[1] - 1)*(x[0] - 1)*(x[0] - 1))/12800 + (x[0]*x[0]*x[1]*exp(t)*(x[0] - 1)*(x[0] - 1))/6400 + (x[0]*x[0]*x[1]*exp(t)*(2*x[1] - 1)*(x[1] - 1))/12800 + (x[1]*exp(t)*(2*x[1] - 1)*(x[0] - 1)*(x[0] - 1)*(x[1] - 1))/12800 + (x[0]*x[1]*exp(t)*(2*x[0] - 2)*(2*x[1] - 1)*(x[1] - 1))/6400 - (x[0]*x[0]*x[1]*exp(t)*(2*x[1] - 1)*(x[0] - 1)*(x[0] - 1)*(x[1] - 1))/256 - (x[0]*x[0]*x[0]*x[1]*x[1]*exp(2*t)*(2*x[0] - 1)*(x[0] - 1)*(x[0] - 1)*(x[0] - 1)*(x[1] - 1)*(x[1] - 1)*(6*x[1]*x[1] - 6*x[1] + 1))/65536 + (x[0]*x[0]*x[0]*x[1]*x[1]*exp(2*t)*(2*x[1] - 1)*(x[0] - 1)*(x[0] - 1)*(x[1] - 1)*(2*x[0]*x[0] - 3*x[0] + 1)*(2*x[1]*x[1] - 3*x[1] + 1))/32768"
fy = "(x[0]*x[1]*x[1]*exp(t)*(2*x[0] - 1)*(x[0] - 1)*(x[1] - 1)*(x[1] - 1))/256 - (x[1]*x[1]*exp(t)*(2*x[0] - 1)*(x[1] - 1)*(x[1] - 1))/12800 - (x[0]*x[1]*x[1]*exp(t)*(x[1] - 1)*(x[1] - 1))/6400 - (x[0]*x[1]*x[1]*exp(t)*(2*x[0] - 1)*(x[0] - 1))/12800 - (x[0]*exp(t)*(2*x[0] - 1)*(x[0] - 1)*(x[1] - 1)*(x[1] - 1))/12800 - (x[0]*x[1]*exp(t)*(2*x[0] - 1)*(2*x[1] - 2)*(x[0] - 1))/6400 - (x[1]*x[1]*exp(t)*(x[0] - 1)*(x[1] - 1)*(x[1] - 1))/6400 - (x[0]*x[0]*x[1]*x[1]*x[1]*exp(2*t)*(2*x[1] - 1)*(x[0] - 1)*(x[0] - 1)*(x[1] - 1)*(x[1] - 1)*(x[1] - 1)*(6*x[0]*x[0] - 6*x[0] + 1))/65536 + (x[0]*x[0]*x[1]*x[1]*x[1]*exp(2*t)*(2*x[0] - 1)*(x[0] - 1)*(x[1] - 1)*(x[1] - 1)*(2*x[0]*x[0] - 3*x[0] + 1)*(2*x[1]*x[1] - 3*x[1] + 1))/32768"

solution_1["ux"] = ux
solution_1["uy"] = uy
solution_1["p"]  = p
solution_1["fx"] = fx
solution_1["fy"] = fy
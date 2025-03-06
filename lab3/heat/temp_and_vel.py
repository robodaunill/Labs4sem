from fenics import *
import matplotlib.pyplot as plt
import numpy as np

# Параметры задачи
rho = 1.0
mu = 0.001     
k = 0.1
c = 1.0
T_final = 2.0
num_steps = 50
dt = T_final / num_steps

# Создаем шестиугольную сетку
mesh = UnitDiscMesh.create(MPI.comm_world, 6, 1, 2)
coords = mesh.coordinates()
coords *= 2.0  # Масштабируем до радиуса 2

# Определяем функциональные пространства
VE = VectorElement("P", mesh.ufl_cell(), 2)
QE = FiniteElement("P", mesh.ufl_cell(), 1)
TE = FiniteElement("P", mesh.ufl_cell(), 1)

W_element = MixedElement([VE, QE])
W = FunctionSpace(mesh, W_element)
V_T = FunctionSpace(mesh, TE)

# Граничные условия для скорости (скольжение по стенкам)
def boundary(x, on_boundary):
    return on_boundary

bc_u = DirichletBC(W.sub(0), Constant((0, 0)), boundary)

# Два источника тепла (гауссовы пики)
class HeatSources(UserExpression):
    def eval(self, value, x):
        # Первый источник в (-0.5, 0)
        d1 = (x[0]+0.5)**2 + x[1]**2
        # Второй источник в (0.5, 0)
        d2 = (x[0]-0.5)**2 + x[1]**2
        value[0] = 3.0*exp(-50*d1) + 2.0*exp(-50*d2)
        
    def value_shape(self):
        return ()

u_D = HeatSources(degree=2)
bc_T = DirichletBC(V_T, u_D, boundary)

# Начальные условия для скорости (вихрь)
class VortexInitialCondition(UserExpression):
    def eval(self, value, x):
        r = np.sqrt(x[0]**2 + x[1]**2)
        theta = np.arctan2(x[1], x[0])  
        value[0] = -5*x[1]*np.exp(-10*r**2)
        value[1] = 5*x[0]*np.exp(-10*r**2)
        
    def value_shape(self):
        return (2,)


# Инициализация функций
w_n = Function(W)
u_init = VortexInitialCondition(degree=2)
assign(w_n.sub(0), interpolate(u_init, W.sub(0).collapse()))

T_n = interpolate(Constant(0.0), V_T)

# Вариационные формы
(u, p) = TrialFunctions(W)
(v, q) = TestFunctions(W)
T = TrialFunction(V_T)
w_T = TestFunction(V_T)

# Навье-Стокс с нелинейной конвекцией
F_u = (rho*dot((u - w_n.sub(0))/dt, v)*dx + 
       rho*dot(dot(w_n.sub(0), nabla_grad(u)), v)*dx +
       mu*inner(grad(u), grad(v))*dx -
       p*div(v)*dx -
       div(u)*q*dx)
a_u, L_u = lhs(F_u), rhs(F_u)

# Уравнение теплопереноса
F_T = (c*(T - T_n)/dt*w_T*dx + 
       c*dot(w_n.sub(0), grad(T_n))*w_T*dx + 
       k*dot(grad(T), grad(w_T))*dx)
a_T, L_T = lhs(F_T), rhs(F_T)

# Решатели
w = Function(W)
T = Function(V_T)

# Файлы для сохранения
file_u = File("velocity.pvd")
file_p = File("pressure.pvd")
file_T = File("temperature.pvd")

for n in range(num_steps):
    t = n*dt
    
    # Шаг Навье-Стокса
    solve(a_u == L_u, w, bcs=bc_u)
    
    # Шаг теплопереноса
    solve(a_T == L_T, T, bcs=bc_T)
    
    # Обновление решений
    w_n.assign(w)
    T_n.assign(T)
    
    # Сохранение результатов
    u, p = w.split()
    file_u << (u, t)
    file_p << (p, t)
    file_T << (T, t)

print("Симуляция завершена!")
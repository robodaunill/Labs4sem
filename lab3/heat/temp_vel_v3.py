from fenics import *
import numpy as np

# Параметры задачи
rho = 1.0
mu = 0.005
k = 0.1
c = 1.0
T_final = 3.0
num_steps = 100
dt = T_final / num_steps

# Создаем шестиугольную сетку
mesh = UnitDiscMesh.create(MPI.comm_world, 6, 1, 2)
mesh.coordinates()[:] *= 1.5  # Масштабируем сетку

# Функциональные пространства
VE = VectorElement("P", mesh.ufl_cell(), 2)
QE = FiniteElement("P", mesh.ufl_cell(), 1)
TE = FiniteElement("P", mesh.ufl_cell(), 1)

W_element = MixedElement([VE, QE])
W = FunctionSpace(mesh, W_element)
V_T = FunctionSpace(mesh, TE)

# Граничные условия для скорости (неподвижные стенки)
bc_u = DirichletBC(W.sub(0), Constant((0, 0)), "on_boundary")

# Граничные условия для температуры (добавлено охлаждение на границах)
bc_T = DirichletBC(V_T, Constant(0.0), "on_boundary")  # Температура 0 на границах

# Начальные условия
class InitialPressure(UserExpression):
    def eval(self, value, x):
        if x[0] > 0.5:
            value[0] = 5.0
        elif x[0] < -0.5:
            value[0] = -2.0
        else:
            value[0] = 1.0*np.sin(4*x[0])
    
    def value_shape(self):
        return () 
            
class InitialTemperature(UserExpression):
    def eval(self, value, x):
        r = np.sqrt(x[0]**2 + x[1]**2)
        value[0] = 10.0*np.exp(-50*r**2)
    
    def value_shape(self):
        return ()


# Инициализация функций
w_n = Function(W)
u_init = Constant((0, 0))
p_init = InitialPressure(degree=2)
T_init = InitialTemperature(degree=2)

# Присваиваем начальные значения
assign(w_n.sub(0), interpolate(u_init, W.sub(0).collapse()))
assign(w_n.sub(1), interpolate(p_init, W.sub(1).collapse()))
T_n = interpolate(T_init, V_T)

# Вариационные формы
(u, p) = TrialFunctions(W)
(v, q) = TestFunctions(W)
T = TrialFunction(V_T)
w_T = TestFunction(V_T)

# Уравнение Навье-Стокса с учетом градиента давления
F_u = (rho*dot((u - w_n.sub(0))/dt, v)*dx +
       rho*dot(dot(w_n.sub(0), nabla_grad(u)), v)*dx +
       mu*inner(grad(u), grad(v))*dx -
       p*div(v)*dx -
       div(u)*q*dx +
       0.01*dot(grad(p), grad(q))*dx) 

a_u, L_u = lhs(F_u), rhs(F_u)

# Уравнение теплопереноса с источником
F_T = (c*(T - T_n)/dt*w_T*dx +
       c*dot(w_n.sub(0), grad(T_n))*w_T*dx +
       k*dot(grad(T), grad(w_T))*dx)

a_T, L_T = lhs(F_T), rhs(F_T)

# Решатели
w = Function(W)
T = Function(V_T)

# Настройки решателя
solver_params = {
    "linear_solver": "gmres",
    "preconditioner": "ilu",
    "krylov_solver": {"relative_tolerance": 1e-6}
}

# Файлы для сохранения
file_u = File("velocity.pvd")
file_p = File("pressure.pvd")
file_T = File("temperature.pvd")

for n in range(num_steps):
    t = n*dt
    
    # Шаг Навье-Стокса
    solve(a_u == L_u, w, bcs=bc_u, solver_parameters=solver_params)
    
    # Шаг теплопереноса
    solve(a_T == L_T, T, bcs=bc_T)  
    
    # Обновление решений
    w_n.assign(w)
    T_n.assign(T)
    
    # Сохранение каждые 5 шагов
    if n % 5 == 0:
        u, p = w.split()
        file_u << (u, t)
        file_p << (p, t)
        file_T << (T, t)

print("Симуляция завершена!")

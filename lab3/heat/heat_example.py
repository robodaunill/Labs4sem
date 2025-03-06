from __future__ import print_function
from fenics import *
import numpy as np
import matplotlib.pyplot as plt

# Параметры задачи
T = 2.0            # Конечный момент времени
num_steps = 10     # Количество шагов по времени
dt = T / num_steps # Шаг по времени
alpha = 1          # Коэффициент теплопроводности

# Создаем сетку (единичный квадрат)
mesh = UnitSquareMesh(16, 16)

# Определяем пространство конечных элементов (P1)
V = FunctionSpace(mesh, 'P', 1)

# Определяем граничные условия
u_D = Expression('1 + x[0]*x[0] + 2*x[1]*x[1]', degree=2) #  Заданная температура на границе
bc = DirichletBC(V, u_D, 'on_boundary') # Условие Дирихле (температура фиксирована на границе)

# Определяем тестовые и пробные функции
u = TrialFunction(V)
v = TestFunction(V)

# Начальное условие (начальная температура)
u_n = interpolate(u_D, V)

# Определяем вариационную форму (схема Crank-Nicolson)
F = (u - u_n) / dt * v * dx + alpha * dot(grad(u), grad(v)) * dx
a, L = lhs(F), rhs(F)

# Решение по времени
u = Function(V)
t = 0
for n in range(num_steps):
    # Обновляем текущий момент времени
    t += dt

    # Решаем уравнение
    solve(a == L, u, bc)

    # Обновляем начальное условие для следующего шага по времени
    u_n.assign(u)


    file = File("heat_solution.pvd")
    file << (u, t)

print("Решение сохранено в файл heat_solution.pvd")
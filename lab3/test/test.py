from dolfin import *

# Создаем сетку (единичный квадрат)
mesh = UnitSquareMesh(8, 8)  # Или другая сетка

# Определяем пространство конечных элементов (P1)
V = FunctionSpace(mesh, 'P', 1)

# Определяем тестовую функцию
u = TrialFunction(V)
v = TestFunction(V)

# Определяем билинейную форму (Laplace оператор)
a = inner(grad(u), grad(v)) * dx

# Определяем линейную форму (функция правой части)
f = Expression("10*exp(-(pow(x[0] - 0.5, 2) + pow(x[1] - 0.5, 2)) / 0.02)", degree=2)  # Гауссиан
L = f * v * dx

# Определяем граничные условия (условия Дирихле - u = 0 на границе)
bc = DirichletBC(V, Constant(0.0), 'on_boundary')

# Решаем задачу
u = Function(V)
solve(a == L, u, bc)

# Сохраняем решение в файл .pvd
file = File("/home/robodaniil/Labs4sem/lab3/solution.pvd")  
file << u

print("Решение сохранено в файл solution.pvd")
'''
g(t) - вынуждающая сила
dq - класс уравнения
Применение к другим уравнениям:
1.Привести уравнение 2-ого порядка к системе уравнений 1-ого порядка
2.Записать в __init__ и F_vec значения self.F = правой части во втором уравнении
3.В строке 67 написать свои начальные условия
4.Подобрать шаг в строке 26 и количество итераций в строке 66
'''

import cmath
import numpy as np
import matplotlib.pyplot as plt


def g(t):
    A = 1
    w = 100
    return A * cmath.sin(w * t)


class dq:
    def __init__(self, y_0, y_dot):
        self.y_0 = y_0    # начальное положение маятника
        self.y_dot = y_dot  # начальная скорость
        self.h = 0.01     # шаг метода
        self.t = 0       # начальное время, которое увеличивается в дальнейшем
        self.w_0_squared = 1    # квадрат частоты колебаний
        self.gamma = 0.01 # коэффициент трения
        self.y = np.array([y_0, y_dot])  # вектор-массив значений функции с начальными значениями
        self.F = np.array([self.y[1],  # вектор-массив значений правой части уравнения с нач.значениями (см.юпитер-файл)
                           (g(self.t) - 2*self.gamma*self.y[1] - self.w_0_squared*self.y[0])])

    def F_vec(self, t, y):
        # значения правой части (см.юпитер-файл) в произвольный момент времени
        self.F = np.array([y[1], (g(t) - 2*self.gamma*y[1] - self.w_0_squared*y[0])])
        return self.F

    def k1(self, t):  # 1-ый вектор метода Рунге-Кутты
        k1 = self.F_vec(t, self.y)
        return k1

    def k2(self, k1, t):  # 2-ой вектор метода Рунге-Кутты
        k2 = self.F_vec(t + self.h/2, self.y + self.h/2*k1)
        return k2

    def k3(self, k2, t):  # 3-ий вектор метода Рунге-Кутты
        k3 = self.F_vec(t + self.h/2, self.y + self.h/2*k2)
        return k3

    def k4(self, k3, t):  # 4-ый вектор метода Рунге-Кутты
        k4 = self.F_vec(t + self.h, self.y + self.h*k3)
        return k4

    def solution(self, t):  # аналитическое решение в произвольный момент времени
        A = 1
        w = 100
        s = cmath.sqrt(self.gamma*self.gamma-self.w_0_squared)
        c = self.w_0_squared*self.w_0_squared-w*w
        c1 = (self.y_0*(s-self.gamma)-self.y_dot)/(2*s)
        c2 = (self.y_dot+self.y_0*(s+self.gamma))/(2*s)
        return c1*cmath.exp(-t*(s+self.gamma)) + c2*cmath.exp(t*(s-self.gamma)) + \
               A*(cmath.sin(w*t)*c + cmath.cos(w*t)*2*self.gamma*w)/(c*c + 4*self.gamma*self.gamma*w*w)


epoch = 4000  # количество итераций
Y = dq(5, -1)  # инстанс уравнения с начальными условиями
func = np.array([0 for _ in range(epoch)])
func_speed = np.array([0 for _ in range(epoch)])
coords = np.array([i*Y.h for i in range(epoch)]) # координатная сетка
sol = np.array([Y.solution(coords[i]) for i in range(epoch)])
error = np.array([0 for _ in range(epoch)])
func[0] = Y.y[0]
for i in range(1, epoch):  # цикл Рунге-Кутты, см. https://ru.wikipedia.org/wiki/Метод_Рунге_—_Кутты
    K1 = Y.k1(coords[i-1])
    K2 = Y.k2(K1, coords[i-1])
    K3 = Y.k3(K2, coords[i-1])
    K4 = Y.k4(K3, coords[i-1])
    Y.y = Y.y + (K1 + 2*K2 + 2*K3 + K4)*Y.h/6
    func[i] = Y.y[0]
    func_speed[i] = Y.y[1]

# Вывод графиков
plt.title("График зависимости координаты маятника от времени")
plt.plot(coords, func, 'r-', label="Численный расчёт")
plt.plot(coords, sol, 'b-', label="Аналитический расчёт")
plt.axis([0, 40, -5, 7])
plt.grid()
plt.text(0.2, 6.5, r"$\gamma = 0.01, w_0 = 1, y_0=5, y'_0=-11$ ", fontsize=11)
plt.text(0.2, 5.5, r"$A = 1, w = 100, h=0.01$ ", fontsize=11)
plt.xlabel("Время, с")
plt.ylabel("Координата, м")
plt.legend()
#plt.savefig('graph1.png')
plt.show()

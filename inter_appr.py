import numpy as np
from sympy import Symbol
import math


def check_distance(coord):
    """
         Проверка равноотстояние точек
        :param coord: Список с координатми
    """
    x_coord = [x[0] for x in coord]
    d = coord[1][0] - coord[0][0]
    s = x_coord[0]
    for x in x_coord:
        if s == x:
            s += d
            continue
        else:
            return False
    return True


# Две вспомогательные функции для дроби
def fraction_u(s=[], x=0, i=0, t=1):
    res = 1
    if t == 1:
        for j in range(len(s)):
            if i != j:
                res *= (x - s[j])
        return res
    elif t == 2:
        for j in range(0, i + 1):
            res *= (x - s[j])
        return res


def fraction_l(s=[], i=0):
    res = 1
    for j in range(len(s)):
        if i != j:
            res *= (s[i] - s[j])
    return res


#########################################

# Интерполяция

def lagrange_interpolation(coord, delta=0):
    """
         Интерполяция методом Лагранжа (для равноотстоящих и неравноотстоящих узлов)
         :param delta: Добавляет построение точек интерполяционным методом -delta к минимальному
         и + delta к маскимальному значению от исходных координат x
        :param coord: Список с координатми
    """
    l_res = dict()
    x_coord = [x[0] for x in coord]
    y_coord = [y[1] for y in coord]
    sort_coord = sorted(x_coord)

    # Рассчёт матриц для нахождения коэффициентов A*X = B
    s1 = []  # A
    s2 = []  # B

    if check_distance(coord):
        k_x = 0
        # Цикл по всем точкам, которые мы хотим найти
        for x in range(int(sort_coord[0]) - delta, int(sort_coord[-1]) + 1 + delta):
            result = 0
            for i in range(len(coord)):
                if k_x == 0:
                    s1.append([x_coord[i] ** j for j in range(len(x_coord))])
                    s2.append([y_coord[i]])
                result += (y_coord[i] * (fraction_u(s=x_coord, x=x, i=i))) / fraction_l(s=x_coord, i=i)
            l_res[x] = result
            k_x += 1

    else:
        k_x = 0
        # Цикл по всем точкам, которые мы хотим найти
        for x in range(int(sort_coord[0]) - delta, int(sort_coord[-1]) + 1 + delta):
            result = 0
            for i in range(len(coord)):
                if k_x == 0:
                    s1.append([x_coord[i] ** j for j in range(len(x_coord))])
                    s2.append([y_coord[i]])
                result += (y_coord[i] * (fraction_u(s=x_coord, x=x, i=i))) / fraction_l(s=x_coord, i=i)
            l_res[x] = result
            k_x += 1

    a = np.array(s1)
    b = np.array(s2)
    cof1 = np.linalg.solve(a, b)
    cof = [i for i in cof1]
    # Массив  точек в формате [x,y]
    ans = [[i[0], i[1]] for i in l_res.items()]
    print(f'Массив  точек в формате [xi,fi] (fi – значение интерполированной функции в точках): {ans}')
    x = Symbol('x')
    expr = 0
    for i in range(len(coord)):
        expr += (x ** i) * cof[i]
    print(expr)
    return ans


def newton_interpolation(coord, delta=0):
    """
         Интерполяция методом Ньютона
        :param delta: Добавляет построение точек интерполяционным методом -delta к минимальному
         и + delta к маскимальному значению от исходных координат x
        :param coord: Список с координатми
    """
    if check_distance(coord):

        # Координаты изначальных точек
        x_coord = [x[0] for x in coord]
        y_coord = [y[1] for y in coord]

        # Некоторые константы
        sort_coord = sorted(x_coord)
        l = len(y_coord)
        y0 = y_coord[0]
        h = x_coord[1] - x_coord[0]

        # Словарь палинома, где ключ = xi, а значение = fi
        l_res = dict()
        # Список конечных разностей это res_dy
        res_dy = []
        # Коэффициенты многочлена Ньютона
        cof = [y0]

        # Цикл по всем точкам, которые мы хотим найти
        for x in range(int(sort_coord[0]) - delta, int(sort_coord[-1]) + 1 + delta):
            result = y0
            new_d = y_coord
            # Цикл по кол-ву изначальных точек
            for i in range(l):
                delta_y = new_d
                new_d = []
                # Нахождение конечных разностей
                if (len(delta_y) - 1) > 0:
                    for j in range(len(delta_y) - 1):
                        dy = delta_y[j + 1] - delta_y[j]
                        new_d.append(dy)
                    # print(new_d, ' ',len(delta_y)-1)
                    qw = new_d[0]
                    res_dy.append(qw)
                    result += (((fraction_u(s=x_coord, x=x, i=i, t=2)) * new_d[0]) / (
                                math.factorial(i + 1) * (h ** (i + 1))))
                    cof.append(((new_d[0]) / (math.factorial(i + 1) * (h ** (i + 1)))))
            l_res[x] = result

        if delta == 0:
            # Массив  точек в формате [xi, yi, fi]
            ans = []
            for i in range(len(y_coord)):
                ans.append([x_coord[i], y_coord[i], l_res[x_coord[i]]])
            print(f'Массив  точек в формате [xi, yi, fi] (fi – значение интерполированной функции в точках): {ans}')
        else:
            # Массив  точек в формате [x,y]
            ans = [[i[0], i[1]] for i in l_res.items()]
            print(f'Массив  точек в формате [xi,fi] (fi – значение интерполированной функции в точках): {ans}')
        # Вывод многочлена Ньютона
        x = Symbol('x')
        expr = cof[0]
        for i in range(len(coord)):
            expr += fraction_u(s=x_coord, x=x, i=i, t=2) * cof[i + 1]
        print(expr)

        return ans

    else:
        return print('Равноотсояние точек не соблюденно, невозможно применить метод!')


# Для теста:

# check_distance(coord = [[0,-1],[1,-3],[2,3],[6,1187]])
# lagrange_interpolation(coord = [[0,-1],[1,-3],[2,3],[6,1187]])
newton_interpolation(coord=[[0, -1], [1, -3], [2, 3], [3, 1187]], delta=0)
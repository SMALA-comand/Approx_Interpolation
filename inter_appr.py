# from input_csv import check_distance
import numpy as np
from sympy import Symbol


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
def fraction_u(s=[], x=0, i=0):
    res = 1
    for j in range(len(s)):
        if i != j:
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

def lagrange_interpolation(coord):
    """
         Интерполяция методом Лагранжа (для равноотстоящих и неравноотстоящих узлов)
        :param coord: Список с координатми
    """
    l_res = dict()
    x_coord = [x[0] for x in coord]
    y_coord = [y[1] for y in coord]
    x_coord = sorted(x_coord)

    # Рассчёт матриц для нахождения коэффициентов A*X = B
    s1 = []  # A
    s2 = []  # B

    if check_distance(coord):
        print('В процессе...')
    else:
        k_x = 0
        for x in range(int(x_coord[0]), int(x_coord[-1]) + 1):
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
    return

# lagrange_interpolation(coord = [[0,-1],[1,-3],[2,3],[6,1187]])
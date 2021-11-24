import numpy as np
import matplotlib.pyplot as plt
from sympy import Symbol
import math
from scipy.signal import argrelextrema
from linear_function import linear_function
from quadratic_function import quadratic_function
from gauss_function import gauss_function
from csv_reader import get_csv_coord
# Нужно прописать установку pip install scipy


def check_distance(coord):
    """
         Проверка равноотстояние точек
        :param coord: Список с координатми
    """
    x_coord = [round(x[0],1) for x in coord]
    #print(x_coord)
    d = x_coord[1] - x_coord[0]
    s = x_coord[0]
    #print(s)
    for x in x_coord:
        if s == x:
            s += d
            s = round(s,1)
            continue
        else:
            #print(False)
            return False
    #print(True)
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
    #print(f'Массив  точек в формате [xi,fi] (fi – значение интерполированной функции в точках): {ans}')
    x = Symbol('x')
    expr = 0
    ans2 = []
    for i in range(len(coord)):
        expr += (x ** i) * cof[i]
    for el in x_coord:
        for i in range(len(cof)):
            ans2.append([el, ((el ** i)*cof[i])])

    #print(expr)
    out_ans = [expr, ans2]
    return out_ans


def newton_interpolation(coord, delta=0):
    """
         Интерполяция методом Ньютона
        :param delta: Добавляет построение точек интерполяционным методом -delta к минимальному
         и + delta к маскимальному значению от исходных координат x
        :param coord: Список с координатми
    """
    if check_distance(coord):

        # Координаты изначальных точек
        x_coord = [round(x[0],1) for x in coord]
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

        out_ans = [[j[0] for j in l_res.items()], [i[1] for i in l_res.items()]]
        if delta == 0:
            # Массив  точек в формате [xi, yi, fi]
            ans = []
            for i in range(len(y_coord)):
                ans.append([x_coord[i], y_coord[i]])
            #print(f'Массив  точек в формате [xi, yi, fi] (fi – значение интерполированной функции в точках): {ans}')
        else:
            # Массив  точек в формате [x,y]
            ans = [[i[0], i[1]] for i in l_res.items()]
            #print(f'Массив  точек в формате [xi,fi] (fi – значение интерполированной функции в точках): {ans}')
        # Вывод многочлена Ньютона
        x = Symbol('x')
        expr = cof[0]
        for i in range(len(coord)):
            expr += fraction_u(s=x_coord, x=x, i=i, t=2) * cof[i + 1]
        #print(expr)
        out_ans = [expr, ans]
        return out_ans

    else:
        return 'Равноотсояние точек не соблюденно, невозможно применить метод!'


def cubic_spline_interpolation(coord, x0=[]):
    """
         Интерполяция методом Кубического сплайна
        :param x0: Точки по x, для которых мы хотим посчитать интерполяцию  (список)
        :param coord: Список с координатми
    """
    # Координаты изначальных точек
    x_coord = [x[0] for x in coord]
    y_coord = [y[1] for y in coord]

    # Некоторые константы
    sort_coord = sorted(x_coord)
    l = len(y_coord)

    # Словарь полинома, где ключ = xi, а значение = fi
    l_res = dict()

    x = np.asfarray(x_coord)
    y = np.asfarray(y_coord)

    # Сортируем данные
    if np.any(np.diff(x) < 0):
        indexes = np.argsort(x)
        x = x[indexes]
        y = y[indexes]

    size = len(x)

    xdiff = np.diff(x)
    ydiff = np.diff(y)

    # Работа с матрицами
    Li = np.empty(size)
    Li_1 = np.empty(size - 1)
    z = np.empty(size)

    # Заполняем диагонали Li и Li-1 и решаем матричное ур-ние L*Y = B
    Li[0] = (2 * xdiff[0]) ** (0.5)
    Li_1[0] = 0.0
    B0 = 0.0
    z[0] = B0 / Li[0]

    for i in range(1, size - 1, 1):
        Li_1[i] = xdiff[i - 1] / Li[i - 1]
        Li[i] = (2 * (xdiff[i - 1] + xdiff[i]) ** (0.5) - Li_1[i - 1] * Li_1[i - 1])
        Bi = 6 * (ydiff[i] / xdiff[i] - ydiff[i - 1] / xdiff[i - 1])
        z[i] = (Bi - Li_1[i - 1] * z[i - 1]) / Li[i]

    i = size - 1
    Li_1[i - 1] = xdiff[-1] / Li[i - 1]
    Li[i] = (2 * xdiff[-1] - Li_1[i - 1] * Li_1[i - 1]) ** (0.5)
    Bi = 0.0
    z[i] = (Bi - Li_1[i - 1] * z[i - 1]) / Li[i]

    # Решаем L(T - транспонированное) * X = Y
    i = size - 1
    z[i] = z[i] / Li[i]
    for i in range(size - 2, -1, -1):
        z[i] = (z[i] - Li_1[i - 1] * z[i + 1]) / Li[i]

    # Поиск индекса
    index = x.searchsorted(x0)
    np.clip(index, 1, size - 1, index)

    xi1, xi0 = x[index], x[index - 1]
    yi1, yi0 = y[index], y[index - 1]
    zi1, zi0 = z[index], z[index - 1]
    hi1 = xi1 - xi0

    ans = zi0 / (6 * hi1) * (xi1 - x0) ** 3 + \
          zi1 / (6 * hi1) * (x0 - xi0) ** 3 + \
          (yi1 / hi1 - zi1 * hi1 / 6) * (x0 - xi0) + \
          (yi0 / hi1 - zi0 * hi1 / 6) * (xi1 - x0)

    return ans


def splitting_appr(coord):
    """
        Алгоритм аппроксимации по подвыборкам
        :param coord: Список с координатми
    """
    # Координаты изначальных точек
    x_coord = [x[0] for x in coord]
    y_coord = np.array([y[1] for y in coord])

    # 1) Находим в выборке локальные экстремумы:
    idx_minimas = argrelextrema(y_coord, np.less)[0]  # Индексы локальных минимумов
    idx_maximas = argrelextrema(y_coord, np.greater)[0]  # Индексы локальных максимумов

    #print('min', idx_minimas)
    #print('max', idx_maximas)
    idx = np.sort(np.concatenate((idx_minimas, idx_maximas))) # Всё вместе и отсортированно

    # 2) Заносим в отдельный массив координаты точек, являющихся соседними противоположными парами:
    # Получееные индексы координат для "сшива"
    dividing_index = []
    prev_min_idx = idx_minimas[0]
    prev_max_idx = idx_maximas[0]

    for extreme in idx:
        if extreme in idx_minimas:
            prev_min_idx = extreme
        elif extreme in idx_maximas:
            prev_max_idx = extreme
        # 3) и 4)
        if abs(prev_min_idx - prev_max_idx) >= 6:
            value = (prev_min_idx + prev_max_idx) // 2
            dividing_index.append(value)

    print(dividing_index)
    def compare(new_coord):
        """
              Алгоритм сравнения и нахождения оптимального варианта для аппроксимации
              :param new_coord: Список с координатми
        """

        y_lin0 = linear_function(new_coord)[1]
        y_quad0 = quadratic_function(new_coord)[1]
        y_gauss0 = gauss_function(new_coord)[1]
        y_lin = [j[2] for j in y_lin0]
        y_quad = [j[3] for j in y_quad0]
        y_gauss = [j[4] for j in y_gauss0]
        list_y = [y_quad, y_gauss]

        y_coord_div = [j[1] for j in new_coord]  # Координаты по y исходной функции
        # Считаем дисперсию и выбираем лучший вариант для аппроксимации
        # sum(y - yi)^2   min ?
        best_sigma = sum([(y_coord_div[j] - y_lin[j]) ** 2 for j in range(len(y_coord_div))])
        best_y = y_lin
        for yi in list_y:
            t = sum([(y_coord_div[j] - yi[j]) ** 2 for j in range(len(y_coord_div))])
            if t <= best_sigma:
                best_sigma = t
                best_y = yi

        y_res1 = best_y

        return y_res1

    index = 0  # Индеск предыдущего начала разделения списка
    y_res = []  # Итоговый результат после кусочной аппроксимации
    len_d_index = len(dividing_index)
    count = 0  # Счётчик
    if len_d_index == 0:
        new_coord = coord
        new_res = compare(new_coord = new_coord)
        y_res += new_res
    else:
        if count != len_d_index:
            for i in dividing_index:
                new_coord = coord[index:i]

                y_res += compare(new_coord=new_coord)
                index = i
                count += 1

            new_coord = coord[index:]
            y_res += compare(new_coord=new_coord)
    #diff = len(y_res) - len(dividing_index)
    #print(diff)
    return y_res

# Для теста:

# check_distance(coord = [[0,-1],[1,-3],[2,3],[6,1187]])
# lagrange_interpolation(coord = [[0,-1],[1,-3],[2,3],[6,1187]])
# newton_interpolation(coord = [[0,-1],[1,-3],[2,3],[3,1187]], delta = 0)[1]
if __name__ == '__main__':
    coord = get_csv_coord("C:\\Users\\Вячеслав\\PycharmProjects\\Approx_Interpolation\\Доллар_США.csv")
    date = []
    x_st = []
    y_st = []
    for i in range(len(coord)):
        date.append(coord[i][0])
        y_st.append(coord[i][1])
        coord[i][0] = i
        x_st.append(coord[i][0])

#    print(coord)
    res = splitting_appr(coord)
    plt.plot(x_st, res, '-', label = 'Апроксимация')
    plt.plot(x_st, y_st, label = 'Исходная выборка')
    plt.legend(loc='best')
    plt.show()


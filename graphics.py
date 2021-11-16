import matplotlib.pyplot as plt
import numpy as np
from csv_reader import get_csv_coord
from linear_function import linear_function
from quadratic_function import quadratic_function
from inter_appr import *
from gauss_function import gauss_function

def get_graphics(path = "cvss.csv"):

    plan = get_csv_coord(path)
    x = [i[0] for i in plan]
    y = [i[1] for i in plan]


    # считаем аппроксимацию линейной функцией
    linear_plan = linear_function(coord=plan)[1]
    y_linear = [i[2] for i in linear_plan]


    # считаем аппроксимацию квадратичной функцией
    quadratic_plan = quadratic_function(coord=plan)[1]
    y_quadratic = [i[3] for i in quadratic_plan]

    # считаем аппроксимацию функцией нормального распределения
    gauss_f, gauss_plan = gauss_function(coords = plan)
    y_normal = [i[4] for i in gauss_plan]
    #y_normal = [1 for i in plan]

    # считаем интерполяцию Лагранжем
    lagrange_plan = lagrange_interpolation(coord=plan)[1]
    y_lagr = [i[1] for i in linear_plan]

    # считаем интерполяцию Ньютоном

    newton_plan = newton_interpolation(coord=plan)
    if newton_plan == 'Равноотсояние точек не соблюденно, невозможно применить метод!':
        newton_f = "Невозможно построить"
        y_newton = [1] * len(x)
    else:
        newton_f = 'Interpolated f(x)'
        newton_plan = newton_plan[1]
        y_newton = [i for i in newton_plan]

    # считаем интерполяцию сплайнами
    y_cube = cubic_spline_interpolation(coord=plan, x0 = np.linspace(x[0], x[-1], 100))

    # считаем аппроксимацию numpy
    c = np.polyfit(x, y, 2)
    p = np.poly1d(c)
    y_num_apr = [p(i) for i in x]
    label_numpy = f'{round(c[0], 3)}x² + {round(c[1], 3)}x + {round(c[2], 3)}'

    # считаем интерполяцию numpy
    y_num_inter = np.interp(x, x, y)




    ##############################################################


    y_plan = [y_linear, y_quadratic, y_normal, y_lagr, y_newton, y_cube, y_num_apr, y_num_inter]

    titles = ['линейная аппроксимация', 'квадратичная аппроксимация', 'аппроксимация норм. распред.',
              'Интерполяция Лагранжем', 'Интерполяция Ньютона', 'Интерполяция куб. сплайном',
              'Аппроксимация NumPy', 'Интерполяция NumPy']

    labels_func = [linear_function(coord=plan)[0], quadratic_function(coord=plan)[0], gauss_f, 'Interpolated f(x)',newton_f , 'Interpolated f(x)', label_numpy, 'Interpolated f(x)']

    fig, axes = plt.subplots(3, 3)

    i = 0
    for ax in axes.ravel():
        if i == 8:
            break
        ax.plot(x, y, linestyle='-.', label='точки')
        ax.plot(x, y_plan[i], '-x', label=labels_func[i])
        ax.legend(loc='best')
        ax.set_title(titles[i])
        i += 1

    fig.set_figheight(20)
    fig.set_figwidth(20)
    plt.show()


if __name__ == '__main__':
    #path = input('Введите пусть:')
    path = "cvss.csv"
    path = "C:\\Users\Вячеслав\Downloads\subject4_Test2.csv"
    coord = get_graphics(path = path)




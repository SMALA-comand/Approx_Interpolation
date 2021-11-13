import matplotlib.pyplot as plt
import numpy as np
from csv_reader import get_csv_coord
from linear_function import linear_function
from quadratic_function import quadratic_function


plan = get_csv_coord("cvss.csv")
x = [i[0] for i in plan]
y = [i[1] for i in plan]

# считаем аппроксимацию линейной функцией
linear_plan = linear_function(coord=plan)[1]
y_linear = [i[2] for i in linear_plan]

# считаем аппроксимацию квадратичной функцией
quadratic_plan = quadratic_function(coord=plan)[1]
y_quadratic = [i[2] for i in quadratic_plan]

# считаем аппроксимацию функцией нормального распределения
y_normal = [1]*50

# считаем интерполяцию Лагранжем
y_lagr = [1]*50

# считаем интерполяцию Ньютоном
y_newton = [1]*50

# считаем аинтерполяцию сплайнами
y_cube = [1]*50

# считаем аппроксимацию numpy
c = np.polyfit(x, y, 2)
p = np.poly1d(c)
y_num_apr = [p(i) for i in x]

# считаем интерполяцию numpy
y_num_inter = [1]*50


##############################################################


y_plan = [y_linear, y_quadratic, y_normal, y_lagr, y_newton, y_cube, y_num_apr, y_num_inter]

titles = ['линейная аппроксимация', 'квадратичная аппроксимация', 'аппроксимация норм. распред.',
          'Интерполяция Лагранжем', 'Интерполяция Ньютона', 'Интерполяция куб. сплайном',
          'Аппроксимация NumPy', 'Интерполяция NumPy']

labels_func = [linear_function(coord=plan)[0], quadratic_function(coord=plan)[0], '', '', '', '', '', '']

fig, axes = plt.subplots(3, 3)

i = 0
for ax in axes.ravel():
    if i == 8:
        break
    ax.plot(x, y, linestyle='-.', label='точки')
    ax.plot(x, y_plan[i], label=labels_func[i])
    ax.legend(loc='upper right')
    ax.set_title(titles[i])
    i += 1

fig.set_figheight(20)
fig.set_figwidth(20)
plt.show()



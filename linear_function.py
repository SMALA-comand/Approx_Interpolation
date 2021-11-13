import numpy as np
import sympy
from csv_reader import get_csv_coord

def linear_function(coord):
    """

    :param coord: координаты точек в форме массива [[x1, y1], [x2, y2], ...]
    """
    ### для первого уравнения
    for_c1_1 = sum([i[0]**2 for i in coord])
    for_c0_1 = sum([i[0] for i in coord])
    ans_1 = sum([i[0]*i[1] for i in coord])

    ### для второго уравнения
    for_c1_2 = for_c0_1
    for_c0_2 = len(coord)
    ans_2 = sum([i[1] for i in coord])

    a = np.array([[for_c1_1, for_c0_1], [for_c1_2, for_c0_2]])
    b = np.array([ans_1, ans_2])
    c = np.linalg.solve(a, b)             # [c1, c0]
    c1 = c[0]
    c0 = c[1]

    answer_b = f'y = {round(c1, 3)}x'
    if c0 < 0:
        answer_b += f'{round(c0, 3)}'
    else:
        answer_b += f'+{round(c0, 3)}'

    answer_a = []

    for i in coord:
        el = i
        fi = c0 + c1*el[0]
        el.append(fi)
        answer_a.append(el)

    return answer_b, answer_a


if __name__ == '__main__':
    coord = get_csv_coord("cvss.csv")
    func, coords = linear_function(coord=coord)
    print(func)
    print(*coords)


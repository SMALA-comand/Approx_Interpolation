import numpy as np
from csv_reader import get_csv_coord


def quadratic_function(coord, is_gauss = 0):
    """
        :param coord: координаты точек в форме массива [[x1, y1], [x2, y2], ...]
        :return: возвращает кортеж (функция в виде строки, [[x1, y1, f1], ...], дисперсия)
    """
    ### для первого уравнения
    for_c2_1 = sum([i[0]**4 for i in coord])
    for_c1_1 = sum([i[0]**3 for i in coord])
    for_c0_1 = sum([i[0]**2 for i in coord])
    ans_1 = sum([(i[0]**2)*i[1] for i in coord])

    ### для второго уравнения
    for_c2_2 = for_c1_1
    for_c1_2 = for_c0_1
    for_c0_2 = sum([i[0] for i in coord])
    ans_2 = sum([i[0]*i[1] for i in coord])

    ### для третьего уравнения
    for_c2_3 = for_c0_1
    for_c1_3 = for_c0_2
    for_c0_3 = len(coord)
    ans_3 = sum([i[1] for i in coord])

    a = np.array([[for_c2_1, for_c1_1, for_c0_1], [for_c2_2, for_c1_2, for_c0_2], [for_c2_3, for_c1_3, for_c0_3]])
    b = np.array([ans_1, ans_2, ans_3])
    c = np.linalg.solve(a, b)             # [c2, c1, c0]
    c2 = c[0]
    c1 = c[1]
    c0 = c[2]

    if is_gauss == 1:
        return(c0, c1, c2)

    answer_b = f'y = {round(c2, 3)}x²'
    if c1 < 0:
        answer_b += f'{round(c1, 3)}x'
    else:
        answer_b += f'+{round(c1, 3)}x'
    if c0 < 0:
        answer_b += f'{round(c0, 3)}'
    else:
        answer_b += f'+{round(c0, 3)}'

    answer_a = []

    for i in coord:
        el = i
        #print(el)
        fi = c0 + c1*el[0] + c2*el[0]**2
        el.append(fi)
        answer_a.append(el)
    #print(answer_a)
    return answer_b, answer_a


if __name__ == '__main__':
    coord = get_csv_coord("cvss.csv")
    func, coords = quadratic_function(coord=coord)
    print(func)
    print(*coords)


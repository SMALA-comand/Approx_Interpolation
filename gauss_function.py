from csv_reader import get_csv_coord
from quadratic_function import quadratic_function
from math import log, exp

def gauss_function(coords):
    for i in range(0, len(coords)):
        if coords[i][1] <= 0:
            coords[i][1] = 0
        else:
            try:
                coords[i][1] = log(coords[i][1])
            except ValueError:
                print(coords[i][1])

    c0, c1, c2 = quadratic_function(coord=coords, is_gauss=1)
    a = exp(c0 - c1 ** 2 / (4 * c2))
    b = (-1) / c2
    c = -c1 / (2 * c2)

    answer_b = f'y = {round(a, 5)} * e^('
    if b > 0:
        answer_b += '-'
    if c > 0:
        answer_b += f'(x - {round(c, 5)}) ^ 2'
    else:
        answer_b += f'(x + {round(-c, 5)}) ^ 2'
    if b < 0:
        answer_b += f' / ({round(-b, 5)}))'
    else:
        answer_b += f' / {round(b, 5)})'

    coord = []

    for i in coords:
        el = i
        fi = a * exp(- (el[0] - c) ** 2 / b)
        el.append(fi)
        el[1] = exp(el[1])
        coord.append(el)

    func = answer_b
    print(coord)
    return(func, coord)

if __name__ == '__main__':
    coords = get_csv_coord("Курс_валюты.csv")
    func, coord = gauss_function(coords)
    print(func)
    print(*coord)
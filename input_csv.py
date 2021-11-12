import csv


def input_csv(csv_name="test_points.csv"):
    """
         Данная функция принимает откуда-то путь к файлу csv с точками для интерполяции и аппроксимации
        :param csv_name: Путь к файлу csv, который нужно ввести пользователю (либо мы заранее должны его знать)
    """
    coord = []

    with open(csv_name, mode='r', encoding='utf-8') as f:
        data = csv.reader(f, delimiter=',')
        for row in data:
            r = []
            for el in row:
                r.append(float(el))
            coord.append(r)
    return coord


# Проверка  на загаловки почему-то работает неправильно, не знаю в чём причина ...
# def checking_headers(csv_name):
#  """
#      Проверка на заголовки
#     :param csv_name: Путь к файлу csv, который нужно ввести пользователю (либо мы заранее должны его знать)
# """
# with open(csv_name, mode='r', encoding='utf-8') as f:
#   sample = f.read(5)
#  header = csv.Sniffer().has_header(sample)
# if header:
#    return True
# else:
#   return False

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

# check_distance([[0,-1],[1,-3],[2,3],[6,1187]])
# input_csv()
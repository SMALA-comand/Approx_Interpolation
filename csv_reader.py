import csv


def get_csv_coord(csv_path):
    """
    Получает на вход csv файл с координатами, разделитель - запятая, десятичная часть числа отделяется точкой
    :param csv_path: путь до csv файла
    :return: массив вида [[x1, y1], [x2, y2], ...]
    """
    coordinates = []
    with open(csv_path, 'r') as f:
        reader = csv.reader(f)
        count = 0
        for line in reader:
            x = line[0]
            y = line[1]
            try:
                x = float(x)
                y = float(y)
            except ValueError:
                if count == 0:
                    count += 1
                    continue
                else:
                    if len(x) == 10 and '.20' in x:
                        count += 1
                        coordinates.append([x, float(y)])
                    else:
                        return 'Лишние символы в координатах!'
            else:
                count += 1
                coordinates.append([x, y])
    return coordinates


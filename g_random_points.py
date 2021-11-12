import random as ra
import csv


def point_generator(k=0):
    """
         Данная функция создаёт csv файл с точками для интерполяции и аппроксимации
        :param k: кол-во точек, которое мы будем генерировать
    """
    if k == 0:
        print('Введите количество точек (k): ')
        k = int(input())

    p = []
    for i in range(k):
        p.append([ra.random() * (10 ** ra.randint(1, 6)), ra.random() * (10 ** ra.randint(1, 6))])
    with open("rand_points.csv", mode='w', encoding='utf-8') as w_file:
        file_writer = csv.writer(w_file, delimiter=',', lineterminator='\r')
        file_writer.writerows(p)
    return p

# print(point_generator(5))


def csv_4_test(p = [[0,-1],[1,-3],[2,3],[6,1187]]):
    with open("test_points.csv", mode='w', encoding='utf-8') as w_file:
        file_writer = csv.writer(w_file, delimiter=',', lineterminator='\r')
        file_writer.writerows(p)
    return 'CSV файл для теста создан успешно!'
# csv_4_test(p = [[0,-1],[1,-3],[2,3],[6,1187]])
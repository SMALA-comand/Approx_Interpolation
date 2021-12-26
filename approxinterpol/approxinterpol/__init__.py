"""
``approxinterpol``
===================

Библиотека, реализующая алгоритмы аппроксимации и интерполяции.
Представлены следующие алгоритмы аппроксимации: линейной и квадратичной функциями + аппрокс. по подвыборкам.
И следующие алгоритмы интерполяции: методом Ньютона, методом Лагранжа и кубическими сплайнами.

"""
from .approxinterpol import (quadratic_function, linear_function, splitting_appr,
                             newton_interpolation, lagrange_interpolation, cubic_spline_interpolation)

__author__ = 'Марк Козлов, Вячеслав Есаков, Артём Радайкин, Александр Савостьянов, Лев Памбухчян'

__version__ = "0.0.3"

__all__ = ['quadratic_function', 'linear_function', 'splitting_appr',
           'newton_interpolation', 'lagrange_interpolation', 'cubic_spline_interpolation']

# -*- coding: utf8 -*-

import math
import time
import multiprocessing
import csv
# import numpy


# TODO: найти более эффективный алгоритм
def razr(x):
    a = 0
    while True:
        if math.modf(x)[0] == 0:
            break
        x *= 10
        a += 1
    return a


# TODO: найти более эффективный алгоритм
def myrange(start, stop, step):
    a = max(razr(start), razr(stop), razr(step))
    start = int(start * math.pow(10, a))
    stop = int(stop * math.pow(10, a))
    step = int(step * math.pow(10, a))
    for i in range(start, stop, step):
        yield float(i) / math.pow(10, a)


# оптимальный по Митрохину параметр m 
def m_opt(psi, mu, fi, betta2, alfa1):
    # промежуточные переменные
    a = math.pow(math.cos(math.radians(betta2)), 2) * (math.pow(psi, 2))
    b = math.pow(math.cos(math.radians(alfa1)), 2) * (math.pow(fi, 2))
    c = math.pow(mu, 2) * (1 - math.pow(fi, 2)) * (1 - a) / (b + math.pow(mu, 2) * (1 - math.pow(fi, 2)))
    return (1 / math.pow(psi, 2)) * (1 - math.pow(c, 0.5))


# оптимальная по Митрохину относительная окружная скорость на наружном диаметре рабочего колеса
def u1_otn_opt(psi, mu, m, betta2, fi, alfa1):
    # промежуточные переменные
    a = math.pow(math.cos(math.radians(betta2)), 2) / math.pow(m, 2) - math.pow(psi, 2)
    b = math.pow(math.cos(math.radians(alfa1)), 2) * math.pow(fi, 2) * (1 - math.pow(m, 2) * math.pow(psi, 4))
    c = math.pow(m, 2) * math.pow(psi, 2) * (1 - math.pow(fi, 2))
    return psi / math.pow(math.pow(mu, 2) * a + b / c, 0.5)


# оптимальная по Митрохину степень реактивности
def ro_opt(m, psi, fi, u1_otn, alfa1):
    # промежуточные переменные
    a = (1 - m * math.pow(psi, 2)) * fi * u1_otn * math.cos(math.radians(alfa1))
    b = m * math.pow(psi, 2) * (1 - math.pow(fi, 2))
    return 1 - math.pow(a / b, 2)


# приведенная скорость
def Lambda(c, k, R, T):
    # промежуточная переменная
    a = math.pow(2 * k * R * T / (k + 1), 0.5)
    return c / a


# газодинамическая функция лямбда от пи
def Lambda_Pi(Pi, k):
    # промежуточная переменная
    a = math.pow(Pi, (k - 1) / k)
    return math.pow((1 - a) * (k + 1) / (k - 1), 0.5)


# газодинамическая функция тау от лямбда
def Tau_Lambda(k, Lambda):
    return 1 - math.pow(Lambda, 2) * (k - 1) / (k + 1)


# газодинамическая функция пи от лямбда
def Pi_Lambda(k, Lambda):
    return math.pow(Tau_Lambda(k, Lambda), k / (k - 1))


# критическая скорость
def a_kr(k, R, T):
    return math.pow(R * T * 2 * k / (k + 1), 0.5)


# скорость звука
def a(k, R, T):
    return math.pow(k * R * T, 0.5)


# угол относительной скорости на входе в рабочее колесо
def get_betta1(alfa1, u1_otn, ro, fi):
    # промежуточная переменная
    a = u1_otn / (math.pow(1 - ro, 0.5) * fi * math.sin(math.radians(alfa1)))
    return math.degrees(math.atan(1 / (1 / math.tan(math.radians(alfa1)) - a)))


# газодинамическая функция М от лямбда
def M_Lambda(k, Lambda):
    # промежуточные переменные
    a = 2 * math.pow(Lambda, 2) / (k + 1)
    b = 1 - (k - 1) * math.pow(Lambda, 2) / (k + 1)
    return math.pow(a / b, 0.5)


# зависимость приведенной скорости от числа Маха
def Lambda_M(k, M):
    # промежуточные переменные
    a = (k + 1) * math.pow(M, 2) / 2
    b = 1 + (k - 1) * math.pow(M, 2) / 2
    return math.pow(a / b, 0.5)


# угол абсолютной скорости за турбиной
def get_alfa2(betta2, u2, w2):
    # промежуточная переменная
    a = u2 / (w2 * math.sin(math.radians(betta2)))
    return math.degrees(math.atan(1 / (1 / math.tan(math.radians(betta2)) - a)))


# расход газа
def get_G(N, c_ad, etta_u):
    return 2 * N * 1000 / (math.pow(c_ad, 2) * etta_u)


# мощность турбины
def get_N(G, c_ad, etta_u):
    return G * math.pow(c_ad, 2) * etta_u / (2 * 1000)


# полное давление на выходе
def get_P2_(P2, c2, R, T2):
    # полное давление на выходе из рабочего колеса
    P2_ = P2 * (1 + math.pow(c2, 2) / (2 * R * T2))

    return P2_


# параметры на входе в рабочее колесо
def inlet_wheel_param(P0_, P2, k, R, T0_, u1_otn, ro, fi, alfa1, n):
    # критическая скорость по параметрам на входе в турбину
    a_kr0 = a_kr(k, R, T0_)

    # отношение давлений
    Pi_ad = P2 / P0_

    # адиабадическая приведенная скорость
    Lambda_ad = Lambda_Pi(Pi_ad, k)

    # адиабатическая скорость истечения
    c_ad = a_kr0 * Lambda_ad

    # теоретическая скорость истечения из соплового аппарата
    c1_t = math.pow((1 - ro), 0.5) * c_ad

    # действительная скорость на выходе из соплового аппарата
    c1 = fi * c1_t

    # окружная скорость на наружном диаметре рабочего колеса
    u1 = u1_otn * c_ad

    # теоретическая приведенная скорость истечения из соплового аппарата
    Lambda_c1_t = math.pow((1 - ro), 0.5) * Lambda_ad

    # действительная приведенная скорость на выходе из соплового аппарата
    Lambda_c1 = fi * Lambda_c1_t

    Tau_c1 = Tau_Lambda(k, Lambda_c1)

    # статическая температура за сопловым аппаратом
    T1 = T0_ * Tau_c1

    Pi_c1_t = Pi_Lambda(k, Lambda_c1_t)

    # статическое давление за сопловым аппаратом
    P1 = P0_ * Pi_Lambda(k, Lambda_c1_t)

    # угол относительной скорости на входе в рабочее колесо
    betta1 = get_betta1(alfa1, u1_otn, ro, fi)

    # число Маха на входе в рабочее колесо
    M_w1 = M_Lambda(k, Lambda_c1) * math.sin(math.radians(alfa1)) / math.sin(math.radians(betta1))

    # приведенная скорость на входе в рабочее колесо
    Lambda_w1 = Lambda_M(k, M_w1)

    Tau_w1 = Tau_Lambda(k, Lambda_w1)

    Pi_w1 = Pi_Lambda(k, Lambda_w1)

    # скорость звука на входе в рабочее колесо, м/с
    a1 = a(k, R, T1)

    # относительная скорость на входе в рабочее колесо м/с
    w1 = M_w1 * a1

    # температура торможения на входе в рабочее колесо в относительном движении
    T_w1_ = T1 / Tau_w1

    # давление торможения на входе в рабочее колесо в относительном движении
    P_w1_ = P1 / Pi_w1

    # приведенная окружная скорость на входе в рабочее колесо
    Lambda_u1 = Lambda(u1, k, R, T0_)

    # средний диаметр на входе в рабочее колесо
    D1 = 60 * u1 / (math.pi * n)

    result = {'a_kr0': a_kr0, 'c_ad': c_ad, 'u1': u1, 'c1_t': c1_t, 'c1': c1,
              'Lambda_c1': Lambda_c1, 'T1': T1, 'P1': P1, 'betta1': betta1,
              'a1': a1, 'w1': w1, 'T_w1_': T_w1_, 'P_w1_': P_w1_, 'D1': D1,
              'Lambda_u1': Lambda_u1, 'Tau_w1': Tau_w1, 'Tau_c1': Tau_c1,
              'Pi_w1': Pi_w1, 'Pi_c1_t': Pi_c1_t}

    return result


# параметры на выходе из рабочего колеса
def outlet_wheel_param(u1, mu, P2, P_w1_, k, R, T_w1_, psi, Lambda_u1, Tau_w1, Tau_c1, betta2, D1):
    # окружная скорость на выходе из турбины
    u2 = u1 * mu

    # отношение давлений в рабочем колесе
    Pi_rk = P2 / P_w1_

    Lambda_rk = Lambda_Pi(Pi_rk, k)

    Tau_rk = Tau_Lambda(k, Lambda_rk)

    # адиабатический теплоперепад, срабатываемый в рабочем колесе
    h_rk = k * R * T_w1_ * (1 - Tau_rk) / (k - 1)

    # теоретическая скорость истечения из рабочего колеса
    w2_t = math.pow(2 * h_rk - math.pow(u1, 2) * (1 - math.pow(mu, 2)), 0.5)

    # действительная скорость на выходе из рабочего колеса
    w2 = psi * w2_t

    # температура торможения на выходе из рабочего колеса
    T_w2_ = (1 - (k - 1) * math.pow(Lambda_u1, 2) * (1 - math.pow(mu, 2)) * Tau_w1 / (Tau_c1 * (k + 1))) * T_w1_

    # критическая скорость на выходе из рабочего колеса в относительном движении
    a_kr_w2 = a_kr(k, R, T_w2_)

    # действительная приведенная скорость на выходе из рабочего колеса в относительном движении
    Lambda_w2 = w2 / a_kr_w2

    Tau_w2 = Tau_Lambda(k, Lambda_w2)

    M_w2 = M_Lambda(k, Lambda_w2)

    # статическая температура за рабочим колесом
    T2 = T_w2_ * Tau_w2

    # скорость звука за рабочим колесом
    a2 = a(k, R, T2)

    # угол абсолютной скорости за турбиной
    alfa2 = get_alfa2(betta2, u2, w2)

    # число Маха по абсолютной скорости на выходе
    M_c2 = M_w2 * math.sin(math.radians(betta2)) / math.sin(math.radians(alfa2))

    # абсолютная скорость на выходе из турбины
    c2 = M_c2 * a2

    # средний диаметр на выходе из рабочего колеса
    D2 = D1 * mu

    P2_ = get_P2_(P2, c2, R, T2)

    result = {'u2': u2, 'h_rk': h_rk, 'w2_t': w2_t, 'w2': w2, 'a_kr_w2': a_kr_w2,
              'Lambda_w2': Lambda_w2, 'alfa2': alfa2, 'T2': T2, 'T_w2_': T_w2_,
              'c2': c2, 'M_c2': M_c2, 'P2_': P2_, 'P2': P2, 'D2': D2, 'Pi_rk': Pi_rk}

    return result


# статическое давление вычисляется перебором
def get_P2(P0_, P2_, k, R, T0_, u1_otn, ro, fi, alfa1, mu, psi, betta2, n):
    P2 = P2_
    while True:
        iwt = inlet_wheel_param(P0_=P0_, P2=P2, k=k, R=R, T0_=T0_, u1_otn=u1_otn, ro=ro, fi=fi, alfa1=alfa1, n=n)
        owt = outlet_wheel_param(u1=iwt['u1'], mu=mu, P2=P2, P_w1_=iwt['P_w1_'], k=k, R=R, T_w1_=iwt['T_w1_'],
                                 psi=psi, Lambda_u1=iwt['Lambda_u1'], Tau_w1=iwt['Tau_w1'],
                                 Tau_c1=iwt['Tau_c1'], betta2=betta2, D1=iwt['D1'])
        _P2_ = get_P2_(P2, owt['c2'], R, owt['T2'])
        if round(_P2_, 4) == round(P2_, 4):
            break
        else:
            P2 = P2 - (_P2_ - P2_)
    return P2


def solution(k, R, T0_, n, alfa1, betta2, u1_otn, mu, ro, fi, psi, G=None, N=None, P0_=None, P2=None, P2_=None):
    global result, result
    if P0_ is not None:
        if P2_ is not None and G is not None:
            P2 = get_P2(P0_=P0_, P2_=P2_, k=k, R=R, T0_=T0_, u1_otn=u1_otn, ro=ro, fi=fi, alfa1=alfa1, mu=mu, psi=psi,
                        betta2=betta2, n=n)
            iwt = inlet_wheel_param(P0_=P0_, P2=P2, k=k, R=R, T0_=T0_, u1_otn=u1_otn, ro=ro, fi=fi, alfa1=alfa1, n=n)
            owt = outlet_wheel_param(u1=iwt['u1'], mu=mu, P2=P2, P_w1_=iwt['P_w1_'], k=k, R=R, T_w1_=iwt['T_w1_'],
                                     psi=psi, Lambda_u1=iwt['Lambda_u1'], Tau_w1=iwt['Tau_w1'],
                                     Tau_c1=iwt['Tau_c1'], betta2=betta2, D1=iwt['D1'])

            # работа ступени на окружности колеса
            h_u = iwt['c1'] * iwt['u1'] * math.cos(math.radians(alfa1)) + owt['c2'] * owt['u2'] * math.cos(
                math.radians(owt['alfa2']))

            # кпд на окружности колеса
            etta_u = 2 * h_u / math.pow(iwt['c_ad'], 2)

            # высота лопаток на выходе из соплового аппарата
            l1 = G * R * iwt['T1'] / (
                math.pi * iwt['D1'] * iwt['c1'] * math.sin(math.radians(alfa1)) * iwt['P1'] * math.pow(10, 3))

            # высота лопаток на выходе из рабочего колеса
            l2 = G * R * owt['T2'] / (
                math.pi * owt['D2'] * owt['w2'] * math.sin(math.radians(betta2)) * P2 * math.pow(10, 3))

            # радиус перефирии на выходе из рабочего колеса
            R2_p = (owt['D2'] + l2) / 2

            # радиус втулки на выходе из рабочего колеса
            R2_v = (owt['D2'] - l2) / 2

            N = get_N(G, iwt['c_ad'], etta_u)

            result = {'G': G, 'P0_': P0_, 'P2_': P2_, 'alfa1': alfa1, 'betta2': betta2,
                      'u1_otn': u1_otn, 'mu': mu, 'fi': fi, 'psi': psi, 'ro': ro,
                      'a_kr0': iwt['a_kr0'], 'c_ad': iwt['c_ad'], 'u1': iwt['u1'], 'c1_t': iwt['c1_t'],
                      'c1': iwt['c1'], 'Lambda_c1': iwt['Lambda_c1'], 'T1': iwt['T1'], 'P1': iwt['P1'],
                      'betta1': iwt['betta1'], 'a1': iwt['a1'], 'w1': iwt['w1'], 'T_w1_': iwt['T_w1_'],
                      'P_w1_': iwt['P_w1_'], 'D1': iwt['D1'],
                      'u2': owt['u2'], 'h_rk': owt['h_rk'], 'w2_t': owt['w2_t'], 'w2': owt['w2'],
                      'a_kr_w2': owt['a_kr_w2'], 'Lambda_w2': owt['Lambda_w2'], 'alfa2': owt['alfa2'],
                      'T2': owt['T2'], 'T_w2_': owt['T_w2_'], 'c2': owt['c2'], 'M_c2': owt['M_c2'],
                      'D2': owt['D2'],
                      'h_u': h_u, 'etta_u': etta_u, 'l1': l1, 'l2': l2, 'R2_p': R2_p, 'R2_v': R2_v,
                      'P2': P2, 'N': N}

        elif P2_ is not None and N is not None:
            P2 = get_P2(P0_=P0_, P2_=P2_, k=k, R=R, T0_=T0_, u1_otn=u1_otn, ro=ro, fi=fi, alfa1=alfa1, mu=mu, psi=psi,
                        betta2=betta2, n=n)
            iwt = inlet_wheel_param(P0_=P0_, P2=P2, k=k, R=R, T0_=T0_, u1_otn=u1_otn, ro=ro, fi=fi, alfa1=alfa1, n=n)
            owt = outlet_wheel_param(u1=iwt['u1'], mu=mu, P2=P2, P_w1_=iwt['P_w1_'], k=k, R=R, T_w1_=iwt['T_w1_'],
                                     psi=psi, Lambda_u1=iwt['Lambda_u1'], Tau_w1=iwt['Tau_w1'],
                                     Tau_c1=iwt['Tau_c1'], betta2=betta2, D1=iwt['D1'])

            # работа ступени на окружности колеса
            h_u = iwt['c1'] * iwt['u1'] * math.cos(math.radians(alfa1)) + owt['c2'] * owt['u2'] * math.cos(
                math.radians(owt['alfa2']))

            # кпд на окружности колеса
            etta_u = 2 * h_u / math.pow(iwt['c_ad'], 2)

            G = get_G(N, iwt['c_ad'], etta_u)

            # высота лопаток на выходе из соплового аппарата
            l1 = G * R * iwt['T1'] / (
                math.pi * iwt['D1'] * iwt['c1'] * math.sin(math.radians(alfa1)) * iwt['P1'] * math.pow(10, 3))

            # высота лопаток на выходе из рабочего колеса
            l2 = G * R * owt['T2'] / (
                math.pi * owt['D2'] * owt['w2'] * math.sin(math.radians(betta2)) * P2 * math.pow(10, 3))

            # радиус перефирии на выходе из рабочего колеса
            R2_p = (owt['D2'] + l2) / 2

            # радиус втулки на выходе из рабочего колеса
            R2_v = (owt['D2'] - l2) / 2

            result = {'N': N, 'P0_': P0_, 'P2_': P2_, 'P2': P2, 'alfa1': alfa1, 'betta2': betta2,
                      'u1_otn': u1_otn, 'mu': mu, 'fi': fi, 'psi': psi, 'ro': ro,
                      'a_kr0': iwt['a_kr0'], 'c_ad': iwt['c_ad'], 'u1': iwt['u1'], 'c1_t': iwt['c1_t'],
                      'c1': iwt['c1'], 'Lambda_c1': iwt['Lambda_c1'], 'T1': iwt['T1'], 'P1': iwt['P1'],
                      'betta1': iwt['betta1'], 'a1': iwt['a1'], 'w1': iwt['w1'], 'T_w1_': iwt['T_w1_'],
                      'P_w1_': iwt['P_w1_'], 'D1': iwt['D1'],
                      'u2': owt['u2'], 'h_rk': owt['h_rk'], 'w2_t': owt['w2_t'], 'w2': owt['w2'],
                      'a_kr_w2': owt['a_kr_w2'], 'Lambda_w2': owt['Lambda_w2'], 'alfa2': owt['alfa2'],
                      'T2': owt['T2'], 'T_w2_': owt['T_w2_'], 'c2': owt['c2'], 'M_c2': owt['M_c2'],
                      'D2': owt['D2'],
                      'h_u': h_u, 'etta_u': etta_u, 'l1': l1, 'l2': l2, 'R2_p': R2_p, 'R2_v': R2_v, 'G': G}

        elif P2 is not None and G is not None:
            iwt = inlet_wheel_param(P0_=P0_, P2=P2, k=k, R=R, T0_=T0_, u1_otn=u1_otn, ro=ro, fi=fi, alfa1=alfa1, n=n)
            owt = outlet_wheel_param(u1=iwt['u1'], mu=mu, P2=P2, P_w1_=iwt['P_w1_'], k=k, R=R, T_w1_=iwt['T_w1_'],
                                     psi=psi, Lambda_u1=iwt['Lambda_u1'], Tau_w1=iwt['Tau_w1'],
                                     Tau_c1=iwt['Tau_c1'], betta2=betta2, D1=iwt['D1'])

            # работа ступени на окружности колеса
            h_u = iwt['c1'] * iwt['u1'] * math.cos(math.radians(alfa1)) + owt['c2'] * owt['u2'] * math.cos(
                math.radians(owt['alfa2']))

            # кпд на окружности колеса
            etta_u = 2 * h_u / math.pow(iwt['c_ad'], 2)

            # высота лопаток на выходе из соплового аппарата
            l1 = G * R * iwt['T1'] / (
                math.pi * iwt['D1'] * iwt['c1'] * math.sin(math.radians(alfa1)) * iwt['P1'] * math.pow(10, 3))

            # высота лопаток на выходе из рабочего колеса
            l2 = G * R * owt['T2'] / (
                math.pi * owt['D2'] * owt['w2'] * math.sin(math.radians(betta2)) * P2 * math.pow(10, 3))

            # радиус перефирии на выходе из рабочего колеса
            R2_p = (owt['D2'] + l2) / 2

            # радиус втулки на выходе из рабочего колеса
            R2_v = (owt['D2'] - l2) / 2

            N = get_N(G, iwt['c_ad'], etta_u)

            result = {'G': G, 'P0_': P0_, 'P2_': owt['P2_'], 'alfa1': alfa1, 'betta2': betta2,
                      'u1_otn': u1_otn, 'mu': mu, 'fi': fi, 'psi': psi, 'ro': ro,
                      'a_kr0': iwt['a_kr0'], 'c_ad': iwt['c_ad'], 'u1': iwt['u1'], 'c1_t': iwt['c1_t'],
                      'c1': iwt['c1'], 'Lambda_c1': iwt['Lambda_c1'], 'T1': iwt['T1'], 'P1': iwt['P1'],
                      'betta1': iwt['betta1'], 'a1': iwt['a1'], 'w1': iwt['w1'], 'T_w1_': iwt['T_w1_'],
                      'P_w1_': iwt['P_w1_'], 'D1': iwt['D1'],
                      'u2': owt['u2'], 'h_rk': owt['h_rk'], 'w2_t': owt['w2_t'], 'w2': owt['w2'],
                      'a_kr_w2': owt['a_kr_w2'], 'Lambda_w2': owt['Lambda_w2'], 'alfa2': owt['alfa2'],
                      'T2': owt['T2'], 'T_w2_': owt['T_w2_'], 'c2': owt['c2'], 'M_c2': owt['M_c2'],
                      'D2': owt['D2'],
                      'h_u': h_u, 'etta_u': etta_u, 'l1': l1, 'l2': l2, 'R2_p': R2_p, 'R2_v': R2_v,
                      'P2': P2, 'N': N}

        elif P2 is not None and N is not None:
            iwt = inlet_wheel_param(P0_=P0_, P2=P2, k=k, R=R, T0_=T0_, u1_otn=u1_otn, ro=ro, fi=fi, alfa1=alfa1, n=n)
            owt = outlet_wheel_param(u1=iwt['u1'], mu=mu, P2=P2, P_w1_=iwt['P_w1_'], k=k, R=R, T_w1_=iwt['T_w1_'],
                                     psi=psi, Lambda_u1=iwt['Lambda_u1'], Tau_w1=iwt['Tau_w1'],
                                     Tau_c1=iwt['Tau_c1'], betta2=betta2, D1=iwt['D1'])

            # работа ступени на окружности колеса
            h_u = iwt['c1'] * iwt['u1'] * math.cos(math.radians(alfa1)) + owt['c2'] * owt['u2'] * math.cos(
                math.radians(owt['alfa2']))

            # кпд на окружности колеса
            etta_u = 2 * h_u / math.pow(iwt['c_ad'], 2)

            G = get_G(N, iwt['c_ad'], etta_u)

            # высота лопаток на выходе из соплового аппарата
            l1 = G * R * iwt['T1'] / (
                math.pi * iwt['D1'] * iwt['c1'] * math.sin(math.radians(alfa1)) * iwt['P1'] * math.pow(10, 3))

            # высота лопаток на выходе из рабочего колеса
            l2 = G * R * owt['T2'] / (
                math.pi * owt['D2'] * owt['w2'] * math.sin(math.radians(betta2)) * P2 * math.pow(10, 3))

            # радиус перефирии на выходе из рабочего колеса
            R2_p = (owt['D2'] + l2) / 2

            # радиус втулки на выходе из рабочего колеса
            R2_v = (owt['D2'] - l2) / 2

            result = {'N': N, 'P0_': P0_, 'P2_': owt['P2_'], 'P2': P2, 'alfa1': alfa1, 'betta2': betta2,
                      'u1_otn': u1_otn, 'mu': mu, 'fi': fi, 'psi': psi, 'ro': ro,
                      'a_kr0': iwt['a_kr0'], 'c_ad': iwt['c_ad'], 'u1': iwt['u1'], 'c1_t': iwt['c1_t'],
                      'c1': iwt['c1'], 'Lambda_c1': iwt['Lambda_c1'], 'T1': iwt['T1'], 'P1': iwt['P1'],
                      'betta1': iwt['betta1'], 'a1': iwt['a1'], 'w1': iwt['w1'], 'T_w1_': iwt['T_w1_'],
                      'P_w1_': iwt['P_w1_'], 'D1': iwt['D1'],
                      'u2': owt['u2'], 'h_rk': owt['h_rk'], 'w2_t': owt['w2_t'], 'w2': owt['w2'],
                      'a_kr_w2': owt['a_kr_w2'], 'Lambda_w2': owt['Lambda_w2'], 'alfa2': owt['alfa2'],
                      'T2': owt['T2'], 'T_w2_': owt['T_w2_'], 'c2': owt['c2'], 'M_c2': owt['M_c2'],
                      'D2': owt['D2'],
                      'h_u': h_u, 'etta_u': etta_u, 'l1': l1, 'l2': l2, 'R2_p': R2_p, 'R2_v': R2_v, 'G': G}

        else:
            raise ValueError(u'Недостаточно параметров')
    elif G is not None and N is not None:
        if P2_ is not None:
            P0_ = 2.1 * P2_
            while True:
                P2 = get_P2(P0_=P0_, P2_=P2_, k=k, R=R, T0_=T0_, u1_otn=u1_otn, ro=ro, fi=fi, alfa1=alfa1, mu=mu,
                            psi=psi, betta2=betta2, n=n)
                result = solution(k, R, T0_, n, alfa1, betta2, u1_otn, mu, ro, fi, psi, G=G, N=None, P0_=P0_, P2=P2,
                                  P2_=None)
                if round(result['N'], 2) != round(N, 2):
                    P0_ = P0_ * N / result['N']
                else:
                    break
        elif P2 is not None:
            P0_ = 2.1 * P2
            while True:
                result = solution(k, R, T0_, n, alfa1, betta2, u1_otn, mu, ro, fi, psi, G=G, N=None, P0_=P0_, P2=P2,
                                  P2_=None)
                if round(result['N'], 2) != round(N, 2):
                    P0_ = P0_ * N / result['N']
                else:
                    break
        else:
            raise ValueError(u'Недостаточно параметров')
    else:
        raise ValueError(u'Недостаточно параметров')
    return result


# обертка, что бы была возможность передать словарь аргументов при использовании map
def solution_wraper(kwargs):
    return solution(**kwargs)


def optimization(pool, k, R, T0_, n, mu, alfa1, betta2, fi, psi, u1_otn, ro, G=None, N=None, P0_=None, P2=None,
                 P2_=None):
    if isinstance(mu, (int, float)):
        mu = [mu, mu+0.01, 0.01]
    if isinstance(alfa1, (int, float)):
        alfa1 = [alfa1, alfa1+1.0, 1.0]
    if isinstance(betta2, (int, float)):
        betta2 = [betta2, betta2+1.0, 1.0]
    if isinstance(fi, (int, float)):
        fi = [fi, fi+0.01, 0.01]
    if isinstance(psi, (int, float)):
        psi = [psi, psi+0.01, 0.01]
    if isinstance(u1_otn, (int, float)):
        u1_otn = [u1_otn, u1_otn+0.01, 0.01]
    if isinstance(ro, (int, float)):
        ro = [ro, ro+0.02, 0.02]
    # print(len(list({'k': k, 'R': R, 'T0_': T0_, 'n': n, 'G': G, 'N': N, 'P0_': P0_, 'P2': P2,
    #                 'P2_': P2_, 'mu': mu_each, 'alfa1': alfa1_each, 'betta2': betta2_each, 'fi': fi_each,
    #                 'psi': psi_each, 'u1_otn': u1_otn_each, 'ro': ro_each}
    #                for mu_each in myrange(mu[0], mu[1], mu[2])
    #                for alfa1_each in myrange(alfa1[0], alfa1[1], alfa1[2])
    #                for betta2_each in myrange(betta2[0], betta2[1], betta2[2])
    #                for fi_each in myrange(fi[0], fi[1], fi[2])
    #                for psi_each in myrange(psi[0], psi[1], psi[2])
    #                for u1_otn_each in myrange(u1_otn[0], u1_otn[1], u1_otn[2])
    #                for ro_each in myrange(ro[0], ro[1], ro[2]))))
    gen_of_init = ({'k': k, 'R': R, 'T0_': T0_, 'n': n, 'G': G, 'N': N, 'P0_': P0_, 'P2': P2,
                    'P2_': P2_, 'mu': mu_each, 'alfa1': alfa1_each, 'betta2': betta2_each, 'fi': fi_each,
                    'psi': psi_each, 'u1_otn': u1_otn_each, 'ro': ro_each}
                   for mu_each in myrange(mu[0], mu[1], mu[2])
                   for alfa1_each in myrange(alfa1[0], alfa1[1], alfa1[2])
                   for betta2_each in myrange(betta2[0], betta2[1], betta2[2])
                   for fi_each in myrange(fi[0], fi[1], fi[2])
                   for psi_each in myrange(psi[0], psi[1], psi[2])
                   for u1_otn_each in myrange(u1_otn[0], u1_otn[1], u1_otn[2])
                   for ro_each in myrange(ro[0], ro[1], ro[2]))
    list_of_solution = []
    list_of_etta_u = []
    count_solution = 0
    err = 0
    math_domain_error = 0
    imap_of_res = pool.imap_unordered(solution_wraper, gen_of_init)
    while True:
        try:
            result = imap_of_res.next()
            if result['R2_p'] <= 0.0723 and result['Lambda_c1'] < 1.05 and result['l1'] <= 0.0147 and \
               result['l2'] <= 0.0175 and (result['D1'] + result['l1']) <= 0.141:
                list_of_solution.append(result)
                list_of_etta_u.append(result['etta_u'])
        except ValueError as error:
            # TODO: выяснить причину возникновения данного исключения и устранить ее
            if error.args[0] == 'math domain error':
                math_domain_error += 1
            else:
                err += 1
        except StopIteration:
            break
        count_solution += 1
        if not count_solution % 1000000:
            print(count_solution)
    return {'list_of_solution': list_of_solution, 'list_of_etta_u': list_of_etta_u, 'count_solution': count_solution,
            'err': err, 'math_domain_error': math_domain_error}


if __name__ == '__main__':
    start = time.time()
    multiprocessing.freeze_support()
    k = 1.34
    R = 287.0
    T0_ = 1223
    n = 60000.0
    mu = [0.98, 1.03, 0.01]
    alfa1 = [10, 31, 1]
    betta2 = [20, 40, 1]
    fi = 0.96
    psi = 0.93
    u1_otn = [0.4, 0.60, 0.01]
    ro = [0.35, 0.5, 0.02]
    G = 0.85
    P2_ = 121.0
    P0_ = 370
    pool = multiprocessing.Pool(multiprocessing.cpu_count())
    res = optimization(pool=pool, k=k, R=R, T0_=T0_, n=n, mu=mu, alfa1=alfa1, betta2=betta2, fi=fi, psi=psi,
                       u1_otn=u1_otn, ro=ro, G=G, P2_=P2_, P0_=P0_)
    pool.close()
    pool.join()
    print('count_solution {0}'.format(res['count_solution']))
    print('math_domain_error {0}'.format(res['math_domain_error']))
    print('err {0}'.format(res['err']))
    print('list_of_solution {0}'.format(len(res['list_of_solution'])))
    etta = max(res['list_of_etta_u'])
    for i in res['list_of_solution']:
        if i['etta_u'] == etta:
            print('etta {0}'.format(etta))
            print('Lambda_c1 {0}'.format(i['Lambda_c1']))
            print('D1 {0}'.format(i['D1']))
            print('R2_p {0}'.format(i['R2_p']))
            print('R2_v {0}'.format(i['R2_v']))
            print('l1 {0}'.format(i['l1']))
            print('G {0}'.format(i['G']))
            print('------------')
            print('mu {0}'.format(i['mu']))
            print('alfa1 {0}'.format(i['alfa1']))
            print('betta2 {0}'.format(i['betta2']))
            print('fi {0}'.format(i['fi']))
            print('psi {0}'.format(i['psi']))
            print('u1_otn {0}'.format(i['u1_otn']))
            print('ro {0}'.format(i['ro']))
            print('P2 {0}'.format(i['P2']))
            print('P0_ {0}'.format(i['P0_']))
            print('N {0}'.format(i['N']))
            print('l2 {0}'.format(i['l2']))
            print('u1 {0}'.format(i['u1']))
            csv_file = open('1.csv', 'wt')
            try:
                csv_dict = csv.DictWriter(csv_file, ['key', 'value'], delimiter=';')
                csv_dict.writerows([{'key': k, 'value': str(v).replace('.', ',')} for k, v in i.items()])
            except Exception as e:
                print(e)
            finally:
                csv_file.close()
    print('{0} sec'.format(time.time() - start))
    # raw_input('')

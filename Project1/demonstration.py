import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

data_by_dim = {}
cordsNSKO = []


def ordinal_to_color(n):
    # Список цветов на английском языке
    colors = [
        "red",  # 0
        "green",  # 1
        "blue",  # 2
        "yellow",  # 3
        "orange",  # 4
        "purple",  # 5
        "cyan",  # 6
        "pink",  # 7
        "brown",  # 8
        "gray",  # 9
        "black",  # 10
        "white"  # 11
    ]

    # Если n больше длины списка, используем остаток от деления
    return colors[n % len(colors)]


with open('obrazi.txt', 'r') as file:
    for line in file:
        line = line.strip()
        if not line:
            continue
        coords = [float(num) for num in line.split()]
        dim = len(coords)
        data_by_dim.setdefault(dim, []).append(coords)

params = np.loadtxt('nsko.txt')     # файл с весами и смещением
weights = params[:2]                # первые два значения - веса
bias = params[2]                    # третье значение - смещение

# Загрузка данных точек
data = np.loadtxt('nskoCords.txt')     # файл с координатами и целевыми значениями
X = data[:, :-1]                        # координаты точек
y = data[:, -1]                   # целевые значения

# Определение границы прямой
x_min, x_max = X[:, 0].min(), X[:, 0].max()
xx = np.linspace(x_min, x_max, 100)
yy = -(weights[0] * xx + bias) / weights[1]  # y = -(w1*x + b) / w2

for dim, points in data_by_dim.items():
    if dim == 3:
        x_vals, y_vals, class_num = zip(*points)
        for n in class_num:
            n = int(n)
            n = ordinal_to_color(n)
        plt.figure()
        plt.scatter(x_vals, y_vals, c=class_num)
        #plt.plot(xx, yy, color='black', label='Decision Boundary')
        plt.title(f'2D Scatter Plot')
        plt.xlabel('X')
        plt.ylabel('Y')
        plt.grid(True)
        plt.show()
    elif dim == 4:
        x_vals, y_vals, z_vals, class_num = zip(*points)
        for n in class_num:
            n = int(n)
            n = ordinal_to_color(n)
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(x_vals, y_vals, z_vals, c=class_num)
        plt.title('3D Scatter Plot')
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        plt.show()
    else:
        print(f'Cannot plot data with dimension {dim}')
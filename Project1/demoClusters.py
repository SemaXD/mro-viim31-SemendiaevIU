import numpy as np
import matplotlib.pyplot as plt


def load_points(filename):
    """
    Загружает координаты точек из файла.
    :param filename: Название файла с точками
    """
    points = []
    with open(filename, 'r') as f:
        for line in f:
            coords = list(map(float, line.strip().split()))
            points.append(coords[:-1])  # игнорируем последнюю координату
    return np.array(points)


def load_centers(filename):
    """
    Загружает координаты центров кластеров из файла.
    :param filename: Название файла с центрами кластеров
    """
    centers = []
    with open(filename, 'r') as f:
        for line in f:
            coords = list(map(float, line.strip().split()))
            centers.append(coords)
    return np.array(centers)


def load_labels(filename):
    """
    Загружает метки кластеров из файла.
    :param filename: Название файла с метками кластеров
    """
    labels = []
    with open(filename, 'r') as f:
        for line in f:
            labels.append(int(line.strip()))
    return np.array(labels)


def plot_clusters(points, centers, labels):
    """
    Визуализирует точки и центры кластеров.
    :param points: Массив с точками
    :param centers: Массив с центрами кластеров
    :param labels: Массив с метками кластеров
    """
    unique_labels = np.unique(labels)

    # Генерация цветов для каждого кластера
    colors = plt.cm.get_cmap('tab10', len(unique_labels))

    # Создание графика
    fig = plt.figure()

    if points.shape[1] == 2:  # 2D
        ax = fig.add_subplot(111)

        # Отображение точек с меткой -1 (шум)
        noise_points = points[labels == -1]
        ax.scatter(noise_points[:, 0], noise_points[:, 1], color='black', label='Noise', marker='x')

        for label in unique_labels:
            if label == -1:
                continue  # Пропускаем шумовые точки
            cluster_points = points[labels == label]
            ax.scatter(cluster_points[:, 0], cluster_points[:, 1],
                       color=colors(label), label=f'Cluster {label}')

        # Отображение центров кластеров
        ax.scatter(centers[:, 0], centers[:, 1],
                   color='black', marker='X', s=100, label='Centers')

        #ax.legend()
        ax.set_title('Кластеризация - 2D')
        ax.set_xlabel('X-ось')
        ax.set_ylabel('Y-ось')

    elif points.shape[1] == 3:  # 3D
        ax = fig.add_subplot(111, projection='3d')

        # Отображение точек с меткой -1 (шум)
        noise_points = points[labels == -1]
        ax.scatter(noise_points[:, 0], noise_points[:, 1], noise_points[:, 2],
                   color='black', label='Noise', marker='x')

        for label in unique_labels:
            if label == -1:
                continue  # Пропускаем шумовые точки
            cluster_points = points[labels == label]
            ax.scatter(cluster_points[:, 0], cluster_points[:, 1], cluster_points[:, 2],
                       color=colors(label), label=f'Cluster {label}')

        # Отображение центров кластеров
        ax.scatter(centers[:, 0], centers[:, 1], centers[:, 2],
                   color='black', marker='X', s=100, label='Centers')

        #ax.legend()
        ax.set_title('Кластеризация - 3D')
        ax.set_xlabel('X-ось')
        ax.set_ylabel('Y-ось')
        ax.set_zlabel('Z-ось')

    plt.show()


# Загрузка данных из файлов
points = load_points('obrazi.txt')

# kMeans
centers1 = load_centers('centers_kMeans.txt')
labels1 = load_labels('labels_kMeans.txt')

# FOREL
centers2 = load_centers('centers_forel.txt')
labels2 = load_labels('labels_forel.txt')

# ISODATA
centers3 = load_centers('centers_isodata.txt')
labels3 = load_labels('labels_isodata.txt')

# Визуализация кластеров
plot_clusters(points, centers1, labels1)
plot_clusters(points, centers2, labels2)
plot_clusters(points, centers3, labels3)
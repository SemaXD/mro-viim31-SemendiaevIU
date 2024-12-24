#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <ctime>
#include <cmath>
#include <algorithm>
//#include <limits>
//#include <cstdlib>


// Задачи:
// 1. Сделать программу более удобной и понятной для пользователя -- На очереди
// 2. Потоковый ввод ----------------------------------------------- Отложено, причина: вечная ошибка 3221225477, решить пока не удалось
// 3. Сделать схему данных для классов, терюсь в иерархии ---------- Отложено, причина: А надо ли?
// 4. Подписать все функции комментариями -------------------------- WIP
// 5. Упростить параметры для функций ------------------------------ WIP
// 6. Размеры массивов перевести в size_t тип ---------------------- На очереди
// 7. Переписать векторы с использованием конструктора ------------- На очереди
// 8. Пора бы поменять vector на array там, где можно -------------- На очереди
// 9. Разбить код на разные файлы, ибо глючит компилятор ----------- На очереди

// ?1 При чтении координат нет обработчика лишних пробелов до/после строки, ВОЗМОЖНО  стоит сделать?
// ?2 Хо-кашьяп работает подозрительно, возможно есть математическая ошибка

//Класс для входных данных генератора случайных образов
class data {
    private:
        int _NKlassov, _NIzmerenii;
        std::vector<int> _NElementov;

    public:
        // get - functions

        int                 getNumberOfClasses()                { return _NKlassov; }
        int                 getNumberOfDimensions()             { return _NIzmerenii; }
        std::vector<int>    getAllNumbersOfElements()           { return _NElementov; }
        int                 getNumberOfElements(int numOfClass) { return _NElementov[numOfClass]; }

        //set - functions

        void                setNumberOfClasses(int numOfClass)                          { _NKlassov = numOfClass; }
        void                setNumberOfDimensions(int numOfDim)                         { _NIzmerenii = numOfDim; }
        void                setAllNumbersOfElements(std::vector<int> listOfElements)    { _NElementov = listOfElements; }
        void                setNumberOfElements(int classID, int numOfElements)         { _NElementov[classID] = numOfElements; }
};

// Класс для точек 
class point {     
    private:
        std::vector<int> _crd;
    public:
        //get - functions

        std::vector<int> getAllCords()          { return _crd; }
        int              getCord(int numOfCord) { return _crd[numOfCord]; }

        //set - functions

        void             setAllCords(std::vector<int> listOfCords) { _crd = listOfCords; }
        void             setCord(int cordID, int cord)             { _crd[cordID] = cord; }
};

// Класс для массива точек
class obrazi {     
    private:
        std::vector<point> _obraz;
    public:
        //get - functions

        std::vector<point>  getAllPoints()           { return _obraz; }
        point               getPoint(int numOfPoint) { return _obraz[numOfPoint]; }

        //set - functions

        void                setAllPoints(std::vector<point> listOfPoints)   { _obraz=listOfPoints; }
};

//Функция рандомизатор
int randoms(int start, int end) {      
    return rand() % (end - start + 1) + start;
};

//Запись точек для nsko в отдельный файл
void nskoCordWriting(std::vector<std::vector<int>>& mat, std::vector<int> trg) { 
    std::ofstream outFile("nskoCords.txt");
    if (outFile.is_open()) {
        int numOfRows = mat.size();
        int numOfColumns = mat[0].size();
         for (int i = 0; i < numOfRows;i++) {
            for (int j = 0; j<numOfColumns;j++) {
                outFile<< mat[i][j] << " ";
            }
        outFile<<trg[i]<<std::endl;
    }
    }
}

//Алгоритм Хо-Кашьяпа
std::vector<double> myNSKO(int numFirst, int numSecond, std::vector<std::vector<int>>& inputData, double& learning_rate, int& max_iterations) {
    std::vector<std::vector<int>> matrix;
    int numOfColumns = inputData[0].size();
    int numOfRows = inputData.size();
    std::vector<int> blunk; blunk.resize(numOfColumns-1);
    std::vector<int> classTargets;

    for (int i = 0;i<numOfRows;i++) { //Заполняем первую половину матрицы элементами 1ого класса
        if (inputData[i][numOfColumns-1]==numFirst) {
            classTargets.push_back(1);
            for(int j = 0; j<numOfColumns-1;j++) blunk[j]=inputData[i][j];
            matrix.push_back(blunk);
        }
    }

     for (int i = 0;i<numOfRows;i++) { //Заполняем вторую половину матрицы элементами 2ого класса
        if (inputData[i][numOfColumns-1]==numSecond) {
            classTargets.push_back(-1);
            for(int j = 0; j<numOfColumns-1;j++) blunk[j]=inputData[i][j];
            matrix.push_back(blunk);
        }
    }

    nskoCordWriting(matrix, classTargets);
    int num_samples = matrix.size();
    int num_features = matrix[0].size();
    
    // Инициализация весов и смещения
    std::vector<double> weights(num_features, 0.0);
    double bias = 0.0;

    for (int iter = 0; iter < max_iterations; ++iter) {
        bool has_updates = false;

        for (int i = 0; i < num_samples; ++i) {
            // Вычисляем предсказание
            double prediction = 0.0;
            for (int j = 0; j < num_features; ++j) {
                prediction += weights[j] * matrix[i][j];
            }
            prediction += bias;

            // Проверяем условие корректности классификации
            if (classTargets[i] * prediction <= 0) {
                // Обновляем веса и смещение
                for (int j = 0; j < num_features; ++j) {
                    weights[j] += learning_rate * classTargets[i] * matrix[i][j];
                }
                bias += learning_rate * classTargets[i];
                has_updates = true;
            }
        }

        // Если не было обновлений, выходим из цикла
        if (!has_updates) {
            break;
        }
    }

    // Возвращаем обученные веса и смещение
    weights.push_back(bias); // Добавляем смещение в конец векторов весов
    return weights;
}

// Чтение данных из txt
data readingFirst(std::string fileName) {    
    std::ifstream in(fileName);
    data tempDate;
    tempDate.setNumberOfClasses(0);
    if (in.is_open()){
        std::string line;

        getline(in, line);
        tempDate.setNumberOfClasses(stoi(line)); 
        tempDate.getAllNumbersOfElements().resize(tempDate.getNumberOfClasses());

        std::vector<int> listOfElements;
        listOfElements.resize(tempDate.getNumberOfClasses());
        for (int i = 0; i<tempDate.getNumberOfClasses(); i++) { 
            getline(in, line);
            listOfElements[i] = stoi(line);
        }
        tempDate.setAllNumbersOfElements(listOfElements);

        getline(in, line);
        tempDate.setNumberOfDimensions(stoi(line));
    in.close(); }
    return tempDate;
} 

//Генератор случайных образов
std::vector<obrazi> generator(data& inputDate, std::vector<obrazi>& mainData, int MIN_COORD, int MAX_COORD, int MAX_RANGE) {  
    srand((time(0)));
    mainData.resize(inputDate.getNumberOfClasses());
    
    for (int i = 0; i < inputDate.getNumberOfClasses(); i++) {
        std::vector<point> points(inputDate.getNumberOfElements(i));
        
        for (int j = 0; j < inputDate.getNumberOfElements(i); j++) {
            points[j].setAllCords(std::vector<int>(inputDate.getNumberOfDimensions()+1));
        }
        
        // Задаём в первых двух точках минимум и максимум соответственно
        for (int izmerenia = 0; izmerenia < inputDate.getNumberOfDimensions(); izmerenia++) {
            points[0].setCord(izmerenia, randoms(MIN_COORD, MAX_COORD - MAX_RANGE)); 
            points[1].setCord(izmerenia, randoms(points[0].getCord(izmerenia), points[0].getCord(izmerenia) + MAX_RANGE));
        }

        // Генерируем последние точки
        for (int j = 2; j < inputDate.getNumberOfElements(i); j++) {
            for (int izmerenia = 0; izmerenia < inputDate.getNumberOfDimensions(); izmerenia++) {
                points[j].setCord(izmerenia, randoms(points[0].getCord(izmerenia), points[1].getCord(izmerenia)));
            }
        }

        //Записываем в каждую точку номер её класса
        for (int t = 0; t< inputDate.getNumberOfElements(i); t++) {
            points[t].setCord(inputDate.getNumberOfDimensions(), i);
        }

        mainData[i].setAllPoints(points);
    }
    return mainData;
}

//Запись образов в отдельный файл
void writtingSubjects(data& inputData, std::vector<obrazi>& mainData, std::string fileName) {
    //Запись в файл
    std::ofstream outFile(fileName);
    if (outFile.is_open()) {
        for (int i = 0; i<inputData.getNumberOfClasses(); i++) {
            for (int j = 0;j<inputData.getNumberOfElements(i);j++) {
                for (int izm = 0; izm<inputData.getNumberOfDimensions()+1;izm++) {
                    outFile << mainData[i].getPoint(j).getCord(izm)<<" ";
                }
                //outFile <<i;
                outFile <<std::endl;
            }
        }
    }
}

//Обобщающая функция генератора классов (чтение input-данных, генерация и запись образов)
void generateSubjects() {
    data inputs = readingFirst("start.txt");
    std::vector<obrazi> allData;
    allData = generator(inputs, allData, 0, 100, 50);
    writtingSubjects(inputs, allData, "obrazi.txt");
}

//Функция подсчёта слов(пробелов) в строке
int wordNumber(std::string& str) {
    int num = 0;
    int lenOfStr = str.length();
    for (int iter = 0;iter<lenOfStr;iter++) if (str[iter]==' ') num++;
    return num;
}

//Чтение матрицы из файла
void readingMatrix(std::vector<std::vector<double>>& matrixOfData, std::string fileName) {
    std::ifstream firstInput(fileName);
    int lenOfRow;
    std::string line;

    //Открываем файл и узнаём кол-во переменных в одной строке
    if (firstInput.is_open()) {              
        getline(firstInput, line);
        lenOfRow = wordNumber(line);
    }
    firstInput.close();

    std::vector<double> blunk; blunk.resize(lenOfRow);
    std::ifstream secondInput(fileName);
    if (secondInput.is_open()) {
        int numOfColumn = 0; int numOfRow = 0;
        matrixOfData.push_back(blunk);
        while (secondInput>>line) {
            if (numOfColumn<lenOfRow) {
                matrixOfData[numOfRow][numOfColumn] = stoi(line);
                numOfColumn++;
            }
            else {
                numOfColumn = 0;
                numOfRow++;
                matrixOfData.push_back(blunk);
                matrixOfData[numOfRow][numOfColumn] = stoi(line);
                numOfColumn++;
            }
        }
    }
}

//Вывод матрицы в консоль для проверки
void checkMatrix(std::vector<std::vector<double>>& matrix) {
    int numOfRows = matrix.size();
    int numOfColumns = matrix[0].size();
     for (int i = 0; i < numOfRows;i++) {
        for (int j = 0; j<numOfColumns;j++) {
            std::cout<< matrix[i][j] << "\t";
        }
        std::cout<<std::endl;
    }
}

//Запись весов и смещения в файл
void writeNSKO(std::vector<double>& res) {
    //Запись в файл
    std::ofstream outFile("nsko.txt");
    if (outFile.is_open()) {
        for (double weight : res) {
            outFile << weight << " ";
        }
    }
}

//Эвклидово рассторяние между двумя точками
double evclidDistance(std::vector<double> firstPoint, std::vector<double> secondPoint) { 
    double distance = 0.0;
    int numOfDim = firstPoint.size()-1;
    for (int i = 0; i < numOfDim; i++) {
        distance += std::pow(firstPoint[i]-secondPoint[i],2);
    }
    return std::sqrt(distance);
}

//Расчёт Манхеттонского смещения между двумя точками
double manhattanDistance(std::vector<double> firstPoint, std::vector<double> secondPoint) {
    double distance = 0.0;
    int numOfDim = firstPoint.size()-1;
    for (int i = 0; i < numOfDim; i++) {
        distance += std::abs(firstPoint[i]-secondPoint[i]);
    }
    return distance;
}

//Обобщающая функция для Хо-кашьяпа (необходимые переменные, сам алгоритм и запись в файл)
void hoKashiap (int fNum, int sNum, std::vector<std::vector<int>>& input) { 
    double learning_rate = 0.1;
    int max_iterations = 100;

    std::vector<double> result = myNSKO(fNum, sNum, input, learning_rate, max_iterations);
    writeNSKO(result);
}

// Функция для выполнения алгоритма K-средних
static std::pair<std::vector<std::vector<double>>, std::vector<int>> kMeans(std::vector<std::vector<double>>& data, std::vector<std::vector<double>> &klasterCenters, int max_iters = 100) {
    size_t n = data.size();
    size_t m = data[0].size()-1;
    int k = klasterCenters.size();

    // Инициализация центров кластеров

    // std::vector<std::vector<double>> centers(k, std::vector<double>(m));
    // std::vector<bool> chosen(n, false);
    // srand((time(0)));

    // for (int i = 0; i < k; i++) {
    //     int index;
    //     do {
    //         index = randoms(0,n-1);
    //     } while (chosen[index]);
    //     chosen[index] = true;
    //     centers[i] = data[index];
    // }

    //Переменная для номеров кластеров
    std::vector<int> labels(n);

    for (int iter = 0; iter < max_iters; iter++) {
        // Присвоение кластеров
        for (size_t i = 0; i < n; i++) {
            double minDist = std::numeric_limits<double>::max(); //Задаём максимально возможное значение double для minDist
            for (int j = 0; j < k; j++) {
                double dist = evclidDistance(data[i], klasterCenters[j]);
                if (dist < minDist) {
                    minDist = dist;
                    labels[i] = j;
                }
            }
        }

        // Новые центры кластеров
        std::vector<std::vector<double>> new_centers(k, std::vector<double>(m, 0.0));
        std::vector<int> counts(k, 0);

        for (size_t i = 0; i < n; i++) {
            int cluster_id = labels[i];
            counts[cluster_id]++;
            for (size_t j = 0; j < m; j++) {
                new_centers[cluster_id][j] += data[i][j];
            }
        }

        for (int j = 0; j < k; j++) {
            if (counts[j] > 0) {
                for (size_t l = 0; l < m; l++) {
                    new_centers[j][l] /= counts[j];
                }
            }
        }

        // Проверка сходимости
        if (klasterCenters == new_centers) {
            break;
        }
        klasterCenters = new_centers;
    }

    return {klasterCenters, labels};
}

//Функция для вывода kMeans
void checking_klusters(std::pair<std::vector<std::vector<double>>, std::vector<int>>& result) {
    std::cout << "Центры кластеров:\n";
    for (const auto& center : result.first) {
        for (const auto& value : center) {
            std::cout << value << "\t";
        }
        std::cout << std::endl; 
    }

    std::cout << "Метки кластеров:\n";
    for (const auto& label : result.second) {
        std::cout << label << " ";
    }
    std::cout << std::endl;
}

//Алгоритм простейшей расстановки центров кластеров 
std::vector<std::vector<double>> centersEasy(const std::vector<std::vector<double>> &points, double h, int k) {

    //Вектор для хранения центров кластеров
    std::vector<std::vector<double>> centers;
    centers.push_back(points[0]);       //Первый центр - первая точка

    for (size_t j = 1; j < points.size() && centers.size() < k; j++) { //Находим остальные центры
        bool isFarEnough = true;

        for (auto &center : centers) {     //Для всех центров проверяем удаление для данной точки 
            if (evclidDistance(points[j], center) <= h) {
                isFarEnough = false;
                break;
            }
        }

        if (isFarEnough) {  // Если точка достаточно удалена от всех центров, добавляем её
            centers.push_back(points[j]);
        }
    }

    return centers;
}

std::vector<std::vector<double>> centersSito(std::vector<std::vector<double>>& points, double h, int k) {

    //Количество точек
    size_t numOfPoints = points.size();

    //Вектор для хранения центров кластеров
    std::vector<std::vector<double>> centers;

    //Плотности точек
    std::vector<std::vector<double>> plotnost;
    plotnost.resize(numOfPoints);

    //Нахождение плотности точек
    for (int i = 0; i<numOfPoints; i++) {
        double sum = 0.0;
        for (int j=0; j<numOfPoints; j++) {
            double distance = evclidDistance(points[i],points[j]);
            if (distance < h) sum += (h * h - distance);
        }
        plotnost[i].push_back(1/(h*h)*sum);
        plotnost[i].push_back(i);
    }

    //Сортируем плотности точек и берём в качестве центров те, что находятся кучнее
    std::sort(plotnost.begin(),plotnost.end());
    for (int i = 0; i<k; i++) centers.push_back(points[plotnost[plotnost.size()-1-i][1]]);
    return centers;
}

double calculateQualityCriterion(const std::vector<std::vector<double>>& points, const std::vector<std::vector<double>>& centers) {
    double totalDistance = 0.0;
    
    for (auto& point : points) {
        double minDistance = std::numeric_limits<double>::max();
        for (auto& center : centers) {
            double distance = evclidDistance(point, center);
            if (distance < minDistance) {
                minDistance = distance;
            }
        }
        totalDistance += minDistance;
    }
    
    return totalDistance;
}

std::vector<std::vector<double>> centersMaximum(const std::vector<std::vector<double>>& points, double h, int k) {
    // Вектор для хранения центров кластеров
    std::vector<std::vector<double>> centers;

    // Количество точек
    size_t numOfPoints = points.size();

    // 1. Назначение первого центра
    centers.push_back(points[0]);

    // 2. Назначение второго центра
    size_t secondCenterIndex = 0;
    double maxDistance = -1.0;

    for (size_t j = 1; j < numOfPoints; j++) {
        double distance = evclidDistance(centers[0], points[j]);
        if (distance > maxDistance) {
            maxDistance = distance;
            secondCenterIndex = j;
        }
    }

    centers.push_back(points[secondCenterIndex]);

    // Критерий качества после выбора первых двух центров
    double quality = calculateQualityCriterion(points, centers);

    // 3. Назначение остальных центров
    for (int i = 2; i < k; i++) {
        size_t newCenterIndex = 0;
        maxDistance = -1.0;

        for (size_t j = 0; j < numOfPoints; j++) {
            double minDistanceToCenters = std::numeric_limits<double>::max();

            // Находим минимальное расстояние до существующих центров
            for (const auto& center : centers) {
                double distance = evclidDistance(center, points[j]);
                if (distance < minDistanceToCenters) {
                    minDistanceToCenters = distance;
                }
            }

            // Если минимальное расстояние больше максимального найденного, обновляем
            if (minDistanceToCenters > maxDistance) {
                maxDistance = minDistanceToCenters;
                newCenterIndex = j;
            }
        }
        centers.push_back(points[newCenterIndex]);

        // Критерий качества после добавления нового центра
        double q_new = calculateQualityCriterion(points, centers);

        // Проверка условия останова
        if (q_new / quality >= h) {
            break; 
        }

        quality = q_new; 
    }

    return centers;
}

std::vector<double> computeCentroid(const std::vector<std::vector<double>> &cluster) {
    std::vector<double> centroid(cluster[0].size()-1, 0.0);
    for (const auto &point : cluster) {
        for (size_t i = 0; i < point.size()-1; i++) {
            centroid[i] += point[i];
        }
    }
    for (size_t i = 0; i < centroid.size(); i++) {
        centroid[i] /= cluster.size();
    }
    return centroid;
}


// std::pair<std::vector<std::vector<double>>, std::vector<int>> forel(
//     std::vector<std::vector<double>> &points,
//     double R) 
// {
//     std::vector<int> clusterAssignments(points.size(), -1); // Вектор принадлежности к кластерам
//     std::vector<std::vector<double>> newCenters; // Вектор новых центров кластеров

//     std::vector<bool> clustered(points.size(), false); // Помечаем, какие точки уже кластеризованы

//     for (size_t i = 0; i < points.size(); i++) {
//         if (clustered[i]) continue; // Пропускаем уже кластеризованные точки

//         std::vector<std::vector<double>> currentCluster; // Текущий кластер
//         std::vector<int> currentClusterAssignments; // Идентификаторы текущего кластера

//         // Начинаем с текущей точки
//         std::vector<double> currentCenter = points[i];
//         currentCluster.push_back(currentCenter);
//         currentClusterAssignments.push_back(i);
//         clustered[i] = true;

//         while (true) {
//             std::vector<std::vector<double>> neighbors;

//             // Находим соседей в радиусе R
//             for (size_t j = 0; j < points.size(); j++) {
//                 if (!clustered[j] && evclidDistance(currentCenter, points[j]) <= R) {
//                     neighbors.push_back(points[j]);
//                     currentClusterAssignments.push_back(j);
//                     clustered[j] = true;
//                 }
//             }

//             if (neighbors.empty()) break; // Если соседей нет, выходим из цикла

//             // Вычисляем новый центр тяжести
//             currentCenter = computeCentroid(neighbors);
//             currentCluster.insert(currentCluster.end(), neighbors.begin(), neighbors.end());
//         }

//         // Добавляем текущий кластер и его идентификаторы
//         newCenters.push_back(computeCentroid(currentCluster));
//         for (auto &id : currentClusterAssignments) {
//             clusterAssignments[id] = newCenters.size() - 1; // Устанавливаем принадлежность к кластеру
//         }
//     }

//     return {newCenters, clusterAssignments}; // Возвращаем новые центры кластеров и принадлежности
// }

std::pair<std::vector<std::vector<double>>, std::vector<int>> forel(
    const std::vector<std::vector<double>>& points,
    double R) 
{
    std::vector<int> clusterAssignments(points.size(), -1); // Вектор принадлежности к кластерам
    std::vector<std::vector<double>> newCenters; // Вектор новых центров кластеров

    std::vector<bool> clustered(points.size(), false); // Помечаем, какие точки уже кластеризованы

    for (size_t i = 0; i < points.size(); i++) {
        if (clustered[i]) continue; // Пропускаем уже кластеризованные точки

        std::vector<std::vector<double>> currentCluster; // Текущий кластер
        std::vector<int> currentClusterAssignments; // Идентификаторы текущего кластера

        // Начинаем с текущей точки
        std::vector<double> currentCenter = points[i];
        currentCluster.push_back(currentCenter);
        currentClusterAssignments.push_back(i);
        clustered[i] = true;

        bool centerChanged = true;

        while (centerChanged) {
            centerChanged = false;
            std::vector<std::vector<double>> neighbors;

            // Находим соседей в радиусе R
            for (size_t j = 0; j < points.size(); j++) {
                if (!clustered[j] && evclidDistance(currentCenter, points[j]) <= R) {
                    neighbors.push_back(points[j]);
                    currentClusterAssignments.push_back(j);
                    clustered[j] = true;
                }
            }

            if (!neighbors.empty()) {
                // Вычисляем новый центр тяжести
                std::vector<double> newCenter = computeCentroid(neighbors);

                // Проверяем изменение центра
                if (evclidDistance(currentCenter, newCenter) > 1e-6) { // Используем небольшой порог для сравнения
                    currentCenter = newCenter;
                    centerChanged = true; // Центр изменился, продолжаем итерацию
                    currentCluster.clear(); // Очищаем текущий кластер для перезаполнения
                    for(int ikr : currentClusterAssignments) {
                        clustered[ikr] = false;
                    }
                    currentClusterAssignments.clear();
                } else {
                    break; // Центр не изменился, выходим из цикла
                }
            } else {
                break; // Если соседей нет, выходим из цикла
            }

            // Добавляем новые соседи к текущему кластеру
            currentCluster.insert(currentCluster.end(), neighbors.begin(), neighbors.end());
        }

        // Добавляем текущий кластер и его идентификаторы только если он не пустой
        if (!currentCluster.empty()) {
            newCenters.push_back(computeCentroid(currentCluster));
            for (auto &id : currentClusterAssignments) {
                clusterAssignments[id] = newCenters.size() - 1; // Устанавливаем принадлежность к кластеру
            }
        }
    }

    return {newCenters, clusterAssignments}; // Возвращаем новые центры кластеров и принадлежности
}

//Запись весов и смещения в файл
void writeCenters(std::pair<std::vector<std::vector<double>>, std::vector<int>>& res, std::string centersName, std::string labelsName) {
    //Запись в файл
    std::ofstream outFile(centersName);
    if (outFile.is_open()) {
        for (std::vector<double> center : res.first) {
            for (double cord : center) {
                outFile << cord << " "; 
            }
            outFile<<std::endl;
        }
    }
    
    std::ofstream outFile_(labelsName);
    if (outFile_.is_open()) {
        for (int label : res.second) {
            outFile_ << label << std::endl;
        }
    }
}

std::pair<std::vector<std::vector<double>>, std::vector<int>> isodata(
    std::vector<std::vector<double>>& points, // Массив точек
    std::vector<std::vector<double>>& centers, // Массив центров
    int max_iter = 100, // Максимальное количество итераций при присвоении точек к центрам кластеров
    int min_cluster_size = 5,  // Минимальное количество точек в кластере
    double distance_threshold = 25.0, // Порог расстояния между кластерами. Если он не будет превышен, то кластера будут объединены в один
    double std_dev_threshold = 15.0) { // Порог отклонения точек от новых центров. Если он будет превышен - будет создан новый класс

    int k = centers.size(); // Количество кластеров, определяется по размеру массива центров
    std::vector<int> labels(points.size(), -1); // Инициализация меток кластеров, изначально значения меток = -1

    for (int iter = 0; iter < max_iter; iter++) {
        // Присвоение точек к ближайшим центрам кластеров
        for (size_t i = 0; i < points.size(); i++) {
            double min_distance = std::numeric_limits<double>::max(); // Переменная минимальной дистанции, изначально хранит максимально возможное double-значение
            for (int j = 0; j < k; j++) {
                double distance = evclidDistance(points[i], centers[j]); // Переменная дистанции точки до центра кластера
                if (distance < min_distance) { 
                    min_distance = distance; 
                    labels[i] = j; 
                }
            }
        }

        // Обновление центров кластеров
        std::vector<std::vector<double>> new_centers(k, std::vector<double>(points[0].size()-1, 0.0)); // Новый массив центров кластеров, ВАЖНО - размерность взята с первой точки, а значит изначально размерность равноа nDimensions + метка класса
        std::vector<int> counts(k, 0); // Счетчики для каждого кластера, изначально = 0

        for (size_t i = 0; i < points.size(); i++) {
            if (labels[i] != -1) { 
                counts[labels[i]]++; 
                for (size_t j = 0; j < points[0].size()-1; j++) {
                    new_centers[labels[i]][j] += points[i][j]; // Суммируем координаты точек в кластере
                }
            }
        }

        for (int i = 0; i < k; i++) {
            if (counts[i] > 0) { 
                for (size_t j = 0; j < new_centers[i].size(); j++) {
                    new_centers[i][j] /= counts[i]; // Вычисляем среднее значение для нового центра кластера
                }
            }
        }

        // Проверка на минимальный размер кластера и обновление меток
        for (int i = 0; i < k; i++) {
            if (counts[i] < min_cluster_size) { // Если размер кластера меньше минимального порога
                for (size_t j = 0; j < labels.size(); j++) {
                    if (labels[j] == i) {
                        labels[j] = -1; // Удаляем метку для малых кластеров
                    }
                }
            }
        }

        // Разделение и объединение кластеров
        for (int i = 0; i < k; i++) {
            if (counts[i] > 0) { // Если в кластере есть точки
                double std_dev = 0.0;
                for (size_t j = 0; j < points.size(); j++) {
                    if (labels[j] == i) {
                        std_dev += evclidDistance(points[j], new_centers[i]); // Суммируем расстояния до нового центра
                    }
                }
                std_dev /= counts[i]; // Вычисляем среднее стандартное отклонение

                if (std_dev > std_dev_threshold) { // Если стандартное отклонение превышает пороговое значение
                    centers.push_back(new_centers[i]); // Добавляем новый центр в список центров кластеров
                    k++; // Увеличиваем количество кластеров
                }
            }

            for (int j = i + 1; j < k; j++) { 
                double distance = evclidDistance(centers[i], centers[j]); 
                if (distance < distance_threshold) { 
                    for (size_t m = 0; m < centers[i].size(); m++) {
                        centers[i][m] += centers[j][m]; 
                    }
                    for (size_t m = 0; m < centers[i].size(); m++) {
                        centers[i][m] /= 2.0; 
                    }
                    centers.erase(centers.begin() + j); 
                    k--; 
                    break;
                }
            }
        }

        centers.swap(new_centers); 

        if (new_centers == centers)
            break;
    }

    return {centers, labels}; 
}












int activation_function(double sum) {
    return (sum >= 0) ? 1 : -1;
}

std::vector<double> train_perceptron(const std::vector<std::vector<double>>& data, const std::vector<int> labels, double learning_rate, int max_iterations) {
    int n = data.size();
    std::vector<double> weights = {0.0, 0.0}; // Инициализация весов
    double bias = 0.0; // Инициализация смещения

    for (int iter = 0; iter < max_iterations; iter++) {
        bool converged = true;
        for (int i = 0; i < n; i++) {
            double sum = weights[0] * data[i][0] + weights[1] * data[i][1] + bias;
            int predicted_class = activation_function(sum);

            if (predicted_class != labels[i]) {
                weights[0] += learning_rate * labels[i] * data[i][0];
                weights[1] += learning_rate * labels[i] * data[i][1];
                bias += learning_rate * labels[i];
                converged = false;
            }
        }
        if (converged) break;
    }
    return weights;
}


// class Perceptron {
// public:
//     Perceptron(size_t inputSize, double learningRate)
//         : weights(inputSize), learningRate(learningRate) {
//         // Инициализируем веса случайными значениями
//         std::srand(static_cast<unsigned>(std::time(0)));
//         for (size_t i = 0; i < inputSize; i++) {
//             weights[i] = static_cast<double>(std::rand()) / RAND_MAX; // Случайные значения от 0 до 1
//         }
//     }

//     // Функция активации (шаговая функция)
//     int activate(const std::vector<double>& inputs) {
//         double sum = 0.0;
//         for (size_t i = 0; i < inputs.size(); i++) {
//             sum += inputs[i] * weights[i];
//         }
//         return sum >= 0.5 ? 1 : 0; // Возвращаем 1, если сумма больше или равна 0.5, иначе 0
//     }

//     // Обучение персептрона
//     void train(const std::vector<std::vector<double>>& trainingData,
//                const std::vector<int>& labels, size_t epochs) {
//         for (size_t epoch = 0; epoch < epochs; epoch++) {
//             for (size_t i = 0; i < trainingData.size(); i++) {
//                 int prediction = activate(trainingData[i]);
//                 int error = labels[i] - prediction;

//                 // Обновляем веса
//                 for (size_t j = 0; j < weights.size(); j++) {
//                     weights[j] += learningRate * error * trainingData[i][j];
//                 }
//             }
//         }
//     }

// private:
//     std::vector<double> weights; // Веса персептрона
//     double learningRate;         // Скорость обучения
// };

void generateSubjectsPerceptron() {
    data inputs = readingFirst("start_pers.txt");
    std::vector<obrazi> allData;
    allData = generator(inputs, allData, 0, 100, 10);
    writtingSubjects(inputs, allData, "obrazi_pers.txt");
}

std::pair<std::vector<std::vector<double>>, std::vector<int>> splitPointsAndLabels(const std::vector<std::vector<double>>& points) {
    std::vector<std::vector<double>> coordinates; // Вектор для хранения точек
    std::vector<int> labels; // Вектор для хранения идентификаторов классов

    for (const auto& point : points) {
        // Добавляем все элементы, кроме последнего, в вектор координат
        std::vector<double> coord(point.begin(), point.end() - 1);
        coordinates.push_back(coord);
        
        // Добавляем последний элемент (идентификатор класса) в вектор меток
        labels.push_back(point.back());
    }
    
    for (int i = 0; i < labels.size(); i++) {
        if (labels[i]==0) {
            labels[i] = 1;
        } else {
            labels[i] = -1;
        }
    }

    return {coordinates, labels}; // Возвращаем пару из координат и меток
}

//Исполнительная функция
int main() { 
    //generateSubjects();
    //generateSubjectsPerceptron();

    //Матрица для данных о точках - их координат и принадлежности по классу
    std::vector<std::vector<double>> inputMatrix;
    readingMatrix(inputMatrix, "obrazi.txt");
    //checkMatrix(inputMatrix);

    //hoKashiap(0, 1, inputMatrix);

    // std::cout<<evclidDistance(inputMatrix[0], inputMatrix[1])<< " Эвклид"<<std::endl;
    // std::cout<<manhattanDistance(inputMatrix[0], inputMatrix[1])<< " Манхеттон";
    
    auto centers = centersSito(inputMatrix, 25, 4);
    auto result = kMeans(inputMatrix, centers);
    checking_klusters(result); writeCenters(result, "centers_kMeans.txt", "labels_kMeans.txt");
    result = forel(inputMatrix, 20);
    checking_klusters(result); writeCenters(result, "centers_forel.txt", "labels_forel.txt");
    result = isodata(inputMatrix, centers);
    checking_klusters(result); writeCenters(result, "centers_isodata.txt", "labels_isodata.txt");

    std::vector<std::vector<double>> inputMatrixPers;
    readingMatrix(inputMatrixPers, "obrazi_pers.txt");
    auto k = splitPointsAndLabels(inputMatrixPers);
    std::vector<std::vector<double>> tochki = k.first;
    std::vector<int> labels = k.second;

    std::vector<double> weigths = train_perceptron(tochki, labels, 0.1, 100);
    double bias = 0.0;
    std::ofstream outfile("results_perc.txt");
    if (outfile.is_open()) {

        for (auto& point : tochki) {
            double sum = weigths[0] * point[0] + weigths[1]*point[1] + bias;
            int predicted_class = activation_function(sum);
            outfile << point[0] << " " << point[1] << " " << predicted_class << std::endl;
        }
        outfile.close();

        std::cout << "Готово, результат в results_perc.txt" << std::endl;

    } else {
        std::cerr << "Ошибка открытия файла!" << std::endl;
        return 1;
    }

    return 0;
}
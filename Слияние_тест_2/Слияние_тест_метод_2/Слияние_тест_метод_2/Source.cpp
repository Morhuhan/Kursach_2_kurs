// Тест слияния

#include <iostream>
#include <stack>
#include <algorithm>
#include <cmath>
#include <fstream>

using namespace std;

class point {

public:
    double x, y;

    // Указатель на соседа по часовой
    point* cw_next ;

    // Указатель на соседа против часовой
    point* ccw_next;

    point() {
        x = 0;
        y = 0;
        cw_next = nullptr;
        ccw_next = nullptr;
    };

    point(double _x, double _y) {
        x = _x;
        y = _y;
    }

    // constructor that takes a pointer to point
    point(point* p) {
        x = p->x;
        y = p->y;
        cw_next = nullptr;
        ccw_next = nullptr;
    }

    point subtract(const point& p) const {
        return point(x - p.x, y - p.y);
    }

    bool operator==(const point& p) const {
        return x == p.x && y == p.y;
    }

    friend std::ostream& operator<<(std::ostream& os, const point& p) {
        os << "(" << p.x << ", " << p.y << ")";
        return os;
    }
};

class ConvexHull {
public:

    // Конструктор, принимающий массив точек и его размер
    ConvexHull(point* points, int size) {
        _size = size;
        _points = new point[size];
        for (int i = 0; i < size; i++) {
            _points[i] = point(points[i].x, points[i].y); // глубокое копирование
        }
    }

    // Конструктор по умолчанию
    ConvexHull() {
        _size = 0;
        _points = nullptr;
    }

    // Метод для получения массива точек
    point* getPoints() const {
        return _points;
    }

    // Метод для получения размера массива точек
    int getSize() const {
        return _size;
    }

    void print_hull() {
        for (int i = 0; i < _size; i++) {
            printf("%f %f NEIGHBOR %f %f\n", _points[i].x, _points[i].y, *_points[i].ccw_next, *_points[i].cw_next);
        }
    }

    void print_in_file(const char* filename) {

        ofstream outfile(filename);

        for (int i = 0; i < _size; i++) {
            outfile << _points[i].x << " " << _points[i].y << endl;
        }
    }

private:
    point* _points; // Указатель на массив точек
    int _size; // Размер массива точек
};

#define N 20
//#define N0 3 // минимальное количество точек для вызова метода Грэхема

// Функция для генерации случайных точек
void generateRandomPoints(int n, point* points_array) {

    // Заполняем массив случайными координатами
    for (int i = 0; i < n; i++) {
        points_array[i].x = (rand() / 500);
        points_array[i].y = (rand() / 500);
    }
}

void set_neighbours(point points[], int n) {

    // найти индекс самой левой точки
    int left_index = 0;
    for (int i = 1; i < n; i++) {
        if (points[i].x < points[left_index].x) {
            left_index = i;
        }
    }

    // установить соседей для каждой точки
    int p = left_index, q;
    do {
        q = (p + 1) % n; // индекс следующей точки по часовой стрелке
        points[p].cw_next = &points[q];
        points[q].ccw_next = &points[p];
        p = q;
    } while (p != left_index);
}


// функция crossProduct возвращает знак псевдоскалярного произведения двух векторов, 
//  которое является положительным, если точка c находится слева от вектора ab, и отрицательным,
//   если точка c находится справа от вектора ab. 

double crossProduct(point a, point b, point c) {
    double ab_x = b.x - a.x;
    double ab_y = b.y - a.y;
    double ac_x = c.x - a.x;
    double ac_y = c.y - a.y;
    return (ab_x * ac_y) - (ab_y * ac_x);
}

// Функция для определения полярного угла между точками
double polar_angle(point p1, point p2) {
    double dx = p2.x - p1.x;
    double dy = p2.y - p1.y;
    return atan2(dy, dx);
}

// Функция для сравнения точек по полярному углу
bool cmp(point p1, point p2, point first_point) {
    double angle1 = polar_angle(first_point, p1);
    double angle2 = polar_angle(first_point, p2);
    if (angle1 == angle2) {
        // Если точки имеют одинаковый угол, то сортируем их по расстоянию до первой точки
        double dist1 = sqrt(pow(p1.x - first_point.x, 2) + pow(p1.y - first_point.y, 2));
        double dist2 = sqrt(pow(p2.x - first_point.x, 2) + pow(p2.y - first_point.y, 2));
        return dist1 < dist2;
    }
    else {
        return angle1 < angle2;
    }
}

// Функция для сортировки массива точек по полярному углу относительно первой точки
void sort_by_polar_angle(point* points, int size, point first_point) {
    std::sort(points, points + size, [first_point](point p1, point p2) { return cmp(p1, p2, first_point); });
}


ConvexHull graham_alg(point* points, int size) {

    // Find the bottom-most point
    int min_idx = 0;

    for (int i = 0; i < size; i++) {
        if (points[i].y < points[min_idx].y) {
            min_idx = i;
        }
        else if (points[i].y == points[min_idx].y && points[i].x < points[min_idx].x) {
            min_idx = i;
        }
    }

    // Swap bottom-most point with the first point
    swap(points[0], points[min_idx]);

    // Sort points by polar angle
    sort_by_polar_angle(points, size, points[0]);

    // Create a stack to store the points of the convex hull
    stack<point> hull;

    // Push the first three points onto the stack
    hull.push(points[0]);
    hull.push(points[1]);
    hull.push(points[2]);

    // Iterate over the remaining points
    for (int i = 3; i < size; i++) {
        // Pop points from the stack while the current point is to the right of the top two points
        while (hull.size() >= 2) {

            // Сохраняем верхнюю точку стека
            point top = hull.top();

            // Удаляем верхнюю точку стека
            hull.pop();

            // Берем следующую точку
            point second = hull.top();

            // Если это условие выполняется, точка отправляется обратно в стек
            if (crossProduct(second, top, points[i]) > 0) {
                hull.push(top);
                break;
            }
        }
        hull.push(points[i]);
    }

    // Create a new array to store the points of the convex hull
    int hull_size = hull.size();
    point* hull_points = new point[hull_size];

    // Pop the points from the stack and store them in the array
    for (int i = hull_size - 1; i >= 0; i--) {
        hull_points[i] = hull.top();
        hull.pop();
    }

    set_neighbours(hull_points, hull_size);

    // Create a new ConvexHull object and return it
    ConvexHull convex_hull(hull_points, hull_size);
    return convex_hull;
}


point find_internal_point(ConvexHull hull) {

    double x_sum = 0;
    double y_sum = 0;

    // Вычисляем сумму координат точек многоугольника
    for (int i = 0; i < hull.getSize(); i++) {
        x_sum += hull.getPoints()[i].x;
        y_sum += hull.getPoints()[i].y;
    }

    // Находим центроид
    point centroid(x_sum / hull.getSize(), y_sum / hull.getSize());

    return centroid;
}

bool isInside(const point p, const ConvexHull P2) {  // test +

    point* points = P2.getPoints();
    int size = P2.getSize();

    for (int i = 0; i < size; i++) {

        point a = points[i];

        // оператор % вычисляет остаток от деления, чтобы обеспечить зацикленность
        //  - когда мы дошли до последней вершины, следующая вершина будет первой вершиной, т.е. индекс 0.
        point b = points[(i + 1) % size];

        double cp = crossProduct(a, b, p);

        if (cp < 0) {
            return false;
        }
    }

    return true;
}

void quicksort_points(point* points, int low, int high) {
    if (low < high) {
        int pivot_index = partition(points, low, high);
        quicksort_points(points, low, pivot_index - 1);
        quicksort_points(points, pivot_index + 1, high);
    }
}

int partition(point* points, int low, int high) {
    point pivot = points[high];
    int i = low - 1;
    for (int j = low; j <= high - 1; j++) {
        if (compare_points(points[j], pivot)) {
            i++;
            std::swap(points[i], points[j]);
        }
    }
    std::swap(points[i + 1], points[high]);
    return i + 1;
}

bool compare_points(const point& a, const point& b) {
    if (a.y == b.y) {
        return a.x < b.x;
    }
    return a.y < b.y;
}

// объединение двух выпуклых многоугольников
//ConvexHull merge_hulls(ConvexHull P1, ConvexHull P2) {
//
//    ConvexHull merged_hull;
//
//    // получить самую правую точку из левого многоугольника
//    point p = most_left_point(P1);
//
//    // получить самую левую точку из правого многоугольника
//    point q = most_left_point(P2);
//
//    // сделать копии p и q  // КОНСТРУКТОР КОПИЙ?
//    point cp_p = p;
//    point cp_q = q;
//
//    // поднять мост pq до верхней касательной
//    point prev_p, prev_q;
//    while (true) {
//        prev_p = p;
//        prev_q = q;
//
//        // Если существует следующая точка по часовой стрелке от текущей точки q
//        if (q.cw_next != nullptr) {
//
//            // двигаться по часовой стрелке, пока происходит левый поворот
//            while (crossProduct(p, q, q.cw_next) < 0) {
//                q = q.cw_next;
//            }
//        }
//
//        // Если существует следующая точка по часовой стрелке от текущей точки p
//        if (p.cw_next != nullptr) {
//
//            // двигаться против часовой стрелки, пока происходит правый поворот
//            while (crossProduct(q, p, p.ccw_next) > 0) {
//                p = p.ccw_next;
//            }
//        }
//
//        //Цикл while продолжается, пока не будет достигнута ситуация, 
//        //когда p и q больше не изменяются после выполнения прохода цикла, 
//        //т.е. пока p и q не будут достигнуты наивысшие точки находящиеся на линии связывающей две выпуклые оболочки.
//        if (p == prev_p && q == prev_q) {
//            break;
//        }
//    }
//
//    // опустить мост cp_p cp_q до нижней касательной
//    prev_p = nullptr;
//    prev_q = nullptr;
//
//    while (true) {
//        prev_p = cp_p;
//        prev_q = cp_q;
//        if (cp_q.ccw_next != nullptr) {
//            // двигаться против часовой стрелки, пока происходит правый поворот
//            while (crossProduct(cp_p, cp_q, cp_q.ccw_next) > 0) {
//                cp_q = cp_q.ccw_next;
//            }
//        }
//        if (cp_p.cw_next != nullptr) {
//            // двигаться по часовой стрелке, пока происходит левый поворот
//            while (crossProduct(cp_q, cp_p, cp_p.cw_next) < 0) {
//                cp_p = cp_p.cw_next;
//            }
//        }
//
//        if (cp_p == prev_p && cp_q == prev_q) {
//            break;
//        }
//    }
//
//    // обновляем связи между точками, чтобы соединить две выпуклые оболочки 
//    p.cw_next = &q;
//    q.ccw_next = &p;
//    cp_p.ccw_next = &cp_q;
//    cp_q.cw_next = &cp_p;
//
//    // конечный результат
//    vector<Point> result;
//    Point *start = p;
//    while (true) {
//        result.push_back(*p);
//        p = p->ccw_next;
//
//        if (p == start) {
//            break;
//        }
//    }
//
//    return result;
//}

//ConvexHull merge_hulls(ConvexHull P1, ConvexHull P2) {
//
//    ConvexHull merged_hull;
//
//    // Найти некоторую внутреннюю точку p многоугольника P1.
//    point p = find_internal_point(P1);
//
//    // Определить, является ли точка p внутренней точкой P2.
//    if (isInside(p, P2)) {
//
//        // Получаем упорядоченный список вершин P1 и P2
//        int size1 = P1.getSize();
//        int size2 = P2.getSize();
//        point* points = new point[size1 + size2];
//
//        for (int i = 0; i < size1; i++) {
//            points[i] = P1.getPoints()[i];
//        }
//        for (int i = 0; i < size2; i++) {
//            points[size1 + i] = P2.getPoints()[i];
//        }
//        quicksort_points(points, 0, size1 + size2 - 1);
//
//        // Применяем обход Грехема к упорядоченному списку вершин
//        merged_hull = graham_alg(points, size1 + size2);
//        delete[] points;
//
//    }
//
//    // Если p не является внутренней точкой P2, применяем описанный шаг
//    else {
//
//
//
//    }
//
//    return merged_hull;
//}








int main() {


    // Инициализируем генератор случайных чисел
    srand(time(NULL));

    point* point_array1 = new point[N];

    generateRandomPoints(N, point_array1);

    point* point_array2 = new point[N];

    generateRandomPoints(N, point_array2);

    //for (int i = 0; i < convex_hull_size; i++) {
    //    printf("%f %f\n", point_array[i].x, point_array[i].y);
    //}

        // 2 выпуклые оболочки
    ConvexHull new_hull1 = graham_alg(point_array1, N);
    ConvexHull new_hull2 = graham_alg(point_array2, N);

    new_hull1.print_in_file("new_hull1.txt");

    new_hull2.print_in_file("new_hull2.txt");

    //printf("\n---------------------\n");
    //new_hull1.print_hull();
    //printf("\n---------------------\n");
    //new_hull2.print_hull();
    //printf("\n---------------------\n");

    //ConvexHull merged_hull = merge_hulls(new_hull1, new_hull1);




    //printf("\nSize: %d", new_hull.getSize());

    return 0;
}
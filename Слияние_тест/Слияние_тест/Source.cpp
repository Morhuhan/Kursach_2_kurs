// ���� �������

#include <iostream>
#include <stack>
#include <algorithm>
#include <cmath>
#include <fstream>

using namespace std;

class point {
public:

    double x, y;

    point() {
        x = 0;
        y = 0;
    };

    point(double _x, double _y) {
        x = _x;
        y = _y;
    }
};

class ConvexHull {
public:
    // �����������, ����������� ������ �����

    ConvexHull(point* points, int size) {
        _size = size;
        _points = new point[size];
        for (int i = 0; i < size; i++) {
            _points[i] = point(points[i].x, points[i].y); // �������� �����������
        }
    }

    // ����� ��� ��������� ������� �����
    point* getPoints() const {
        return _points;
    }

    // ����� ��� ��������� ������� ������� �����
    int getSize() const {
        return _size;
    }

    void print_hull() {
        for (int i = 0; i < _size; i++) {
            printf("%f %f\n", _points[i].x, _points[i].y);
        }
    }

    void print_in_file(const char* filename) {

        ofstream outfile(filename);

        for (int i = 0; i < _size; i++) {
            outfile << _points[i].x << " " << _points[i].y << endl;
        }
    }

private:
    point* _points; // ��������� �� ������ �����
    int _size; // ������ ������� �����
};

#define N 20
//#define N0 3 // ����������� ���������� ����� ��� ������ ������ �������

// ������� ��� ��������� ��������� �����
void generateRandomPoints(int n, point* points_array) {
    // �������������� ��������� ��������� �����
    srand(time(NULL));

    // ��������� ������ ���������� ������������
    for (int i = 0; i < n; i++) {
        points_array[i].x = (rand() / 100);
        points_array[i].y = (rand() / 100);
    }
}




// ������� crossProduct ���������� ���� ���������������� ������������ ���� ��������, 
//  ������� �������� �������������, ���� ����� c ��������� ����� �� ������� ab, � �������������,
//   ���� ����� c ��������� ������ �� ������� ab. 

double crossProduct(point a, point b, point c) {
    double ab_x = b.x - a.x;
    double ab_y = b.y - a.y;
    double ac_x = c.x - a.x;
    double ac_y = c.y - a.y;
    return (ab_x * ac_y) - (ab_y * ac_x);
}

// ������� ��� ����������� ��������� ���� ����� �������
double polar_angle(point p1, point p2) {
    double dx = p2.x - p1.x;
    double dy = p2.y - p1.y;
    return atan2(dy, dx);
}

// ������� ��� ��������� ����� �� ��������� ����
bool cmp(point p1, point p2, point first_point) {
    double angle1 = polar_angle(first_point, p1);
    double angle2 = polar_angle(first_point, p2);
    if (angle1 == angle2) {
        // ���� ����� ����� ���������� ����, �� ��������� �� �� ���������� �� ������ �����
        double dist1 = sqrt(pow(p1.x - first_point.x, 2) + pow(p1.y - first_point.y, 2));
        double dist2 = sqrt(pow(p2.x - first_point.x, 2) + pow(p2.y - first_point.y, 2));
        return dist1 < dist2;
    }
    else {
        return angle1 < angle2;
    }
}

// ������� ��� ���������� ������� ����� �� ��������� ���� ������������ ������ �����
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

            // ��������� ������� ����� �����
            point top = hull.top();

            // ������� ������� ����� �����
            hull.pop();

            // ����� ��������� �����
            point second = hull.top();

            // ���� ��� ������� �����������, ����� ������������ ������� � ����
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

    // Create a new ConvexHull object and return it
    ConvexHull convex_hull(hull_points, hull_size);
    return convex_hull;
}


point find_internal_point(ConvexHull hull) {

    double x_sum = 0;
    double y_sum = 0;

    // ��������� ����� ��������� ����� ��������������
    for (int i = 0; i < hull.getSize(); i++) {
        x_sum += hull.getPoints()[i].x;
        y_sum += hull.getPoints()[i].y;
    }

    // ������� ��������
    point centroid(x_sum / hull.getSize(), y_sum / hull.getSize());

    return centroid;
}

bool isInside(const point p, const ConvexHull P2) {

    point* points = P2.getPoints();
    int size = P2.getSize();

    for (int i = 0; i < size; i++) {

        point a = points[i];

        // �������� % ��������� ������� �� �������, ����� ���������� �������������
        //  - ����� �� ����� �� ��������� �������, ��������� ������� ����� ������ ��������, �.�. ������ 0.
        point b = points[(i + 1) % size];

        double cp = crossProduct(a, b, p);

        if (cp < 0) {
            return false;
        }
    }

    return true;
}


ConvexHull merge_hulls(ConvexHull P1, ConvexHull P2) {

    // ����� ��������� ���������� ����� p �������������� P1. � �������� ����� ����� ����� ����� �������� ���� ����� ������� P1. 
    point p = find_internal_point(P1);

    
    cout << p.x << " " << p.y << endl;

    // ����������, �������� �� ����� p ���������� ������ P2.
    cout << isInside(p, P2);



    return P1;


}








int main() {

    point* point_array1 = new point[N];

    generateRandomPoints(N, point_array1);

    point* point_array2 = new point[N];

    generateRandomPoints(N, point_array2);

    //for (int i = 0; i < convex_hull_size; i++) {
    //    printf("%f %f\n", point_array[i].x, point_array[i].y);
    //}

        // 2 �������� ��������
    ConvexHull new_hull1 = graham_alg(point_array1, N);
    ConvexHull new_hull2 = graham_alg(point_array2, N);

    new_hull1.print_in_file("new_hull1.txt");

    new_hull2.print_in_file("new_hull2.txt");

    //printf("\n---------------------\n");
    //new_hull1.print_hull();
    //printf("\n---------------------\n");
    //new_hull2.print_hull();
    //printf("\n---------------------\n");

    ConvexHull merged_hull = merge_hulls(new_hull1, new_hull1);

    


    //printf("\nSize: %d", new_hull.getSize());

    return 0;
}
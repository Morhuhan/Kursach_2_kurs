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

    point operator-(const point& other) const {
        return point(x - other.x, y - other.y);
    }

    point operator+(const point& other) const {
        return point(x + other.x, y + other.y);
    }

    point operator*(const double k) const {
        return point(x * k, y * k);
    }

    double operator^(const point& other) const {
        return x * other.y - y * other.x;
    }

    bool operator==(const point& other) const {
        if (x == other.x && y == other.y) return true;
        else { return false; }
    }

};

class ConvexHull {
public:

    // �����������, ����������� ������ ����� � ��� ������
    ConvexHull(point* points, int size) {
        _size = size;
        _points = new point[size];
        for (int i = 0; i < size; i++) {
            _points[i] = point(points[i].x, points[i].y); // �������� �����������
        }
    }

    // ����������� �� ���������
    ConvexHull() {
        _size = 0;
        _points = nullptr;
    }

    // ����� ��� ��������� ������� �����
    point* getPoints() const {
        return _points;
    }

    // ����� ��� ��������� ������� ������� �����
    int getSize() const {
        return _size;
    }

    // ����� ��� ���������� ������� �����
    void add_point(point p) {
        // ������� ����� ������ ����� � �������� �� 1 ������ ��������
        point* newPoints = new point[_size + 1];
        // �������� ������� ����� � ����� ������
        for (int i = 0; i < _size; i++) {
            newPoints[i] = _points[i];
        }
        // ��������� ����� ����� � ����� �������
        newPoints[_size] = p;
        // ����������� ������ ������� �������
        delete[] _points;
        // ����������� ��������� �� ����� ������
        _points = newPoints;
        // ����������� ������ �������
        _size++;
    }

    int indexOf(point p) {
        for (int i = 0; i < _size; i++) {
            if (_points[i] == p) return i;
        }
    }

    void print_hull() {
        for (int i = 0; i < _size; i++) {
            printf("%f %f NEIGHBOR %f %f\n", _points[i].x, _points[i].y);
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

    // ��������� ������ ���������� ������������
    for (int i = 0; i < n; i++) {
        points_array[i].x = (rand() / 500);
        points_array[i].y = (rand() / 500);
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

// ������� ��� ���������� ������� ����� �� ��������� ���� ������������ ������ ����� (����� ����������� � ������� ���������� ���� -> ��� ����� ������������ � ������� ������ �� ������� �������)
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

    int x_sum = 0;
    int y_sum = 0;

    // ��������� ����� ��������� ����� ��������������
    for (int i = 0; i < hull.getSize(); i++) {
        x_sum += hull.getPoints()[i].x;
        y_sum += hull.getPoints()[i].y;
    }

    // ������� ��������
    point centroid(x_sum / hull.getSize(), y_sum / hull.getSize());

    //cout <<""<< centroid.x<< endl<< "" << centroid.y << endl;

    return centroid;
}

bool isInside(const point p, const ConvexHull P2) {  // test +

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

bool compare_points(const point& a, const point& b) {
    if (a.y == b.y) {
        return a.x < b.x;
    }
    return a.y < b.y;
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

void quicksort_points(point* points, int low, int high) {
    if (low < high) {
        int pivot_index = partition(points, low, high);
        quicksort_points(points, low, pivot_index - 1);
        quicksort_points(points, pivot_index + 1, high);
    }
}


// isLeft(): tests if a point is Left|On|Right of an infinite line.
//    Input:  three points P0, P1, and P2
//    Return: >0 for P2 left of the line through P0 and P1
//            =0 for P2 on the line
//            <0 for P2 right of the line
inline float isLeft(const point& p0, const point& p1, const point& p2) {
    return (p1.x - p0.x) * (p2.y - p0.y) - (p2.x - p0.x) * (p1.y - p0.y);
}


// tangent_PointPoly(): find any polygon's exterior tangents
//    Input:  p = a 2D point (exterior to the polygon)
//            P = vertices for any 2D polygon
//    Output: rtan = pointer to the rightmost tangent point
//            ltan = pointer of leftmost tangent point
void tangent_PointPoly(const point& p, point* P, int n, point** rtan, point** ltan) {

    double  eprev, enext;
    int    i;

    *rtan = &P[0];  // initially assume P[0] = both tangents
    *ltan = &P[0];

    eprev = isLeft(P[0], P[1], p);

    for (i = 1; i <= n; ++i) {
        enext = isLeft(P[i % n], P[(i + 1) % n], p);
        if ((eprev <= 0) && (enext > 0)) {
            if (isLeft(p, P[i % n], **rtan) >= 0)
                *rtan = &P[i % n];
        }
        else if ((eprev > 0) && (enext <= 0)) {
            if (isLeft(p, P[i % n], **ltan) <= 0)
                *ltan = &P[i % n];
        }
        eprev = enext;
    }
}

// ��������: 
// �������� hull, ������� ����� ������� �� 2 ����
// ����� rtan, ltan ������� ��������� hull
// ������ �� ������� ����� chain1 � chain2, �� ������� ����� ������� hull

// � �������� ������ ��������� ������� �������� chain1 � chain2 � ��������� �� ���������������� �������. 
void split_hull(ConvexHull& hull, const point& rtan, const point& ltan, point*& chain1, point*& chain2) {
    int index_rtan = hull.indexOf(rtan);
    int index_ltan = hull.indexOf(ltan);

    // ��������� ������ ������ ����
    int chain1_size = 0;
    if (index_rtan < index_ltan) {
        chain1_size = index_ltan - index_rtan + 1;
    }
    else {
        chain1_size = hull.getSize() - index_rtan + index_ltan + 1;
    }
    chain1 = new point[chain1_size];

    // ��������� ������ ����
    int j = 0;
    for (int i = index_rtan; i <= hull.getSize() - 1; i++) {
        chain1[j++] = hull.getPoints()[i];
        if (i == index_ltan) break;
    }
    if (index_ltan < index_rtan) {
        for (int i = 0; i <= index_ltan; i++) {
            chain1[j++] = hull.getPoints()[i];
        }
    }

    // ��������� ������ ������ ����
    int chain2_size = hull.getSize() - chain1_size;
    chain2 = new point[chain2_size];

    // ��������� ������ ����
    j = 0;
    for (int i = index_ltan; i <= hull.getSize() - 1; i++) {
        chain2[j++] = hull.getPoints()[i];
        if (i == index_rtan) break;
    }
    if (index_rtan < index_ltan) {
        for (int i = 0; i <= index_rtan; i++) {
            chain2[j++] = hull.getPoints()[i];
        }
    }
}

// ��� ������� ���������, ��������� �� ����� p ����� ��� ������ �� ����������� �� ������ ����� ���� � ������ ����� ����, ��������� ������� isLeft(). ��� ���������� true, ���� ����� ��������� ����� ��� ������ �� ���� ����������� � false, ���� ����� ����� ������ �� ������ ����������� � ����� �� �������, �� ���� ��������� ������� ����.
inline bool isInsideChain(const point& p, point* chain, int size) {
    float isLeftVal = isLeft(chain[0], chain[1], p);
    for (int i = 1; i < size - 1; i++) {
        float nextIsLeftVal = isLeft(chain[i], chain[i + 1], p);
        if ((isLeftVal > 0 && nextIsLeftVal < 0) || (isLeftVal < 0 && nextIsLeftVal > 0))
            return false;
    }
    return true;
}

ConvexHull merge_hulls(ConvexHull P1, ConvexHull P2) {

    ConvexHull merged_hull;

    // ����� ��������� ���������� ����� p �������������� P1.
    point p = find_internal_point(P1);

    int size1 = P1.getSize();
    int size2 = P2.getSize();

    p.x = 100, p.y = 100;

    // ����������, �������� �� ����� p ���������� ������ P2.
    if (isInside(p, P2)) {

        // �������� ������������� ������ ������ P1 � P2
        point* points = new point[size1 + size2];

        for (int i = 0; i < size1; i++) {
            points[i] = P1.getPoints()[i];
        }
        for (int i = 0; i < size2; i++) {
            points[size1 + i] = P2.getPoints()[i];
        }
        quicksort_points(points, 0, size1 + size2 - 1);

        // ��������� ����� ������� � �������������� ������ ������
        merged_hull = graham_alg(points, size1 + size2);
        delete[] points;

    }

    // ���� p �� �������� ���������� ������ P2, ��������� ��������� ���
    else {

        // ������� �������, ����������� P2 �� ��� ���������� ����
        point* rtan;
        point* ltan;
        tangent_PointPoly(p, P2.getPoints(), P2.getSize(), &rtan, &ltan);

        point* chain1 = new point;
        point* chain2 = new point;

        // ��������� P2 �� 2 ���� ������������ ����� rtan � ltan
        split_hull(P2, *rtan, *ltan, chain1, chain2);


        //���� ������� isInsideChain ���������� true ��� ����� �� ���� �����, �� ��� ��������, ��� ��� ���� �������� �������� �� ����������� � ����� p, � �� ����� �������.
        if (isInsideChain(p, chain1, size1) != true) {
            delete chain2;
        }

        // �������� ������������� ������ ������ P1 � P2
        point* points = new point[size1 + size2];

        for (int i = 0; i < size1; i++) {
            points[i] = P1.getPoints()[i];
        }
        for (int i = 0; i < size2; i++) {
            points[size1 + i] = P2.getPoints()[i];
        }

        quicksort_points(points, 0, size1 + size2 - 1);

        // ��������� ����� ������� � �������������� ������ ������
        merged_hull = graham_alg(points, size1 + size2);
        delete[] points;
    }

    merged_hull.print_in_file("merged_hull");

    return merged_hull;
}


int main() {


    // �������������� ��������� ��������� �����
    srand(time(NULL));

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

    //point p = find_internal_point(new_hull1);

    //point p(100, 100);

    merge_hulls(new_hull1, new_hull2);





    //cout << "" << res.first << endl;
    //cout << "" << res.second << endl;

    //printf("\n---------------------\n");
    //new_hull1.print_hull();
    //printf("\n---------------------\n");
    //new_hull2.print_hull();
    //printf("\n---------------------\n");

    //ConvexHull merged_hull = merge_hulls(new_hull1, new_hull1);




    //printf("\nSize: %d", new_hull.getSize());

    return 0;
}
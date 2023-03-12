#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define N 20
#define N0 6 // ����������� ���������� ����� ��� ������ ������ �������

typedef struct {
    double x;
    double y;
} point;

// ������� ��� ��������� ��������� �����
point* generateRandomPoints(int n) {

    point* points_array = (point*)malloc(n * sizeof(point));

    srand(time(NULL)); // ��������� seed ��� ��������� ��������� �����

    // ���������� ������� ����� ���������� ������������
    for (int i = 0; i < n; i++) {
        points_array[i].x = (rand() % 21) - 10; // ��������� ��������� ���������� �� x � �������� [-10, 10]
        points_array[i].y = (rand() % 21) - 10; // ��������� ��������� ���������� �� y � �������� [-10, 10]
    }
    
    return points_array;
}

// Returns the orientation of three points (clockwise, counterclockwise, or collinear).
int orientation(point p1, point p2, point p3) {
    double val = (p2.y - p1.y) * (p3.x - p2.x) - (p2.x - p1.x) * (p3.y - p2.y);
    if (val == 0) {
        return 0; // collinear
    }
    else if (val > 0) {
        return 1; // counterclockwise
    }
    else {
        return 2; // clockwise
    }
}



void merge_hulls(point* left_hull, int left_size, point* right_hull, int right_size, point* result_hull, int* result_size) {

    // ������� ������ � ������� �����������
    int left_tangent_idx = 0;
    int right_tangent_idx = 0;
    for (int i = 0; i < left_size; i++) {
        if (left_hull[i].y < left_hull[left_tangent_idx].y) {
            left_tangent_idx = i;
        }
    }
    for (int i = 0; i < right_size; i++) {
        if (right_hull[i].y > right_hull[right_tangent_idx].y) {
            right_tangent_idx = i;
        }
    }
    int left_idx = left_tangent_idx;
    int right_idx = right_tangent_idx;
    int done = 0;

    // ��������� ����� ����������� ����������� � ��������� ����������
    while (!done) {
        result_hull[(*result_size)++] = left_hull[left_idx];
        result_hull[(*result_size)++] = right_hull[right_idx];
        int next_left_idx = (left_idx + 1) % left_size;
        int next_right_idx = (right_idx + 1) % right_size;

        int orientation_val = orientation(left_hull[left_idx], right_hull[right_idx], left_hull[next_left_idx]);
        if (orientation_val < 0) {
            left_idx = next_left_idx;
        }
        else if (orientation_val > 0) {
            right_idx = next_right_idx;
        }
        else {
            left_idx = next_left_idx;
            right_idx = next_right_idx;
        }

        // ���� �� ��������� � ��������� ��������, �� ��������� �����������
        if (left_idx == left_tangent_idx && right_idx == right_tangent_idx) {
            done = 1;
        }
    }
}


void sort_points(point* points, int n, point center) {
    int i, j;
    point temp;
    double angle1, angle2;

    for (i = 1; i < n; i++) {
        temp = points[i];
        angle1 = atan2(points[i].y - center.y, points[i].x - center.x);

        for (j = i - 1; j >= 0; j--) {
            angle2 = atan2(points[j].y - center.y, points[j].x - center.x);
            if (angle1 < angle2) {
                points[j + 1] = points[j];
            }
            else {
                break;
            }
        }
        points[j + 1] = temp;
    }
}

// �������, ����������� �������� ������� 
point* graham_alg(point* point_array, int size_of_point_array) {

    // ������� ����� � ���������� y-����������� (��� ��������� �������� ��, � ������� ���������� x-����������).
    int min_idx = 0;
    for (int i = 1; i < size_of_point_array; i++) {
        if (point_array[i].y < point_array[min_idx].y || (point_array[i].y == point_array[min_idx].y && point_array[i].x < point_array[min_idx].x)) {
            min_idx = i;
        }
    }

    // ���������� ����� � ���������� x, y-����������� � ������ �������.
    point temp = point_array[0];
    point_array[0] = point_array[min_idx];
    point_array[min_idx] = temp;

    // ��������� ����� �� ��������� ���� ������������ ������ �����.
    sort_points(point_array + 1, size_of_point_array - 1, point_array[0]);

    // ���������� ���� ��� ������������ ����� �������� ��������.
    point* stack = (point*)malloc(size_of_point_array * sizeof(point));
    int stack_size = 0;

    // �������� ������ ��� ����� � ����.
    stack[stack_size++] = point_array[0];
    stack[stack_size++] = point_array[1];

    // ���������� ���������� �����.
    for (int i = 2; i < size_of_point_array; i++) {

        // ������� ����� �� �����, ���� ������� ����� ��������� ������ �� �����, ������������ ����� �������� �������.
        while (stack_size > 1 && orientation(stack[stack_size - 2], stack[stack_size - 1], point_array[i]) <= 0) {
            stack_size--;
        }

        // ��������� ������� ����� � ����.
        stack[stack_size++] = point_array[i];
    }

    // ����� �� ����� �������� ��������� �������� ��������.
    // �������� �� � �������� ������ � ����������� ����.
    for (int i = 0; i < stack_size; i++) {
        point_array[i] = stack[i];
    }

    free(stack);

    return point_array;
}





point* find_convex_hull(point* point_array, int size_of_point_array, int* size_of_hull) {
    if (size_of_point_array < N0) {
        return graham_alg(point_array, size_of_point_array);
    }

    int mid = size_of_point_array / 2;
    point* left_hull = find_convex_hull(point_array, mid, size_of_hull);
    point* right_hull = find_convex_hull(point_array + mid, size_of_point_array - mid, size_of_hull);

    point* result_hull = (point*)malloc(size_of_point_array * sizeof(point));
    int result_size = 0;
    merge_hulls(left_hull, mid, right_hull, size_of_point_array - mid, result_hull, &result_size);

    free(left_hull);
    free(right_hull);

    *size_of_hull = result_size;
    return result_hull;
}



//// �������, ���������� ��������� �������� �������� ��������� �����
//point* find_convex_hull(point* point_array, int size_of_point_array) {
//
//    // ���������, �������� �� ������ ��������� ������ ���������� �������� N0.
//    if (size_of_point_array <= N0) {
//
//        // ���� ��, �� ������� �������� �������� ������� �������
//        point_array = graham_alg(point_array, size_of_point_array);
//        return point_array;
//    }
//
//    // ���� ������ ��������� ������ N0, �� ��������� ��� �� ��� �������� ������ �����
//    int midpoint = size_of_point_array / 2;
//
//    // ���������� ������� �������� �������� ��� ���� �����������
//    point* hull1 = find_convex_hull(point_array, midpoint);
//    point* hull2 = find_convex_hull(point_array + midpoint, size_of_point_array - midpoint);
//
//    // ������� ��� �������� �������� ������� "����������� ������"
//    return merge_hulls(hull1, midpoint, hull2, size_of_point_array - midpoint);
//}

int main() {

    point* points_array = generateRandomPoints(N);
    int convex_hull_size;
    point* convex_hull = find_convex_hull(points_array, N, &convex_hull_size);

    printf("Convex hull:\n");

    for (int i = 0; i < convex_hull_size; i++) {
        printf("(%lf, %lf)\n", convex_hull[i].x, convex_hull[i].y);
    }

    free(points_array);
    free(convex_hull);

    return 0;
}
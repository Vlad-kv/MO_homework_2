#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <cmath>
#include <map>
#include <functional>
#include <cassert>
#include <Windows.h>
#include <iomanip>
using namespace std;

template <typename T>
T sqr(T val) {
    return val * val;
}

typedef long double ld;
typedef function<ld (const vector<ld >& )> multi_var_func;


vector<vector<ld>> points;

const ld  F = 1.61803398874989484;

long long iterations_count = 0;

ostream& operator << (ostream& out, const vector<ld> &v) {
    for (auto val : v) {
        out << setprecision(10) << val << " ";
    }
    return out;
}

void print_m(const vector<vector<ld>>& m) {
    for (auto row : m) {
        for (ld val : row) {
            cout << val << " ";
        }
        cout << "\n";
    }
}

ld dichotomy_method(function<ld (ld)> f, ld lim_1, ld lim_2, ld eps) {
    assert(lim_1 <= lim_2);

    while (abs(lim_2 - lim_1) > eps) {
        ld p1 = (lim_2 + lim_1) / 2 - eps / 4;
        ld p2 = (lim_2 + lim_1) / 2 + eps / 4;

        ld p1_val = f(p1);
        ld p2_val = f(p2);

        if (p1_val < p2_val) {
            lim_2 = p2;
        } else {
            lim_1 = p1;
        }
        iterations_count++;
    }
    return (lim_1 + lim_2) / 2;
}

ld golden_ratio_method(function<ld (ld)> f, ld lim_1, ld lim_2, ld eps) {
    assert(lim_1 <= lim_2);

    ld p1_val = INFINITE;
    ld p2_val = INFINITE;

    while (abs(lim_2 - lim_1) >= eps) {
//        cout << "abs(lim_2 - lim_1) : " << abs(lim_2 - lim_1) << "\n";

        ld p1 = lim_2 - (lim_2 - lim_1) / F;
        ld p2 = lim_1 + (lim_2 - lim_1) / F;

        if (p1_val == INFINITE) {
            p1_val = f(p1);
        }
        if (p2_val == INFINITE) {
            p2_val = f(p2);
        }

        if (p1_val < p2_val) {
            lim_2 = p2;
            p2_val = p1_val;
            p1_val = INFINITE;
        } else {
            lim_1 = p1;
            p1_val = p2_val;
            p2_val = INFINITE;
        }
        iterations_count++;
    }
    return (lim_1 + lim_2) / 2;
}

ld fibonacci_method(function<ld (ld)> f, ld lim_1, ld lim_2, ld eps) {
    assert(lim_1 <= lim_2);

    ld f1 = 0, f2 = 1;
    int n = 1;

    while (f2 <= (lim_2 - lim_1) / eps) {
        ld f3 = f1 + f2;
        f1 = f2;
        f2 = f3;
        n++;
    }

    ld p1 = lim_1 + (f2 - f1) / f2 * (lim_2 - lim_1);
    ld p2 = lim_1 +        f1 / f2 * (lim_2 - lim_1);

    ld p1_val = f(p1);
    ld p2_val = f(p2);

    for (int k = 1; k <= n - 2; k++) {
        ld f0 = f2 - f1;
        f2 = f1;
        f1 = f0;

        if (p1_val > p2_val) {
            lim_1 = p1;
            p1 = p2;
            p2 = lim_1 + f1 / f2 * (lim_2 - lim_1);
            p1_val = p2_val;
            p2_val = f(p2);
        } else {
            lim_2 = p2;
            p2 = p1;
            p1 = lim_1 + (f2 - f1) / f2 * (lim_2 - lim_1);
            p2_val = p1_val;
            p1_val = f(p1);
        }
        iterations_count++;
    }

    p2 = p1 + eps / 4;
    p2_val = f(p2);

    if (p1_val < p2_val) {
        lim_2 = p2;
    } else {
        lim_1 = p1;
    }

    return (lim_1 + lim_2) / 2;
}

vector<vector<ld>> invert_matrix(const vector<vector<ld>>& m) {
    vector<vector<ld>> m1 = m;
    vector<vector<ld>> m2(m.size(), vector<ld>(m.size(), 0.0));
    for (int i = 0; i < m.size(); i++) {
        m2[i][i] = 1;
    }

    for (int i = 0; i < m.size(); i++) {
//        cout << "i : " << i << "\n";
        int max_row = i;
        for (int j = i; j < m.size(); j++) {
            if (abs(m1[max_row][i]) < abs(m1[j][i])) {
                max_row = j;
            }
        }
        swap(m1[i], m1[max_row]);
        swap(m2[i], m2[max_row]);

        ld c = m1[i][i];
        assert(c != 0);

        for (int j = 0; j < m.size(); j++) {
            m1[i][j] /= c;
            m2[i][j] /= c;
        }

//        print_m(m1);
//        cout << "\n";
//        print_m(m2);
//        cout << "\n\n\n";

        for (int j = 0; j < m.size(); j++) {
            if (j != i) {
                c = m1[j][i];
                for (int k = 0; k < m.size(); k++) {
                    m1[j][k] -= m1[i][k] * c;
                    m2[j][k] -= m2[i][k] * c;
                }
            }
        }
//        print_m(m1);
//        cout << "\n";
//        print_m(m2);
//        cout << "-----\n";
    }
    return m2;
}

vector<ld> newton_method(vector<ld> p,
                             multi_var_func f,
                             vector<multi_var_func> f_grad,
                             vector<vector<multi_var_func>> hessian,
                             ld eps) {
    points.push_back(p);
    int n_step = 10000;
    while (true) {
        n_step--;
        if (n_step == 0) {
            break;
        }
        vector<vector<ld>> H(p.size(), vector<ld>(p.size(), 0.0));
        for (int i = 0; i < p.size(); i++) {
            for (int j = 0; j < p.size(); j++) {
                H[i][j] = hessian[i][j](p);
            }
        }
        H = invert_matrix(H);
        vector<ld> grad;
        for (auto f : f_grad) {
            grad.push_back(f(p));
        }
        ld grad_norm = 0;
        for (int i = 0; i < p.size(); i++) {
            for (int j = 0; j < p.size(); j++) {
                p[i] -= H[i][j] * grad[j];
            }
            grad_norm += abs(grad[i]);
        }
        cout << p << "  " << f(p) << "  \n";
        points.push_back(p);

        if (grad_norm < eps) {
            break;
        }
        {
            ld p_norm = 0;
            for (ld val : p) {
                p_norm = max(p_norm, abs(val));
            }
            if (p_norm > 100000) {
                cout << 10000 - n_step << "\n";
                return p;
            }
        }
    }
    cout << 10000 - n_step << "\n";
    return p;
}

vector<ld> conjugate_gradient_method(
    vector<ld> p,
    multi_var_func f,
    vector<multi_var_func> f_grad,
    function<ld ((function<ld (ld)> f, ld lim_1, ld lim_2, ld eps))> one_dimensional_search_method,
    ld eps) {

        points.push_back(p);
        int n_step = 10000;
        while (n_step > 0) {
            vector<ld> dir;
            for (auto f_g : f_grad) {
                dir.push_back(-f_g(p));
            }

            for (int k = 0; k < p.size(); k++) {
//                cout << "k : " << k << "\n";
                n_step--;
                auto one_dim_task = [&dir, &p, &f] (ld c) {
//                    cout << "in one_dim_task : " << c << "\n";

                    vector<ld> p2;
                    for (int i = 0; i < dir.size(); i++) {
                        p2.push_back(p[i] + dir[i] * c);
                    }
                    return f(p2);
                };
                ld c = one_dimensional_search_method(one_dim_task, -10, 10, eps);

//                cout << "?\n";

                vector<ld> p2;
                for (int i = 0; i < p.size(); i++) {
                    p2.push_back(p[i] + dir[i] * c);
                }

                vector<ld> g1;
                for (auto f_g : f_grad) {
                    g1.push_back(f_g(p));
                }
                vector<ld> g2;
                for (auto f_g : f_grad) {
                    g2.push_back(f_g(p2));
                }

                ld w1 = 0, w2 = 0;

                for (ld val : g1) {
                    w1 += sqr(val);
                }
                for (ld val : g2) {
                    w2 += sqr(val);
                }
                ld w = w2 / w1;

                ld dir_norm = 0;

                for (int i = 0; i < dir.size(); i++) {
                    dir[i] = dir[i] * w - g2[i];
                    dir_norm += abs(dir[i]);
                }
                swap(p, p2);

                cout << p << "  " << f(p) << "\n";
                points.push_back(p);

                if (dir_norm < eps) {
                    cout << 10000 - n_step << "\n";
                    return p;
                }
                {
                    ld p_norm = 0;
                    for (ld val : p) {
                        p_norm = max(p_norm, abs(val));
                    }
                    if (p_norm > 100000) {
                        cout << 10000 - n_step << "\n";
                        return p;
                    }
                }
            }
        }
    cout << 10000 - n_step << "\n";
    return p;
}

vector<ld> simple_gradient_method(
    vector<ld> p,
    multi_var_func f,
    vector<multi_var_func> f_grad,
    function<ld ((function<ld (ld)> f, ld lim_1, ld lim_2, ld eps))> one_dimensional_search_method,
    ld eps) {
        points.push_back(p);
        int n_step = 10000;
        while (n_step > 0) {
            n_step--;
            vector<ld> dir;
            for (auto f_g : f_grad) {
                dir.push_back(-f_g(p));
            }

            auto one_dim_task = [&dir, &p, &f] (ld c) {
                vector<ld> p2;
                for (int i = 0; i < dir.size(); i++) {
                    p2.push_back(p[i] + dir[i] * c);
                }
                ld val = f(p2);

//                cout << p2 << "   " << c << "    " << val << "\n";

                return val;
            };
            ld c = one_dimensional_search_method(one_dim_task, -10, 10, eps);

            vector<ld> p2;
            ld grad_norm = 0;
            for (int i = 0; i < p.size(); i++) {
                p2.push_back(p[i] + dir[i] * c);
                grad_norm += abs(dir[i]);
            }
            swap(p, p2);

//            cout << "grad_norm : " << grad_norm << "\n";

            if (grad_norm < eps) {
                cout << 10000 - n_step << "\n";
                return p;
            }
            {
                ld p_norm = 0;
                for (ld val : p) {
                    p_norm = max(p_norm, abs(val));
                }
                if (p_norm > 100000) {
                    cout << 10000 - n_step << "\n";
                    return p;
                }
            }
            cout << p << "  " << f(p) << "\n";
            points.push_back(p);
        }
    cout << 10000 - n_step << "\n";
    return p;
}

long long calc_counter = 0;

void test_func(
    multi_var_func f,
    const vector<multi_var_func>& f_grad,
    const vector<vector<multi_var_func>>& f_hessian,
    vector<ld> p,
    ld eps) {

    points.clear();

    vector<ld> res_1_1 = simple_gradient_method(p, f, f_grad, golden_ratio_method, eps);
    cout << res_1_1 << "\n";

    cout << "(" << points[0][0];
    for (int i = 1; i < points.size(); i++) {
        cout << ", " << points[i][0];
    }
    cout << ")\n";

    cout << "(" << points[0][1];
    for (int i = 1; i < points.size(); i++) {
        cout << ", " << points[i][1];
    }
    cout << ")\n";

    cout << "(" << f({points[0]}) + abs(f({points[0]})) * 0.01;
    for (int i = 1; i < points.size(); i++) {
        cout << ", " << f({points[i]}) + abs(f({points[i]})) * 0.01;
    }
    cout << ")\n";

    cout << "--------\n";


    points.clear();

    vector<ld> res_1_2 = conjugate_gradient_method(p, f, f_grad, golden_ratio_method, eps);

    cout << "?\n";

    cout << res_1_2 << "\n";

        cout << "(" << points[0][0];
    for (int i = 1; i < points.size(); i++) {
        cout << ", " << points[i][0];
    }
    cout << ")\n";

    cout << "(" << points[0][1];
    for (int i = 1; i < points.size(); i++) {
        cout << ", " << points[i][1];
    }
    cout << ")\n";

    cout << "(" << f({points[0]}) + abs(f({points[0]})) * 0.01;
    for (int i = 1; i < points.size(); i++) {
        cout << ", " << f({points[i]}) + abs(f({points[i]})) * 0.01;
    }
    cout << ")\n";

    cout << "--------\n";

    points.clear();

    vector<ld> res_1_3 = newton_method(p, f, f_grad, f_hessian, eps);
    cout << res_1_3 << "\n";

    cout << "(" << points[0][0];
    for (int i = 1; i < points.size(); i++) {
        cout << ", " << points[i][0];
    }
    cout << ")\n";

    cout << "(" << points[0][1];
    for (int i = 1; i < points.size(); i++) {
        cout << ", " << points[i][1];
    }
    cout << ")\n";

    cout << "(" << f({points[0]}) + abs(f({points[0]})) * 0.01;
    for (int i = 1; i < points.size(); i++) {
        cout << ", " << f({points[i]}) + abs(f({points[i]})) * 0.01;
    }
    cout << ")\n";

    cout << "--------\n";
}

void test_multidim_methods() {
    multi_var_func f1 = [] (const vector<ld>& v) {
        return 100 * sqr(v[1] - v[0]) + sqr(1 - v[0]);
    };
    vector<multi_var_func> f1_grad = {
        [] (const vector<ld>& v) {
            return 202 * v[0] - 200 * v[1] - 2;
        },
        [] (const vector<ld>& v) {
            return 200 * (v[1] - v[0]);
        },
    };
    vector<vector<multi_var_func>> f1_hessian = {
        {
            [] (const vector<ld>& v) {
                return 202;
            },
            [] (const vector<ld>& v) {
                return -200;
            }
        },
        {
            [] (const vector<ld>& v) {
                return -200;
            },
            [] (const vector<ld>& v) {
                return 200;
            }
        }
    };

    multi_var_func f2 = [](const vector<ld>& v) {
        return 100 * sqr(v[1] - sqr(v[0])) + sqr(1 - v[0]);
    };
    vector<multi_var_func> f2_grad = {
        [] (const vector<ld>& v) {
            return 2 * (200 * v[0] * v[0] * v[0] - 200 * v[0] * v[1] + v[0] - 1);
        },
        [] (const vector<ld>& v) {
            return 200 * (v[1] - sqr(v[0]));
        },
    };
    vector<vector<multi_var_func>> f2_hessian = {
        {
            [] (const vector<ld>& v) {
                return 1200 * sqr(v[0]) - 400 * v[1] + 2;
            },
            [] (const vector<ld>& v) {
                return -400 * v[0];
            }
        },
        {
            [] (const vector<ld>& v) {
                return -400 * v[0];
            },
            [] (const vector<ld>& v) {
                return 200;
            }
        }
    };

    multi_var_func f3 = [](const vector<ld>& v) {
        return -2 * exp(-sqr((v[0] - 1) / 2) - sqr((v[1] - 1) / 1))
              - 3 * exp(-sqr((v[0] - 2) / 3) - sqr((v[1] - 3) / 2));
    };
    vector<multi_var_func> f3_grad = {
        [] (const vector<ld>& v) {
            return 2 / 3.0 * (v[0] - 2) * exp(-sqr((v[0] - 2) / 3) - sqr((v[1] - 3) / 2))
                 + (v[0] - 2)           * exp(-sqr((v[0] - 1) / 2) - sqr((v[1] - 1) / 1));
        },
        [] (const vector<ld>& v) {
            return 3 / 2.0 * (v[1] - 3) * exp(-sqr((v[0] - 2) / 3) - sqr((v[1] - 3) / 2))
                 + 4 * (v[1] - 1)       * exp(-sqr((v[0] - 1) / 2) - sqr((v[1] - 1) / 1));
        },
    };
    vector<vector<multi_var_func>> f3_hessian = {
        {
            [] (const vector<ld>& v) {
                return exp(-sqr((v[0] - 2) / 3) - sqr((v[1] - 3) / 2)) * (-4 / 27.0 * sqr(v[0] - 2) + 2 / 3.0)
                     + exp(-sqr((v[0] - 1) / 2) - sqr((v[1] - 1) / 1)) * (1 - 1 / 2.0 * sqr(v[0] - 2));
            },
            [] (const vector<ld>& v) {
                return exp(-sqr((v[0] - 2) / 3) - sqr((v[1] - 3) / 2)) * (-1 / 3.0 * (v[0] - 2) * (v[1] - 3))
                     + exp(-sqr((v[0] - 1) / 2) - sqr((v[1] - 1) / 1)) * (-2 * (v[0] - 2) * (v[1] - 1));
            }
        },
        {
            [] (const vector<ld>& v) {
                return exp(-sqr((v[0] - 2) / 3) - sqr((v[1] - 3) / 2)) * (-1 / 3.0 * (v[0] - 2) * (v[1] - 3))
                     + exp(-sqr((v[0] - 1) / 2) - sqr((v[1] - 1) / 1)) * (-2 * (v[0] - 2) * (v[1] - 1));
            },
            [] (const vector<ld>& v) {
                return exp(-sqr((v[0] - 2) / 3) - sqr((v[1] - 3) / 2)) * (-3 / 4.0 * sqr(v[1] - 3) + 3 / 2.0)
                     + exp(-sqr((v[0] - 1) / 2) - sqr((v[1] - 1) / 1)) * (4 - 8 * sqr(v[1] - 1));
            }
        }
    };

    test_func(f1, f1_grad, f1_hessian, {2, 3}, 1e-11);
//    test_func(f1, f1_grad, f1_hessian, {-1, -1}, 1e-9);
//    test_func(f1, f1_grad, f1_hessian, {0, 1}, 1e-9);
//    test_func(f1, f1_grad, f1_hessian, {10, 10}, 1e-9);
//    test_func(f1, f1_grad, f1_hessian, {-100, 100}, 1e-9);
//
//
//    test_func(f2, f2_grad, f2_hessian, {0, 0}, 1e-9);
//    test_func(f2, f2_grad, f2_hessian, {2, 2}, 1e-9);
//    test_func(f2, f2_grad, f2_hessian, {1.1, 1.1}, 1e-9);
//    test_func(f2, f2_grad, f2_hessian, {10, -10}, 1e-9);
//    test_func(f2, f2_grad, f2_hessian, {-100, 100}, 1e-9);
//
//    test_func(f3, f3_grad, f3_hessian, {1, 1}, 1e-6);
//    test_func(f3, f3_grad, f3_hessian, {0, 0}, 1e-9);
//    test_func(f3, f3_grad, f3_hessian, {2, 2}, 1e-9);
//    test_func(f3, f3_grad, f3_hessian, {-2, 4}, 1e-9);
//    test_func(f3, f3_grad, f3_hessian, {4, 2}, 1e-9);
//    test_func(f3, f3_grad, f3_hessian, {2, 1.4}, 1e-9);
}


int main() {
    #ifdef Vlad_kv
        freopen("input.txt", "r", stdin);
        freopen("output.txt", "w", stdout);
    #endif // Vlad_kv

    test_multidim_methods();
    return 0;
}

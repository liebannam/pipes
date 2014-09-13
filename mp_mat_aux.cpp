// this file is included repeatedly at end of mp_mat.cpp

template class mp_mat<T>;
template void dump_vec(const vector<T> &x, const char *name);
template T operator*(const vector<T> &u, const vector<T> &v);
template vector<T> operator*(const vector<T> &u, T a);
template vector<T> operator-(const vector<T> &u);
template vector<T>
    operator+(const vector<T> &u, const vector<T> &v);
template vector<T>
    operator-(const vector<T> &u, const vector<T> &v);
template mp_mat<T>
    operator*(const mp_mat<T> &A, const mp_mat<T> &B);
template mp_mat<T>
    operator+(const mp_mat<T> &A, const mp_mat<T> &B);
template mp_mat<T>
    operator-(const mp_mat<T> &A, const mp_mat<T> &B);
template mp_mat<T>
    operator*(const mp_mat<T> &A, T a);
template mp_mat<T>
    operator*(T a, const mp_mat<T> &A);
template mp_mat<T>
    trans(const mp_mat<T> &A);
template vector<T>
    operator*(const mp_mat<T> &A, const vector<T> &x);
template vector<T>
    operator*(const vector<T> &x, const mp_mat<T> &A);
template void transpose(const mp_mat<T>& A, mp_mat<T>& B);
template void
sandwichL(char trans, mp_mat<T>& X, mp_mat<T>& K, mp_mat<T>& W);
template T compute_norm(vector<T> &x);
template void scale_vec(vector<T> &x, const vector<T> &s);
template void scale_mat(mp_mat<T> &A, const vector<T> &s);
template void scale_mat(const vector<T> &s, mp_mat<T> &A);

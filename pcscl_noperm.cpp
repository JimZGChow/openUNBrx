#include "pcscl_noperm.h"

using namespace std;

#define MINSUM


float phi(float x)
{
    static const float lim = 31.;
    static const float inv_lim = log((exp((double)lim) + 1) / (exp((double)lim) - 1));

    if (x > lim)
        return(0);
    else if (x < inv_lim)
        return(lim);
    else
    {
        double t = exp((double)x);
        return (float)log((t + 1) / (t - 1));
    }
}

inline float sign(float x)
{
    if (x < 0)
    {
        return -1;
    }
    return 1;
}


inline float jacoblog(float x)
{
    if (x < -17)
    {
        return 0;
    }
    if (x > 15)
    {
        return x;
    }
    else
    {
        return log(1 + exp(x));
    }
}

inline float cnop(float a, float b)
{
    // return maxS(a,b) - jacoblog(a+b);
#ifdef MINSUM
    return sign(a) * sign(b) * min(abs(a), abs(b));
#else
    return sign(a) * sign(b) * phi(phi(abs(l1)) + phi(abs(l2)));
#endif
}

inline float cnop(bool a, float b)
{
    if (a)
    {
        return -b;
    }
    else
    {
        return b;
    }
}

inline float vnop(float a, float b)
{
    return a + b;
}

vector<vector<float>> prob;
vector<vector<uint8_t>> u;
unsigned idx;
unsigned L;
vector<bool> used;
vector<tuple<float, int, bool> > sorted;


vector<pcscl_list> pcscl(const vector<vector<float> >& y, vector<int>::const_iterator f_it) {

    const unsigned N = y[0].size();
    const unsigned L0 = y.size();

    if (N == 1)
    { // End recursion
        vector<pcscl_list> res;
        if (*f_it == 2)
        { // Information bit
            sorted.clear();
            sorted.reserve(2 * L0);
            for (unsigned i = 0; i < L0; i++)
            {
                float logy = jacoblog(y[i][0]);
                sorted.emplace_back(prob[i][idx] + y[i][0] - logy, i, false);
                sorted.emplace_back(prob[i][idx] - logy, i, true);
            }
            sort(sorted.begin(), sorted.end(), greater<decltype(sorted[0])>());
            if (L < 2 * L0)
            {
                sorted.resize(L);
            }
            res.resize(sorted.size());

            vector<int8_t> last(L0, -1);
            for (unsigned i = 0; i < sorted.size(); i++)
            {
                last[get<1>(sorted[i])] = i; // Find the last use
            }
            for (unsigned i = 0; i < L0; i++)
            {
                if (last[i] < 0) used[i] = false;
            }
            for (unsigned i = 0; i < sorted.size(); i++)
            {
                if (last[get<1>(sorted[i])] == (int)i)
                {
                    unsigned j = get<1>(sorted[i]);
                    u[j][idx] = get<2>(sorted[i]);
                    prob[j][idx + 1] = get<0>(sorted[i]);
                    res[j].x.push_back(get<2>(sorted[i]));
                    res[j].idx = j;
                }
                else if (last[get<1>(sorted[i])] >= 0)
                {
                    unsigned j = 0;
                    while (j < L && used[j]) j++; // Find empty space
                    copy_n(u[get<1>(sorted[i])].begin(), idx, u[j].begin()); // CLONE PATH
                    u[j][idx] = get<2>(sorted[i]);
                    prob[j][idx + 1] = get<0>(sorted[i]);
                    res[j].x.push_back(get<2>(sorted[i]));
                    res[j].idx = get<1>(sorted[i]);
                    used[j] = true;
                }
            }
        }
        else
        { // Frozen
            res.resize(L0);
            for (size_t l = 0; l < L0; l++)
            {
                u[l][idx] = *f_it;
                prob[l][idx + 1] = prob[l][idx] - jacoblog(y[l][0]);
                if (*f_it == 0)
                {
                    prob[l][idx + 1] += y[l][0];
                }
                res[l].x.push_back(*f_it);
                res[l].idx = l;
            }
        }

        ++idx;
        return res;
    }
    else
    {   // Recursion
        assert(N % 2 == 0);
        vector<vector<float>> u_est(y.size(), vector<float>(N / 2));
        for (unsigned i = 0; i < y.size(); i++)
        {
            for (unsigned j = 0; j < N / 2; j++)
            {
                u_est[i][j] = cnop(y[i][j], y[i][j + N / 2]);
            }
        }


        vector<pcscl_list> res1 = pcscl(u_est, f_it);

        u_est.resize(res1.size(), vector<float>(N / 2));
        for (size_t i = 0; i < res1.size(); i++)
        {
            const vector<bool>& u1hardprev = res1[i].x;
            const vector<float>& yy = y[res1[i].idx];
            for (unsigned j = 0; j < N / 2; j++)
            {
                u_est[i][j] = vnop(cnop(u1hardprev[j], yy[j]), yy[j + N / 2]);
            }
        }
        vector<pcscl_list> res2 = pcscl(u_est, f_it + N / 2);

        for (size_t i = 0; i < res2.size(); i++)
        {
            vector<bool> x(N);
            const vector<bool>& x1 = res1[res2[i].idx].x;
            const vector<bool>& x2 = res2[i].x;
            for (unsigned j = 0; j < N / 2; j++)
            {
                x[j] = x1[j] ^ x2[j];
                x[j + N / 2] = x2[j];
            }
            res2[i].x = move(x);
            res2[i].idx = res1[res2[i].idx].idx;
        }
        return res2;
    }
}

vector<vector<uint8_t>> pcscl_prep(int m, int _L, const std::vector<float>& data, vector<vector<float>>* _prob, uint8_t* info_bit_pattern) {
    /* Read arguments */
    int N = 1 << m;
    L = _L;

    vector<int> f_matrix(N);
    for (int i = 0; i < N; ++i) {
        f_matrix[i] = info_bit_pattern[i] * 2;
    }


    // LLR
    vector<float> in_llr(N);
    for (int i = 0; i < N; ++i) {
        in_llr[i] = data[i];
    }

    /* Initialization */
    used.assign(L, false);
    used[0] = true;
    u.assign(L, vector<uint8_t>(N));
    prob.assign(L, vector<float>(N + 1));
    idx = 0;

    /* Decode */
    //vector<vector<float>> tmp(1, vector<float>());
    //tmp[0] = in_llr;
    vector<pcscl_list> li = pcscl({ in_llr }, f_matrix.begin());
    //vector<pcscl_list> li = pcscl(tmp, f_matrix.begin());

    /* Return results to MATLAB */
    // iwd list

    vector<vector<uint8_t>> ret = u;
    if (_prob) {
        *_prob = prob;
    }

    prob.clear();
    u.clear();
    used.clear();
    sorted.clear();

    return ret;

}

vector<vector<uint8_t>> pcscl_prep(int m, int _L,const std::vector<uint8_t>& data, vector<vector<float>>* _prob, uint8_t* info_bit_pattern) {
    /* Read arguments */
    int N = 1 << m;
    L = _L;

    vector<int> f_matrix(N);
    for (int i = 0; i < N; ++i) {
        f_matrix[i] = info_bit_pattern[i] * 2;
    }


    // LLR
    vector<float> in_llr(N);
    for (int i = 0; i < N; ++i) {
        in_llr[i] = data[i] ? -1 : 1;
    }

    /* Initialization */
    used.assign(L, false);
    used[0] = true;
    u.assign(L, vector<uint8_t>(N));
    prob.assign(L, vector<float>(N + 1));
    idx = 0;

    /* Decode */
    vector<pcscl_list> li = pcscl({ in_llr }, f_matrix.begin());

    /* Return results to MATLAB */
    // iwd list

    vector<vector<uint8_t>> ret = u;
    if (_prob) {
        *_prob = prob;
    }

    prob.clear();
    u.clear();
    used.clear();
    sorted.clear();

    return ret;

}

uint16_t crc(uint32_t polynom, std::vector<uint8_t> a) {
    uint32_t reg = 0;
    for (auto b : a) {
        reg ^= b;
        if (reg & 1) reg = (reg >> 1) ^ polynom;
        else reg >>= 1;
    }
    return reg & 0x3FF;
}

uint8_t crc_ok(uint32_t polynom, std::vector<uint8_t> a) {
    uint32_t reg = 0;
    for (auto b : a) {
        reg ^= b;
        if (reg & 1) reg = (reg >> 1) ^ polynom;
        else reg >>= 1;
    }
    return reg == 0;
}

vector<uint8_t> crc_ok_array(uint32_t polynom, vector<vector<uint8_t>> a)
{
    vector<uint8_t> res(a.size());
    for (unsigned i = 0; i < a.size(); i++) {
        res[i] = crc_ok(polynom, a[i]);
    }
    return res;
}

vector<uint8_t> remove_crc(vector<uint8_t> a)
{
    vector<uint8_t> res = a;
    res.erase(res.end() - 10, res.end());
    return res;
}

#ifdef MEX
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    /* Read arguments */

    int m = (int)*(mxGetPr(prhs[0])); // Length
    int N = 1 << m;
    L = (int)*(mxGetPr(prhs[1])); // List size

    // f matrix
    double* input = mxGetPr(prhs[2]);

    vector<int> f_matrix(N);
    for (int i = 0; i < N; ++i)
    {
        f_matrix[i] = (int)input[i];
    }


    // LLR
    input = mxGetPr(prhs[3]);

    vector<float> in_llr(N);
    for (int i = 0; i < N; ++i)
    {
        in_llr[i] = input[i];
    }

    /* Initialization */
    used.assign(L, false);
    used[0] = true;
    u.assign(L, vector<bool>(N));
    prob.assign(L, vector<float>(N + 1));
    idx = 0;

    /* Decode */
    pcscl({ in_llr }, f_matrix.begin());

    /* Return results to MATLAB */
    // iwd list
    plhs[0] = mxCreateDoubleMatrix(L, N, mxREAL);
    double* output = mxGetPr(plhs[0]);
    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < L; ++j)
        {
            output[i * L + j] = u[j][i];
        }
    }


    // probabilities
    plhs[1] = mxCreateDoubleMatrix(L, N + 1, mxREAL);
    output = mxGetPr(plhs[1]);
    for (int i = 0; i < N + 1; ++i)
    {
        for (int j = 0; j < L; ++j)
        {
            output[i * L + j] = prob[j][i];
        }
    }
}
#endif

std::vector<std::vector<uint8_t>> polar_transform_noperm(std::vector<std::vector<uint8_t>> data) {
    if (data[0].size() == 1) {
        return data;
    }
    else {
        int N = data[0].size();
        std::vector<std::vector<uint8_t>> u1u2(data.size(), std::vector<uint8_t>(N/2));

        for (int i = 0; i < data.size(); i++) {
            for (int j = 0; j < data[0].size()/2; j++) {
                u1u2[i][j] = data[i][j] ^ data[i][j + N / 2];
            }
        }

        std::vector<std::vector<uint8_t>> u2(data.size(), std::vector<uint8_t>(N / 2));

        for (int i = 0; i < data.size(); i++) {
            for (int j = 0; j < data[0].size()/2; j++) {
                u2[i][j] = data[i][j + N / 2];
            }
        }

        std::vector<std::vector<uint8_t>> x(data.size(), std::vector<uint8_t>(N));
        std::vector<std::vector<uint8_t>> ret1 = polar_transform_noperm(u1u2);
        std::vector<std::vector<uint8_t>> ret2 = polar_transform_noperm(u2);

        for (int i = 0; i < x.size(); i++) {
            for (int j = 0; j < x[0].size(); j++) {
                if (j < x[0].size() / 2) {
                    x[i][j] = ret1[i][j];
                }
                else {
                    x[i][j] = ret2[i][j - x[0].size() / 2];
                }
            }
        }

        return x;
    }

}

std::vector<std::vector<uint8_t>> extract_with_filter(const std::vector<std::vector<uint8_t>>& data, const uint16_t* filtr, int K) {
    std::vector<std::vector<uint8_t>> ret(16, std::vector<uint8_t>(K));
    for (int i = 0; i < K; i++) {
        for (int x = 0; x < data.size(); x++) {

            ret[x][i] = data[x][filtr[i]];
        }
    }

    return ret;
}

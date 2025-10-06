#include <vector>
#include <iostream>
#include <string>
#include <algorithm>
#include <cmath>
#include <map>
#include <numeric>

#pragma once

// Simplifies sqrt(i) as much as possible into a * sqrt(b), returns {a, b}
std::pair<int, int> simplify_sqrt(int i)
{
    int out_sqrt = 1;
    int in_sqrt = i;

    int d = 2;
    while (d * d <= in_sqrt)
    {
        if (in_sqrt % (d * d) == 0)
        {
            in_sqrt /= d * d;
            out_sqrt *= d;
        }
        else
        {
            d++;
        }
    }
    return { out_sqrt, in_sqrt };
}

//A (constant) expression, just a integer or fraction right now.
struct expr
{
    std::map<int, int[2]> values; //values[n] gives the numerator & denominator of the part of the expression multiplied by sqrt(n).

    void simplify()
    {
        std::pair<int, int> s;
        std::vector<int> indicesToErase;
        for (auto const& p : values)
        {
            s = simplify_sqrt(p.first);
            if (p.first != s.second)
            {
                int val[2];
                val[0] = values[s.second][0];
                val[1] = values[s.second][1];

                if (val[0] == 0)
                {
                    values[s.second][0] = p.second[0] * s.first;
                    values[s.second][1] = p.second[1];
                }
                else
                {
                    values[s.second][0] = val[0] * p.second[1] + val[1] * p.second[0] * s.first;
                    values[s.second][1] = val[1] * p.second[1];
                }

                indicesToErase.push_back(p.first); //don't change the size of a map while iterating through it
            }
        }

        for (int i : indicesToErase)
        {
            values.erase(i);
        }

        for (auto const& p : values)
        {
            int gcd = std::__gcd(p.second[0], p.second[1]);

            if (gcd != 0)
            {
                values[p.first][0] /= gcd;
                values[p.first][1] /= gcd;

                if (p.second[1] < 0) //always want negative in numerator
                {
                    values[p.first][0] *= -1;
                    values[p.first][1] *= -1;
                }
            }
        }
    }

    expr() {}

    expr(int i)
    {
        if (i != 0)
        {
            values[1][0] = i;
            values[1][1] = 1; // i = sqrt(1) * i/1.
        }
    }

    expr(int i, int j)
    {
        values[1][0] = i;
        values[1][1] = j; // i/j = sqrt(1) * i/j.
        simplify();
    }

    expr(int sqrt, int i, int j)
    {
        values[sqrt][0] = i;
        values[sqrt][1] = j;
        simplify();
    }

    expr(std::map<int, int[2]> _values)
    {
        values = _values;
        simplify();
    }

    double eval()
    {
        double total = 0;
        for (auto const& p : values)
        {
            total += std::sqrt(p.first) * (double)(p.second[0]) / (double)(p.second[1]);
        }

        return total;
    }

    //Operator overloading
    expr operator+(const expr& ex)
    {
        expr sum(values);

        for (auto const& p : ex.values)
        {
            int val[2];
            val[0] = sum.values[p.first][0];
            val[1] = sum.values[p.first][1];

            if (sum.values[p.first][0] == 0)
            {
                sum.values[p.first][0] = p.second[0];
                sum.values[p.first][1] = p.second[1];
            }
            else
            {
                sum.values[p.first][0] = val[0] * p.second[1] + val[1] * p.second[0];
                sum.values[p.first][1] = val[1] * p.second[1];
            }

            if (sum.values[p.first][0] == 0)
            {
                sum.values.erase(p.first);
            }
        }
        sum.simplify();
        return sum;
    }

    expr operator-()
    {
        expr negate(values);
        for (auto const& p : negate.values)
        {
            negate.values[p.first][0] *= -1;
        }
        return negate;
    }

    expr operator-(expr& ex)
    {
        return (*this) + -ex;
    }

    expr operator*(const expr& ex)
    {
        expr prod;
        for (auto const& p : values) //Foil the two expressions when multiplying
        {
            for (auto const& q : ex.values)
            {
                std::pair<int, int> v = simplify_sqrt(p.first * q.first);

                int val1[2];
                val1[0] = prod.values[v.second][0];
                val1[1] = prod.values[v.second][1];
                if (val1[1] == 0)
                {
                    val1[1]++;
                }

                int val2[2];
                val2[0] = v.first * p.second[0] * q.second[0];
                val2[1] = p.second[1] * q.second[1];

                prod.values[v.second][0] = val1[0] * val2[1] + val1[1] * val2[0];
                prod.values[v.second][1] = val1[1] * val2[1];

                if (prod.values[v.second][0] == 0)
                {
                    prod.values.erase(v.second);
                }
            }
        }

        prod.simplify();
        return prod;
    }

    expr operator/(const int& i)
    {
        expr div;
        for (auto const& p : values)
        {
            values[p.first][0] = p.second[0];
            values[p.first][1] = p.second[1] * i;
        }
        div.simplify();
        return div;
    }

    bool operator==(const expr& ex)
    {
        for (auto& p : values)
        {
            if (ex.values.count(p.first) == 0)
            {
                return false;
            }

            else if (p.second[0] != ex.values.at(p.first)[0])
            {
                return false;
            }

            else if (p.second[1] != ex.values.at(p.first)[1])
            {
                return false;
            }
        }

        for (auto& p : ex.values)
        {
            if (values.count(p.first) == 0)
            {
                return false;
            }
        }

        return true;
    }

    bool operator==(const int& i)
    {
        if ((values.size() == 0) && (i == 0))
        {
            return true;
        }

        if (values.size() != 1)
        {
            return false;
        }

        if (values.count(1))
        {
            return (values[1][0] == i) && (values[1][1] == 1); //expression equals sqrt(1)*i/1 = i
        }
        
        return false;
    }

    friend std::ostream& operator<<(std::ostream& os, const expr& ex)
    {
        if (ex.values.size() == 0)
        {
            os << 0;
            return os;
        }

        int frac[2];

        bool firstCoeff = true;
        for (auto const& p : ex.values)
        {
            if (firstCoeff)
            {
                firstCoeff = false;
            }
            else
            {
                os << " + ";
            }

            if (p.first != 1)
            {
                os << "sqrt(" << p.first << ")*";
            }
            frac[0] = ex.values.at(p.first)[0];
            frac[1] = ex.values.at(p.first)[1];

            if (frac[1] == 1)
            {
                os << frac[0];
            }
            else
            {
                os << frac[0] << "/" << frac[1];
            }
        }

        return os;
    }
};

//A (possibly multivariate) monomial.
static std::string variableNames = "xyzwvutsr";
struct monom
{
    std::vector<int> powers; //Powers of variables in the polynomial.
    expr coeff;

    monom(std::vector<int> pow, expr c)
    {
        powers = pow;
        coeff = c;
    }

    double eval(std::vector<double> subst)
    {
        if (powers.size() > subst.size())
        {
            throw std::invalid_argument("Invalid substitution into vector.");
        }

        double a = coeff.eval();

        for (int i = 0; i < powers.size(); i++)
        {
            a *= std::pow(subst[i], powers[i]);
        }

        return a;
    }

    int degree()
    {
        return std::accumulate(powers.begin(), powers.end(), 0);
    }

    //Operator overloading
    monom operator-()
    {
        return monom(powers, -coeff);
    }

    monom operator*(const monom& m)
    {
        std::vector<int> productPowers = {};

        std::vector<int> pow1 = powers;
        std::vector<int> pow2 = m.powers;

        //make lengths of vectors equal
        if (pow1.size() < pow2.size())
        {
            pow1.insert(pow1.end(), pow2.size() - pow1.size(), 0);
        }
        else if (pow2.size() < pow1.size())
        {
            pow2.insert(pow2.end(), pow1.size() - pow2.size(), 0);
        }

        for (int i = 0; i < pow1.size(); i++)
        {
            productPowers.push_back(pow1[i] + pow2[i]);
        }

        return monom(productPowers, coeff * m.coeff);
    }

    monom operator*(expr ex)
    {
        return monom(powers, coeff * ex);
    }

    monom operator/(int i)
    {
        return monom(powers, coeff / i);
    }

    bool operator==(monom f)
    {
        return (coeff == f.coeff) && (powers == f.powers);
    }

    friend std::ostream& operator<<(std::ostream& os, const monom& m)
    {
        if (m.coeff.values.size() <= 1)
        {
            os << m.coeff;
        }
        else
        {
            os << "(" << m.coeff << ")";
        }

        for (int i = 0; i < m.powers.size(); i++)
        {
            int p = m.powers[i];

            if (p == 0)
            {
                continue;
            }
            else if (p == 1)
            {
                os << "*" << variableNames[i];
            }
            else
            {
                os << "*" << variableNames[i] << "^" << p;
            }
        }

        return os;
    }
};

//Order two monomials lexicographically based on their powers - returns whether m1 <= m2 in the ordering.
bool lexicographic(monom &m1, monom &m2)
{
    if (m1.powers.size() == 0)
    {
        return true;
    }
    if (m2.powers.size() == 0)
    {
        return false;
    }

    if (m1.powers.size() <= m2.powers.size())
    {
        for (int i = 0; i < m1.powers.size(); i++)
        {
            if (m1.powers[i] < m2.powers[i]) //m1 < m2
            {
                return true;
            }
            else if (m1.powers[i] > m2.powers[i]) //m1 > m2
            {
                return false;
            }
        }
    }
    else
    {
        for (int i = 0; i < m2.powers.size(); i++)
        {
            if (m1.powers[i] < m2.powers[i])
            {
                return true;
            }
            else if (m1.powers[i] > m2.powers[i])
            {
                return false;
            }
        }
    }

    return true; //m1 == m2
}

//A (possibly multivariate) polynomial.
struct poly
{
    std::vector<monom> coeffs;

    poly(std::vector<monom> c)
    {
        coeffs = c;
        simplify();
    }

    poly(expr ex)
    {
        if (ex == 0)
        {
            coeffs = {};
        }
        else
        {
            coeffs = { monom({}, ex) };
        }
    }

    poly(char varName)
    {
        int i = variableNames.find(varName);

        std::vector<int> pow(i, 0);
        pow.push_back(1);
        coeffs = { monom(pow, expr(1)) };
    }

    void simplify()
    {
        std::sort(coeffs.begin(), coeffs.end(), lexicographic);

        //Remove coefficients
        for (int i = 0; i < coeffs.size(); i++)
        {
            if (coeffs[i].coeff == 0)
            {
                coeffs.erase(coeffs.begin() + i);

                i--; //decrement i because we just deleted a member of coeffs
            }
        }

        //Add together coefficients with the same powers
        for (int i = 1; i < coeffs.size(); i++)
        {
            if (coeffs[i].powers == coeffs[i - 1].powers)
            {
                coeffs[i - 1].coeff = coeffs[i].coeff + coeffs[i - 1].coeff;

                coeffs.erase(coeffs.begin() + i);

                if (coeffs[i-1].coeff == 0)
                {
                    coeffs.erase(coeffs.begin() + i - 1);
                    i--; //see above
                }

                i--; //see above
            }
        }
    }

    poly normalize()
    {
        if (coeffs.size() == 0)
        {
            return poly(expr(0));
        }

        expr ex = coeffs.back().coeff;
        int gcd = (*(ex.values.begin())).second[0];
        for (auto& p : ex.values)
        {
            gcd = std::__gcd(gcd, p.second[0]);
        }

        if (gcd == 0)
        {
            throw std::invalid_argument( "Bad gcd encountered during normalization." );
        }

        return (*this) * expr(1,gcd);
    }

    double eval(const std::vector<double> subst)
    {
        double a = 0;
        for (monom m : coeffs)
        {
            a += m.eval(subst);
        }
        return a;
    }

    expr get_coeff(const std::vector<int> _powers)
    {
        for (int i = 0; i < coeffs.size(); i++)
        {
            if (coeffs[i].powers == _powers)
            {
                return coeffs[i].coeff;
            }
        }

        return expr(0);
    }

    int degree()
    {
        int deg = 0;
        for (int i = 0; i < coeffs.size(); i++)
        {
            deg = std::max(deg, coeffs[i].degree());
        }
        return deg;
    }

    int num_vars()
    {
        int num = 0;
        for (int i = 0; i < coeffs.size(); i++)
        {
            num = std::max(num, (int)(coeffs[i].powers.size()));
        }
        return num;
    }

    //Operator overloading
    poly operator+(const poly& f)
    {
        std::vector<monom> join;
        join.reserve(coeffs.size() + f.coeffs.size());
        join.insert(join.end(), coeffs.begin(), coeffs.end());
        join.insert(join.end(), f.coeffs.begin(), f.coeffs.end());
        return poly(join);
    }

    poly operator-()
    {
        return (*this)*expr(-1);
    }

    poly operator-(poly f)
    {
        return (*this) + f * expr(-1);
    }

    poly operator*(const poly& f)
    {
        std::vector<monom> c;

        c.reserve(coeffs.size() * f.coeffs.size());
        for (int i = 0; i < coeffs.size(); i++)
        {
            for (int j = 0; j < f.coeffs.size(); j++)
            {
                c.push_back(coeffs[i] * f.coeffs[j]);
            }
        }

        return poly(c);
    }

    poly operator*(expr ex)
    {
        std::vector<monom> c;
        c.reserve(coeffs.size());
        for (int i = 0; i < coeffs.size(); i++)
        {
            c.push_back(coeffs[i] * ex);
        }

        return poly(c);
    }

    bool operator==(poly f) const
    {
        if (coeffs.size() != f.coeffs.size())
        {
            return false;
        }
        
        bool eq = true;
        for (int i = 0; i < coeffs.size(); i++)
        {
            monom m1 = coeffs[i];
            monom m2 = f.coeffs[i];

            eq = eq && (m1 == m2);
        }

        return eq;
    }

    friend std::ostream& operator<<(std::ostream& os, const poly& f)
    {
        os << f.coeffs[0];

        for (int i = 1; i < f.coeffs.size(); i++)
        {
            os << " + " << f.coeffs[i];
        }

        return os;
    }
};

//A vector with polynomial entries. Class name is 'vec' to differentiate from std::vector.
struct vec
{
    int dim;
    std::vector<poly> entries;

    vec(std::vector<poly> _entries)
    {
        entries = _entries;
        dim = entries.size();
    }

    vec(int _dim)
    {
        dim = _dim;
        entries = std::vector<poly>(dim, poly(expr(0)));
    }

    vec(std::vector<expr> _entries)
    {
        entries.reserve(_entries.size());
        dim = _entries.size();
        for (int i = 0; i < dim; i++)
        {
            entries[i] = poly(expr(_entries[i]));
        }
    }

    //Dot product.
    poly dot(const vec& v) const
    {
        if (v.dim != dim)
        {
            throw std::invalid_argument("Vector dimensions do not match.");
        }

        poly f(expr(0));
        for (int i = 0; i < dim; i++)
        {
            f = f + v[i] * entries[i];
        }

        return f;
    }

    std::vector<double> eval(const std::vector<double> subst)
    {
        std::vector<double> v = {};
        v.reserve(dim);

        for (int i = 0; i < dim; i++)
        {
            v.push_back(entries[i].eval(subst));
        }

        return v;
    }

    int num_vars()
    {
        int num = 0;
        for (poly f : entries)
        {
            num = std::max(num, f.num_vars());
        }
        return num;
    }

    bool is_invalid() // Used for debugging - tests if any coefficients of a vector have a value of 0 that has not been simplified
    {
        for (poly f : entries)
        {
            for (monom m : f.coeffs)
            {
                for (auto& p : m.coeff.values)
                {
                    if ((p.second[0] == 0) || (p.second[1] == 0))
                    {
                        return true;
                    }
                }
            }
        }
        return false;
    }

    //Operator overloading
    vec operator+(const vec& v)
    {
        std::vector<poly> sum = {};
        for (int i = 0; i < dim; i++)
        {
            sum.push_back(entries[i] + v.entries[i]);
        }
        return vec(sum);
    }

    vec operator*(const poly& f)
    {
        std::vector<poly> prod = {};
        for (int i = 0; i < dim; i++)
        {
            prod.push_back(entries[i]*f);
        }
        return vec(prod);
    }

    vec operator*(const expr& f)
    {
        std::vector<poly> prod = {};
        for (int i = 0; i < dim; i++)
        {
            prod.push_back(entries[i] * f);
        }
        return vec(prod);
    }

    poly operator[](size_t index) const
    {
        if (index >= entries.size())
        {
            throw std::out_of_range("Index out of bounds");
        }
        return entries[index];
    }

    bool operator==(vec v)
    {
        if (dim != v.dim)
        {
            return false;
        }

        bool eq = true;
        for (int i = 0; i < dim; i++)
        {
            eq = eq && (entries[i] == v.entries[i]);
        }

        return eq;
    }

    friend std::ostream& operator<<(std::ostream& os, const vec& v)
    {
        for (poly f : v.entries)
        {
            os << "[ " << f << " ]\n";
        }

        return os;
    }
};

//A matrix containing entries that are expr objects.
struct matrix
{
    std::vector<std::vector<expr>> entries;
    int dim[2];

    matrix(std::vector<std::vector<expr>> _entries)
    {
        entries = _entries;
        dim[0] = entries[0].size();
        dim[1] = entries.size();
    }

    matrix(int dimX, int dimY)
    {
        dim[0] = dimX;
        dim[1] = dimY;
        entries = std::vector<std::vector<expr>>(dimY, std::vector<expr>(dimX, 0));
    }

    //Operator overloading
    expr& operator()(int indexX, int indexY)
    {
        if (indexX >= dim[0] || indexY >= dim[1])
        {
            throw std::out_of_range("Index out of bounds");
        }
        return entries[indexY][indexX];
    }

    //Returns row y as a vector.
    vec rowvector(int y) const
    {
        std::vector<poly> u = {};
        for (int x = 0; x < dim[0]; x++)
        {
            u.push_back(poly(entries[y][x]));
        }
        return vec(u);
    }
    
    //Returns row x as a vector.
    vec columnvector(int x) const
    {
        std::vector<poly> u = {};
        for (int y = 0; y < dim[1]; y++)
        {
            u.push_back(poly(entries[y][x]));
        }
        return vec(u);
    }

    vec operator*(const vec& v)
    {
        if (v.entries.size() != dim[0])
        {
            throw std::invalid_argument("Dimensions do not match for matrix multiplication.");
        }

        std::vector<poly> u;
        u.reserve(dim[1]);

        for (int i = 0; i < dim[1]; i++)
        {
            u.push_back(v.dot(rowvector(i)));
        }
        
        return vec(u);
    }

    matrix operator*(const matrix& M)
    {
        if (M.dim[1] != dim[0])
        {
            throw std::invalid_argument("Dimensions do not match for matrix multiplication.");
        }

        matrix N = matrix(M.dim[0], dim[1]);

        for (int i = 0; i < dim[1]; i++)
        {
            for (int j = 0; j < M.dim[0]; j++)
            {
                N(j,i) = (M.columnvector(j).dot(rowvector(i))).coeffs[0].coeff;
            }
        }

        return N;
    }

    friend std::ostream& operator<<(std::ostream& os, matrix M)
    {
        for (int i = 0; i < M.dim[1]; i++)
        {
            os << "[";
            for (int j = 0; j < M.dim[0]; j++)
            {
                os << " " << M(j, i) << " ";
            }
            os << "]\n";
        }

        return os;
    }
};

poly operator*(expr ex, poly p)
{
    return p * ex;
}

vec operator*(poly f, vec v)
{
    return v * f;
}

vec operator*(expr f, vec v)
{
    return v * f;
}

//Hashing for the various types
template<typename T>
struct std::hash<std::vector<T>>
{
    std::size_t operator()(std::vector<T> const& v) const {
        std::size_t seed = v.size();
        for (auto& i : v) {
            seed ^= std::hash<T>{}(i) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};

template<>
struct std::hash<expr>
{
    std::size_t operator()(const expr& ex) const noexcept
    {
        std::size_t seed = ex.values.size();
        for (auto& p : ex.values)
        {
            seed ^= std::hash<int>{}(p.first)     + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            seed ^= std::hash<int>{}(p.second[0]) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            seed ^= std::hash<int>{}(p.second[1]) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};

template<>
struct std::hash<monom>
{
    std::size_t operator()(const monom& m) const noexcept
    {
        return std::hash<std::vector<int>>{}(m.powers) ^ (std::hash<expr>{}(m.coeff) << 1);
    }
};

template<>
struct std::hash<poly>
{
    std::size_t operator()(const poly& f) const noexcept
    {
        return std::hash<std::vector<monom>>{}(f.coeffs);
    }
};

template<>
struct std::hash<vec>
{
    std::size_t operator()(const vec& m) const noexcept
    {
        return std::hash<std::vector<poly>>{}(m.entries) ^ (std::hash<int>{}(m.dim) << 1);
    }
};
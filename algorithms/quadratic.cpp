#include "operations.cpp"
#include <unordered_map>

namespace quadratic
{
    bool check(std::vector<vec> orbit)
    {
        for (vec v : orbit)
        {
            for (poly f : v.entries)
            {
                if (f.degree() > 2)
                {
                    return false;
                }
                if (f.num_vars() > 1)
                {
                    return false;
                }
            }
        }
        return true;
    }

    std::map<double, int> search(std::vector<vec> orbit)
    {
        std::cout << "Finding polynomials...\n";

        std::unordered_map<poly, int> counter;
        std::unordered_set<poly> polynomials = {};
        for (int i = 1; i < orbit.size(); i++)
        {
            std::cout << i << "/" << orbit.size() << "\n";

            for (int j = i + 1; j < orbit.size(); j++)
            {
                poly f = sqdist(orbit[0], orbit[i]) - sqdist(orbit[0], orbit[j]);

                poly g = f.normalize();

                if (!all_positive(g))
                {
                    polynomials.insert(g);
                    counter[g]++;
                }
            }
        }

        std::cout << "Found " << polynomials.size() << " polynomials.\nSolving polynomials...\n";

        std::map<double, int> roots;
        for (poly f : polynomials)
        {
            for (double root : solve_quadratic(f))
            {
                if (root > 0)
                {
                    if (roots.find(root) == roots.end())
                    {
                        roots.insert({ root,counter[f] });
                    }
                    else
                    {
                        roots[root] += counter[f];
                    }
                }
            }
        }

        std::cout << "Found " << roots.size() << " roots.\n";

        return roots;
    }
}
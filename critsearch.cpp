#include <vector>
#include <typeinfo>
#include <unordered_set>
#include <iomanip>
#include <map>
#include <sys/stat.h>
#include "operations.cpp"
#include "commands.cpp"

#include "algorithms/quadratic.cpp"

jsondata importData;

struct stat sb;

int numAlgorithms = 1;
std::vector<bool (*)(std::vector<vec>)>                 algorithmChecks =   { quadratic::check };
std::vector<std::map<double, int>(*)(std::vector<vec>)> algorithmSearches = { quadratic::search };
std::vector<std::string>                                algorithmNames =    { "quaduni" };

std::vector<vec> orbit;
std::map<double, int> results;

int main()
{
    std::cout << std::setprecision(13);

    std::vector<std::string> input = initialize();

    while (true)
    {
        if (input[0] == "help")
        {
            help();
        }
        else if (input[0] == "import")
        {
            std::string dir = "library/" + input[1];
            if (stat(dir.c_str(), &sb) == 0)
            {
                importData = importjson(dir);
                std::cout << "\tSuccessfully imported " << importData.name << ".\n";

                std::cout << "\tObtaining orbit... \n";

                vec v(importData.dim);
                for (int i = 0; i < importData.vectorGens.size(); i++)
                {
                    v = v + importData.vectorGens[i] * importData.nodeActivations[i];
                }

                orbit = generate_orbit(v, importData.matrixGens);

                std::cout << "Success! Orbit has " << orbit.size() << " points.\n";
            }
            else
            {
                std::cout << "\tThat is not a valid directory. Files are imported through the library folder.\n";
            }
        }
        else if (input[0] == "search")
        {
            for (int i = 0; i < numAlgorithms; i++)
            {
                std::cout << "\tApplicability of implemented search algorithms:\n";
                std::vector<std::string> validAlgorithms;

                if (algorithmChecks[i](orbit))
                {
                    std::cout << "\t\t- Algorithm '" + algorithmNames[i] + "': yes\n";
                    validAlgorithms.push_back(algorithmNames[i]);
                }
                else
                {
                    std::cout << "\t\t- Algorithm '" + algorithmNames[i] + "': no\n";
                }

                std::cout << "\tPlease choose an algorithm:\n>";
                std::string alg = get_command()[0];

                if (std::find(validAlgorithms.begin(), validAlgorithms.end(), alg) != validAlgorithms.end())
                {
                    std::cout << "\tStarting search with algorithm '" + alg + "'...\n";
                    int index = std::find(algorithmNames.begin(), algorithmNames.end(), alg) - algorithmNames.begin();
                    results = algorithmSearches[index](orbit);
                    std::cout << "\tSearch finished.\n";
                }
                else
                {
                    std::cout << "That is not a valid algorithm.\n";
                }
            }
        }
        else if (input[0] == "results")
        {
            leaderboard(results);
        }
        else if (input[0] == "orbit")
        {
            print_example_orbit(orbit);
        }
        else if (input[0] == "example")
        {
            float exval = 1;

            if (input.size() > 1)
            {
                exval = std::stod(input[1]);
            }

            exportoff(importData, eval_orbit(orbit, std::vector<double>(importData.num_vars, exval)), "orbit-example.off");
            std::cout << "Exported example orbit.\n";
        }
        else if (input[0] == "export")
        {

            std::vector<double> roots;
            if (input.size() == 1)
            {
                for (std::pair<double, int> p : results)
                {
                    roots.push_back(p.first);
                }
            }
            else if (input.size() == 2)
            {
                for (std::pair<double, int> p : results)
                {
                    if ( p.second >= std::stoi(input[1]) )
                    {
                        roots.push_back(p.first);
                    }
                }
            }
            else
            {
                for (std::pair<double, int> p : results)
                {
                    if ( (std::stoi(input[2]) >= p.second) && (p.second >= std::stoi(input[1])) )
                    {
                        roots.push_back(p.first);
                    }
                }
            }

            std::cout << "Exporting " << roots.size() << " orbits." << "\n";

            for (int i = 0; i < roots.size(); i++)
            {
                std::vector<std::vector<double>> pts = eval_orbit(orbit, { roots[i] });
                exportoff(importData, pts, "orbit-" + std::to_string(i) + ".off");
            }
        }
        else if (input[0] == "quit")
        {
            break;
        }
        else
        {
            std::cout << "\tThat is not a recognized command. Type help for a list of commands.\n";
        }

        input = get_command();
    }
}
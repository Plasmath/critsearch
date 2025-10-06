#include <fstream>
#include <iostream>

#include "symbolic.cpp"
#include "dependencies/json.hpp"

using json = nlohmann::json;

struct token //A single token, used for processing string expressions
{
    bool isNumeric;
    int value;
    std::string operatorName;

    token(int i)
    {
        isNumeric = true;
        value = i;
    }

    token(std::string s)
    {
        isNumeric = false;
        operatorName = s;
    }
};

struct jsondata
{
    std::string name;

    std::vector<matrix> matrixGens;
    std::vector<vec> vectorGens;
    int dim;
    int num_vars;

    std::string offheader;
    std::string offbody;

    std::vector<poly> nodeActivations;

    jsondata() {};
};

expr parseString(std::string str)
{
    //Tokenization
    std::vector<token> tokens = {};

    int k = 0;
    std::string s = "";
    bool parsingNumber = false;
    bool parsingString = false;
    for (int i = 0; i < str.size(); i++)
    {
        char current = str[i];

        if (current >= '0' && current <= '9') //current character is a number
        {
            parsingNumber = true;
            k *= 10;
            k += current - '0';

            if (parsingString)
            {
                parsingString = false;
                tokens.push_back(token(s));
                s = "";
            }
        }
        else
        {
            if (parsingNumber) //if we were parsing a number, stop and add it to the list of tokens
            {
                parsingNumber = false;
                tokens.push_back(token(k));
                k = 0;
            }

            if (current == '-')
            {
                if (parsingString)
                {
                    parsingString = false;
                    tokens.push_back(token(s));
                    s = "";
                }
                tokens.push_back(token("-"));
            }
            else if (current == '/')
            {
                if (parsingString)
                {
                    parsingString = false;
                    tokens.push_back(token(s));
                    s = "";
                }
                tokens.push_back(token("/"));
            }
            else if (current == ' ')
            {
                if (parsingString)
                {
                    parsingString = false;
                    tokens.push_back(token(s));
                    s = "";
                }
            }
            else
            {
                parsingString = true;
                s += current;
            }
        }
    }

    //check what was being parsed when the string ended
    if (parsingNumber)
    {
        tokens.push_back(token(k));
    }
    else if (parsingString)
    {
        tokens.push_back(token(s));
    }
    
    /*
    // for debugging
    for (token t : tokens)
    {
        if (t.isNumeric)
        {
            std::cout << t.value << " ";
        }
        else
        {
            std::cout << t.operatorName << " ";
        }
    }
    std::cout << "\n";*/

    //feeling lazy right now
    if (tokens[0].operatorName == "-")
    {
        // -sqrt(a)/b
        if (tokens.size() == 6 && tokens[1].operatorName == "sqrt(")
        {
            return expr(tokens[2].value, -1, tokens[5].value);
        }

        // -sqrt(a)
        if (tokens.size() == 4 && tokens[1].operatorName == "sqrt(")
        {
            return expr(tokens[2].value, -1, 1);
        }

        if (tokens.size() == 4 && tokens[2].operatorName == "/") // -a/b
        {
            return -expr(tokens[1].value, tokens[3].value);
        }

        if (tokens.size() == 2) // -a
        {
            return expr(-tokens[1].value);
        }
    }
    else
    {
        // sqrt(a)/b
        if (tokens.size() == 5 && tokens[0].operatorName == "sqrt(")
        {
            return expr(tokens[1].value, 1, tokens[4].value);
        }

        // sqrt(a)
        if (tokens.size() == 3 && tokens[0].operatorName == "sqrt(")
        {
            return expr(tokens[1].value, 1, 1);
        }

        if (tokens.size() == 3 && tokens[1].operatorName == "/") // a/b
        {
            return expr(tokens[0].value, tokens[2].value);
        }

        if (tokens.size() == 1) // a
        {
            return expr(tokens[0].value);
        }
    }

    throw std::invalid_argument("Invalid string expression encountered during import.");
}

jsondata importjson(std::string dir)
{
    std::ifstream f(dir);
    json data = json::parse(f);

    std::string name = data["name"];
    int dim = data["dim"];
    int num_vars = data["num_variables"];

    std::vector<std::vector<std::vector<std::string>>> generators = data["generating_matrices"];
    std::vector<std::vector<std::string>> vectors = data["generating_vectors"];

    std::vector<matrix> matrixGens = {};
    std::vector<vec> vectorGens = {};

    for (std::vector<std::vector<std::string>> mat : generators)
    {
        matrix M(dim, dim);

        for (int i = 0; i < dim; i++)
        {
            for (int j = 0; j < dim; j++)
            {
                M(i, j) = parseString(mat[i][j]);
            }
        }

        matrixGens.push_back(M);
    }

    for (std::vector<std::string> vector : vectors)
    {
        std::vector<poly> v;

        for (int i = 0; i < dim; i++)
        {
            v.push_back(poly(parseString(vector[i])));
        }

        vectorGens.push_back(vec(v));
    }

    std::vector<poly> nodeActivations = {};
    for (int i = 0; i < vectorGens.size(); i++)
    {
        std::string str = data["node_activations"][i];
        if (variableNames.find(str[0]) != std::string::npos)
        {
            nodeActivations.push_back(poly(str[0]));
        }
        else
        {
            nodeActivations.push_back(poly(parseString(str)));
        }
    }

    jsondata jdat;

    jdat.name = name;
    jdat.matrixGens = matrixGens;
    jdat.vectorGens = vectorGens;
    jdat.dim = dim;
    jdat.num_vars = num_vars;
    jdat.nodeActivations = nodeActivations;
    jdat.offheader = data["off_header"];
    jdat.offbody = data["off_body"];

    return jdat;
}

std::vector<vec> generate_orbit(vec &basePoint, std::vector<matrix> &matrixGens)
{
    std::vector<vec> orbit = { basePoint };
    std::vector<vec> untestedPts = { basePoint };
    std::vector<vec> newUntestedPts;
    vec u = vec({});

    int prevSize = 0;
    while (orbit.size() != prevSize)
    {
        std::cout << "Found " << orbit.size() << " points so far in orbit.\n";

        prevSize = orbit.size();
        newUntestedPts = {};

        for (matrix M : matrixGens)
        {
            for (vec v : untestedPts)
            {
                u = M * v;
                //the new vector u is not yet in the orbit
                if (std::find(orbit.begin(), orbit.end(), u) == orbit.end())
                {
                    if (u.is_invalid())
                    {
                        //std::cout << "mat:\n" << M << "\nvec:\n" << v << "\noutput:\n" << u << "\n";
                    }
                    
                    orbit.push_back(u);
                    newUntestedPts.push_back(u);
                }
            }
        }

        untestedPts = newUntestedPts;
    }

    return orbit;
}

std::vector<std::vector<double>> eval_orbit(std::vector<vec>& orbit, const std::vector<double> subst)
{
    std::vector<std::vector<double>> vertices = {};
    vertices.reserve(orbit.size());

    for (vec v : orbit)
    {
        vertices.push_back(v.eval(subst));
    }

    return vertices;
}

void exportoff(jsondata jdat, std::vector<std::vector<double>> pts, std::string dir)
{
    std::ofstream off;
    off.open(dir);
    off << std::setprecision(13);
    off << jdat.offheader << "\n\n";

    for (int i = 0; i < pts.size(); i++)
    {
        for (int j = 0; j < jdat.dim; j++)
        {
            off << pts[i][j] << " ";
        }

        off << "\n";
    }

    off << "\n" << jdat.offbody << "\n";
    off.close();
}
#include <fstream>
#include <iostream>

#include "symbolic.cpp"
#include "dependencies/json.hpp"

using json = nlohmann::json;

struct token;
expr parseTokensConstant(std::vector<token> tokens);

struct token //A single token, used for processing string expressions
{
    bool isNumeric = false;
    bool hasExprValue = false;
    int value;
    std::string operatorName = "";
    std::vector<token> subtokens;
    expr expr_value;

    token(int i)
    {
        isNumeric = true;
        value = i;

        hasExprValue = true;
        expr_value = expr(i);
    }

    token(std::string s)
    {
        operatorName = s;
    }

    token(expr ex)
    {
        hasExprValue = true;
        expr_value = ex;
    }

    token(std::vector<token> tokens)
    {
        if (tokens.size() == 1)
        {
            this->isNumeric = tokens[0].isNumeric;
            this->hasExprValue = tokens[0].hasExprValue;
            this->value = tokens[0].value;
            this->operatorName = tokens[0].operatorName;
            this->subtokens = tokens[0].subtokens;
            this->expr_value = tokens[0].expr_value;
        }
        else
        {
            hasExprValue = true;
            expr_value = parseTokensConstant(tokens);
            subtokens = tokens;
        }
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

std::vector<token> organizeTokensConstant(std::vector<token> tokens)
{
    std::vector<token> new_tokens;
    std::vector<token> paren_tokens;
    int paren_depth = 0;
    bool in_sqrt = false;

    //Parse parentheses
    for (int i = 0; i < tokens.size(); i++)
    {
        token t = tokens[i];

        if (t.operatorName == "sqrt" && (paren_depth == 0))
        { 
            in_sqrt = true;
        }
        else if ((t.operatorName != "(") && (paren_depth == 0))
        {
            new_tokens.push_back(t);
        }
        else if (paren_depth > 0)
        {
            paren_tokens.push_back(t);
        }

        if (t.operatorName == "(")
        {
            paren_depth++;
        }
        else if (t.operatorName == ")")
        {
            paren_depth--;

            if (paren_depth == 0)
            {
                paren_tokens.pop_back(); //remove closing parenthesis ')';
                token u = token(paren_tokens);

                if (in_sqrt)
                {
                    if (u.isNumeric)
                    {
                        u = token(expr(u.value, 1, 1));
                    }
                    else
                    {
                        throw std::invalid_argument("Error when parsing square root. Only square roots of integers are currently supported.");
                    }
                    in_sqrt = false;
                }

                new_tokens.push_back(u);
                paren_tokens = {};
            }
        }
    }

    if (paren_depth != 0)
    {
        throw std::invalid_argument("Invalid string expression encountered during import.");
    }

    return new_tokens;
}

expr parseTokensConstant(std::vector<token> tokens)
{
    tokens = organizeTokensConstant(tokens);

    // Parse multiplication and division (PE[MD]AS)
    std::vector<token> md_tokens = { tokens[0] };

    for (int i = 1; i < tokens.size(); i++)
    {
        if (tokens[i].operatorName == "*")
        {
            md_tokens.back().expr_value = md_tokens.back().expr_value * tokens[i + 1].expr_value;
            i++; //skip the next token because we just processed it
        }
        else if (tokens[i].operatorName == "/")
        {
            if (tokens[i + 1].isNumeric)
            {
                md_tokens.back().expr_value = md_tokens.back().expr_value / tokens[i + 1].value;
                i++; //skip the next token because we just processed it
            }
            else
            {
                throw std::invalid_argument("Error while attempting to parse division by non-integral value.");
            }
        }
        else
        {
            md_tokens.push_back(tokens[i]);
        }
    }

    // Parse addition and subtraction (PEMD[AS])
    std::vector<token> as_tokens = { md_tokens[0] };

    for (int i = 1; i < md_tokens.size(); i++)
    {
        if (md_tokens[i].operatorName == "+")
        {
            as_tokens.back().expr_value = as_tokens.back().expr_value + md_tokens[i + 1].expr_value;
            i++; //skip the next token because we just processed it
        }
        else if (md_tokens[i].operatorName == "-")
        {
            as_tokens.back().expr_value = as_tokens.back().expr_value - md_tokens[i + 1].expr_value;
            i++; //skip the next token because we just processed it
        }
        else
        {
            as_tokens.push_back(md_tokens[i]);
        }
    }

    
    if (as_tokens.size() == 1)
    {
        return as_tokens[0].expr_value;
    }
    else if ((as_tokens[0].operatorName == "-") && (as_tokens.size() == 2)) //It's possible the expression starts with a minus sign, which wouldn't have been processed yet, so we take that into account
    {
        return -as_tokens[1].expr_value;
    }
    else
    {
        throw std::invalid_argument("Error when parsing expression."); //By this point our expression should have been fully simplified.
    }

}

//Parses a string into a constant expression.
expr parseStringConstant(std::string str)
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
            else if (current == '(')
            {
                if (parsingString)
                {
                    parsingString = false;
                    tokens.push_back(token(s));
                    s = "";
                }
                tokens.push_back(token("("));
            }
            else if (current == ')')
            {
                if (parsingString)
                {
                    parsingString = false;
                    tokens.push_back(token(s));
                    s = "";
                }
                tokens.push_back(token(")"));
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

    return parseTokensConstant(tokens);
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
                M(i, j) = parseStringConstant(mat[i][j]);
            }
        }

        matrixGens.push_back(M);
    }

    for (std::vector<std::string> vector : vectors)
    {
        std::vector<poly> v;

        for (int i = 0; i < dim; i++)
        {
            v.push_back(poly(parseStringConstant(vector[i])));
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
            nodeActivations.push_back(poly(parseStringConstant(str)));
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
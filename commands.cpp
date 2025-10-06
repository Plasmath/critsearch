#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include "processing.cpp"

bool sort_pairs(std::pair<double, int> a, std::pair<double, int> b)
{
	return a.second >= b.second;
}

void leaderboard(std::map<double, int> results)
{
	std::cout << " value    | hits\n";
	std::cout << "----------+------\n";

	std::vector<std::pair<double, int>> roots;
	for (auto& p : results)
	{
		roots.push_back(p);
	}
	std::sort(roots.begin(), roots.end(), sort_pairs);

	for (int i = 0; i < std::min( 25, (int)(roots.size()) ); i++)
	{
		std::pair<double, int> root = roots[i];

		std::string val = std::to_string(root.first);
		val = val.substr( 0, std::min(8, (int)(val.size()) ) ) +
			  std::string( std::max( 0, 8 - (int)(val.size()) ) , '0' );

		std::cout << " " << val << " | " << root.second << "\n";
	}
}

std::vector<std::string> get_command()
{
	std::string str;
	
	std::cout << "> ";
	std::getline(std::cin, str);

	std::vector<std::string> input;

	std::string s = "";
	for (int i = 0; i < str.length(); i++)
	{
		if (str[i] == ' ')
		{
			input.push_back(s);
			s = "";
		}
		else
		{
			s += str[i];
		}
	}
	input.push_back(s);

	return input;
}

std::vector<std::string> initialize()
{
	std::cout << "\n------------------\n";
	std::cout <<   "  critsearch.cpp  \n";
	std::cout <<   "------------------\n\n";
	std::cout << "Welcome to critsearch. Please input a command (type help for a list of commands).\n";
	return get_command();
}

void help()
{
	std::cout << "\thelp - Lists commands implemented in the program.\n";
	std::cout << "\timport (filename) - Used to import a JSON file from the library folder into the program.\n";
	std::cout << "\tsearch - Used to search the imported JSON with a search algorithm.\n";
	std::cout << "\tresults - Display a table of the results of a search.\n";
	std::cout << "\texport - Export all solutions to a search program.\n";
	std::cout << "\texport (number) - Export all solutions to a search program with more than (number) hits.\n";
	std::cout << "\texport (num1) (num2) - Export all solutions to a search program with hits greater than (num1) and less than (num2).\n";
	std::cout << "\torbit - Gives the coordinates of an example orbit (helpful when creating .json files).\n";
	std::cout << "\texample - Exports an .off file of an example orbit (helpful when creating .json files).\n";
	std::cout << "\tquit - Quits the program.\n";
}

void print_example_orbit(std::vector<vec> &orbit)
{
	int num = 0;
	for (vec v : orbit)
	{
		num = std::max(num, v.num_vars());
	}

	std::vector<double> subst(num, 2);

	for (std::vector<double> v : eval_orbit(orbit, subst))
	{
		for (int i = 0; i < v.size(); i++)
		{
			std::cout << v[i] << " ";
		}
		std::cout << "\n";
	}
}
#include "symbolic.cpp"
#include <algorithm>
#include <cmath>

#pragma once

//The squared distance between two vectors v and u. Assumes that the vectors are the same length.
poly sqdist(vec v, vec u)
{
	poly sum(expr(0));

	for (int i = 0; i < v.dim; i++)
	{
		poly subt = v[i] - u[i];
		sum = sum + subt * subt;
	}

	return sum;
}

//Returns whether all coefficients of a poly object are positive.
bool all_positive(poly &f)
{
	for (monom m : f.coeffs)
	{
		if (m.coeff.eval() < 0)
		{
			return false;
		}
	}
	return true;
}

//Partial derivative with respect to var.
poly derivative(poly f, char var)
{
	int vIndex = variableNames.find(var);
	
	std::vector<monom> dCoeffs = {};

	for (monom m : f.coeffs)
	{
		if (m.powers.size() <= vIndex) //remove factors that are constant with respect to var
		{
			continue;
		}
		else if (m.powers[vIndex] == 0)
		{
			continue;
		}

		monom n = monom(m.powers, m.coeff);
		n.coeff = n.coeff * expr(m.powers[vIndex]);
		n.powers[vIndex]--;
		dCoeffs.push_back(n);
	}

	return poly(dCoeffs);
}

matrix jacobian(vec v, std::vector<double> subst)
{
	matrix M(v.dim, v.dim);

	for (int i = 0; i < v.dim; i++)
	{
		for (int j = 0; j < v.dim; j++)
		{
			poly f = derivative(v.entries[i], variableNames[j]);
			M(i, j) = f.eval(subst);
		}
	}
	return M;
}

//Find all real solutions to a univariate quadratic polynomial.
std::vector<double> solve_quadratic(poly f)
{
	double a = f.get_coeff({ 2 }).eval();
	double b = f.get_coeff({ 1 }).eval();
	double c = f.get_coeff({ }).eval();

	if (a == 0 && b == 0)
	{
		return {}; //constant, cannot solve for x
	}

	if (a == 0)
	{
		//polynomial is linear & of the form b*x + c
		return { - c / b };
	}
	else
	{
		//polynomial is actually quadratic & of the form a*x^2 + b*x + c
		double discriminant = b*b - 4*a*c;

		if (discriminant < 0) //roots are both complex
		{
			return {};
		}
		else
		{
			//quadratic formula
			double sq = std::sqrt(discriminant);
			return { (-b + sq) / (2 * a), (-b - sq) / (2 * a) };
		}
	}
}

//Perform one iteration of Newton's method on a given univariate polynomial & approximation.
double newton_iter(double approx, poly f)
{
	return approx - f.eval({ approx }) / derivative(f, 'x').eval({ approx });
}
#include "solver.hpp"
#include <math.h>
#include <iostream>
#include <stdexcept>
#include <complex>

namespace solver
{
	bool doublesEqual(double d1, double d2)
	{
		return abs(d1 - d2) <= EPSILON;
	}

	bool complexAndDoublEqual(const std::complex<double>& complex, double real)
	{
		return doublesEqual(std::imag(complex), 0.0) && doublesEqual(std::real(complex), real);
	}

	bool isReal(const std::complex<double>& complex)
	{
		return doublesEqual(std::imag(complex), 0.0);
	}

	bool complexEqual(const std::complex<double>& complex1, const std::complex<double>& complex2)
	{
		return doublesEqual(std::real(complex1), std::real(complex2)) && doublesEqual(std::imag(complex1), std::imag(complex2));
	}

	// assumes real is a quadratic polynomial
	void findQuadraticEquationRoots(const RealVariable& real, double roots[2], int& numRoots)
	{
		double a = real.getA();
		double b = real.getB();
		double c = real.getC();

		double delta = b * b - 4 * a * c;

		if (delta < 0) // no roots
		{
			numRoots = 0;
		}
		else
		{
			roots[0] = (-b + sqrt(delta)) / (2 * a);
			roots[1] = (-b - sqrt(delta)) / (2 * a);

			numRoots = doublesEqual(roots[0], roots[1]) ? 1 : 2;
		}
	}

	void findQuadraticEquationRoots(const ComplexVariable& complex, std::complex<double> roots[2], int& numRoots)
	{
		const std::complex<double>& a = complex.getA();
		const std::complex<double>& b = complex.getB();
		const std::complex<double>& c = complex.getC();

		std::complex<double> delta = b * b - 4.0 * a * c;

		roots[0] = (-b + sqrt(delta)) / (2.0 * a);
		roots[1] = (-b - sqrt(delta)) / (2.0 * a);

		numRoots = complexEqual(roots[0], roots[1]) ? 1 : 2;
	}

	double solve(const RealVariable& real)
	{
		double a = real.getA();
		double b = real.getB();
		double c = real.getC();

		if (doublesEqual(a, 0.0))
		{
			if (doublesEqual(b, 0.0))
			{
				throw std::invalid_argument("There is no solution");
			}
			else
			{
				return -c / b;
			}
		}
		else
		{
			double roots[2];
			int numRoots;

			findQuadraticEquationRoots(real, roots, numRoots);

			if (numRoots == 0)
			{
				throw std::invalid_argument("Negative root is not allowed");
			}
			else
			{
				return roots[0]; // return first root
			}
		}
	}

    RealVariable::RealVariable(double c, double b, double a)
        : c(c), b(b), a(a)
    {

    }

	RealVariable::RealVariable() : c(0), b(1), a(0)
	{

	}

	RealVariable::RealVariable(double c) : c(c), b(0), a(0)
	{

	}

	RealVariable RealVariable::fromComplexVariable(const ComplexVariable& complex)
	{
		return RealVariable(std::real(complex.getC()), std::real(complex.getB()), std::real(complex.getA()));
	}

	RealVariable::operator ComplexVariable() const
	{
		return ComplexVariable(c, b, a);
	}

    RealVariable operator +(double num, const RealVariable& real)
    {
        return real + num;
    }

    RealVariable RealVariable::operator +(const RealVariable& other) const
    {
         ComplexVariable complexRes = (ComplexVariable)(*this) + (ComplexVariable)other;
		 return RealVariable::fromComplexVariable(complexRes);
	}
    RealVariable operator -(double num, const RealVariable& real)
    {
        return RealVariable(-real.c + num, -real.b, -real.a);
    }

    RealVariable RealVariable::operator -(const RealVariable& real) const
    {
        return RealVariable(c - real.c, b - real.b, a - real.a);
    }

    RealVariable operator *(double num, const RealVariable& real)
    {
        return RealVariable(num, 0.0, 0.0) * real;
    }

    RealVariable RealVariable::operator *(const RealVariable& real) const
    {
        if(doublesEqual(real.a, 0.0) && doublesEqual(real.b, 0.0)) // (ax^2 + bx + c) * k
        {
            return RealVariable(c * real.c, b * real.c, a * real.c);
        }
        else if(doublesEqual(real.a, 0.0)) // ( bx + c) * (dx + k)
        {
            return RealVariable(c * real.c, b * real.c + c * real.b,
                b * real.b);
        }
        else // (c)* (ax^2 + bx + c)
        {
            return real * (*this);
        }
    }

    RealVariable RealVariable::operator /(const RealVariable& real) const
    {
        if(real.a == 0 && real.b == 0) // (a1x^2 +b1x + c1) / (c2)
        {
            if(real.c == 0)
            {
                // throw ...
                throw std::invalid_argument("Divide zero");
            }
            else
            {
                return RealVariable(c / real.c, b / real.c, a / real.c);
            }
        }
        else if(doublesEqual(real.a,0) && !doublesEqual(real.b, 0)) 
        {
			if (doublesEqual(a, 0.0)) // (b1x + c1) / (b2x + c2) where (b1x + c1) = k(b2x + c2)
			{
				return RealVariable(b / real.b, 0.0, 0.0);
			}
			else // (a1x^2 +b1x + c1) / (b2x + c2)
			{
				// for division to be valid, the solution of (b2x + c2) = 0, must be a root
				// of (a1x^2 +b1x + c1)
				// then, the result of the division will be (a1 / b2) * (x - p) where p
				// is the other root of (a1x^2 +b1x + c1)
				double roots[2];
				int numRoots;

				double realSolution = -real.c / real.b;
				findQuadraticEquationRoots(*this, roots, numRoots);

				double rootLeft = doublesEqual(roots[0], realSolution) ? roots[1] : roots[0];

				RealVariable x;
				return (a / real.b) * (x - rootLeft);
			}
        }
         
        else // (a1x^2 + b1x + c1) = k(a2x^2 + b2x + c2)
        {
            return RealVariable(a / real.a, 0.0, 0.0);
        }
    }

    RealVariable operator /(double num, const RealVariable& real)
    {
        return RealVariable(num, 0.0, 0.0) / real;
    }

    RealVariable RealVariable::operator==(const RealVariable& real) const
    {
        return RealVariable(c - real.c, b - real.b, a - real.a);
    }

    RealVariable operator==(double num, const RealVariable& real)
    {
        return RealVariable(num, 0.0, 0.0) == real;
    }

    // 0 ^ (k) | k <= 0  : error
    // 0 ^ (k) | k > 0 : 0
    // 0 ^ (ax^2 + bx + c) : illegal
    // 1 ^ (ax^2 + bx + c)) = 1

    // (ax^2 +bx + c) ^ (0) = 1
    // (ax^2 +bx + c) ^ (1) = (ax^2 +bx + c)

    // (c1) ^ (c2) | c2 != 0. 1
    // (3x + 5) ^ (c2) | c2 != 0. 1
    
    // (bx + c) ^ 2 = b^2 + 2bcx + c^2
    RealVariable RealVariable::operator ^(const RealVariable& real) const
    {
		if(real.c > 2 || !doublesEqual(real.b, 0) || !doublesEqual(real.a, 0))
		{
			throw std::invalid_argument("Power invalid");
		}                                                                                                                                                                                                                                                                                                                                                                    
        // 0 ^ c2
        if(doublesEqual(a, 0) && doublesEqual(b, 0) && doublesEqual(c, 0))
        {
            if(real.c <= 0)
            {
                // throw error
                throw std::invalid_argument("Nagative power or zero power");
            }
            else // c2 > 0
            {
                return RealVariable(0.0, 0.0, 0.0);
            }
        }
        // 1 ^ (a2x^2 + b2x + c2)
        else if(doublesEqual(a, 0) && doublesEqual(b, 0) && doublesEqual(c, 1))
        {
            return RealVariable(1.0, 0.0, 0.0);
        }
        // (a1x^2 +b1x + c1) ^ 0
        else if(doublesEqual(real.a, 0) && doublesEqual(real.b, 0) && doublesEqual(real.c, 0))
        {
            return RealVariable(1.0, 0.0, 0.0);
        }
          // (a1x^2 +b1x + c1) ^ 1
        else if(doublesEqual(real.a, 0) && doublesEqual(real.b, 0) && doublesEqual(real.c, 1))
        {
            return *this;
        }
        //(c1) ^ (c2) | c2 != 0. 1
        else if(doublesEqual(b, 0))
        {
            return pow(c, real.c);
        }
        // (b1x + c1) ^ 2
        else
        {
            // b1^2x^2 + 2b1c1x + c1^2
            return RealVariable(c * c, 2 * b * c, b * b);
        }   
    }
   
	RealVariable operator ^(double num, const RealVariable& real)
	{
		return RealVariable(num, 0.0, 0.0) ^ real;
	}

// -----------------------------------------------

    ComplexVariable::ComplexVariable(const std::complex<double>& c, 
	const std::complex<double>& b, const std::complex<double>& a)
        : c(c), b(b), a(a)
    {

    }

	ComplexVariable::ComplexVariable() : c(), b(1.0), a()
	{

	}

	ComplexVariable::ComplexVariable(double c) : c(c), b(), a()
	{

	}

	ComplexVariable::ComplexVariable(const std::complex<double>& c)
		: c(c), b(), a()
	{

	}

	ComplexVariable ComplexVariable::operator +(const ComplexVariable& complex) const
    {
        return ComplexVariable(c + complex.c, b + complex.b, a + complex.a);
    }


	ComplexVariable ComplexVariable::operator -(const ComplexVariable& complex) const
	{
		return ComplexVariable(c - complex.c, b - complex.b, a - complex.a);
	}

	ComplexVariable ComplexVariable::operator *(const ComplexVariable& complex) const
	{
		if(complexAndDoublEqual(complex.a, 0.0) && complexAndDoublEqual(complex.b, 0.0)) // (ax^2 + bx + c) * k
        {
            return ComplexVariable(c * complex.c, b * complex.c, a * complex.c);
        }
        else if(complexAndDoublEqual(complex.a, 0.0)) // ( bx + c) * (dx + k)
        {
            return ComplexVariable(c * complex.c, b * complex.c + c * complex.b,
                b * complex.b);
        }
        else // (c)* (ax^2 + bx + c)
        {
            return complex * (*this);
        }
	}

	ComplexVariable ComplexVariable::operator /(const ComplexVariable& complex) const
	{
		if(complexAndDoublEqual(complex.a, 0.0) && complexAndDoublEqual(complex.b, 0.0)) // (a1x^2 +b1x + c1) / (c2)
        {
            if(complexAndDoublEqual(complex.c, 0.0))
            {
                // throw ...
                throw std::invalid_argument("Divide zero");
            }
            else
            {
                return ComplexVariable(c / complex.c, b / complex.c, a / complex.c);
            }
        }
        else if(complexAndDoublEqual(complex.a, 0.0) && !complexAndDoublEqual(complex.b, 0.0)) 
        {
			if (complexAndDoublEqual(a, 0.0)) // (b1x + c1) / (b2x + c2) where (b1x + c1) = k(b2x + c2)
			{
				return ComplexVariable(b / complex.b, 0.0, 0.0);
			}
			else // (a1x^2 +b1x + c1) / (b2x + c2)
			{
				// for division to be valid, the solution of (b2x + c2) = 0, must be a root
				// of (a1x^2 +b1x + c1)
				// then, the result of the division will be (a1 / b2) * (x - p) where p
				// is the other root of (a1x^2 +b1x + c1)
				std::complex<double> roots[2];
				int numRoots;

				std::complex<double> compSolution = -complex.c / complex.b;
				findQuadraticEquationRoots(*this, roots, numRoots);

				std::complex<double>& rootLeft = complexEqual(roots[0], compSolution) ? roots[1] : roots[0];

				ComplexVariable x;
				return (a / complex.b) * (x - rootLeft);
			}
        }
         
        else // (a1x^2 + b1x + c1) = k(a2x^2 + b2x + c2)
        {
            return ComplexVariable(a / complex.a, 0.0, 0.0);
        }
	}

	ComplexVariable ComplexVariable::operator ^(const ComplexVariable& complex) const
	{
		if(complex.getC().real() > 2.0 || complex.getC().real() < 0.0 || !isReal(complex.c) || !complexAndDoublEqual(complex.b, 0) || !isReal(complex.b) || !complexAndDoublEqual(complex.a, 0) || !isReal(complex.a))
		{
			throw std::invalid_argument("Power invalid");
		}     
		// 0 ^ c2
        if(complexAndDoublEqual(a, 0.0) && complexAndDoublEqual(b, 0.0) && complexAndDoublEqual(c, 0.0))
        {
            if(isReal(complex.c) && std::real(complex.c) <= 0.0)
            {
                // throw error
                throw std::invalid_argument("Nagative power or zero power");
            }
            else // c2 > 0
            {
                return ComplexVariable(0.0, 0.0, 0.0);
            }
        }
        // 1 ^ (a2x^2 + b2x + c2)
        else if(complexAndDoublEqual(a, 0.0) && complexAndDoublEqual(b, 0.0) && complexAndDoublEqual(c, 1.0))
        {
            return ComplexVariable(1.0, 0.0, 0.0);
        }
        // (a1x^2 +b1x + c1) ^ 0
        else if(complexAndDoublEqual(complex.a, 0.0) && complexAndDoublEqual(complex.b, 0.0) && complexAndDoublEqual(complex.c, 0.0))
        {
            return ComplexVariable(1.0, 0.0, 0.0);
        }
          // (a1x^2 +b1x + c1) ^ 1
        else if(complexAndDoublEqual(complex.a, 0.0) && complexAndDoublEqual(complex.b, 0.0) && complexAndDoublEqual(complex.c, 1.0))
        {
            return *this;
        }
        //(c1) ^ (c2) | c2 != 0. 1
        else if(complexAndDoublEqual(b, 0.0))
        {
            return pow(c, complex.c);
        }
        // (b1x + c1) ^ 2
        else
        {
            // b1^2x^2 + 2b1c1x + c1^2
            return ComplexVariable(c * c, 2.0 * b * c, b * b);
        }   
	}

	ComplexVariable ComplexVariable::operator==(const ComplexVariable& complex) const
	{
        return ComplexVariable(c - complex.c, b - complex.b, a - complex.a);
	}

	ComplexVariable operator *(double num, const ComplexVariable& complex)
	{
		return complex * num;
	}

	ComplexVariable operator *(const std::complex<double>& compNum, const ComplexVariable& complex)
	{
		return complex * compNum;
	}

	ComplexVariable operator /(double num, const ComplexVariable& complex)
	{
		return ComplexVariable(num) / complex;
	}

	ComplexVariable operator /(const std::complex<double>& compNum, const ComplexVariable& complex)
	{
		return ComplexVariable(compNum) / complex;
	}

	ComplexVariable operator ^(double num, const ComplexVariable& complex)
	{
		return ComplexVariable(num) ^ complex;
	}

	ComplexVariable operator ^(const std::complex<double>& compNum, const ComplexVariable& complex)
	{
		return ComplexVariable(compNum) ^ complex;
	}

	ComplexVariable operator +(double num, const ComplexVariable& complex)
	{
		return complex + num;
	}

	ComplexVariable operator +(const std::complex<double>& compNum, const ComplexVariable& complex)
	{
		return complex + compNum;
	}

	ComplexVariable operator -(double num, const ComplexVariable& complex)
	{
		return ComplexVariable(num) - complex;
	}

	ComplexVariable operator -(const std::complex<double>& compNum, const ComplexVariable& complex)
	{
		return ComplexVariable(compNum) - complex;
	}

	ComplexVariable operator ==(double num, const ComplexVariable& complex)
	{
		return ComplexVariable(num) == complex;
	}

	ComplexVariable operator ==(const std::complex<double>& compNum, const ComplexVariable& complex)
	{
		return ComplexVariable(compNum) == complex;
	}

    std::complex<double> solve(const ComplexVariable& complex)
    {
        const std::complex<double>& a = complex.getA();
		const std::complex<double>& b = complex.getB();
		const std::complex<double>& c = complex.getC();

		if (complexAndDoublEqual(a, 0.0))
		{
			if (complexAndDoublEqual(b, 0.0))
			{
				throw std::invalid_argument("There is no solution");
			}
			else // bx = c
			{
				return -c / b;
			}
		}
		else
		{
			std::complex<double> roots[2];
			int numRoots;

			findQuadraticEquationRoots(complex, roots, numRoots);

			return roots[0]; // return first root
		}
    }
};

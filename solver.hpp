#define EPSILON 0.0001
#include <vector>
#include <complex>
namespace solver
{

        class ComplexVariable
    {
        private:   
            std::complex<double> a, b, c;

        public:
		 		ComplexVariable();
			    ComplexVariable(double c);
                ComplexVariable(const std::complex<double>& c);
                ComplexVariable(const std::complex<double>& c,
                 const std::complex<double>& b, const std::complex<double>& a);

                const std::complex<double>& getA() const
                {
                    return a;
                }
                const std::complex<double>& getB() const
                {
                    return b;
                }
                const std::complex<double>& getC() const
                {
                    return c;
                }

                ComplexVariable operator +(const ComplexVariable& complex) const;
                ComplexVariable operator -(const ComplexVariable& complex) const;
                ComplexVariable operator *(const ComplexVariable& complex) const;
                ComplexVariable operator /(const ComplexVariable& complex) const;
                ComplexVariable operator ^(const ComplexVariable& complex) const;
                ComplexVariable operator==(const ComplexVariable& complex) const;

                friend ComplexVariable operator *(double, const ComplexVariable&);
                friend ComplexVariable operator *(const std::complex<double>&, const ComplexVariable&);
                friend ComplexVariable operator /(double, const ComplexVariable&);
                friend ComplexVariable operator /(const std::complex<double>&, const ComplexVariable&);                                                                                                                                                                                                                    
                friend ComplexVariable operator ^(double, const ComplexVariable&);
                friend ComplexVariable operator ^(const std::complex<double>&, const ComplexVariable&);                
                friend ComplexVariable operator +(double, const ComplexVariable&);
                friend ComplexVariable operator +(const std::complex<double>&, const ComplexVariable&);
                friend ComplexVariable operator -(double, const ComplexVariable&);
                friend ComplexVariable operator -(const std::complex<double>&, const ComplexVariable&);
		        friend ComplexVariable operator ==(double, const ComplexVariable&);
                friend ComplexVariable operator ==(const std::complex<double>&, const ComplexVariable&);
                

    };

    std::complex<double> solve(const ComplexVariable& complex);

    
    class RealVariable;
    double solve(const RealVariable& real);

    class RealVariable
    {
        private:   
            double a, b, c;

        public:
			    RealVariable();
			    RealVariable(double c);
                RealVariable(double c, double b, double a);

                double getA() const
                {
                    return a;
                }
                double getB() const
                {
                    return b;
                }
                double getC() const
                {
                    return c;
                }

                static RealVariable fromComplexVariable(const ComplexVariable& complex);

                operator ComplexVariable() const;

                RealVariable operator +(const RealVariable& real) const;
                RealVariable operator -(const RealVariable& real) const;
                RealVariable operator *(const RealVariable& real) const;
                RealVariable operator /(const RealVariable& real) const;
                RealVariable operator ^(const RealVariable& real) const;
                RealVariable operator==(const RealVariable& real) const;                                                                                                                                                                                                                                              


                friend RealVariable operator *(double, const RealVariable&);
                friend RealVariable operator /(double, const RealVariable&);                                                                                                                                                                                                                    
                friend RealVariable operator ^(double, const RealVariable&);                
                friend RealVariable operator +(double, const RealVariable&);
                friend RealVariable operator -(double, const RealVariable&);
		        friend RealVariable operator ==(double, const RealVariable&);
                

    };
    
    bool isReal(const std::complex<double>& complex);
    bool complexAndDoublEqual(const std::complex<double>& complex, double real);
    bool doublesEqual(double d1, double d2);
    bool complexEqual(const std::complex<double>& complex1, const std::complex<double>& complex2);
	void findQuadraticEquationRoots(const RealVariable& real, double roots[2], int& numRoots);
    void findQuadraticEquationRoots(const ComplexVariable& complex, std::complex<double> roots[2], int& numRoots);
};
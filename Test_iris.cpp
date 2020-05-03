#include "doctest.h"
#include "solver.hpp"

using namespace std;
using solver::solve, solver::RealVariable, solver::ComplexVariable;
TEST_CASE("Test") 
{
 RealVariable x;


     CHECK(solve(2*x-4 == 10) ==7);
     CHECK(solve(3*x-5 == 1) ==2);
     CHECK(solve(4*x-4 == 16) ==5);
     CHECK(solve(2*x-4 == 10) ==7);
     CHECK(solve(5*x-6 == 14) ==4);
     CHECK(solve(1*x-4 == 3) ==7);
     CHECK(solve(5*x-4 == 36) ==8);
     CHECK(solve(2*x-3 == 15) ==9);
     CHECK(solve(5*x-2 == 8) ==2);
     CHECK(solve(5*x-3 == 2) ==1);

     CHECK(solve(3*x-1 == 2) ==1);
     CHECK(solve(3*x-2 == 4) ==2);
     CHECK(solve(4*x-1 == 3) ==1);
     CHECK(solve(4*x-2 == 6) ==2);
     CHECK(solve(5*x-1 == 4) ==1);
     CHECK(solve(5*x-2 == 8) ==2);
     CHECK(solve(6*x-1 == 5) ==1);
     CHECK(solve(6*x-2 == 10) ==2);
     CHECK(solve(7*x-1 == 6) ==1);
     CHECK(solve(7*x-2 == 12) ==2);

     CHECK(solve(3*x-3 == 6) ==3);
     CHECK(solve(3*x-4 == 5) ==3);
     CHECK(solve(4*x-3 == 9) ==3);
     CHECK(solve(4*x-4 == 8) ==3);
     CHECK(solve(5*x-3 == 12) ==3);
     CHECK(solve(5*x-4 == 11) ==3);
     CHECK(solve(6*x-3 == 15) ==3);
     CHECK(solve(6*x-4 == 14) ==3);
     CHECK(solve(7*x-3 == 18) ==3);
     CHECK(solve(7*x-4 == 17) ==3);

    
     CHECK(solve((x^2) == 1) ==1); 
     CHECK(solve((x^2) == 4) ==2);
     CHECK(solve((x^2) == 9) ==3);
     CHECK(solve((x^2) == 16) ==4);
     CHECK(solve((x^2) == 25) ==5);
     CHECK(solve((x^2) == 36) ==6);
     CHECK(solve((x^2) == 49) ==7);
     CHECK(solve((x^2) == 64) ==8);
     CHECK(solve((x^2) == 81) ==9);
     CHECK(solve((x^2) == 100) ==10); 

     CHECK(solve((x^2) == 121) ==11); 
     CHECK(solve((x^2) == 1936) ==44);
     CHECK(solve((x^2) == 9801) ==99);
     CHECK(solve((x^2) == 1089) ==33);
     CHECK(solve((x^2) == 484) ==22);
     CHECK(solve((x^2) == 3025) ==55);
     CHECK(solve((x^2) == 4356) ==66);
     CHECK(solve((x^2) == 5929) ==77);
     CHECK(solve((x^2) == 7744) ==88);
     CHECK(solve((x^2) == 12321) ==111);  

     CHECK(solve((x^2) + 2*x + 4.0 == 20 + 6.0*x/2 - x) ==4);   
     CHECK(solve((x^2) + 3*x + 6.0 == 12 + 6.0*x/2 - x) ==2);//16
     CHECK(solve((x^2) + 2*x + 5.0 == 27 + 3.0*x/2 - x) ==4);//29
     CHECK(solve((x^2) + 5*x + 7.0 == 18 + 5.0*x/2 - x) ==2);//21
     CHECK(solve((x^2) + 5*x + 1.0 == 22 + 4.0*x/2 - x) ==3);//25
     CHECK(solve((x^2) + 9*x + 3.0 == 51 + 4.0*x/2 - x) ==4);//55
     CHECK(solve((x^2) + 1*x + 1.0 == 1 + 6.0*x/2 - x) ==1);//3
     CHECK(solve((x^2) + 2*x + 8.0 == 33 + 6.0*x/2 - x) ==5);//43
     CHECK(solve((x^2) + 2*x + 5.0 == 41 + 6.0*x/2 - x) ==6);//53
     CHECK(solve((x^2) + 2*x + 10.0 == 16 + 4.0*x/2 - x) ==2);//18
     
     CHECK(solve((x^2) + 1*x + 1.0 == 4 + 0.0*x/2 - x) ==1);//3   
     CHECK(solve((x^2) + 1*x + 2.0 == 10 + 0.0*x/2 - x) ==2);//8
     CHECK(solve((x^2) + 3*x + 3.0 == 24 + 0.0*x/2 - x) ==3);//21
     CHECK(solve((x^2) + 4*x + 4.0 == 40 + 0.0*x/2 - x) ==4);//36
     CHECK(solve((x^2) + 5*x + 5.0 == 60 + 0.0*x/2 - x) ==5);//55
     CHECK(solve((x^2) + 6*x + 6.0 == 84 + 0.0*x/2 - x) ==6);//78
     CHECK(solve((x^2) + 7*x + 7.0 == 112 + 0.0*x/2 - x) ==7);//105
     CHECK(solve((x^2) + 8*x + 8.0 == 144 + 0.0*x/2 - x) ==8);//136
     CHECK(solve((x^2) + 9*x + 9.0 == 180 + 0.0*x/2 - x) ==9);//171
     CHECK(solve((x^2) + 10*x + 10.0 == 220 + 0.0*x/2 - x) ==10);//210
     
     
     CHECK(solve(2*x-4.0 == 10.0) ==7);
     CHECK(solve(3*x-5.0 == 1.0) ==2);
     CHECK(solve(4*x-4.0 == 16.0) ==5);
     CHECK(solve(2*x-4.0 == 10.0) ==7);
     CHECK(solve(5*x-6.0 == 14.0) ==4);
     CHECK(solve(1*x-4.0 == 3.0) ==7);
     CHECK(solve(5*x-4.0 == 36.0) ==8);
     CHECK(solve(2*x-3.0 == 15.0) ==9);
     CHECK(solve(5*x-2.0 == 8.0) ==2);
     CHECK(solve(5*x-3.0 == 2.0) ==1);

     CHECK(solve(3*x-1.0 == 2.0) ==1);
     CHECK(solve(3*x-2.0 == 4.0) ==2);
     CHECK(solve(4*x-1.0 == 3.0) ==1);
     CHECK(solve(4*x-2.0 == 6.0) ==2);
     CHECK(solve(5*x-1.0 == 4.0) ==1);
     CHECK(solve(5*x-2.0 == 8.0) ==2);
     CHECK(solve(6*x-1.0 == 5.0) ==1);
     CHECK(solve(6*x-2.0 == 10.0) ==2);
     CHECK(solve(7*x-1.0 == 6.0) ==1);
     CHECK(solve(7*x-2.0 == 12.0) ==2);

     CHECK_THROWS(solve(0*x==20));
     CHECK_THROWS(solve(2*x-2*x == 10));
     CHECK_THROWS(solve((x^2)==-4));
     CHECK_THROWS(solve((x^2) + 1*x + 1.0 == 4 + 6.0*x/0.0 - x));
     CHECK_THROWS(solve(3*0*x==45));
     CHECK_THROWS(solve(3*x-3*x+1*x-1*x==20));
     CHECK_THROWS(solve(x/0-0+0==20));
     CHECK_THROWS(solve(0*x+5.0==10));
     CHECK_THROWS(solve(0*x-7.0==30));
     CHECK_THROWS(solve(0*x+11==40));

     
}

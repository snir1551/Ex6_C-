//
// Created by shahar on 24/04/2020.
//
#include <string>
#include <iostream>
#include "doctest.h"
#include "solver.hpp"
using namespace std;
using namespace solver;


TEST_CASE("simple equation") {
     RealVariable x;
     CHECK(solve(2*x==10)==5);
     CHECK(solve(-2*x==-10)==5);
     CHECK(solve(2*x==0)==0);
     CHECK(solve(8*x==0)==0);
     CHECK(solve(8*x/8==1)==1);
     CHECK(solve(3*x==90)==30);
     CHECK(solve(4*x+10==90)==20);
     CHECK(solve(4*x-10==90)==25);
     CHECK(solve(2*x+5==10)==2.5);
     CHECK(solve(2*x==9)==4.5);
     CHECK(solve(2*x==18)==9);
     CHECK(solve(4*x==36)==9);
     CHECK(solve(2*x==100)==50);
     CHECK(solve(x+5==20 )==15);
     CHECK(solve(x-20==-20 )==0);
     CHECK(solve(x-20==-5)==15);
     CHECK(solve(x+20==5)==-15);
     CHECK(solve(x+25==5)==-20);
     CHECK(solve(2*x== -4 )==-2);
     CHECK(solve(3*x+5==11)==2);
     CHECK(solve(2*x-4.0 == 10.0)==7);
     CHECK(solve(4*x/2-4.0 == 10.0)==7);
     CHECK(solve(4*x-4.0 == 36.0)==10);
     CHECK(solve(19*x==0)==0);
     CHECK(solve(8*x/2-4.0 == 36.0)==10);
     CHECK(solve(x==5)==5);
     CHECK(solve(x==20-10)==10);
     CHECK(solve(3*x+2==14)==4);
     CHECK(solve(3*x+4==16)==4);
     CHECK(solve(x+15==-8)==-23);
     CHECK(solve(6*x/2+2==14)==4);


}
TEST_CASE("test 2"){
    RealVariable x;
            CHECK(solve((x^2)==9)==3);
            CHECK(solve(2*(x^2)==18)==3);
            CHECK(solve((x^2)+5==30)==5);
            CHECK(solve(2*x-2==4)==3);
            CHECK(((solve((x^2)==16)==4)||(solve((x^2)==16)==-4)));
            CHECK((solve((x^2)==25)==5 ||solve((x^2)==25)==-5));
            CHECK(solve((x^2)+2==18)==4);
            CHECK((solve((x^2) + 2*x + 4.0 == 20 + 6.0*x/2 - x)==4||solve((x^2) + 2*x + 4.0 == 20 + 6.0*x/2 - x)==-4));
            CHECK((solve((x^2) + 2*x == 16 + 6.0*x/2 - x)==4||solve((x^2) + 2*x == 16 + 6.0*x/2 - x)==-4));
            CHECK((solve(2*(x^2)  == 32 )==4||solve(2*(x^2)  == 32 )==-4));
            CHECK(solve(3*(x^2)-15 == 60 )==5);
            CHECK((solve((x^2)==81)==9||solve((x^2)==81)==-9));
            CHECK((solve(2*(x^2)==162)==9||solve(2*(x^2)==162)==-9));
            CHECK((solve((x^2) + 2*x+5 == 21 + 6.0*x/2 - x)==4||solve((x^2) + 2*x+5 == 21 + 6.0*x/2 - x)==-4));
            CHECK(solve(3*(x^2)-30 == 45 )==5);
            CHECK((solve((x^2)==1)==1||solve((x^2)==1)==-1));// also -1
            CHECK((solve((x^2)==100)==10||solve((x^2)==100)==-10));
}

TEST_CASE("CHECK_THROWS"){
    RealVariable x;
    CHECK_THROWS(solve((x^2)==-100));
    CHECK_THROWS( solve((x^2)==-16));
    CHECK_THROWS( solve((x^2)==-1));
    CHECK_THROWS( solve(2*(x^2)==-16));
    CHECK_THROWS( solve(2*(x^2)==-100));
    CHECK_THROWS(solve((x^2)==-25));
    CHECK_THROWS(solve(2*(x^2)==-50));
    CHECK_THROWS(solve(2*(x^2)==-7));
    CHECK_THROWS(solve((x^2)==-9));
    CHECK_THROWS(solve((x^2)==-12));

}
TEST_CASE("test ComplexVariable"){

    ComplexVariable y;
    CHECK(solve(2*y-4 == 10)==std::complex<double>(7,0));
    CHECK(solve(y-4 == 10)==std::complex<double>(14,0));
    CHECK((solve((y^2) == 16)==std::complex<double>(4,0)||solve((y^2) == 16)==std::complex<double>(-4,0)));
    CHECK((solve((y^2) == -16)==std::complex<double>(0,4)||solve((y^2) == -16)==std::complex<double>(0,-4)));
    CHECK((solve((y^2) + 2*y + 4 == 20 + 6*y/2 - y)==std::complex<double>(4,0)||solve((y^2) + 2*y + 4 == 20 + 6*y/2 - y)==std::complex<double>(-4,0))); //(-4,0)
    CHECK((solve(y+5i == 2*y+3i)==std::complex<double>(0,2)||solve(y+5i == 2*y+3i)==std::complex<double>(0,-2)));
    CHECK((solve((y^2)==-100)==std::complex<double>(0,10)||solve((y^2)==-100)==std::complex<double>(0,-10)));
    CHECK((solve(2*(y^2)==-200)==std::complex<double>(0,10)||solve(2*(y^2)==-200)==std::complex<double>(0,-10)));
    CHECK((solve((y^2)==100)==std::complex<double>(10,0)||solve((y^2)==100)==std::complex<double>(-10,0)));
    CHECK((solve((y^2)==-25)==std::complex<double>(0,5)||solve((y^2)==-25)==std::complex<double>(0,-5)));
    CHECK((solve(2*(y^2)==-50)==std::complex<double>(0,5)||solve(2*(y^2)==-50)==std::complex<double>(0,-5)));
    CHECK((solve((y^2)==-9)==std::complex<double>(0,3)||solve((y^2)==-9)==std::complex<double>(0,-3)));
    CHECK((solve((y^2)==-81)==std::complex<double>(0,9)||solve((y^2)==-81)==std::complex<double>(0,-9)));
    CHECK((solve((y^2) + 3*y + 4 == 20 + 6*y/2 )==std::complex<double>(4,0)||solve((y^2) + 3*y + 4 == 20 + 6*y/2 )==std::complex<double>(-4,0))); //(-4,0)
    CHECK((solve((y^2) + 3*y == 16 + 3*y )==std::complex<double>(4,0)||solve((y^2) + 3*y == 16 + 3*y )==std::complex<double>(-4,0))); //(-4,0)
    CHECK((solve((y^2) + 6*y/2 == 16 + 9*y/3 )==std::complex<double>(4,0)||solve((y^2) + 6*y/2 == 16 + 9*y/3 )==std::complex<double>(-4,0))); //(-4,0)
    CHECK((solve((y^2)==9)==std::complex<double>(3,0)||solve((y^2)==9)==std::complex<double>(-3,0)));  //-3
    CHECK((solve((y^2)+5==30)==std::complex<double>(5,0)||solve((y^2)+5==30)==std::complex<double>(-5,0))); //-5
    CHECK((solve(2*y-2==4)==std::complex<double>(3,0)||solve(2*y-2==4)==std::complex<double>(-3,0)));   //-3
    CHECK((solve((y^2) == 16)==std::complex<double>(4,0)||solve((y^2) == 16)==std::complex<double>(-4,0))); //(-4,0)
    CHECK((solve((y^2)==25)==std::complex<double>(5,0)||solve((y^2)==25)==std::complex<double>(-5,0)));  //-5
    CHECK((solve((y^2)+2==18)==std::complex<double>(4,0)||solve((y^2)+2==18)==std::complex<double>(-4,0)));  //-4



}
TEST_CASE("CHECK_THROWS ComplexVariable ") {
    ComplexVariable y;
    CHECK_THROWS(solve((y ^ 0)==-16));
}
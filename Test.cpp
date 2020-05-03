#include "doctest.h"
#include <string.h>
#include "solver.hpp"
#include "math.h"
#include <complex>
#define EP 0.0001
using namespace solver;
using namespace std;
TEST_CASE("Check the correctness of addition (operator +)") {

    SUBCASE("on RealVariable")
    {
        RealVariable x;
        CHECK_LE(abs(solve(35+x == 5)+30), EP);//-30
        CHECK_LE(abs(solve(x+35 == 5)+30), EP);//-30
        CHECK_LE(abs(solve(5+x == 35)-30), EP);//30
        CHECK_LE(abs(solve(x+5 == 35)-30), EP);//30
        CHECK_LE(abs(solve(x+5+3+8+10 == 35)-9), EP);//9
        CHECK_LE(abs(solve(5+x+3+8+10 == 35)-9), EP);//9
        CHECK_LE(abs(solve(5+3+8+10+x == 35)-9), EP);//9
        CHECK_LE(abs(solve(5+3+8+10+x == 35+x+x+x+x)+3), EP);//-3
        CHECK_LE(abs(solve(5+3+8+10+x == x+35+x+x+x)+3), EP);//-3
        CHECK_LE(abs(solve(5+3+x+8+10 == x+x+x+x+35)+3), EP);//-3
        CHECK_LE(abs(solve(5+5+5+5+5+5 == x+x+x+x+x+x+x+x+x+x+x+x+x+x+x)-2), EP);//2
        CHECK_LE(abs(solve(x+x+x+x+x+x+x+x+x+x+x+x+x+x+x == 5+5+5+5+5+5)-2), EP);//2
    }

    SUBCASE("on ComplexVariable")
    {
        ComplexVariable x;

        CHECK_LE(abs(solve(35+x == 5).real()+30), EP);//-30+0i
        CHECK_LE(abs(solve(x+35 == 5).real()+30), EP);//-30+0i
        CHECK_LE(abs(solve(5+x == 35).real()-30), EP);//30+0i
        CHECK_LE(abs(solve(x+5 == 35).real()-30), EP);//30+0i
        CHECK_LE(abs(solve(x+5+3+8+10 == 35).real()-9), EP);//9+0i
        CHECK_LE(abs(solve(5+x+3+8+10 == 35).real()-9), EP);//9+0i
        CHECK_LE(abs(solve(5+3+8+10+x == 35).real()-9), EP);//9+0i
        CHECK_LE(abs(solve(5+3+8+10+x == 35+x+x+x+x).real()+3), EP);//-3+0i
        CHECK_LE(abs(solve(5+3+8+10+x == x+35+x+x+x).real()+3), EP);//-3+0i
        CHECK_LE(abs(solve(5+3+x+8+10 == x+x+x+x+35).real()+3), EP);//-3+0i
        CHECK_LE(abs(solve(5+5+5+5+5+5 == x+x+x+x+x+x+x+x+x+x+x+x+x+x+x).real()-2), EP);//2+0i
        CHECK_LE(abs(solve(x+x+x+x+x+x+x+x+x+x+x+x+x+x+x == 5+5+5+5+5+5).real()-2), EP);//2+0i

        CHECK_LE(abs(solve(35+x == 5).real()+30), EP);//-30+0i
        CHECK_LE(abs(solve(35+x == 5).imag()+0), EP);//-30+0i
        CHECK_LE(abs(solve(x+3i+5 == 5).real()+0), EP);//0-3i
        CHECK_LE(abs(solve(x+3i+5 == 5).imag()+3), EP);//0-3i
        CHECK_LE(abs(solve(2i+3i == x).real()+0), EP);//5i
        CHECK_LE(abs(solve(2i+3i == x).imag()-5), EP);//5i
        CHECK_LE(abs(solve(2i+3i+2i+x+x+x+x+x+10i+x+3i == x+x+19i+x+0i+x).real()+0), EP);//x = 0-i
        CHECK_LE(abs(solve(0i+1i+5i+x+3+2i+61+3i == 7i+3i+x+x+3i+x+1).imag()+1), EP);//30-i
        CHECK_LE(abs(solve(2i+3i+3.0 == x).imag()-5), EP);//5i
    }

}

TEST_CASE("Check the correctness of subtraction (operator -)") {

    SUBCASE("on RealVariable")
    {
        RealVariable x;
        CHECK_LE(abs(solve(35-x == 5)-30), EP);//30
        CHECK_LE(abs(solve(x-35 == 5)-40), EP);//40
        CHECK_LE(abs(solve(5-x == 35)+30), EP);//-30
        CHECK_LE(abs(solve(x-5 == 35)-40), EP);//40
        CHECK_LE(abs(solve(x-5-3-8-10 == 35)-61), EP);//61
        CHECK_LE(abs(solve(5-x-3-8-10 == 35)+51), EP);//-51
        CHECK_LE(abs(solve(5-3-8-10-x == 35)+51), EP);//-51
        CHECK_LE(abs(solve(5-3-8-10-x == 35-x-x-x-x)-17), EP);//17
        CHECK_LE(abs(solve(5-3-8-10-x == x-35-x-x-x)+19), EP);//-19
        CHECK_LE(abs(solve(5-3-x-8-10 == x-x-x-x-35)+19), EP);//-19
        CHECK_LE(abs(solve(5-5-5-5-5-5-5-5 == x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x)-2), EP);//2
        CHECK_LE(abs(solve(x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x == 5-5-5-5-5-5-5-5)-2), EP);//2

    }

    SUBCASE("on ComplexVariable")
    {
        ComplexVariable y;

        CHECK_LE(abs(solve(35-y == 5).real()-30), EP);//30+0i
        CHECK_LE(abs(solve(y-35 == 5).real()-40), EP);//40+0i
        CHECK_LE(abs(solve(5-y == 35).real()+30), EP);//-30+0i
        CHECK_LE(abs(solve(y-5 == 35).real()-40), EP);//40+0i
        CHECK_LE(abs(solve(y-5-3-8-10 == 35).real()-61), EP);//61+0i
        CHECK_LE(abs(solve(5-y-3-8-10 == 35).real()+51), EP);//-51+0i
        CHECK_LE(abs(solve(5-3-8-10-y == 35).real()+51), EP);//-51+0i
        CHECK_LE(abs(solve(5-3-8-10-y == 35-y-y-y-y).real()-17), EP);//17+0i
        CHECK_LE(abs(solve(5-3-8-10-y == y-35-y-y-y).real()+19), EP);//-19+0i
        CHECK_LE(abs(solve(5-3-y-8-10 == y-y-y-y-35).real()+19), EP);//-19+0i
        CHECK_LE(abs(solve(5-5-5-5-5-5-5-5 == y-y-y-y-y-y-y-y-y-y-y-y-y-y-y-y-y).real()-2), EP);//2+0i
        CHECK_LE(abs(solve(y-y-y-y-y-y-y-y-y-y-y-y-y-y-y-y-y == 5-5-5-5-5-5-5-5).real()-2), EP);//2+0i

        CHECK_LE(abs(solve(35-y == 5).real()-30), EP);//30+0i
        CHECK_LE(abs(solve(35-y == 5).imag()-0), EP);// 30+0i
        CHECK_LE(abs(solve(y-3i-5 == 5).real()-10), EP);//10+3i
        CHECK_LE(abs(solve(y-3i-5 == 5).imag()-3), EP);// 10+3i
        CHECK_LE(abs(solve(2i-3i == y).real()+0), EP);//0-i
        CHECK_LE(abs(solve(2i-3i == y).imag()+1), EP);//0-i
        CHECK_LE(abs(solve(2i-3i-2i-y-y-11-y-y-y-10i-y-3i == 7-y-y-19i-y-0i-y).real()+9), EP);//-9+(3/2)i

    }

}

TEST_CASE("Check the correctness of multiplication (operator *)") {

    SUBCASE("on RealVariable")
    {
        RealVariable y;
        CHECK_LE(abs(solve(2*y == 6)-3), EP);//3
        CHECK_LE(abs(solve(y*2 == 6)-3), EP);//3
        CHECK_LE(abs(solve(6 == 2*y)-3), EP);//3
        CHECK_LE(abs(solve(6 == y*2)-3), EP);//3
        CHECK_LE(abs(solve(6*y == 6)-1), EP);//1
        CHECK_LE(abs(solve(y*6 == 6)-1), EP);//1
        CHECK_LE(abs(solve(6 == 6*y)-1), EP);//1
        CHECK_LE(abs(solve(6 == y*6)-1), EP);//1
        CHECK_LE(abs(solve(6*y == 0)-0), EP);//0
        CHECK_LE(abs(solve(0 == 6*y)-0), EP);//0
        CHECK_LE(abs(solve(3*y*2 == 6)-1), EP);//1
        CHECK_LE(abs(solve(3*y*2 == 6*4)-4), EP);//4
        CHECK_LE(abs(solve(3*y*2*15 == 6*4*90)-24), EP);//24
        CHECK_LE(abs(solve(3*y*2*15 == 6*4*90)-24), EP);//24
        CHECK_LE(abs(solve(y*1*2*3*4*5*6 == 4320)-6), EP);//6
        CHECK_LE(abs(solve(1*y*2*3*4*5*6 == 4320)-6), EP);//6
        CHECK_LE(abs(solve(1*2*y*3*4*5*6 == 4320)-6), EP);//6
        CHECK_LE(abs(solve(1*2*3*4*5*6*y == 4320)-6), EP);//6
    }

    SUBCASE("on ComplexVariable")
    {
        ComplexVariable y;
        CHECK_LE(abs(solve(2*y == 6).real()-3), EP);//3+0i
        CHECK_LE(abs(solve(y*2 == 6).real()-3), EP);//3+0i
        CHECK_LE(abs(solve(6 == 2*y).real()-3), EP);//3+0i
        CHECK_LE(abs(solve(6 == y*2).real()-3), EP);//3+0i
        CHECK_LE(abs(solve(6*y == 6).real()-1), EP);//1+0i
        CHECK_LE(abs(solve(y*6 == 6).real()-1), EP);//1+0i
        CHECK_LE(abs(solve(6 == 6*y).real()-1), EP);//1+0i
        CHECK_LE(abs(solve(6 == y*6).real()-1), EP);//1+0i
        CHECK_LE(abs(solve(6*y == 0).real()-0), EP);//0+0i
        CHECK_LE(abs(solve(0 == 6*y).real()-0), EP);//0+0i
        CHECK_LE(abs(solve(3*y*2 == 6).real()-1), EP);//1+0i
        CHECK_LE(abs(solve(3*y*2 == 6*4).real()-4), EP);//4+0i
        CHECK_LE(abs(solve(3*y*2*15 == 6*4*90).real()-24), EP);//24+0i
        CHECK_LE(abs(solve(3*y*2*15 == 6*4*90).real()-24), EP);//24+0i
        CHECK_LE(abs(solve(y*1*2*3*4*5*6 == 4320).real()-6), EP);//6+0i
        CHECK_LE(abs(solve(1*y*2*3*4*5*6 == 4320).real()-6), EP);//6+0i
        CHECK_LE(abs(solve(1*2*y*3*4*5*6 == 4320).real()-6), EP);//6+0i

        CHECK_LE(abs(solve(4*y == 2i).real()-0), EP);//0+(1/2)i
        CHECK_LE(abs(solve(4*y*0i*4*5i == y).real()-0), EP);//0+0i
        CHECK_LE(abs(solve(4*y*0i*4*5i == y).real()-0), EP);//0+0i

    }

}

TEST_CASE("Check the correctness of divition (operator /)") {

    SUBCASE("on RealVariable")
    {
        RealVariable z;
        CHECK_LE(abs(solve(z/2 == 6)-12), EP);//12
        CHECK_LE(abs(solve(z/10 == 6)-60), EP);//60
        CHECK_LE(abs(solve(6 == z/10)-60), EP);//60
        CHECK_LE(abs(solve(z/z == z)-1), EP);//1
        CHECK_LE(abs(solve(z/2/5 == 6)-60), EP);//60
        CHECK_LE(abs(solve(z/2/5/3 == 6)-180), EP);//180
        CHECK_LE(abs(solve(6/2 == z/2/5/3)-90), EP);//90
        CHECK_LE(abs(solve(z/2/5/3 == 12/2/2/1/1/1/1)-90), EP);//90
        CHECK_LE(abs(solve((z^2)/z == 12/2/2/1/1/1/1)-3), EP);//3
    }

    SUBCASE("on ComplexVariable")
    {
        ComplexVariable z;
        CHECK_LE(abs(solve(z/2 == 6).real()-12), EP);//12+0i
        CHECK_LE(abs(solve(z/10 == 6).real()-60), EP);//60+0i
        CHECK_LE(abs(solve(6 == z/10).real()-60), EP);//60+0i
        CHECK_LE(abs(solve(z/z == z).real()-1), EP);//1+0i
        CHECK_LE(abs(solve(z/2/5 == 6).real()-60), EP);//60+0i
        CHECK_LE(abs(solve(z/2/5/3 == 6).real()-180), EP);//180+0i
        CHECK_LE(abs(solve(6/2 == z/2/5/3).real()-90), EP);//90+0i
        CHECK_LE(abs(solve(z/2/5/3 == 12/2/2/1/1/1/1).real()-90), EP);//90+0i
        CHECK_LE(abs(solve((z^2)/z == 12/2/2/1/1/1/1).real()-3), EP);//3+0i

        CHECK_LE(abs(solve(1i/1i == z).real()-1), EP);//1+0i
        CHECK_LE(abs(solve(1i/1i == z).imag()-0), EP);//1+0i

    }

}

TEST_CASE("Check the correctness of power (operator ^)") {

    SUBCASE("on RealVariable")
    {
        RealVariable x;
        CHECK_LE(abs(solve((x^1) == 6)-6), EP);//6
        CHECK_LE(abs(solve((x^0) == x)-1), EP);//1
        CHECK(((abs(solve((x^2) == 64)-8) <= EP) || (abs(solve((x^2) == 64)+8) <= EP)));//8 or -8
        CHECK_LE(abs(solve((x^2) == 0)-0), EP);//0
        CHECK_LE(abs(solve(0 == (x^2))-0), EP);//0
        CHECK(((abs(solve((((x^1)^1)^2) == 4)-2) <= EP) || (abs(solve((((x^1)^1)^2) == 4)+2) <= EP)));//2 or -2
    }

    SUBCASE("on ComplexVariable")
    {
        ComplexVariable x;
        CHECK_LE(abs(solve((x^1) == 6).real()-6), EP);//6+0i
        CHECK_LE(abs(solve((x^0) == x).real()-1), EP);//1+0i
        CHECK(((abs(solve((x^2) == 64).real()-8) <= EP) || (abs(solve((x^2) == 64).real()+8) <= EP)));//8+0i or -8+0i
        CHECK_LE(abs(solve((x^2) == 0).real()-0), EP);//0+0i
        CHECK_LE(abs(solve(0 == (x^2)).real()-0), EP);//0+0i
        CHECK(((abs(solve((((x^1)^1)^2) == 4).real()-2) <= EP) || (abs(solve((((x^1)^1)^2) == 4).real()+2) <= EP)));//2+0i or -2+0i
    }

}

TEST_CASE("Check the correctness of equation 1") {

    SUBCASE("on RealVariable")
    {
        RealVariable x;
        CHECK_LE(abs(solve(x-x+x-x+2*x+6-x*2+x == x-x+x-x+2*x+6-x*2)-0), EP);//0
        CHECK_LE(abs(solve(6-x == x-x+x-x+2*x+6-x*2)-0), EP);//0
    }

    SUBCASE("on ComplexVariable")
    {
        ComplexVariable x;
        CHECK_LE(abs(solve(x-x+x-x+2*x+6-x*2+x == x-x+x-x+2*x+6-x*2).real()-0), EP);//0+0i
        CHECK_LE(abs(solve(6-x == x-x+x-x+2*x+6-x*2).real()-0), EP);//0+0i
    }
    
}

TEST_CASE("Check the correctness of equation 2") {

    SUBCASE("on RealVariable")
    {
        RealVariable x;
        CHECK(((abs(solve(x*x-9 == 0)-3) <= EP) || (abs(solve(x*x-9 == 0)+3) <= EP)));//3 or -3
        CHECK_LE(abs(solve(x*x/x+27-4*x == 0/2*3*x*x)-9), EP);//9
        CHECK(((abs(solve(3+2*x*x+x*x+10*x+12*x+3+10+2+39 == x*x-3)+6) <= EP) || (abs(solve(3+2*x*x+x*x+10*x+12*x+3+10+2+39 == x*x-3)+5) <= EP)));//-5 or -6
    }

    SUBCASE("on ComplexVariable")
    {
        ComplexVariable x;
        CHECK(((abs(solve(x*x-9 == 0).real()-3) <= EP) || (abs(solve(x*x-9 == 0).real()+3) <= EP)));//3+0i or -3+0i
        CHECK_LE(abs(solve(x*x/x+27-4*x == 0/2*3*x*x).real()-9), EP);//9
        CHECK(((abs(solve(3+2*x*x+x*x+10*x+12*x+3+10+2+39 == x*x-3).real()+6) <= EP) || (abs(solve(3+2*x*x+x*x+10*x+12*x+3+10+2+39 == x*x-3).real()+5) <= EP)));//-5+0i or -6+0i
        CHECK_LE(abs(solve(x+5i == 2*x+3i+10).real()+10), EP);//-10+2i
    }
}

TEST_CASE("Check the correctness of trow exception") {
    RealVariable x;
    CHECK_THROWS(solve(((x + 2 - x - 2)^0) == x));//0 power 0 not illegal
    CHECK_THROWS(solve(((x + 2 - x - 2)^-1) == x));//0 power minus not illegal
    CHECK_THROWS(solve(((x + 2 - x - 2)^-10) == x));//0 power minus not illegal
    CHECK_THROWS(solve((x^2) == -16));
    CHECK_THROWS(solve(x-x == x-x+2));
}

TEST_CASE("Check the correctness of trow exception") {
    ComplexVariable x;
    CHECK_THROWS(solve(((x + 2 - x - 2)^0) == x));//0 power 0 not illegal
    CHECK_THROWS(solve(((x + 2 - x - 2)^-1) == x));//0 power minus not illegal
    CHECK_THROWS(solve(((x + 2 - x - 2)^-10) == x));//0 power minus not illegal
    CHECK_THROWS(solve(x-x == x-x+2));//Not true
}










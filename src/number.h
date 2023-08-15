#ifndef NUMBER_H
#define NUMBER_H

#include <cmath>
#include <iostream>

// Number = sector + precistion
// sector is integer part
// preciistion is floating point part limited by treshold
class Number{
public:
    Number() {};
    Number(float num) : m_precision(num) { normalize(); }
    Number(double num){
        double sector = std::floor( num / (double) c_treshold );
        m_precision = (num - sector);
        m_sector = sector;
    }
    Number(long long sector, float precistion) : m_sector(sector), m_precision(precistion) { normalize(); }
    Number(const Number& num): m_sector(num.m_sector), m_precision(num.m_precision) {};
    ~Number() {};

    long long m_sector{0};
    float m_precision{0.0};
    const float c_treshold{1000.0};

    inline void normalize(){
        float q = std::floor( m_precision / c_treshold );
        m_sector += q;
        m_precision -= q * c_treshold;
        m_sector -= 1 * ( m_precision < 0 );
        m_precision += 1.0 * ( m_precision < 0);
    }

    Number operator=(Number& n){
        m_sector = n.m_sector;
        m_precision = n.m_precision;
        return *this;
    }

    Number operator+(Number& n){
        return Number(m_sector + n.m_sector,
            m_precision + n.m_precision);
    }

    Number operator-(Number& n){
        return Number(m_sector - n.m_sector,
            m_precision - n.m_precision);
    }

    Number operator+=(Number& n){
        m_sector += n.m_sector;
        m_precision += n.m_precision;
        normalize();
        return *this;
    }

    Number operator-=(Number& n){
        m_sector -= n.m_sector;
        m_precision -= n.m_precision;
        normalize();
        return *this;
    }

    Number operator*(double d){
        return Number((double)m_sector * d,
                    (double)m_precision * d);
    }

    Number operator*(Number& n){
        long long sector = m_sector * n.m_sector;
        Number n1((double)m_sector * (double)n.m_precision);
        Number n2((double)n.m_sector * (double)m_precision);
        float precision = m_precision * n.m_precision;

        return Number(sector + n1.m_sector + n2.m_sector,
                    precision + n1.m_precision + n2.m_precision);
    }

     Number operator*=(Number& n){
        Number n1 = *this * n;
        m_sector = n1.m_sector;
        m_precision = n.m_precision;
        return *this;
    }

    Number operator/(Number& n){
        double div = ((double) n.m_sector + (double) n.m_precision);
        return Number(m_sector / div, (double)m_precision / div); 
    }

    Number operator/=(Number& n){
        Number n1 = *this/n;
        m_sector = n1.m_sector;
        m_precision = n1.m_precision;
        return *this;
    }
   
    // Number operator%(Number& n);
    // Number operator%=(Number& n);
    
    Number operator==(Number& n);
    Number operator!=(Number& n);
    Number operator<(Number& n);
    Number operator>(Number& n);
    Number operator<=(Number& n);
    Number operator>=(Number& n); 

    std::string pretty_string();    
};

std::ostream & operator<< (std::ostream& os, Number & number){
    os << "Number( "<<number.m_sector<<" + "<<number.m_precision<<" )";
    return os;
}

#endif // NUMBER_H
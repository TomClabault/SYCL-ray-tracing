
#include <algorithm>

#include "color.h"

Color& Color::operator=(const Vector& vec)
{
    r = vec.x;
    g = vec.y;
    b = vec.z;

    return *this;
}

Color& Color::operator+=(const Color &other)
{
    r += other.r;
    g += other.g;
    b += other.b;

    return *this;
}

Color& Color::operator*=(const Color &other)
{
    r *= other.r;
    g *= other.g;
    b *= other.b;

    return *this;
}

Color& Color::operator*=(float k)
{
    r *= k;
    g *= k;
    b *= k;

    return *this;
}

Color& Color::operator/=(const float k)
{
    r /= k;
    g /= k;
    b /= k;

    return *this;
}

float Color::power( ) const
{
    return (r+g+b) / 3;
}

float Color::max( ) const
{
    return std::max(r, std::max(g, std::max(b, float(0))));
}

Color Black( )
{
    return Color(0, 0, 0);
}

Color White( )
{
    return Color(1, 1, 1);
}

Color Red( )
{
    return Color(1, 0, 0);
}

Color Green( )
{
    return Color(0, 1, 0);
}

Color Blue( )
{
    return Color(0, 0, 1);
}

Color Yellow( )
{
    return Color(1, 1, 0);
}


Color operator+ ( const Color& a, const Color& b )
{
    return Color(a.r + b.r, a.g + b.g, a.b + b.b, a.a + b.a);
}

Color operator- ( const Color& c )
{
    return Color(-c.r, -c.g, -c.b, -c.a);
}

Color operator- ( const Color& a, const Color& b )
{
    return a + (-b);
}

Color operator* ( const Color& a, const Color& b )
{
    return Color(a.r * b.r, a.g * b.g, a.b * b.b, a.a * b.a);
}

Color operator* ( const float k, const Color& c )
{
    return Color(c.r * k, c.g * k, c.b * k, c.a * k);
}

Color operator* ( const Color& c, const float k )
{
    return k * c;
}

Color operator/ ( const Color& a, const Color& b )
{
    return Color(a.r / b.r, a.g / b.g, a.b / b.b, a.a / b.a);
}

Color operator/ ( const float k, const Color& c )
{
    return Color(k / c.r, k / c.g, k / c.b, k / c.a);
}

Color operator/ ( const Color& c, const float k )
{
    float kk= 1 / k;
    return kk * c;
}

std::ostream& operator << (std::ostream& os, const Color& color)
{
    os << "Color[" << color.r << ", " << color.g << ", " << color.b << "]";
    return os;
}

Color exp(const Color &col)
{
    return Color(std::exp(col.r), std::exp(col.g), std::exp(col.b), col.a);
}

Color pow(const Color &col, float k)
{
    return Color(std::pow(col.r, k), std::pow(col.g, k), std::pow(col.b, k), col.a);
}


#ifndef _MAT_H
#define _MAT_H

#include "vec.h"

//! \addtogroup math manipulations de points, vecteur, matrices, transformations
///@{

//! \file
//! transformation de points et vecteurs

//! conversion en radians.
SYCL_EXTERNAL float radians(const float deg);
//! conversion en degres.
SYCL_EXTERNAL float degrees(const float rad);

//! representation d'une transformation, une matrice 4x4, organisee par ligne / row major.
struct Transform
{
    //! constructeur.
    SYCL_EXTERNAL Transform(
        const float t00 = 1.0f, const float t01 = 0.0f, const float t02 = 0.0f, const float t03 = 0.0f,
        const float t10 = 0.0f, const float t11 = 1.0f, const float t12 = 0.0f, const float t13 = 0.0f,
        const float t20 = 0.0f, const float t21 = 0.0f, const float t22 = 1.0f, const float t23 = 0.0f,
        const float t30 = 0.0f, const float t31 = 0.0f, const float t32 = 0.0f, const float t33 = 1.0f)
    {
        m[0][0]= t00; m[0][1]= t01; m[0][2]= t02; m[0][3]= t03;
        m[1][0]= t10; m[1][1]= t11; m[1][2]= t12; m[1][3]= t13;
        m[2][0]= t20; m[2][1]= t21; m[2][2]= t22; m[2][3]= t23;
        m[3][0]= t30; m[3][1]= t31; m[3][2]= t32; m[3][3]= t33;
    }

    //! constructeur a partir de 4 Vector colonnes, met (0, 0, 0, 1) dans la derniere ligne.
    SYCL_EXTERNAL Transform(const Vector& x, const Vector& y, const Vector& z, const Vector& w)
    {
        m[0][0] = x.x;	m[0][1] = y.x;	m[0][2] = z.x;	m[0][3] = w.x;
        m[1][0] = x.y;	m[1][1] = y.y;	m[1][2] = z.y;	m[1][3] = w.y;
        m[2][0] = x.z;	m[2][1] = y.z;	m[2][2] = z.z;	m[2][3] = w.z;
        m[3][0] = 0;	m[3][1] = 0;	m[3][2] = 0;	m[3][3] = 1;
    }
    //! constructeur a partir de 4 colonnes
    Transform(const vec4& x, const vec4& y, const vec4& z, const vec4& w);

    //! initialise une colonne de la matrice a partir de 4 floats.
    SYCL_EXTERNAL Transform& column(const unsigned id, const float t0, const float t1, const float t2, const float t3);
    //!renvoie une colonne.
    SYCL_EXTERNAL vec4 column(const unsigned id) const;
    //!renvoie une colonne.
    SYCL_EXTERNAL vec4 column(const unsigned id);

    //! initialise une ligne de la matrice.
    SYCL_EXTERNAL Transform& row(const unsigned id, const float t0, const float t1, const float t2, const float t3);
    //!renvoie une ligne.
    SYCL_EXTERNAL vec4 row(const unsigned id) const;
    //!renvoie une ligne.
    SYCL_EXTERNAL vec4 row(const unsigned id);

    //! initialise la matrice avec 16 floats organises par colonne.
    SYCL_EXTERNAL Transform& column_major(const float matrix[16]);

    //! initialise la matrice avec 16 floats organises par ligne.
    SYCL_EXTERNAL Transform& row_major(const float matrix[16]);

    //! renvoie le Vector colonne c de la matrice
    SYCL_EXTERNAL Vector operator[] (const unsigned c) const;

    //! renvoie le point transforme.
    SYCL_EXTERNAL Point operator() (const Point& p) const;
    //! renvoie le vecteur transforme.
    SYCL_EXTERNAL Vector operator() (const Vector& v) const;
    //! renvoie le point/vecteur homogene transforme.
    SYCL_EXTERNAL vec4 operator() (const vec4& v) const;

    //! renvoie la composition de la transformation this et b, t = this * b. permet de transformer un point sans "ambiguite" Point q= a(b(c(p)));
    SYCL_EXTERNAL Transform operator() (const Transform& b) const;

    //! renvoie la transposee de la matrice.
    SYCL_EXTERNAL Transform transpose() const;
    //! renvoie l'inverse de la matrice.
    SYCL_EXTERNAL Transform inverse() const;
    //! renvoie la transformation a appliquer aux normales d'un objet transforme par la matrice m.
    SYCL_EXTERNAL Transform normal() const;

    //! renvoie l'adresse de la premiere valeur de la matrice.
    SYCL_EXTERNAL const float* data() const { return &m[0][0]; }

    float m[4][4];
};

//! construit la transformation identite.
SYCL_EXTERNAL Transform Identity();

//! renvoie la transposee de la matrice.
SYCL_EXTERNAL Transform Transpose(const Transform& m);
//! renvoie l'inverse de la matrice.
SYCL_EXTERNAL Transform Inverse(const Transform& m);
//! renvoie la transformation a appliquer aux normales d'un objet transforme par la matrice m.
SYCL_EXTERNAL Transform Normal(const Transform& m);

//! renvoie la matrice representant une mise a l'echelle / etirement.
SYCL_EXTERNAL Transform Scale(const float x, const float y, const float z);
SYCL_EXTERNAL inline Transform Scale(const float s) { return Scale(s, s, s); }

//! renvoie la matrice representant une translation par un vecteur.
SYCL_EXTERNAL Transform Translation(const Vector& v);
//! renvoie la matrice representant une translation par un vecteur x y z.
SYCL_EXTERNAL Transform Translation(const float x, const float y, const float z);

//! renvoie la matrice representation une rotation de angle degree autour de l'axe X.
SYCL_EXTERNAL Transform RotationX(const float angle);
//! renvoie la matrice representation une rotation de a degree autour de l'axe Y.
SYCL_EXTERNAL Transform RotationY(const float angle);
//! renvoie la matrice representation une rotation de angle degree autour de l'axe Z.
SYCL_EXTERNAL Transform RotationZ(const float angle);
//! renvoie la matrice representation une rotation de angle degree autour de l'axe axis.
SYCL_EXTERNAL Transform Rotation(const Vector& axis, const float angle);

//! renvoie la matrice de rotation entre u et v.
SYCL_EXTERNAL Transform Rotation(const Vector& u, const Vector& v);

//! renvoie la matrice representant une transformation viewport.
SYCL_EXTERNAL Transform Viewport(const float width, const float height);
//! renvoie la matrice representant une transformation projection perspective.
SYCL_EXTERNAL Transform Perspective(const float fov, const float aspect, const float znear, const float zfar);
//! renvoie la matrice representant une transformation orthographique, passage d'un cube []x[]x[] vers [-1 1]x[-1 1]x[-1 1].
SYCL_EXTERNAL Transform Ortho(const float left, const float right, const float bottom, const float top, const float znear, const float zfar);
//! renvoie la matrice representant le placement et l'orientation d'une camera pour observer le point to.
SYCL_EXTERNAL Transform Lookat(const Point& from, const Point& to, const Vector& up);

//! renvoie la composition des transformations a et b, t= a * b.
SYCL_EXTERNAL Transform compose_transform(const Transform& a, const Transform& b);
//! renvoie la composition des transformations a et b, t = a * b.
SYCL_EXTERNAL Transform operator* (const Transform& a, const Transform& b);

#include <iostream>

inline std::ostream& operator<<(std::ostream& o, const Transform& t)
{
    o << t.m[0][0] << " " << t.m[0][1] << " " << t.m[0][2] << " " << t.m[0][3] << " " << std::endl;
    o << t.m[1][0] << " " << t.m[1][1] << " " << t.m[1][2] << " " << t.m[1][3] << " " << std::endl;
    o << t.m[2][0] << " " << t.m[2][1] << " " << t.m[2][2] << " " << t.m[2][3] << " " << std::endl;
    o << t.m[3][0] << " " << t.m[3][1] << " " << t.m[3][2] << " " << t.m[3][3] << " " << std::endl;
    return o;
}


///@}
#endif

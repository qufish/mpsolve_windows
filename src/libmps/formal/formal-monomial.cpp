#include <mps/mps.h>
#include <iostream>
#include <sstream>

using namespace mps::formal;

mps_formal_monomial*
mps_formal_monomial_new_with_string(const char* coeff_string, long degree)
{
    Monomial temp(coeff_string, degree);
    mps_new_obj(Monomial, m, sizeof(Monomial));
    // that may be a memory leak
    *m = temp;
    return (mps_formal_monomial*)m;
}

mps_formal_monomial*
mps_formal_monomial_new_with_strings(const char* real, const char* imag,
    long degree)
{
    Monomial temp(real, imag, degree);
    mps_new_obj(Monomial, m, sizeof(Monomial));
    *m = temp;
    return (mps_formal_monomial*)m;
}

void
mps_formal_monomial_free(mps_formal_monomial* m)
{
    mps_del_obj(m);
}

void
mps_formal_monomial_print(mps_formal_monomial* m)
{
    char* buf;
    buf = get_string(*((Monomial*)m));
    std::cout << buf;
    mps_free(buf);
}

mps_formal_monomial*
mps_formal_monomial_neg(mps_formal_monomial* m)
{
    mps_new_obj(Monomial, m2, sizeof(Monomial));
    *m2 = -*((Monomial*)m);
    return (mps_formal_monomial*)m2;
}

mps_formal_monomial*
mps_formal_monomial_mul_eq(mps_formal_monomial* m,
    mps_formal_monomial* other)
{
    Monomial* mm = (Monomial*)m;
    *mm *= *((Monomial*)other);
    return (mps_formal_monomial*)mm;
}

mps_formal_monomial*
mps_formal_monomial_mul(mps_formal_monomial* m,
    mps_formal_monomial* other)
{
    mps_new_obj(Monomial, result, sizeof(Monomial));
    *result = *((Monomial*)m) * *((Monomial*)other);
    return (mps_formal_monomial*)(result);
}

char*
mps_formal_monomial_get_str(mps_formal_monomial* m)
{
    char* buf;
    buf = get_string(*((Monomial*)m));
    return buf;
}

int
mps_formal_monomial_degree(mps_formal_monomial* m)
{
    return ((Monomial*)m)->degree();
}

Monomial::Monomial()
{
    mCoeffR = 0;
    mCoeffI = 0;
    mDegree = 0;
}

Monomial::Monomial(const char* coeff_string, long degree)
{
    char* er = mps_utils_build_equivalent_rational_string(NULL, coeff_string);

    mCoeffR = er;
    mDegree = degree;
    mps_free(er);

    mCoeffR.canonicalize();
}

Monomial::Monomial(const char* real_part, const char* imag_part, long degree)
{
    char* er = mps_utils_build_equivalent_rational_string(NULL, real_part);
    char* ei = mps_utils_build_equivalent_rational_string(NULL, imag_part);

    mDegree = degree;

    mCoeffR = er;
    mCoeffI = ei;

    mps_free(er);
    mps_free(ei);

    mCoeffR.canonicalize();
    mCoeffI.canonicalize();
}

Monomial::Monomial(const mpq_class coeff, long degree)
{
    mCoeffR = coeff;
    mCoeffR.canonicalize();
    mDegree = degree;
}

Monomial::Monomial(const mpq_class realpart, const mpq_class imagpart, long degree)
{
    mCoeffR = realpart;
    mCoeffI = imagpart;

    mCoeffR.canonicalize();
    mCoeffI.canonicalize();

    mDegree = degree;
}

Monomial::Monomial(const Monomial & rhs)
{
    mCoeffR = rhs.coefficientReal();
    mCoeffI = rhs.coefficientImag();
    mDegree = rhs.degree();
}

Monomial::~Monomial()
{
}

bool
    Monomial::isZero() const
{
    return mCoeffR == 0 && mCoeffI == 0;
}

bool
    Monomial::isReal() const
{
    return mCoeffI == 0;
}

bool
    Monomial::isImag() const
{
    return mCoeffR == 0;
}



Monomial
    Monomial::operator-()
{
    return Monomial(-mCoeffR, -mCoeffI, mDegree);
}

Monomial&
    Monomial::operator*=(const Monomial & other)
{
    mpq_class tmp;

    tmp = mCoeffR * other.mCoeffR - mCoeffI * other.mCoeffI;
    mCoeffI = mCoeffI * other.mCoeffR + mCoeffR * other.mCoeffI;
    mCoeffR = tmp;
    mDegree += other.mDegree;

    return *this;
}

Monomial
    Monomial::operator*(const Monomial & other) const
{
    Monomial result = *this;
    result *= other;
    return result;
}

namespace mps {
    namespace formal
    {
        char* get_string(const mps::formal::Monomial& m)
        {
            char* bufR, * bufI, * buf, * newbuf;
            int lenR, lenI;
            std::string strR, strI;

            if (m.coefficientReal() == 0)
            {
                if (m.coefficientImag() == 0)
                {
                    buf = (char*)mps_malloc(1);
                    buf[0] = 0;
                }
                else
                {
                    strI = m.coefficientImag().get_str();
                    lenI = (int)strI.length();
                    buf = (char*)mps_malloc(lenI + 1 + 1);
                    strcpy(buf, strI.c_str());
                    strcat(buf, "i");
                }
            }
            else        // there is a real coefficient
            {
                if (m.coefficientImag() == 0)
                {
                    strR = m.coefficientReal().get_str();
                    lenR = (int)strR.length();
                    buf = (char*)mps_malloc(lenR + 1);
                    strcpy(buf, strR.c_str());
                }
                else
                {
                    strR = m.coefficientReal().get_str();
                    lenR = (int)strR.length();
                    bufR = (char*)mps_malloc(lenR + 1);
                    strI = m.coefficientImag().get_str();
                    lenI = (int)strI.length();
                    bufI = (char*)mps_malloc(lenI + 1 + 1);
                    buf = (char*)mps_malloc(lenR + lenI + 2 + 1 + 1);
                    sprintf(buf, "(%s%s%s)", bufR, m.coefficientImag() > 0 ? "+" : "", bufI);
                    mps_free(bufR);
                    mps_free(bufI);
                }
            }

            if (buf[0])
            {
                switch (m.degree())
                {
                case 0:
                    break;
                case 1:
                {
                    newbuf = (char*)mps_malloc((int)strlen(buf) + 1 + 1);
                    strcpy(newbuf, buf);
                    strcat(newbuf, "x");
                    mps_free(buf);
                    buf = newbuf;
                }
                break;
                default:
                    newbuf = (char*)mps_malloc((int)strlen(buf) + 2 + 20 + 1);
                    strcpy(newbuf, buf);
                    strcat(newbuf, "x^");
                    sprintf(&newbuf[(int)strlen(newbuf)], "%d", m.degree());
                    mps_free(buf);
                    buf = newbuf;
                    break;
                }
            }

            return buf;
        }
    }
}

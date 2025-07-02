def gf_mult_no_table(x, y, prim=0x11d):
    """Multiply two elements in GF(2^8) without lookup tables.
    
    Multiplication is performed as polynomial multiplication 
    modulo the primitive polynomial `prim`.
    """

    def bit_mult(a, b):
        """Bitwise multiply integers without carry.
        
        Equivalent to polynomial multiplication over GF(2).
        """
        result = 0
        i = 0
        while (b >> i) > 0:
            if b & (1 << i):
                result ^= a << i
            i += 1
        return result

    def bit_length(n):
        """Return the number of bits required to represent n."""
        bits = 0
        while n >> bits:
            bits += 1
        return bits

    def bit_div(dividend, divisor):
        """Polynomial division over GF(2).
        
        Returns the remainder of dividend divided by divisor.
        """
        len_dividend = bit_length(dividend)
        len_divisor = bit_length(divisor)
        if len_dividend < len_divisor:
            return dividend
        
        for i in range(len_dividend - len_divisor, -1, -1):
            if dividend & (1 << (i + len_divisor - 1)):
                dividend ^= divisor << i
        return dividend

    # Multiply as polynomials over GF(2)
    product = bit_mult(x, y)
    # Reduce modulo primitive polynomial if specified
    if prim > 0:
        product = bit_div(product, prim)

    return product


# Global lookup tables for GF(2^8) logarithm and exponentiation
gf_exp = [0] * 512
gf_log = [0] * 256


def create_tables(prim=0x11d):
    """Create logarithm and exponentiation tables for GF(2^8)."""

    global gf_exp, gf_log
    gf_exp = [0] * 512
    gf_log = [0] * 256

    x = 1
    for i in range(255):
        gf_exp[i] = x
        gf_log[x] = i
        x = gf_mult_no_table(x, 2, prim)

    # Extend gf_exp table to avoid modulus operations in multiplication
    for i in range(255, 512):
        gf_exp[i] = gf_exp[i - 255]

    return gf_log, gf_exp


def gf_mul(x, y):
    """Multiply two GF(2^8) elements using log/exp tables."""
    if x == 0 or y == 0:
        return 0
    return gf_exp[gf_log[x] + gf_log[y]]


def gf_div(x, y):
    """Divide two GF(2^8) elements using log/exp tables."""
    if y == 0:
        raise ZeroDivisionError("Division by zero in GF(2^8)")
    if x == 0:
        return 0
    return gf_exp[(gf_log[x] + 255 - gf_log[y]) % 255]


def gf_pow(x, power):
    """Raise GF(2^8) element x to a power."""
    return gf_exp[(gf_log[x] * power) % 255]


def gf_inverse(x):
    """Calculate multiplicative inverse of x in GF(2^8)."""
    return gf_div(1, x)


def gf_poly_num_mul(poly, x):
    """Multiply polynomial by scalar x in GF(2^8)."""
    return [gf_mul(coeff, x) for coeff in poly]


def gf_poly_add(p, q):
    """Add two polynomials in GF(2^8)."""
    max_len = max(len(p), len(q))
    result = [0] * max_len
    # Copy p into result aligned to the right
    for i in range(len(p)):
        result[i + max_len - len(p)] = p[i]
    # XOR q into result aligned to the right
    for i in range(len(q)):
        result[i + max_len - len(q)] ^= q[i]
    return result


def gf_poly_mul(p, q):
    """Multiply two polynomials in GF(2^8)."""
    result = [0] * (len(p) + len(q) - 1)
    for j in range(len(q)):
        for i in range(len(p)):
            result[i + j] ^= gf_mul(p[i], q[j])
    return result


def gf_poly_point_eval(poly, x):
    """Evaluate polynomial at point x in GF(2^8)."""
    y = poly[0]
    for i in range(1, len(poly)):
        y = gf_mul(y, x) ^ poly[i]
    return y


def gf_poly_div(dividend, divisor):
    """Divide two polynomials in GF(2^8).
    
    Returns tuple (quotient, remainder).
    """
    msg_out = list(dividend)
    for i in range(len(dividend) - len(divisor) + 1):
        coef = msg_out[i]
        if coef != 0:
            for j in range(1, len(divisor)):
                if divisor[j] != 0:
                    msg_out[i + j] ^= gf_mul(divisor[j], coef)

    separator = -(len(divisor) - 1)
    return msg_out[:separator], msg_out[separator:]


def gf_add(x, y):
    """Addition in GF(2^8) is just XOR."""
    return x ^ y

# Reed-Solomon Error Correction Code Implementation

This repository contains an implementation of Reed-Solomon (RS) error correction codes over GF(2^8) with support for error and erasure correction. The code includes polynomial arithmetic, syndrome calculation, Forney algorithm for erasure handling, and Berlekamp-Massey algorithm for error locator polynomial calculation.

## Features

- Generate Reed-Solomon generator polynomial
- Encode messages with RS code
- Calculate syndromes and Forney syndromes for erasures
- Find error locator polynomial using Berlekamp-Massey algorithm
- Locate and correct errors and erasures
- Works with arbitrary byte sequences (supports UTF-8 encoding)
- Can correct errors and erasures as long as `2*errors + erasures <= r`

## Usage

### Requirements

- Python 3.x
- `GFmath.py` module (contains finite field arithmetic over GF(2^8))

# GF(2^8) Arithmetic Library in Python

This library provides functions for arithmetic operations in the finite field GF(2^8) without using any external libraries.

## Features

- Multiplication in GF(2^8) without lookup tables (bitwise polynomial multiplication)
- Creation of logarithm and exponentiation tables for fast multiplication and division
- Multiplication, division, exponentiation, and inversion in GF(2^8) using tables
- Polynomial operations: addition, multiplication, division, evaluation
- All operations implemented from scratch, no third-party dependencies

## Usage

1. Initialize tables with the primitive polynomial (default `0x11d`):

```python
    n = 20 # n is the total number of characters in the code word (message + verification characters)
    message = "BÏÕHÄzÃRÐ"
    k = len(message)
    r = n - k

    create_tables()

    print(f"Original message:       \033[95m{message}\033[0m")
    ascii_msg = [ord(c) for c in message]
    print(f"Message in ASCII:       \033[92m{ascii_msg}\033[0m")

    encoded = encode_message(ascii_msg, r)
    print(f"Encoded message:        \033[92m{encoded}\033[0m")

    # Introduce some errors and erasures
    encoded[0] = 12
    encoded[1] = 101
    encoded[6] = 15
    encoded[3] = 8
    print(f"Corrupted message:      \033[91m{encoded}\033[0m")

    corrected_msg, corrected_ecc = correct_message(encoded, r, erase_pos=[7])
    print(f"Corrected ASCII:        \033[92m{corrected_msg + corrected_ecc}\033[0m")
    print(f"Recovered message:      \033[96m{''.join([chr(c) for c in corrected_msg])}\033[0m")


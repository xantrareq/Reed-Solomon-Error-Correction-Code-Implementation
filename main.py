from GFmath import *

# Generate the generator polynomial of degree r
def generate_generator_poly(r):
    g = [1]
    for i in range(r):
        g = gf_poly_mul(g, [1, gf_pow(2, i)])
    return g

# Encode a message with r parity symbols
def encode_message(msg_in, r):
    gen = generate_generator_poly(r)
    _, remainder = gf_poly_div(msg_in + [0] * (len(gen) - 1), gen)
    return msg_in + remainder

# Calculate syndrome vector
def calculate_syndromes(msg, r):
    syndromes = [gf_poly_point_eval(msg, gf_pow(2, i)) for i in range(r)]
    return [0] + syndromes  # Add a leading zero for alignment

# Forney's syndrome calculation (used with erasures)
def calculate_forney_syndromes(syndromes, erase_pos, n):
    erase_pos = [n - 1 - p for p in erase_pos]  # reverse indices
    forney_synd = list(syndromes[1:])
    for i in range(len(erase_pos)):
        x = gf_pow(2, erase_pos[i])
        for j in range(len(forney_synd) - 1):
            forney_synd[j] = gf_mul(forney_synd[j], x) ^ forney_synd[j + 1]
    return forney_synd

# Berlekamp-Massey algorithm for error locator polynomial
def berlekamp_massey(synd, r, erase_loc=None, erase_count=0):
    err_loc = list(erase_loc) if erase_loc else [1]
    old_loc = list(erase_loc) if erase_loc else [1]

    for i in range(r - erase_count):
        K = i + (len(synd) - r) + (erase_count if erase_loc else 0)
        delta = synd[K]
        for j in range(1, len(err_loc)):
            delta ^= gf_mul(err_loc[-(j + 1)], synd[K - j])

        old_loc.append(0)

        if delta != 0:
            if len(old_loc) > len(err_loc):
                new_loc = gf_poly_num_mul(old_loc, delta)
                old_loc = gf_poly_num_mul(err_loc, gf_inverse(delta))
                err_loc = new_loc
            err_loc = gf_poly_add(err_loc, gf_poly_num_mul(old_loc, delta))

    while err_loc and err_loc[0] == 0:
        del err_loc[0]

    errs = len(err_loc) - 1
    if (errs - erase_count) * 2 + erase_count > r:
        print("Too many errors to correct")
        exit(0)

    return err_loc

# Find error positions from error locator polynomial
def find_error_positions(err_loc, n):
    err_count = len(err_loc) - 1
    positions = [n - 1 - i for i in range(n) if gf_poly_point_eval(err_loc, gf_pow(2, i)) == 0]
    if len(positions) != err_count:
        print("Mismatch in error count")
        exit(0)
    return positions

# Compute error locator polynomial from erasure positions
def compute_erasure_locator_poly(erase_positions):
    loc = [1]
    for i in erase_positions:
        loc = gf_poly_mul(loc, gf_poly_add([1], [gf_pow(2, i), 0]))
    return loc

# Compute error evaluator polynomial
def compute_error_evaluator(synd, err_loc, r):
    _, remainder = gf_poly_div(gf_poly_mul(synd, err_loc), [1] + [0] * (r + 1))
    return remainder

# Forney's algorithm to correct errors
def forney_correct_errors(msg, synd, error_positions):
    n = len(msg)
    coef_pos = [n - 1 - p for p in error_positions]
    err_loc = compute_erasure_locator_poly(coef_pos)
    err_eval = compute_error_evaluator(synd[::-1], err_loc, len(err_loc) - 1)[::-1]

    X = [gf_pow(2, -(255 - i)) for i in coef_pos]
    E = [0] * n

    for i, Xi in enumerate(X):
        Xi_inv = gf_inverse(Xi)
        err_loc_derivative = 1
        for j in range(len(X)):
            if j != i:
                err_loc_derivative = gf_mul(err_loc_derivative, gf_add(1, gf_mul(Xi_inv, X[j])))

        y = gf_poly_point_eval(err_eval[::-1], Xi_inv)
        y = gf_mul(Xi, y)

        if err_loc_derivative == 0:
            return None

        magnitude = gf_div(y, err_loc_derivative)
        E[error_positions[i]] = magnitude

    return gf_poly_add(msg, E)

# High-level message correction function
def correct_message(msg, r, erase_pos=None):
    if len(msg) > 255:
        print(f"Message too long ({len(msg)} > 255)")
        exit(0)

    msg_copy = list(msg)
    erase_pos = erase_pos or []

    for e in erase_pos:
        msg_copy[e] = 0

    if len(erase_pos) > r:
        print("Too many erasures")
        exit(0)

    synd = calculate_syndromes(msg_copy, r)
    if max(synd) == 0:
        return msg_copy[:-r], msg_copy[-r:]

    forney_synd = calculate_forney_syndromes(synd, erase_pos, len(msg_copy))
    err_loc = berlekamp_massey(forney_synd, r, erase_count=len(erase_pos))
    err_pos = find_error_positions(err_loc[::-1], len(msg_copy))

    if err_pos is None:
        print("Unable to locate errors")
        exit(0)

    msg_corrected = forney_correct_errors(msg_copy, synd, erase_pos + err_pos)

    synd = calculate_syndromes(msg_corrected, r)
    if max(synd) > 0:
        print("Correction failed")
        exit(0)

    return msg_corrected[:-r], msg_corrected[-r:]

# Example usage
if __name__ == "__main__":
    
    # Original message and parameters
    message = "BÏÕHÄzÃRÐ"
    n = 20 # Total codeword length (message + parity)
    k = len(message) # Message length
    r = n - k  # Number of parity symbols (corrects up to 4 errors)


    # Initialize Galois Field tables
    create_tables()

    # Convert string message to ASCII codes
    ascii_msg = [ord(c) for c in message]

    # Encode the message by appending parity symbols
    encoded = encode_message(ascii_msg, r)
    encoded_no_corr = encoded.copy()

    # Introduce errors in the encoded message
    encoded[0] = 12
    encoded[1] = 101
    encoded[6] = 15
    encoded[3] = 8

    # Specify known erasure positions (indices where errors are known)
    erase_pos = [7]

    # Attempt error correction with erasures
    corrected_msg, corrected_ecc = correct_message(encoded, r, erase_pos=erase_pos)

    # Output results
    print(f"Original message:       \033[95m{message}\033[0m")
    print(f"Encoded message:        \033[92m{encoded_no_corr}\033[0m")
    print(f"Corrupted message:      \033[91m{encoded}\033[0m")
    print(f"Corrected ASCII:        \033[92m{corrected_msg + corrected_ecc}\033[0m")
    print(f"Recovered message:      \033[96m{''.join([chr(c) for c in corrected_msg])}\033[0m")



